/*
 * Fixed-point PI with sigma-delta low-load handling; see pi_sigma_delta.h.
 */
#include "pi_sigma_delta.h"

int64_t fw_rshift_round(int64_t x, int32_t k)
{
    /* x / 2^k rounded to nearest, ties away from zero (k >= 1) */
    int64_t half = (int64_t)1 << (k - 1);
    return x >= 0 ? (x + half) >> k : -((-x + half) >> k);
}

int32_t fw_mul_q24(int32_t gain_q24, int32_t x_q16)
{
    /* Q8.24 x Q16.16 -> Q16.16 */
    return (int32_t)fw_rshift_round((int64_t)gain_q24 * (int64_t)x_q16, 24);
}

int32_t fw_size_of_state(void)
{
    return (int32_t)sizeof(fw_state);
}

void fw_init(fw_state *st, const fw_params *p)
{
    st->p = *p;
    st->integrator_q32 = 0;
    st->acc_q16 = 0;
    st->s_prev_q16 = 0;
    st->on = 0;
    st->switched = 0;
    st->t_switch_s = 0;
    st->starts = 0;
    st->demand_q16 = 0;
}

static int can_start(const fw_state *st, int64_t t_s)
{
    return !st->on && (!st->switched || t_s - st->t_switch_s >= st->p.min_off_s);
}

static int can_stop(const fw_state *st, int64_t t_s)
{
    return st->on && (!st->switched || t_s - st->t_switch_s >= st->p.min_on_s);
}

static void set_on(fw_state *st, int64_t t_s, int on)
{
    if (on != st->on) {
        st->on = on;
        st->switched = 1;
        st->t_switch_s = t_s;
        st->starts += on;
    }
}

int32_t fw_step(fw_state *st, int64_t t_s, int32_t t_air_q16, int32_t t_set_q16)
{
    const fw_params *p = &st->p;
    int32_t e = t_air_q16 - t_set_q16;                       /* cooling: warm room, demand */
    int32_t v = fw_mul_q24(p->kp_q24, e) + (int32_t)fw_rshift_round(st->integrator_q32, 16);
    int32_t s, u;
    int64_t thr = p->threshold_q16;

    if (v >= p->s_min_q16) {
        /* modulating range: follow the demand */
        st->acc_q16 = 0;
        if (st->on || can_start(st, t_s))
            set_on(st, t_s, 1);
        s = st->on ? (v < p->s_max_q16 ? v : p->s_max_q16) : 0;
        u = s;
    } else if (v <= 0) {
        /* no demand: stop as soon as the guard allows */
        st->acc_q16 = 0;
        if (st->on && can_stop(st, t_s))
            set_on(st, t_s, 0);
        s = st->on ? p->s_min_q16 : 0;
        u = 0;
    } else {
        /* low load: sigma-delta modulation between 0 and s_min */
        st->acc_q16 += (int64_t)(v - st->s_prev_q16) * p->ts_s;
        if (st->acc_q16 > thr)
            st->acc_q16 = thr;
        if (st->acc_q16 < -thr)
            st->acc_q16 = -thr;
        if (!st->on && st->acc_q16 >= thr && can_start(st, t_s))
            set_on(st, t_s, 1);
        else if (st->on && st->acc_q16 <= -thr && can_stop(st, t_s))
            set_on(st, t_s, 0);
        s = st->on ? p->s_min_q16 : 0;
        u = v;                                               /* met on average */
    }
    st->s_prev_q16 = s;
    st->demand_q16 = v;

    /* integrator, Q32.32: Q8.24 x Q16.16 = Q40, shifted by 8 */
    st->integrator_q32 += fw_rshift_round((int64_t)p->ki_q24 * e, 8)
                        + fw_rshift_round((int64_t)p->kt_q24 * (int64_t)(u - v), 8);
    return s;
}
