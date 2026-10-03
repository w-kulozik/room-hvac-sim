/*
 * Temperature PI with back-calculation anti-windup, sigma-delta low-load handling and a
 * compressor guard, in fixed point (decision D59 in docs/decisions.md).
 *
 * The algorithm is that of roomsim.controllers (PIUnitController with SigmaDeltaLowLoad and
 * CompressorGuard, cooling): every Ts seconds
 *     e = T - T_set,   v = Kp e + I,
 *     s = sigma-delta(v)               compressor command, 0 or in [s_min, s_max]
 *     I += Kp Ts/Ti e + Ts/Tt (u - v)  u = demand met (s above s_min, v while cycling)
 *
 * Number formats (no floating point):
 *     Q16.16 in int32_t   temperatures (degC), errors (K), demands and speeds (-)
 *     Q8.24  in int32_t   gains Kp (1/K), Kp Ts/Ti (1/K), Ts/Tt (-)
 *     Q32.32 in int64_t   PI integrator
 *     Q16.16 in int64_t   sigma-delta accumulator (seconds)
 *     int64_t             times, whole seconds
 * Products are formed in 64 bits and rounded to nearest, ties away from zero.
 */
#ifndef PI_SIGMA_DELTA_H
#define PI_SIGMA_DELTA_H

#include <stdint.h>

#ifdef _WIN32
#define FW_API __declspec(dllexport)
#else
#define FW_API
#endif

#define FW_Q16_ONE 65536

typedef struct {
    int32_t kp_q24;       /* Kp, demand per K */
    int32_t ki_q24;       /* Kp Ts / Ti, demand per K per sample */
    int32_t kt_q24;       /* Ts / Tt, back-calculation gain per sample */
    int32_t s_min_q16;    /* minimum compressor speed */
    int32_t s_max_q16;    /* maximum compressor speed */
    int32_t threshold_q16; /* sigma-delta switching threshold, s */
    int32_t ts_s;         /* sampling period, s */
    int32_t min_on_s;     /* minimum run time, s */
    int32_t min_off_s;    /* minimum stop time, s */
} fw_params;

typedef struct {
    fw_params p;
    int64_t integrator_q32; /* PI integrator I */
    int64_t acc_q16;        /* sigma-delta accumulator, s */
    int32_t s_prev_q16;     /* command of the previous sample */
    int32_t on;             /* compressor running */
    int32_t switched;       /* a start or stop has happened (t_switch is valid) */
    int64_t t_switch_s;     /* time of the last start or stop, s */
    int32_t starts;         /* number of starts */
    int32_t demand_q16;     /* unsaturated demand v of the last sample */
} fw_state;

/* Set the parameters and reset the state (compressor off, integrator zero). */
FW_API void fw_init(fw_state *st, const fw_params *p);

/* One controller sample at time t_s: returns the compressor command (Q16.16; 0 = off). */
FW_API int32_t fw_step(fw_state *st, int64_t t_s, int32_t t_air_q16, int32_t t_set_q16);

/* Primitives, exported for testing. */
FW_API int64_t fw_rshift_round(int64_t x, int32_t k);
FW_API int32_t fw_mul_q24(int32_t gain_q24, int32_t x_q16);
FW_API int32_t fw_size_of_state(void);

#endif
