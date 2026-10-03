"""The C firmware version of the PI with sigma-delta low-load handling, run in the loop (D59).

``firmware/pi_sigma_delta.c`` implements the compressor part of ``PIUnitController`` with
``SigmaDeltaLowLoad`` in fixed point. ``build_library`` compiles it into a shared library with
any C99 compiler (``$CC``, or the zig toolchain from PyPI); ``FirmwareController`` loads it with
``ctypes`` and has the interface of the Python controllers, so the simulator runs the firmware
in place of the Python controller (software in the loop). The set point schedule, the fan and
the ventilation stay in Python.

The firmware sources are found relative to a source checkout of the repository.
"""

from __future__ import annotations

import ctypes
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

from .._typing import FloatArray
from ..plant import IDX_Y_TA
from .base import fan_command, set_point_at
from .lowload import SigmaDeltaLowLoad
from .pi_unit import PIUnitController

FIRMWARE_DIR = Path(__file__).resolve().parents[3] / "firmware"
SOURCE = FIRMWARE_DIR / "pi_sigma_delta.c"
Q16, Q24 = 1 << 16, 1 << 24


class FwParams(ctypes.Structure):
    """``fw_params`` of pi_sigma_delta.h."""

    _fields_ = [("kp_q24", ctypes.c_int32), ("ki_q24", ctypes.c_int32),  # noqa: RUF012  ctypes
                ("kt_q24", ctypes.c_int32), ("s_min_q16", ctypes.c_int32),
                ("s_max_q16", ctypes.c_int32), ("threshold_q16", ctypes.c_int32),
                ("ts_s", ctypes.c_int32), ("min_on_s", ctypes.c_int32),
                ("min_off_s", ctypes.c_int32)]


class FwState(ctypes.Structure):
    """``fw_state`` of pi_sigma_delta.h."""

    _fields_ = [("p", FwParams), ("integrator_q32", ctypes.c_int64),  # noqa: RUF012  ctypes
                ("acc_q16", ctypes.c_int64), ("s_prev_q16", ctypes.c_int32),
                ("on", ctypes.c_int32), ("switched", ctypes.c_int32),
                ("t_switch_s", ctypes.c_int64), ("starts", ctypes.c_int32),
                ("demand_q16", ctypes.c_int32)]


def build_library(out_dir: Path | None = None) -> Path:
    """Compile the firmware into a shared library, unless an up-to-date one exists.

    Parameters
    ----------
    out_dir : Path, optional
        Output directory; ``firmware/build`` by default.

    Returns
    -------
    Path
        The shared library.

    Raises
    ------
    RuntimeError
        If the compiler fails.
    """
    out_dir = FIRMWARE_DIR / "build" if out_dir is None else out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    lib = out_dir / ("pi_sigma_delta.dll" if sys.platform == "win32" else "pi_sigma_delta.so")
    sources = [SOURCE, SOURCE.with_suffix(".h")]
    if lib.exists() and lib.stat().st_mtime >= max(s.stat().st_mtime for s in sources):
        return lib
    cc = os.environ.get("CC")
    cmd = [cc] if cc else [sys.executable, "-m", "ziglang", "cc"]
    flags = ["-std=c99", "-O2", "-Wall", "-Wextra", "-Werror", "-shared"]
    if sys.platform != "win32":
        flags.append("-fPIC")
    proc = subprocess.run([*cmd, *flags, "-o", str(lib), str(SOURCE)], capture_output=True,
                          text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"building the firmware failed:\n{proc.stdout}{proc.stderr}")
    return lib


def load_library(path: Path | None = None) -> ctypes.CDLL:
    """Load the firmware library and declare its function signatures.

    Parameters
    ----------
    path : Path, optional
        Shared library; built with ``build_library`` if not given.

    Returns
    -------
    ctypes.CDLL
        The library, with argument and result types set.
    """
    lib = ctypes.CDLL(str(build_library() if path is None else path))
    lib.fw_init.argtypes = [ctypes.POINTER(FwState), ctypes.POINTER(FwParams)]
    lib.fw_init.restype = None
    lib.fw_step.argtypes = [ctypes.POINTER(FwState), ctypes.c_int64, ctypes.c_int32,
                            ctypes.c_int32]
    lib.fw_step.restype = ctypes.c_int32
    lib.fw_rshift_round.argtypes = [ctypes.c_int64, ctypes.c_int32]
    lib.fw_rshift_round.restype = ctypes.c_int64
    lib.fw_mul_q24.argtypes = [ctypes.c_int32, ctypes.c_int32]
    lib.fw_mul_q24.restype = ctypes.c_int32
    lib.fw_size_of_state.argtypes = []
    lib.fw_size_of_state.restype = ctypes.c_int32
    if lib.fw_size_of_state() != ctypes.sizeof(FwState):
        raise RuntimeError("the firmware state layout differs from FwState")
    return lib


def to_q16(x: float) -> int:
    """Return x in Q16.16, rounded to nearest."""
    return round(x * Q16)


def to_q24(x: float) -> int:
    """Return x in Q8.24, rounded to nearest."""
    return round(x * Q24)


def firmware_params(controller: PIUnitController) -> FwParams:
    """Convert the tuning of a Python PI with sigma-delta handling to firmware parameters.

    Parameters
    ----------
    controller : PIUnitController
        Cooling PI with a ``SigmaDeltaLowLoad`` handler.

    Returns
    -------
    FwParams
        The same gains, limits and times in the firmware's number formats.

    Raises
    ------
    ValueError
        If the controller is in heating mode or has another low-load handler.
    """
    ll = controller.lowload
    if controller.heating or not isinstance(ll, SigmaDeltaLowLoad):
        raise ValueError("the firmware implements the cooling PI with sigma-delta handling only")
    pi = controller.pi
    Tt = 0.5 * pi.Ti if pi.Tt is None else pi.Tt
    return FwParams(kp_q24=to_q24(pi.Kp), ki_q24=to_q24(pi.Kp * pi.Ts / pi.Ti),
                    kt_q24=to_q24(pi.Ts / Tt), s_min_q16=to_q16(ll.s_min),
                    s_max_q16=to_q16(ll.s_max), threshold_q16=to_q16(ll.threshold),
                    ts_s=round(pi.Ts), min_on_s=round(ll.guard.min_on),
                    min_off_s=round(ll.guard.min_off))


class FirmwareController:
    """The firmware in the loop, with the interface of ``PIUnitController``.

    Parameters
    ----------
    template : PIUnitController
        Python controller whose tuning, set point, fan mode and ventilation are used; its own
        state is not touched.
    lib : ctypes.CDLL, optional
        Loaded firmware library (``load_library``); loaded on demand if not given.

    Attributes
    ----------
    state : FwState
        Firmware state, readable for tests.
    """

    def __init__(self, template: PIUnitController, lib: ctypes.CDLL | None = None) -> None:
        self.lib = load_library() if lib is None else lib
        self.template = template
        self.state = FwState()
        self.lib.fw_init(ctypes.byref(self.state), ctypes.byref(firmware_params(template)))

    @property
    def Ts(self) -> float:
        """Sampling period, s."""
        return float(self.state.p.ts_s)

    def step(self, t: float, y: FloatArray) -> FloatArray:
        """Return the command for one controller sample, computed by the firmware.

        Parameters
        ----------
        t : float
            Current time, s (a whole number of seconds).
        y : ndarray
            Measurements; uses the air temperature y[IDX_Y_TA], degC.

        Returns
        -------
        ndarray, shape (3,)
            Command [s (-), f (-), V_vent (m3/s)].
        """
        c = self.template
        T_set = set_point_at(c.T_set, t)
        s_q16 = self.lib.fw_step(ctypes.byref(self.state), round(t), to_q16(float(y[IDX_Y_TA])),
                                 to_q16(T_set))
        s = s_q16 / Q16
        return np.array([s, fan_command(s, c.fan, c.fan_mode), c.ventilation(t)])
