# Firmware

The temperature PI with back-calculation anti-windup, the sigma-delta low-load handler and the
compressor guard (handler C of the studies, cooling) in fixed-point C99, as it would run on the
unit's microcontroller. Design and number formats: decision D59 in
[docs/decisions.md](../docs/decisions.md); results: D60 and `figures/firmware_sil.md`.

- `pi_sigma_delta.h`, `pi_sigma_delta.c`: no floating point, no dynamic memory, one state
  structure.
- Any C99 compiler builds it. `roomsim.controllers.firmware.build_library` uses `$CC` if set,
  otherwise the zig toolchain from PyPI (`pip install ziglang`, part of the `dev` extras), and
  writes the shared library to `firmware/build/`.
- `roomsim.controllers.firmware.FirmwareController` loads it with `ctypes`, so the simulator
  runs the firmware in the loop in place of the Python controller.
- Checks: `tests/test_firmware.py`; measured comparison: `python scripts/compare_firmware.py`.
