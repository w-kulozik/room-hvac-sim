"""Type aliases shared by the package."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]
"""Array of float64 values (temperatures, flows, powers, ...)."""

BoolArray = NDArray[np.bool_]
"""Array of booleans (masks)."""

FloatOrArray = float | FloatArray
"""A scalar or an array; functions taking it apply element-wise."""
