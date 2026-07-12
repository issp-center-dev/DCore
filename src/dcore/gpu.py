#
# DCore -- Integrated DMFT software for correlated electrons
# Copyright (C) 2017 The University of Tokyo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
"""Array-backend selection (numpy / cupy) for optional GPU acceleration.

CuPy is an optional dependency imported lazily; GPU use is strictly opt-in and
falls back to numpy (with a warning) when CuPy or a CUDA device is unavailable.
"""
import warnings
import numpy as np

_INSTALL_HINT = ("install the CuPy wheel matching your CUDA "
                 "(e.g. 'pip install cupy-cuda12x'); see "
                 "https://docs.cupy.dev/en/stable/install.html")


def _import_cupy():
    import cupy
    return cupy


def get_backend(use_gpu):
    """Return (xp, gpu_active): (cupy, True) if use_gpu and a usable CUDA device
    exists, else (numpy, False) with a one-time fallback warning."""
    if not use_gpu:
        return np, False
    try:
        cupy = _import_cupy()
    except ImportError:
        warnings.warn("gpu=true requested but CuPy is not installed; using the "
                      "numpy (CPU) backend. To enable the GPU, " + _INSTALL_HINT)
        return np, False
    try:
        if cupy.cuda.runtime.getDeviceCount() < 1:
            raise RuntimeError("no CUDA device")
    except Exception as exc:
        warnings.warn("gpu=true requested but no usable CUDA device was found "
                      "({}); using the numpy (CPU) backend.".format(exc))
        return np, False
    return cupy, True


def array_module_of(arr):
    """Return numpy or cupy, whichever owns arr."""
    if type(arr).__module__.split(".")[0] == "cupy":
        return _import_cupy()
    return np


def to_host(arr):
    """Return arr as a numpy array (device->host copy for cupy; identity for numpy)."""
    xp = array_module_of(arr)
    return arr if xp is np else xp.asnumpy(arr)
