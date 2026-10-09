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
import numpy as np
from dcore import gpu


def test_get_backend_cpu_default():
    xp, active = gpu.get_backend(False)
    assert xp is np
    assert active is False


def test_get_backend_gpu_missing_falls_back(monkeypatch):
    # Force the cupy import to fail -> graceful numpy fallback + warning.
    monkeypatch.setattr(gpu, "_import_cupy",
                        lambda: (_ for _ in ()).throw(ImportError("no cupy")))
    import warnings
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        xp, active = gpu.get_backend(True)
    assert xp is np and active is False
    assert any("cupy" in str(x.message).lower() for x in w)


def test_to_host_identity_for_numpy():
    a = np.arange(3)
    assert gpu.to_host(a) is a



def test_unusable_cuda_context_falls_back(monkeypatch):
    import types
    import pytest

    def allocate(*args):
        raise RuntimeError("CUDA context unavailable")

    cupy = types.SimpleNamespace(
        cuda=types.SimpleNamespace(runtime=types.SimpleNamespace(getDeviceCount=lambda: 1)),
        zeros=allocate)
    monkeypatch.setattr(gpu, '_import_cupy', lambda: cupy)
    with pytest.warns(UserWarning, match="CUDA context unavailable"):
        xp, active = gpu.get_backend(True)
    assert xp is np and active is False


def test_missing_cusolver_falls_back(monkeypatch):
    import types
    import pytest

    def eigh(*args):
        raise ImportError("libcusolver.so.11: cannot open shared object file")

    cupy = types.SimpleNamespace(
        cuda=types.SimpleNamespace(runtime=types.SimpleNamespace(getDeviceCount=lambda: 1)),
        zeros=np.zeros, eye=np.eye, linalg=types.SimpleNamespace(eigh=eigh))
    monkeypatch.setattr(gpu, '_import_cupy', lambda: cupy)
    with pytest.warns(UserWarning, match="libcusolver"):
        xp, active = gpu.get_backend(True)
    assert xp is np and active is False
