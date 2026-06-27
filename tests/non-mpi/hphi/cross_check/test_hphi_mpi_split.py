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
"""
Integration test for the HPhi Gf two-level parallelism (needs MPI).

It checks *split invariance*: running each HPhi of the Gf step under
``mpirun -np 4`` (n_procs_per_hphi = 4, n_inner = 4) must give the same
self-energy as the plain serial layout (n_procs_per_hphi = 1). The physics is
identical, only the parallel decomposition of the Gf step differs.

Because it compares two layouts of the *same* computation, a small ``exct`` is
enough (we are not checking thermal convergence here), which keeps it fast.

Requires three things, otherwise the test is skipped (safe in CI):
  * DCORE_HPHI_MPI_EXEC : an MPI build of HPhi (used with the launcher),
  * DCORE_HPHI_EXEC     : a (serial-capable) HPhi for the reference run,
  * mpirun on PATH.
"""

import os
import shutil

import numpy
import pytest


def _exe(env_var):
    p = os.environ.get(env_var)
    if p and os.path.isfile(p) and os.access(p, os.X_OK):
        return p
    return None


HPHI_MPI = _exe("DCORE_HPHI_MPI_EXEC")
HPHI_SERIAL = _exe("DCORE_HPHI_EXEC") or shutil.which("HPhi")
MPIRUN = shutil.which("mpirun")

pytestmark = pytest.mark.skipif(
    not (HPHI_MPI and HPHI_SERIAL and MPIRUN),
    reason="needs DCORE_HPHI_MPI_EXEC + DCORE_HPHI_EXEC + mpirun",
)


# Complex off-diagonal crystal field: a genuinely anti-symmetric off-diagonal Gf,
# so the parallel path is exercised on the full reconstruction, not just diagonals.
_CRYSTAL_FIELD = """\
# sp o1 o2 re im
0 0 1 0.3  0.2
0 1 0 0.3 -0.2
1 0 1 0.3  0.2
1 1 0 0.3 -0.2
"""

_MODEL = """\
[model]
seedname = {seed}
lattice = square
norb = 2
nelec = 1.0
t = -1.0
kanamori = [(4.0, 0.8, 0.0)]
nk = 8
local_potential_matrix = {{0: 'cf.in'}}
local_potential_factor = 1.0

[mpi]
command = mpirun -np #

[system]
T = 0.1
n_iw = 100
fix_mu = True
mu = 0.5

[impurity_solver]
name = HPhi
exec_path{{str}} = {exec_path}
n_bath{{int}} = 0
exct{{int}} = 4
n_procs_per_hphi{{int}} = {n_procs_per_hphi}

[control]
max_step = 1
sigma_mix = 1.0
"""


def _run(work_dir, seed, exec_path, n_procs_per_hphi, np_total):
    from dcore.dcore_pre import dcore_pre
    from dcore.dcore import dcore
    import h5py

    with open(os.path.join(work_dir, "cf.in"), "w") as f:
        f.write(_CRYSTAL_FIELD)
    ini = os.path.join(work_dir, seed + ".ini")
    with open(ini, "w") as f:
        f.write(_MODEL.format(seed=seed, exec_path=exec_path,
                              n_procs_per_hphi=n_procs_per_hphi))

    cwd = os.getcwd()
    os.chdir(work_dir)
    try:
        dcore_pre(ini)
        dcore(ini, np_total)
        with h5py.File(seed + ".out.h5", "r") as h:
            d = h["dmft_out"]["Sigma_iw"]["ite1"]["sh0"]["up"]["data"][()]
    finally:
        os.chdir(cwd)
    return d[..., 0] + 1j * d[..., 1]


def test_gf_split_is_invariant(tmp_path):
    # serial layout (n_inner = 1): the reference
    sigma_serial = _run(str(tmp_path), "serial_ref", HPHI_SERIAL,
                        n_procs_per_hphi=1, np_total=1)
    # two-level layout: each HPhi of the Gf step runs under mpirun -np 4
    sigma_mpi = _run(str(tmp_path), "mpi_split", HPHI_MPI,
                     n_procs_per_hphi=4, np_total=4)

    assert sigma_mpi.shape == sigma_serial.shape
    # genuine off-diagonal, so the parallel path really exercises the reconstruction
    assert numpy.abs(sigma_serial[:, 0, 1]).max() > 0.1
    # the parallel decomposition must not change the answer
    diff = numpy.abs(sigma_mpi - sigma_serial).max()
    assert diff < 5e-3, f"Gf split changed the self-energy: {diff}"
