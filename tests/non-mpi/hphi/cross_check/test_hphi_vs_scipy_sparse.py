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
Cross-validation of the HPhi impurity solver against the scipy/sparse solver.

Both are exact-diagonalization (ED) solvers, so for the same impurity model
they must yield the same self-energy.  This test catches a class of bugs that
single-solver unit tests cannot, e.g. a wrong off-diagonal Green's-function
reconstruction *or* an insufficient number of computed eigenstates (``exct``).

Physics being checked (2-orbital Kanamori atom, Hubbard-I / ``n_bath = 0``):
  * the two orbitals are degenerate, so Sigma must be orbital-symmetric
    (Sigma[0,0] == Sigma[1,1]) with vanishing off-diagonal (Sigma[0,1] == 0);
  * the 1-electron ground state is 4-fold degenerate, so HPhi needs
    ``exct`` large enough to cover the whole low-energy multiplet at finite T
    (here the full 16-dimensional space) -- otherwise the thermal trace is
    incomplete and the orbital symmetry is spuriously broken.

The HPhi part runs only when an HPhi executable is available; set the path via
the ``DCORE_HPHI_EXEC`` environment variable (or have ``HPhi`` on ``PATH``).
Without it the whole module is skipped, so this is safe in CI.
"""

import os
import shutil

import numpy
import pytest


def _hphi_exec():
    """Resolve an HPhi executable from the environment, else None."""
    path = os.environ.get("DCORE_HPHI_EXEC")
    if path and os.path.isfile(path) and os.access(path, os.X_OK):
        return path
    return shutil.which("HPhi")


HPHI_EXEC = _hphi_exec()

pytestmark = pytest.mark.skipif(
    HPHI_EXEC is None,
    reason="HPhi executable not found; set DCORE_HPHI_EXEC to enable the cross-check",
)


# Shared 2-orbital model.  Only the [impurity_solver] block differs between runs.
_MODEL_TEMPLATE = """\
[model]
seedname = {seed}
lattice = square
norb = 2
nelec = 1.0
t = -1.0
kanamori = [(4.0, 0.8, 0.0)]
nk = 8

[system]
T = 0.1
n_iw = 500
fix_mu = True
mu = 0.5

[impurity_solver]
{solver_block}

[control]
max_step = 1
sigma_mix = 1.0
"""


def _run_solver(work_dir, seed, solver_block):
    """Run dcore_pre + dcore for one solver and return Sigma_iw as (n_iw, norb, norb)."""
    from dcore.dcore_pre import dcore_pre
    from dcore.dcore import dcore
    import h5py

    ini = os.path.join(work_dir, seed + ".ini")
    with open(ini, "w") as f:
        f.write(_MODEL_TEMPLATE.format(seed=seed, solver_block=solver_block))

    cwd = os.getcwd()
    os.chdir(work_dir)
    try:
        dcore_pre(ini)
        dcore(ini)
        with h5py.File(seed + ".out.h5", "r") as h:
            data = h["dmft_out"]["Sigma_iw"]["ite1"]["sh0"]["up"]["data"][()]
    finally:
        os.chdir(cwd)
    return data[..., 0] + 1j * data[..., 1]


@pytest.fixture(scope="module")
def sigma_scipy(tmp_path_factory):
    """Self-energy from the scipy/sparse ED solver -- the ground-truth reference."""
    pytest.importorskip("scipy")
    work_dir = str(tmp_path_factory.mktemp("scipy_sparse"))
    solver_block = "name = scipy/sparse\nn_bath{int} = 0"
    return _run_solver(work_dir, "scipy_ref", solver_block)


def test_scipy_sparse_is_orbital_symmetric(sigma_scipy):
    """Sanity check on the ED reference: degenerate orbitals -> symmetric, diagonal Sigma."""
    n2 = sigma_scipy.shape[0] // 2
    low = slice(n2 - 50, n2 + 50)
    diag_asym = numpy.abs(sigma_scipy[low, 0, 0] - sigma_scipy[low, 1, 1]).max()
    offdiag = numpy.abs(sigma_scipy[:, 0, 1]).max()
    assert diag_asym < 1e-4, f"reference diagonal not symmetric: {diag_asym}"
    assert offdiag < 1e-4, f"reference off-diagonal not zero: {offdiag}"


def test_hphi_matches_scipy_sparse(sigma_scipy, tmp_path):
    """
    HPhi with enough eigenstates must reproduce the scipy/sparse ED self-energy.

    exct = 16 spans the full Hilbert space of the 2-site (2-orbital, n_bath=0)
    problem, so the finite-T thermal trace is complete and the 4-fold degenerate
    ground multiplet is fully included.
    """
    solver_block = (
        "name = HPhi\n"
        f"exec_path{{str}} = {HPHI_EXEC}\n"
        "n_bath{int} = 0\n"
        "exct{int} = 16\n"
        "np{int} = 1"
    )
    sigma_hphi = _run_solver(str(tmp_path), "hphi_run", solver_block)

    assert sigma_hphi.shape == sigma_scipy.shape

    n2 = sigma_hphi.shape[0] // 2
    low = slice(n2 - 50, n2 + 50)

    # (1) HPhi must keep the orbital symmetry it should have (broken when exct=1).
    diag_asym = numpy.abs(sigma_hphi[low, 0, 0] - sigma_hphi[low, 1, 1]).max()
    assert diag_asym < 5e-3, f"HPhi diagonal not orbital-symmetric: {diag_asym}"

    # (2) HPhi off-diagonal must vanish for this model.
    offdiag = numpy.abs(sigma_hphi[:, 0, 1]).max()
    assert offdiag < 5e-3, f"HPhi spurious off-diagonal: {offdiag}"

    # (3) HPhi diagonal must match the scipy/sparse ED reference.
    diff = numpy.abs(sigma_hphi[low, 0, 0] - sigma_scipy[low, 0, 0]).max()
    assert diff < 5e-3, f"HPhi vs scipy/sparse diagonal mismatch: {diff}"
