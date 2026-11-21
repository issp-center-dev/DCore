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
HPhi-level golden test for the multi-operator finite-T spectrum loop
(``SpectrumNumOp`` / op-inner eigenvector reuse).

When ``SpectrumNumOp = N`` the spectrum calculation reads each eigenvector ONCE
and, for that eigenstate, builds and solves the spectrum for every operator set
op = 0 .. N-1 in turn (op-inner), writing ``zvo_DynamicalGreen_<idx>_<op>.dat``.
This must be numerically identical to running each operator separately with
``SpectrumNumOp = 1`` (the Stage-1a per-operator path): the batching only avoids
redundant eigenvector reads, it must not change the result.

The test drives HPhi directly (standard mode for the eigenvalue step, expert mode
for the spectrum) on a 2-site Hubbard dimer, so it does not depend on the DCore
HPhi driver and validates the HPhi feature in isolation.  It runs only when an
HPhi executable is available (``DCORE_HPHI_EXEC`` or ``HPhi`` on ``PATH``).
"""

import os
import shutil
import subprocess

import numpy
import pytest


def _hphi_exec():
    path = os.environ.get("DCORE_HPHI_EXEC")
    if path and os.path.isfile(path) and os.access(path, os.X_OK):
        return path
    return shutil.which("HPhi")


HPHI_EXEC = _hphi_exec()

pytestmark = pytest.mark.skipif(
    HPHI_EXEC is None,
    reason="HPhi executable not found; set DCORE_HPHI_EXEC to enable the multi-op test",
)

NEXCT = 4  # full Hilbert space of the 2-site, 2Sz=0, nelec=2 Hubbard dimer

_STAN = """\
L = 2
model = "Hubbard"
method = "CG"
lattice = "chain"
t = 1.0
U = 4.0
2Sz = 0
nelec = 2
exct = {nexct}
EigenVecIO = "out"
"""

_CALCMOD_SPEC = """\
CalcType   3
CalcModel   0
ReStart   0
CalcSpec   1
CalcEigenVec   0
InitialVecType   0
InputEigenVec   0
OutputEigenVec   0
InputHam   0
OutputHam   0
OutputExVec   0
"""

# single-excitation operator sets (HPhi single_ex format): 5 header lines, the
# 2nd being "NSingle <n>", then "<site> <spin> <type> <re> <im>" lines.
# op0 = c_{0,up}; op1 = c_{0,up} + c_{1,up}.  Both annihilation (type 0), spin
# up, so they map to the SAME Hilbert sector but give DIFFERENT spectra.
_OP0 = "=====\nNSingle 1\n=====\n=====\n=====\n0 0 0 1.0 0.0\n"
_OP1 = "=====\nNSingle 2\n=====\n=====\n=====\n0 0 0 1.0 0.0\n1 0 0 1.0 0.0\n"

_NAMELIST = """\
         ModPara  modpara_spec.def
         LocSpin  locspn.def
           Trans  trans.def
    CoulombIntra  coulombintra.def
         CalcMod  calcmod_spec.def
     SpectrumVec  zvo_eigenvec
  SingleExcitation  {single_ex}
"""


def _write(path, content):
    with open(path, "w") as f:
        f.write(content)


def _make_modpara(work, num_op):
    """Spectrum modpara = the eigenvalue modpara + SpectrumLoopExct/SpectrumNumOp."""
    with open(os.path.join(work, "modpara.def")) as f:
        lines = [ln for ln in f.readlines() if not ln.startswith("PreCG")]
    lines += ["SpectrumLoopExct  {}\n".format(NEXCT),
              "SpectrumNumOp  {}\n".format(num_op),
              "PreCG          1\n"]
    _write(os.path.join(work, "modpara_spec.def"), "".join(lines))


def _run_hphi(work, namelist):
    subprocess.run([HPHI_EXEC, "-e", namelist], cwd=work, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def _green(work, name):
    return numpy.loadtxt(os.path.join(work, "output", name))


def _eigenvalue_run(work):
    """Standard-mode eigenvalue step -> eigenvectors + zvo_energy.dat + def files."""
    _write(os.path.join(work, "stan1.in"), _STAN.format(nexct=NEXCT))
    subprocess.run([HPHI_EXEC, "-s", "stan1.in"], cwd=work, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    assert os.path.exists(os.path.join(work, "output", "zvo_energy.dat"))
    assert os.path.exists(os.path.join(work, "output", "zvo_eigenvec_0_rank_0.dat"))
    _write(os.path.join(work, "calcmod_spec.def"), _CALCMOD_SPEC)


def test_multiop_loop_matches_per_operator(tmp_path):
    work = str(tmp_path)
    _eigenvalue_run(work)
    _write(os.path.join(work, "single_ex_0.def"), _OP0)
    _write(os.path.join(work, "single_ex_1.def"), _OP1)

    # --- combined run: SpectrumNumOp=2, one launch -> _<idx>_<op>.dat ---
    _make_modpara(work, 2)
    _write(os.path.join(work, "namelist_spec.def"),
           _NAMELIST.format(single_ex="single_ex_0.def"))
    _run_hphi(work, "namelist_spec.def")
    combined = {(i, op): _green(work, "zvo_DynamicalGreen_{}_{}.dat".format(i, op))
                for i in range(NEXCT) for op in (0, 1)}

    # --- reference runs: SpectrumNumOp=1, one launch per operator -> _<idx>.dat ---
    ref = {}
    for op, single_ex in ((0, "single_ex_0.def"), (1, "single_ex_1.def")):
        _make_modpara(work, 1)
        _write(os.path.join(work, "namelist_spec.def"),
               _NAMELIST.format(single_ex=single_ex))
        _run_hphi(work, "namelist_spec.def")
        for i in range(NEXCT):
            ref[(i, op)] = _green(work, "zvo_DynamicalGreen_{}.dat".format(i))

    # (sanity) the two operators must give genuinely different spectra, otherwise
    # the op->slot mapping would not be exercised.
    op_spread = numpy.abs(ref[(0, 0)] - ref[(0, 1)]).max()
    assert op_spread > 1e-2, f"operators not distinct enough ({op_spread}); test is vacuous"

    # (main) batched op-inner result == per-operator result, for every (idx, op).
    worst = 0.0
    for key in combined:
        worst = max(worst, numpy.abs(combined[key] - ref[key]).max())
    assert worst < 1e-10, f"multi-op loop differs from per-operator runs: {worst}"


def _run_multiop_expect_failure(work, single_ex_1):
    """Run a SpectrumNumOp=2 launch whose set 1 is invalid; expect a clean abort."""
    _write(os.path.join(work, "single_ex_0.def"), _OP0)
    _write(os.path.join(work, "single_ex_1.def"), single_ex_1)
    _make_modpara(work, 2)
    _write(os.path.join(work, "namelist_spec.def"),
           _NAMELIST.format(single_ex="single_ex_0.def"))
    for f in os.listdir(os.path.join(work, "output")):
        if "DynamicalGreen" in f:
            os.remove(os.path.join(work, "output", f))
    proc = subprocess.run([HPHI_EXEC, "-e", "namelist_spec.def"], cwd=work,
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    assert proc.returncode != 0, "HPhi should abort on the invalid operator set"
    # and it must not leave any partial per-state output behind.
    leftover = [f for f in os.listdir(os.path.join(work, "output")) if "DynamicalGreen" in f]
    assert not leftover, f"partial output left behind: {leftover}"


def test_multiop_rejects_sector_mismatch(tmp_path):
    """Set 1 in a different Hilbert sector (creation vs annihilation) must be rejected."""
    work = str(tmp_path)
    _eigenvalue_run(work)
    # op0 = c_{0,up} (annihilation, type 0); set 1 = c_{1,up}^dag (creation, type 1).
    _run_multiop_expect_failure(work, "=====\nNSingle 1\n=====\n=====\n=====\n1 0 1 1.0 0.0\n")


def test_multiop_rejects_truncated_operator_file(tmp_path):
    """A single_ex_<op>.def truncated before its operator line must be rejected."""
    work = str(tmp_path)
    _eigenvalue_run(work)
    # full 5-line header announcing one operator, but the operator line is missing.
    _run_multiop_expect_failure(work, "=====\nNSingle 1\n=====\n=====\n=====\n")
