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
HPhi-level golden test for the multi-bra finite-T spectrum loop
(``SpectrumNumBra`` / Stage-3 bra-column reuse).

When ``SpectrumNumBra = L`` the spectrum calculation solves the resolvent for ONE
ket A|phi> a single time and projects it onto every bra B_b|phi>, b = 0..L-1
(Komega's ``nl`` left vectors), writing ``zvo_DynamicalGreen_<idx>_<op>_<b>.dat``.
Bra set 0 is the namelist ``SingleExcitationBra``; sets 1.. come from
``single_ex_bra_<b>.def``.  This must be numerically identical to running each
bra separately with ``SpectrumNumBra = 1`` (the single-bra off-diagonal path,
writing ``zvo_DynamicalGreen_<idx>.dat``): reusing one ket solve for all bras
only avoids redundant BiCG solves, it must not change the result.

This is the core enabler that lets the DCore HPhi driver cut the BiCG count for
the one-body Green's function from n_orb**2 to n_orb.  The test drives HPhi
directly on a 2-site Hubbard dimer, independent of the DCore driver, and runs
only when an HPhi executable is available (``DCORE_HPHI_EXEC`` or ``HPhi`` on
``PATH``).
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
    reason="HPhi executable not found; set DCORE_HPHI_EXEC to enable the multi-bra test",
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
# All annihilation (type 0), spin up, so ket and both bras map to the SAME
# Hilbert sector but give DIFFERENT (diagonal vs off-diagonal) spectra.
_KET = "=====\nNSingle 1\n=====\n=====\n=====\n0 0 0 1.0 0.0\n"      # c_{0,up}
_BRA0 = "=====\nNSingle 1\n=====\n=====\n=====\n0 0 0 1.0 0.0\n"     # c_{0,up}  -> diagonal G_00
_BRA1 = "=====\nNSingle 1\n=====\n=====\n=====\n1 0 0 1.0 0.0\n"     # c_{1,up}  -> off-diagonal G_10

_NAMELIST = """\
         ModPara  modpara_spec.def
         LocSpin  locspn.def
           Trans  trans.def
    CoulombIntra  coulombintra.def
         CalcMod  calcmod_spec.def
     SpectrumVec  zvo_eigenvec
  SingleExcitation  {ket}
  SingleExcitationBra  {bra}
"""


def _write(path, content):
    with open(path, "w") as f:
        f.write(content)


def _make_modpara(work, num_bra):
    """Spectrum modpara = the eigenvalue modpara + SpectrumLoopExct/SpectrumNumBra."""
    with open(os.path.join(work, "modpara.def")) as f:
        lines = [ln for ln in f.readlines() if not ln.startswith("PreCG")]
    lines += ["SpectrumLoopExct  {}\n".format(NEXCT),
              "SpectrumNumBra  {}\n".format(num_bra),
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


def test_multibra_loop_matches_per_bra(tmp_path):
    work = str(tmp_path)
    _eigenvalue_run(work)
    _write(os.path.join(work, "single_ex_0.def"), _KET)
    _write(os.path.join(work, "single_ex_bra_0.def"), _BRA0)
    _write(os.path.join(work, "single_ex_bra_1.def"), _BRA1)

    # --- combined run: SpectrumNumBra=2, one ket solve -> _<idx>_<op>_<bra>.dat ---
    # bra set 0 = namelist SingleExcitationBra (= single_ex_bra_0.def);
    # bra set 1 = single_ex_bra_1.def.  op field is forced on (op=0).
    _make_modpara(work, 2)
    _write(os.path.join(work, "namelist_spec.def"),
           _NAMELIST.format(ket="single_ex_0.def", bra="single_ex_bra_0.def"))
    _run_hphi(work, "namelist_spec.def")
    combined = {(i, b): _green(work, "zvo_DynamicalGreen_{}_0_{}.dat".format(i, b))
                for i in range(NEXCT) for b in (0, 1)}

    # --- reference runs: SpectrumNumBra=1, one launch per bra -> _<idx>.dat ---
    ref = {}
    for b, bra in ((0, "single_ex_bra_0.def"), (1, "single_ex_bra_1.def")):
        _make_modpara(work, 1)
        _write(os.path.join(work, "namelist_spec.def"),
               _NAMELIST.format(ket="single_ex_0.def", bra=bra))
        _run_hphi(work, "namelist_spec.def")
        for i in range(NEXCT):
            ref[(i, b)] = _green(work, "zvo_DynamicalGreen_{}.dat".format(i))

    # (sanity) the diagonal (bra 0) and off-diagonal (bra 1) spectra must genuinely
    # differ, otherwise the bra->slot mapping would not be exercised.
    bra_spread = numpy.abs(ref[(0, 0)] - ref[(0, 1)]).max()
    assert bra_spread > 1e-2, f"bras not distinct enough ({bra_spread}); test is vacuous"

    # (main) one-ket-solve multi-bra result == per-bra result, for every (idx, bra).
    worst = 0.0
    for key in combined:
        worst = max(worst, numpy.abs(combined[key] - ref[key]).max())
    assert worst < 1e-10, f"multi-bra loop differs from per-bra runs: {worst}"


def test_multibra_rejects_sector_mismatch(tmp_path):
    """A bra set in a different Hilbert sector (creation vs annihilation) must be rejected,
    leaving no partial output behind."""
    work = str(tmp_path)
    _eigenvalue_run(work)
    _write(os.path.join(work, "single_ex_0.def"), _KET)
    _write(os.path.join(work, "single_ex_bra_0.def"), _BRA0)
    # bra set 1 = c_{1,up}^dag (creation, type 1) -> opposite sector shift to the ket.
    _write(os.path.join(work, "single_ex_bra_1.def"),
           "=====\nNSingle 1\n=====\n=====\n=====\n1 0 1 1.0 0.0\n")
    _make_modpara(work, 2)
    _write(os.path.join(work, "namelist_spec.def"),
           _NAMELIST.format(ket="single_ex_0.def", bra="single_ex_bra_0.def"))
    for f in os.listdir(os.path.join(work, "output")):
        if "DynamicalGreen" in f:
            os.remove(os.path.join(work, "output", f))
    proc = subprocess.run([HPHI_EXEC, "-e", "namelist_spec.def"], cwd=work,
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    assert proc.returncode != 0, "HPhi should abort on the mismatched bra set"
    leftover = [f for f in os.listdir(os.path.join(work, "output")) if "DynamicalGreen" in f]
    assert not leftover, f"partial output left behind: {leftover}"
