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
Unit tests for SumkDFT_opt.calc_mu_newton (safeguarded Newton chemical-potential
search). The root-finder is driven by a hand-built n(mu)/n'(mu) so the test needs
no real lattice model. In particular it checks the charge-gap case, where
dn/dmu ~ 0 and a naive Newton step would blow up.
"""
import math

import pytest

SumkDFT_opt = pytest.importorskip("dcore.sumkdft_opt").SumkDFT_opt


class _Fake:
    pass


def _make_self(n_func, dn_func, density_required, mu0):
    s = _Fake()
    s.density_required = density_required
    s.charge_below = 0.0
    s.chemical_potential = mu0
    s.calls = [0]

    def tdad(mu, **kwargs):
        s.calls[0] += 1
        return n_func(mu), dn_func(mu)

    s.total_density_and_derivative = tdad
    return s


def test_newton_smooth_monotonic():
    # n(mu) = 2 + tanh(10 (mu - 0.7)),  target 2.5 -> mu* = 0.7 + atanh(0.5)/10
    n = lambda mu: 2.0 + math.tanh(10.0 * (mu - 0.7))
    dn = lambda mu: 10.0 * (1.0 - math.tanh(10.0 * (mu - 0.7)) ** 2)
    mu_star = 0.7 + math.atanh(0.5) / 10.0

    # warm start near the root (as in DMFT, mu carries over between iterations)
    s = _make_self(n, dn, density_required=2.5, mu0=0.6)
    mu = SumkDFT_opt.calc_mu_newton(s, precision=1e-8, max_loops=200)

    assert abs(n(mu) - 2.5) < 1e-8
    assert abs(mu - mu_star) < 1e-5
    assert s.calls[0] < 12          # quadratic convergence from a warm start


def test_newton_charge_gap():
    # Two levels at -1 and +1 (sharp): n is ~flat (=1) inside the gap so
    # dn/dmu ~ 0 there. A safeguarded search must still converge (via bisection).
    def n(mu):
        return 0.5 * (1.0 + math.tanh(20.0 * (mu + 1.0))) \
            + 0.5 * (1.0 + math.tanh(20.0 * (mu - 1.0)))

    def dn(mu):
        return 10.0 * (1.0 - math.tanh(20.0 * (mu + 1.0)) ** 2) \
            + 10.0 * (1.0 - math.tanh(20.0 * (mu - 1.0)) ** 2)

    # start far below the gap; target sits on the in-gap plateau (n = 1)
    s = _make_self(n, dn, density_required=1.0, mu0=-3.0)
    mu = SumkDFT_opt.calc_mu_newton(s, precision=1e-6, max_loops=200)

    assert abs(n(mu) - 1.0) < 1e-6          # converged despite dn/dmu ~ 0
    assert -1.0 <= mu <= 1.0                # landed inside the gap
    assert math.isfinite(mu)
    assert s.calls[0] < 80                  # no runaway


def test_newton_already_converged():
    # if the starting mu already matches the target, no extra work is done
    n = lambda mu: 1.0 + mu
    dn = lambda mu: 1.0
    s = _make_self(n, dn, density_required=1.0, mu0=0.0)
    mu = SumkDFT_opt.calc_mu_newton(s, precision=1e-8)
    assert abs(mu) < 1e-12
    assert s.calls[0] == 1


def test_newton_wrong_signed_initial_derivative():
    # n(mu) is monotonic with a valid root, but the derivative reported at the
    # starting point is given the WRONG sign (as can happen from numerical noise
    # near a plateau). The two-directional bracketing must still find the root.
    n = lambda mu: 1.0 + math.tanh(mu - 1.0)
    def dn(mu):
        d = 1.0 - math.tanh(mu - 1.0) ** 2
        return -d if abs(mu) < 1e-9 else d   # wrong sign only at the start mu=0
    mu_star = 1.0 + math.atanh(0.4)          # n = 1.4
    s = _make_self(n, dn, density_required=1.4, mu0=0.0)
    mu = SumkDFT_opt.calc_mu_newton(s, precision=1e-8, max_loops=200)
    assert mu is not None
    assert abs(n(mu) - 1.4) < 1e-8
    assert abs(mu - mu_star) < 1e-5


def test_newton_failure_returns_none():
    # a constant density never reaches a different target -> cannot bracket
    n = lambda mu: 1.0
    dn = lambda mu: 0.0
    s = _make_self(n, dn, density_required=2.0, mu0=0.0)
    mu = SumkDFT_opt.calc_mu_newton(s, precision=1e-8, max_loops=20)
    assert mu is None
    assert s.chemical_potential is None


# --- real total_density_and_derivative against the explicit formula / FD -------

import numpy as _np


class _Mesh:
    def __init__(self, beta):
        self.beta = beta


class _Gf2:
    def __init__(self, data, beta):
        self.data = data
        self.mesh = _Mesh(beta)


class _BlockGf2:
    def __init__(self, blocks):
        self._b = blocks

    def __iter__(self):
        return iter(self._b.items())


class _FakeSumk:
    """Minimal SumkDFT_opt stand-in: lattice_gf(mu) = [(iw+mu) I - H]^{-1}."""
    def __init__(self, H, iw, beta):
        self.H = H
        self.iw = iw
        self.beta = beta
        self.n_orb = H.shape[0]
        self.n_k = 1
        self.bz_weights = _np.array([1.0])

    def lattice_gf(self, ik, mu, **kwargs):
        M = (self.iw[:, None, None] + mu) * _np.eye(self.n_orb)[None] - self.H[None]
        return _BlockGf2({"up": _Gf2(_np.linalg.inv(M), self.beta)})


def test_total_density_and_derivative_rejects_real_freq():
    s = _FakeSumk(_np.eye(2, dtype=complex), _np.array([1j]), 40.0)
    with pytest.raises(ValueError):
        SumkDFT_opt.total_density_and_derivative(s, 0.0, iw_or_w='w')


def test_mu_search_parameter_validation():
    from dcore.program_options import parse_parameters

    def _p(mu_search, no_tail_fit=False):
        return {'system': {'T': -1.0, 'mu_search': mu_search, 'no_tail_fit': no_tail_fit}}

    # parse_parameters only validates the mu_search value; the newton/no_tail_fit
    # consistency is enforced where the search runs (SumkDFTWorkerGloc), so these
    # all parse without raising regardless of no_tail_fit.
    parse_parameters(_p('brent'))
    parse_parameters(_p('newton', no_tail_fit=True))
    parse_parameters(_p('newton', no_tail_fit=False))
    # an unknown value is rejected
    with pytest.raises(SystemExit):
        parse_parameters(_p('foo'))

    # a legacy/pruned [system] without mu_search defaults to 'brent'
    legacy = {'system': {'T': -1.0, 'no_tail_fit': False}}
    parse_parameters(legacy)
    assert legacy['system']['mu_search'] == 'brent'


def test_total_density_and_derivative_formula_and_fd():
    rng = _np.random.RandomState(0)
    n_orb, n_iw, beta = 4, 400, 40.0
    iw = 1j * (2 * _np.arange(-n_iw, n_iw) + 1) * _np.pi / beta
    H = rng.standard_normal((n_orb, n_orb)) + 1j * rng.standard_normal((n_orb, n_orb))
    H = 0.5 * (H + H.conj().T)
    s = _FakeSumk(H, iw, beta)

    mu = 0.3
    n_val, dn_val = SumkDFT_opt.total_density_and_derivative(s, mu)

    # n via the explicit Matsubara formula (same as total_density_matsubara)
    G = s.lattice_gf(0, mu)._b["up"].data
    n_ref = (_np.sum(_np.trace(G, axis1=1, axis2=2)) / beta + 0.5 * n_orb).real
    assert abs(n_val - n_ref) < 1e-10

    # derivative against a central finite difference
    h = 1e-5
    n_plus = SumkDFT_opt.total_density_and_derivative(s, mu + h)[0]
    n_minus = SumkDFT_opt.total_density_and_derivative(s, mu - h)[0]
    assert abs(dn_val - (n_plus - n_minus) / (2 * h)) < 1e-5
