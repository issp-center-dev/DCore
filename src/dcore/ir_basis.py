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
"""Cached sparse-ir (IR basis) helpers.

Optional dependency: sparse-ir is imported lazily, so the default dense-Matsubara
code paths never require it.
"""

_INSTALL_HINT = (
    "sparse-ir is required for the IR basis option; "
    "install it with 'pip install dcore[ir]' (or 'pip install sparse-ir')."
)

# Cache of built bases keyed by (beta, wmax, eps, statistics). Building the SVE is
# expensive and the same basis is reused across all orbital pairs / blocks.
_basis_cache = {}


def _import_sparse_ir():
    try:
        import sparse_ir
    except ImportError as exc:  # pragma: no cover - exercised only without sparse-ir
        raise ImportError(_INSTALL_HINT) from exc
    return sparse_ir


def get_basis(beta, wmax, eps=1e-10, statistics='F'):
    """Return a cached sparse_ir.FiniteTempBasis.

    Args:
        beta (float): Inverse temperature.
        wmax (float): Real-frequency cutoff.
        eps (float): Basis truncation tolerance.
        statistics (str): 'F' (fermion) or 'B' (boson).

    Returns:
        sparse_ir.FiniteTempBasis
    """
    key = (float(beta), float(wmax), float(eps), statistics)
    if key not in _basis_cache:
        sparse_ir = _import_sparse_ir()
        _basis_cache[key] = sparse_ir.FiniteTempBasis(
            statistics, beta=key[0], wmax=key[1], eps=key[2])
    return _basis_cache[key]
