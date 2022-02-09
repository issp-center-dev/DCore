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


import os
from dcore._dispatcher import *
from dcore.impurity_solvers.alps_cthyb import *
from dcore.tools import to_spin_full_U_matrix

beta = 10.0
n_iw = 1000
n_l = 20

def test_copy_between_numpy_blockgf():
    block_names = ['up', 'down']
    g1 = GfLegendre(indices=[0, 1], beta=beta, n_points=n_l, name="up")
    g1.data[:,:,:] = numpy.random.rand(g1.data.shape[0], g1.data.shape[1], g1.data.shape[2])
    g2 = GfLegendre(indices=[0, 1], beta=beta, n_points=n_l, name="down")
    g2.data[:,:,:] = numpy.random.rand(g2.data.shape[0], g2.data.shape[1], g2.data.shape[2])
    G = BlockGf(name_list=('up', 'down'), block_list=(g1, g2), make_copies=False)

    data = to_numpy_array(G, block_names)

    G_reconst = G.copy()
    # data is dimension of (l, spin_orb, spin_orb)
    assign_from_numpy_array_legendre(G_reconst, data.transpose((1,2,0)), block_names)

    for name, g in G:
        numpy.allclose(G[name].data, G_reconst[name].data, 1e-10)

def test_solver_dry_run(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    # Half band width
    D = 2.0
    mu = 1.0
    U = 4.0
    J = 1.0

    l = 1
    n_orbs = 2 * l +1

    # Inter-orbital hybridization
    off_diag = True

    spin_names = ['ud']
    orb_names = [i for i in range(2*n_orbs)]

    # Block structure of Green's functions
    gf_struct = set_operator_structure(spin_names, orb_names, off_diag=off_diag)
    # Convert to dict
    if isinstance(gf_struct, list):
        gf_struct = {x[0]: x[1] for x in gf_struct}

    # Local interaction Hamiltonian
    U_mat = to_spin_full_U_matrix(U_matrix(l=l, U_int=U, J_hund=J, basis='spherical'))

    s = ALPSCTHYBSolver(beta, gf_struct, U_mat, n_iw=n_iw)

    Delta_iw = make_block_gf(GfImFreq, gf_struct, beta, n_iw)

    for name, b in Delta_iw:
        b << (D/2.0)**2 * SemiCircular(D)

    # Random local transfer matrix
    H0 = [numpy.random.rand(2*n_orbs, 2*n_orbs) for name, b in Delta_iw]
    H0 = [h + h.conjugate().transpose() for h in H0]

    G0_iw = make_block_gf(GfImFreq, gf_struct, beta, n_iw)
    G0_iw.zero()
    for ib, (name, g0) in enumerate(G0_iw):
        g0 << inverse(iOmega_n - H0[ib] - Delta_iw[name])
    s.set_G0_iw(G0_iw)

    #rot = compute_diag_basis(G0_iw)
    #s.rotate_basis(rot, direction='forward')

    p = {
        'exec_path'     : './dummy/hybmat',
        'work_dir'      : './',
        'dry_run'       : True,
        'random_seed_offset' : 20,
        'timelimit'     : 60,
    }
    s.solve(None, '', p)

    #s.rotate_basis(rot, direction='backward')

    # Check Delta(iwn)
    diff = s.get_Delta_iw() - Delta_iw

    for name, b in diff:
        # FIXME: The precision is worse with sparse sampling. Why?
        assert numpy.all(numpy.abs(b.data) < 1e-4)

    os.chdir(org_dir)
