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

from dcore.fourier import _fft_fermion_w2t, _fft_fermion_t2w, _matsubara_freq_fermion
import numpy
import os


def test_ftt_fermion(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    beta = 10.0
    # nw = 2048
    nw = 1024
    # nw = 512
    nt = 2 * nw
    e0 = 2.0

    print("nw =", nw)
    print("nt =", nt)

    iw = _matsubara_freq_fermion(beta, nw)
    g0_w = 1 / (iw - e0)

    t = numpy.linspace(0, beta, nt+1)
    g0_t = -numpy.exp(-t*e0) / (numpy.exp(-beta*e0) + 1)
    a = - g0_t[0] - g0_t[-1]

    # test for t2w

    g0_w_fft = _fft_fermion_t2w(g0_t, beta)

    # numpy.savetxt('g0_w.dat', g0_w.view(float).reshape(-1,2))
    # numpy.savetxt('g0_w_fft.dat', g0_w_fft.view(float).reshape(-1,2))
    # numpy.savetxt('g0_w_diff.dat', (g0_w_fft - g0_w).view(float).reshape(-1,2))

    assert numpy.allclose(g0_w_fft, g0_w, atol=1e-3)

    # test for w2t

    g0_t_fft = _fft_fermion_w2t(g0_w, beta, a=a)

    # numpy.savetxt('g0_t.dat', g0_t)
    # numpy.savetxt('g0_t_fft.dat', g0_t_fft)
    # numpy.savetxt('g0_t_diff.dat', g0_t_fft - g0_t)

    assert numpy.allclose(g0_t_fft, g0_t, atol=1e-3)
