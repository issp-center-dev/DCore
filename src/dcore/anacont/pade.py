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

import numpy

from dcore._dispatcher import MeshImFreq, GfReFreq, GfImFreq
from dcore.program_options import create_parser, parse_parameters


def _set_n_pade(omega_cutoff, beta, n_min, n_max):
    """
    Return (int)n_pade: the number of Matsubara frequencies below the cutoff frequency.
    n_pade is bounded between n_min and n_max
    """
    n_pade = int((beta * omega_cutoff + numpy.pi) / (2.0 * numpy.pi))
    print("n_pade = {} (evaluated from omega_pade)".format(n_pade))
    n_pade = max(n_pade, n_min)
    n_pade = min(n_pade, n_max)
    print("n_pade = {}".format(n_pade))
    return n_pade


def anacont(sigma_iw_npz, beta, mesh_w, params_pade):
    iw_max = params_pade["iomega_max"]
    n_min = params_pade["n_min"]
    n_max = params_pade["n_max"]
    eta = params_pade["eta"]
    n_pade = _set_n_pade(iw_max, beta, n_min=n_min, n_max=n_max)

    data_w = {}
    num_data = numpy.sum([key.startswith("data") for key in sigma_iw_npz.keys()])
    for idata in range(num_data):
        key = f"data{idata}"
        data = sigma_iw_npz[key]
        mesh_iw = MeshImFreq(beta, "Fermion", data.shape[0] // 2)
        sigma_iw = GfImFreq(data=data, beta=beta, mesh=mesh_iw)
        sigma_w = GfReFreq(mesh=mesh_w, target_shape=data.shape[1:])
        sigma_w.set_from_pade(sigma_iw, n_points=n_pade, freq_offset=eta)
        data_w[key] = sigma_w.data
    return data_w


def parameters_from_ini(inifile):
    parser = create_parser(["post.anacont.pade"])
    parser.read(inifile)
    params = parser.as_dict()
    return params["post.anacont.pade"]
