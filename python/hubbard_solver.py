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
import pytriqs.applications.impurity_solvers.hubbard_I

class Solver(pytriqs.applications.impurity_solvers.hubbard_I.Solver, object):
    """
       Hubbard I Solver
    """

    # initialisation:
    def __init__(self, beta, l, n_msb=1025, use_spin_orbit=False, Nmoments=5):
        super(Solver, self).__init__(beta, l, n_msb, use_spin_orbit, Nmoments)

    # Make read-only getter
    @property
    def Sigma_iw(self):
        return self.Sigma

    # Make read-only getter
    @property
    def G_iw(self):
        return self.G

    # Make read-only getter
    @property
    def G0_iw(self):
        return self.G0
