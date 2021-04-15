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


from .base import PytriqsMPISolver

class TRIQSHubbardISolver(PytriqsMPISolver):

    def __init__(self, beta, gf_struct, u_mat, n_iw=1025):
        """

        Initialize the solver.

        """

        super(TRIQSHubbardISolver, self).__init__(beta, gf_struct, u_mat, n_iw)

    def _impl_module_name(self):
        return "dcore.impurity_solvers.triqs_hubbard_I_impl"

    def name(self):
        return "TRIQS/hubbard-I"

    @classmethod
    def is_gf_realomega_available(cls):
        return True
