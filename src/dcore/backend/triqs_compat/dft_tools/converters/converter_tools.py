
##########################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Aichhorn
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################
from ... import mpi

class ConverterTools:

    def __init__(self):
        pass

    def read_fortran_file(self, filename, to_replace):
        """
        Returns a generator that yields all numbers in the Fortran file as float, with possible replacements.

        Parameters
        ----------
        filename : string
                   Name of Fortran-produced file.
        to_replace : dict of str:str
                     Dictionary defining old_char:new_char.

        Yields
        ------
        string
            The next number in file.

        """
        import os.path
        import string
        if not(os.path.exists(filename)):
            raise IOError("File %s does not exist." % filename)
        for line in open(filename, 'r'):
            for old, new in to_replace.items():
                line = line.replace(old, new)
            for x in line.split():
                yield float(x)

    def repack(self):
        """
        Calls the h5repack routine in order to reduce the file size of the hdf5 archive.

        Note
        ----
        Should only be used before the first invokation of HDFArchive in the program, 
        otherwise the hdf5 linking will be broken.

        """

        import subprocess

        if not (mpi.is_master_node()):
            return
        mpi.report("Repacking the file %s" % self.hdf_file)

        retcode = subprocess.call(
            ["h5repack", "-i%s" % self.hdf_file, "-otemphgfrt.h5"])
        if retcode != 0:
            mpi.report("h5repack failed!")
        else:
            subprocess.call(["mv", "-f", "temphgfrt.h5", "%s" % self.hdf_file])

    def det_shell_equivalence(self, corr_shells):
        """
        Determine the equivalence of correlated shells.

        Parameters
        ----------
        corr_shells : list of dicts
                      See documentation of necessary hdf5 elements.

        Returns
        -------
        n_inequiv_shells : integer
                           Number of inequivalent shells.
        corr_to_inequiv : list
                          Mapping between correlated shell index and inequivalent shell index.
                          corr_to_inequiv(i_corr_shells) = i_inequiv_shells 
        inequiv_to_corr : list
                          Mapping between inequivalent shell index and correlated shell index.
                          inequiv_to_corr(i_inequiv_shells) = i_corr_shells

        Note
        ----
        This is needed to set the self energies of all equivalent shells and to extract G_loc.

        """
        corr_to_inequiv = [0 for i in range(len(corr_shells))]
        inequiv_to_corr = [0]
        n_inequiv_shells = 1

        if len(corr_shells) > 1:
            inequiv_sort = [corr_shells[0]['sort']]
            inequiv_l = [corr_shells[0]['l']]
            for i in range(len(corr_shells) - 1):
                is_equiv = False
                for j in range(n_inequiv_shells):
                    if (inequiv_sort[j] == corr_shells[i + 1]['sort']) and (inequiv_l[j] == corr_shells[i + 1]['l']):
                        is_equiv = True
                        corr_to_inequiv[i + 1] = j
                if is_equiv == False:
                    corr_to_inequiv[i + 1] = n_inequiv_shells
                    n_inequiv_shells += 1
                    inequiv_sort.append(corr_shells[i + 1]['sort'])
                    inequiv_l.append(corr_shells[i + 1]['l'])
                    inequiv_to_corr.append(i + 1)

        return n_inequiv_shells, corr_to_inequiv, inequiv_to_corr
