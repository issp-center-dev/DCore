
##########################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2018 by G. J. Kraberger
# Copyright (C) 2018 by Simons Foundation
# Authors: G. J. Kraberger, O. Parcollet
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


import copy
import numpy as np
from ..gf import GfImFreq, BlockGf
from ast import literal_eval
from .. import mpi
from warnings import warn
from collections import defaultdict

class BlockStructure(object):
    """ Contains information about the Green function structure.

    This class contains information about the structure of the solver
    and sumk Green functions and the mapping between them.

    Parameters
    ----------
    gf_struct_sumk : list of list of tuple
        gf_struct_sumk[ish][idx] = (block_name,list of indices in block)

        for correlated shell ish; idx is just a counter in the list
    gf_struct_solver : list of dict
        gf_struct_solver[ish][block] = list of indices in that block

        for *inequivalent* correlated shell ish
    solver_to_sumk : list of dict
        solver_to_sumk[ish][(from_block,from_idx)] = (to_block,to_idx)

        maps from the solver block and index to the sumk block and index
        for *inequivalent* correlated shell ish
    sumk_to_solver : list of dict
        sumk_to_solver[ish][(from_block,from_idx)] = (to_block,to_idx)

        maps from the sumk block and index to the solver block and index
        for *inequivalent* correlated shell ish
    solver_to_sumk_block : list of dict
        solver_to_sumk_block[ish][from_block] = to_block

        maps from the solver block to the sumk block
        for *inequivalent* correlated shell ish
    deg_shells : list of lists of lists OR list of lists of dicts
        In the simple format, ``deg_shells[ish][grp]`` is a list of
        block names; ``ish`` is the index of the inequivalent correlated shell,
        ``grp`` is the index of the group of symmetry-related blocks.
        When the name of two blocks is in the same group, that means that
        these two blocks are the same by symmetry.

        In the more general format, ``deg_shells[ish][grp]`` is a
        dictionary whose keys are the block names that are part of the
        group. The values of the dict for each key are tuples ``(v, conj)``,
        where ``v`` is a transformation matrix with the same matrix dimensions
        as the block and ``conj`` is a bool (whether or not to conjugate).
        Two blocks with ``v_1, conj_1`` and ``v_2, conj_2`` being in the same
        symmetry group means that

        .. math::

            C_1(v_1^\dagger G_1 v_1) = C_2(v_2^\dagger G_2 v_2),

        where the :math:`G_i` are the Green's functions of the block,
        and the functions :math:`C_i` conjugate their argument if the bool
        ``conj_i`` is ``True``.
    corr_to_inequiv : list
        a list where, for each correlated shell, the index of the corresponding
        inequivalent correlated shell is given
    transformation : list of numpy.array or list of dict
        a list with entries for each ``ish`` giving transformation matrices
        that are used on the Green's function in ``sumk`` space before
        converting to the ``solver`` space
        Up to the change in block structure,

        .. math::

            G_{solver} = T G_{sumk} T^\dagger

        if :math:`T` is the ``transformation`` of that particular shell.

        Note that for each shell this can either be a numpy array which
        applies to all blocks or a dict with a transformation matrix for
        each block.
    """

    def __init__(self, gf_struct_sumk=None,
                 gf_struct_solver=None,
                 solver_to_sumk=None,
                 sumk_to_solver=None,
                 solver_to_sumk_block=None,
                 deg_shells=None,
                 corr_to_inequiv = None,
                 transformation=None):
        self.gf_struct_sumk = gf_struct_sumk
        self.gf_struct_solver = gf_struct_solver
        self.solver_to_sumk = solver_to_sumk
        self.sumk_to_solver = sumk_to_solver
        self.solver_to_sumk_block = solver_to_sumk_block
        self.deg_shells = deg_shells
        self.corr_to_inequiv = corr_to_inequiv
        self.transformation = transformation

    @property
    def gf_struct_solver_list(self):
        """ The structure of the solver Green's function

        This is returned as a
        list (for each shell)
        of lists (for each block)
        of tuples (block_name, block_indices).

        That is,
        ``gf_struct_solver_list[ish][b][0]``
        is the name of the block number ``b`` of shell ``ish``, and
        ``gf_struct_solver_list[ish][b][1]``
        is a list of its indices.

        The list for each shell is sorted alphabetically by block name.
        """
        if self.gf_struct_solver is None:
            return None
        # we sort by block name in order to get a reproducible result
        return [sorted([(k, v) for k, v in list(gfs.items())], key=lambda x: x[0])
                for gfs in self.gf_struct_solver]

    @property
    def gf_struct_sumk_list(self):
        """ The structure of the sumk Green's function

        This is returned as a
        list (for each shell)
        of lists (for each block)
        of tuples (block_name, block_indices)

        That is,
        ``gf_struct_sumk_list[ish][b][0]``
        is the name of the block number ``b`` of shell ``ish``, and
        ``gf_struct_sumk_list[ish][b][1]``
        is a list of its indices.
        """
        return self.gf_struct_sumk

    @property
    def gf_struct_solver_dict(self):
        """ The structure of the solver Green's function

        This is returned as a
        list (for each shell)
        of dictionaries.

        That is,
        ``gf_struct_solver_dict[ish][bname]``
        is a list of the indices of block ``bname`` of shell ``ish``.
        """
        return self.gf_struct_solver

    @property
    def gf_struct_sumk_dict(self):
        """ The structure of the sumk Green's function

        This is returned as a
        list (for each shell)
        of dictionaries.

        That is,
        ``gf_struct_sumk_dict[ish][bname]``
        is a list of the indices of block ``bname`` of shell ``ish``.
        """
        if self.gf_struct_sumk is None:
            return None
        return [{block: indices for block, indices in gfs}
                for gfs in self.gf_struct_sumk]

    @property
    def inequiv_to_corr(self):
        """ A list mapping an inequivalent correlated shell to a correlated shell
        """

        if self.corr_to_inequiv is None:
            return None
        N_solver = len(np.unique(self.corr_to_inequiv))
        if self.gf_struct_solver is not None:
            assert N_solver == len(self.gf_struct_solver)
        assert sorted(np.unique(self.corr_to_inequiv)) == list(range(N_solver)),\
            "an inequivalent shell is missing in corr_to_inequiv"
        return [self.corr_to_inequiv.index(icrsh)
                for icrsh in list(range(N_solver))]

    @inequiv_to_corr.setter
    def inequiv_to_corr(self, value):
        # a check for backward compatibility
        if value is None:
            return
        assert self.inequiv_to_corr == value, "Trying to set incompatible inequiv_to_corr"

    @property
    def sumk_to_solver_block(self):
        if self.inequiv_to_corr is None:
            return None
        ret = []
        for ish, icrsh in enumerate(self.inequiv_to_corr):
            d = defaultdict(list)
            for block_solver, block_sumk in list(self.solver_to_sumk_block[ish].items()):
                d[block_sumk].append(block_solver)
            ret.append(d)
        return ret

    @property
    def effective_transformation_sumk(self):
        """ Return the effective transformation matrix

        A list of dicts, one for every correlated shell. In the dict,
        there is a transformation matrix (as numpy array) for each
        block in sumk space, that is used to transform the block.
        """
        trans = copy.deepcopy(self.transformation)
        if self.gf_struct_sumk is None:
            raise Exception('gf_struct_sumk not set.')
        if self.gf_struct_solver is None:
            raise Exception('gf_struct_solver not set.')

        if trans is None:
            trans = [{block: np.eye(len(indices)) for block, indices in gfs}
                     for gfs in self.gf_struct_sumk]

        assert isinstance(trans, list),\
            "transformation has to be a list"

        assert len(trans) == len(self.gf_struct_sumk),\
            "give one transformation per correlated shell"

        for icrsh in list(range(len(trans))):
            ish = self.corr_to_inequiv[icrsh]
            if trans[icrsh] is None:
                trans[icrsh] = {block: np.eye(len(indices))
                                for block, indices in self.gf_struct_sumk[icrsh]}

            if not isinstance(trans[icrsh], dict):
                trans[icrsh] = {block: copy.deepcopy(trans[icrsh])
                                for block, indices in self.gf_struct_sumk[icrsh]}

            assert list(trans[icrsh].keys()) == list(self.gf_struct_sumk_dict[icrsh].keys()),\
                "wrong block names used in transformation (icrsh = {})".format(icrsh)

            for block in trans[icrsh]:
                assert trans[icrsh][block].shape[0] == trans[icrsh][block].shape[1],\
                    "Transformation has to be quadratic; throwing away orbitals can be achieved on the level of the mapping. (icrsh = {}, block = {})".format(icrsh, block)

                assert trans[icrsh][block].shape[0] == len(self.gf_struct_sumk_dict[icrsh][block]),\
                    "Transformation block shape does not match gf_struct_sumk. (icrsh = {}, block = {})".format(icrsh, block)

                # zero out all the lines of the transformation that are
                # not included in gf_struct_solver
                for iorb, norb in enumerate(self.gf_struct_sumk_dict[icrsh][block]):
                    if self.sumk_to_solver[ish][(block, norb)][0] is None:
                        trans[icrsh][block][iorb, :] = 0.0
        return trans

    @property
    def effective_transformation_solver(self):
        """ Return the effective transformation matrix

        A list of dicts, one for every inequivalent correlated shell.
        In the dict, there is a transformation matrix (as numpy array)
        for each block in solver space, that is used to transform from
        the sumk block (see :py:meth:`.solver_to_sumk_block`) to the
        solver block.


        For a solver block ``b`` for inequivalent correlated shell ``ish``,
        the corresponding block of the solver Green's function is::

            # the effective transformation matrix for the block
            T = block_structure.effective_transformation_solver[ish][b]
            # the index of the correlated shell
            icrsh = block_structure.inequiv_to_corr[ish]
            # the name of the corresponding sumk block
            block_sumk = block_structure.solver_to_sumk_block[icrsh][b]
            # transform the Green's function
            G_solver[ish][b].from_L_G_R(T, G_sumk[icrsh][block_sumk], T.conjugate().transpose())

        The functionality of that code block is implemented in
        :py:meth:`.convert_gf` (i.e., you don't need to use this directly).
        """

        eff_trans_sumk = self.effective_transformation_sumk

        ets = []
        for ish in range(len(self.gf_struct_solver)):
            icrsh = self.inequiv_to_corr[ish]
            ets.append(dict())
            for block in self.gf_struct_solver[ish]:
                block_sumk = self.solver_to_sumk_block[ish][block]
                T = eff_trans_sumk[icrsh][block_sumk]
                ets[ish][block] = np.zeros((len(self.gf_struct_solver[ish][block]),
                                            len(T)),
                                           dtype=T.dtype)
                for i in self.gf_struct_solver[ish][block]:
                    i_sumk = self.solver_to_sumk[ish][block, i]
                    assert i_sumk[0] == block_sumk,\
                        "Wrong block in solver_to_sumk"
                    i_sumk = i_sumk[1]
                    ets[ish][block][i, :] = T[i_sumk, :]
        return ets


    @classmethod
    def full_structure(cls,gf_struct,corr_to_inequiv):
        """ Construct structure that maps to itself.

        This has the same structure for sumk and solver, and the
        mapping solver_to_sumk and sumk_to_solver is one-to-one.

        Parameters
        ----------
        gf_struct : list of dict
            gf_struct[ish][block] = list of indices in that block

            for (inequivalent) correlated shell ish
        corr_to_inequiv : list
            gives the mapping from correlated shell csh to inequivalent
            correlated shell icsh, so that corr_to_inequiv[csh]=icsh
            e.g. SumkDFT.corr_to_inequiv

            if None, each inequivalent correlated shell is supposed to
            be correspond to just one correlated shell with the same
            index; there is not default, None has to be set explicitly!
        """

        solver_to_sumk = []
        s2sblock = []
        gs_sumk = []
        for ish in range(len(gf_struct)):
            so2su = {}
            so2sublock = {}
            gss = []
            for block in gf_struct[ish]:
                so2sublock[block]=block
                for ind in gf_struct[ish][block]:
                    so2su[(block,ind)]=(block,ind)
                gss.append((block,gf_struct[ish][block]))
            solver_to_sumk.append(so2su)
            s2sblock.append(so2sublock)
            gs_sumk.append(gss)

        # gf_struct_sumk is not given for each inequivalent correlated
        # shell, but for every correlated shell!
        if corr_to_inequiv is not None:
            gs_sumk_all = [None]*len(corr_to_inequiv)
            for i in range(len(corr_to_inequiv)):
                gs_sumk_all[i] = gs_sumk[corr_to_inequiv[i]]
        else:
            gs_sumk_all = gs_sumk

        return cls(gf_struct_solver=copy.deepcopy(gf_struct),
                gf_struct_sumk = gs_sumk_all,
                solver_to_sumk = copy.deepcopy(solver_to_sumk),
                sumk_to_solver = solver_to_sumk,
                solver_to_sumk_block = s2sblock,
                deg_shells = [[] for ish in range(len(gf_struct))],
                corr_to_inequiv = corr_to_inequiv)

    def pick_gf_struct_solver(self, new_gf_struct):
        """ Pick selected orbitals within blocks.

        This throws away parts of the Green's function that (for some
        reason - be sure that you know what you're doing) shouldn't be
        included in the calculation.

        To drop an entire block, just don't include it.
        To drop a certain index within a block, just don't include it.

        If it was before:

        [{'up':[0,1],'down':[0,1],'left':[0,1]}]

        to choose the 0th index of the up block and the 1st index of
        the down block and drop the left block, the new_gf_struct would
        have to be

        [{'up':[0],'down':[1]}]

        Note that the indices will be renamed to be a 0-based
        sequence of integers, i.e. the new structure will actually
        be  [{'up':[0],'down':[0]}].

        For dropped indices, sumk_to_solver will map to (None,None).

        Parameters
        ----------
        new_gf_struct : list of dict
            formatted the same as gf_struct_solver:

            new_gf_struct[ish][block]=list of indices in that block.
        """

        for ish in range(len(self.gf_struct_solver)):
            gf_struct = new_gf_struct[ish]

            # create new solver_to_sumk
            so2su = {}
            so2su_block = {}
            for blk,idxs in list(gf_struct.items()):
                for i in range(len(idxs)):
                    so2su[(blk, i)] = self.solver_to_sumk[ish][(blk, idxs[i])]
                    so2su_block[blk] = so2su[(blk, i)][0]
            self.solver_to_sumk[ish] = so2su
            self.solver_to_sumk_block[ish] = so2su_block
            # create new sumk_to_solver
            for k,v in list(self.sumk_to_solver[ish].items()):
                blk,ind=v
                if blk in gf_struct and ind in gf_struct[blk]:
                    new_ind = gf_struct[blk].index(ind)
                    self.sumk_to_solver[ish][k] = (blk, new_ind)
                else:
                    self.sumk_to_solver[ish][k] = (None, None)
            # adapt deg_shells

            self.adapt_deg_shells(gf_struct, ish)


            # reindexing gf_struct so that it starts with 0
            for k in gf_struct:
                gf_struct[k]=list(range(len(gf_struct[k])))
            self.gf_struct_solver[ish]=gf_struct

    def adapt_deg_shells(self, gf_struct, ish=0):
        """ Adapts the deg_shells to a new gf_struct
            Internally called when using pick_gf_struct and map_gf_struct
        """
        if self.deg_shells is not None:
            for degsh in self.deg_shells[ish]:
                if isinstance(degsh, dict):
                    for key in list(degsh.keys()):
                        if not key in gf_struct:
                            del degsh[key]
                            continue
                        if gf_struct[key] != self.gf_struct_solver[ish][key]:
                            v, C = degsh[key]
                            degsh[key][0] = \
                                v[gf_struct[key], :][:, gf_struct[key]]
                            warn(
                                'Removing shells from degenerate shell {}. Cannot guarantee that they continue to be equivalent.')
                else: # degshell is list
                    degsh1 = copy.deepcopy(degsh) # in order to not remove a key while iterating
                    for key in degsh1:
                        if not key in gf_struct:
                            warn('Removing shells from degenerate shell {}.')
                            degsh.remove(key)

    def pick_gf_struct_sumk(self, new_gf_struct):
        """ Pick selected orbitals within blocks.

        This throws away parts of the Green's function that (for some
        reason - be sure that you know what you're doing) shouldn't be
        included in the calculation.

        To drop an entire block, just don't include it.
        To drop a certain index within a block, just don't include it.

        If it was before:

        [{'up':[0,1],'down':[0,1],'left':[0,1]}]

        to choose the 0th index of the up block and the 1st index of
        the down block and drop the left block, the new_gf_struct would
        have to be

        [{'up':[0],'down':[1]}]

        Note that the indices will be renamed to be a 0-based
        sequence of integers.

        For dropped indices, sumk_to_solver will map to (None,None).

        Parameters
        ----------
        new_gf_struct : list of dict
            formatted the same as gf_struct_solver:

            new_gf_struct[ish][block]=list of indices in that block.

            However, the indices are not according to the solver Gf
            but the sumk Gf.
        """
        eff_trans_sumk = self.effective_transformation_sumk

        assert len(eff_trans_sumk) == len(new_gf_struct),\
            "new_gf_struct has the wrong length"

        new_gf_struct_transformed = copy.deepcopy(new_gf_struct)

        # when there is a transformation matrix, this first zeroes out
        # the corresponding rows of (a copy of) T and then applies
        # pick_gf_struct_solver for all lines of  T that have at least
        # one non-zero entry

        for icrsh in range(len(new_gf_struct)):
            for block, indices in self.gf_struct_sumk[icrsh]:
                if not block in new_gf_struct[icrsh]:
                    #del new_gf_struct_transformed[block] # this MUST be wrong, as new_gf_struct_transformed needs to have a integer index for icrsh... # error when index is not kept at all
                    continue
                T = eff_trans_sumk[icrsh][block]
                for index in indices:
                    if not index in new_gf_struct[icrsh][block]:
                        T[:, index] = 0.0
                new_indices = []
                for index in indices:
                    if np.any(np.abs(T[index, :]) > 1.e-15):
                        new_indices.append(index)
                new_gf_struct_transformed[icrsh][block] = new_indices

        gfs = []
        # construct gfs, which is the equivalent of new_gf_struct
        # but according to the solver Gf, by using the sumk_to_solver
        # mapping
        for icrsh in range(len(new_gf_struct_transformed)):
            ish = self.corr_to_inequiv[icrsh]
            gfs.append({})
            for block in list(new_gf_struct_transformed[icrsh].keys()):
                for ind in new_gf_struct_transformed[icrsh][block]:
                    ind_sol = self.sumk_to_solver[ish][(block,ind)]
                    if not ind_sol[0] in gfs[icrsh]:
                        gfs[icrsh][ind_sol[0]]=[]
                    gfs[icrsh][ind_sol[0]].append(ind_sol[1])
        self.pick_gf_struct_solver(gfs)

    def map_gf_struct_solver(self, mapping):
        r""" Map the Green function structure from one struct to another.

        Parameters
        ----------
        mapping : list of dict
            the dict consists of elements
            (from_block,from_index) : (to_block,to_index)
            that maps from one structure to the other
            (one for each shell; use a mapping ``None`` for a shell
            you want to leave unchanged)

        Examples
        --------

        Consider a `gf_struct_solver` consisting of two :math:`1 \times 1`
        blocks, `block_1` and `block_2`. Say you want to have a new block
        structure where you merge them into one block because you want to
        introduce an off-diagonal element. You could perform the mapping
        via::

            map_gf_struct_solver([{('block_1',0) : ('block', 0)
                                   ('block_2',0) : ('block', 1)}])
        """

        for ish in range(len(mapping)):
            if mapping[ish] is None:
                continue
            gf_struct = {}
            so2su = {}
            su2so = {}
            so2su_block = {}
            for frm,to in list(mapping[ish].items()):
                if not to[0] in gf_struct:
                    gf_struct[to[0]] = []
                gf_struct[to[0]].append(to[1])

                so2su[to] = self.solver_to_sumk[ish][frm]
                su2so[self.solver_to_sumk[ish][frm]] = to
                if to[0] in so2su_block:
                    if so2su_block[to[0]] != \
                            self.solver_to_sumk_block[ish][frm[0]]:
                            warn("solver block '{}' maps to more than one sumk block: '{}', '{}'".format(
                                to[0], so2su_block[to[0]], self.solver_to_sumk_block[ish][frm[0]]))
                else:
                    so2su_block[to[0]] =\
                        self.solver_to_sumk_block[ish][frm[0]]
            for k in list(self.sumk_to_solver[ish].keys()):
                if not k in su2so:
                    su2so[k] = (None, None)

            self.adapt_deg_shells(gf_struct, ish)

            self.gf_struct_solver[ish] = gf_struct
            self.solver_to_sumk[ish] = so2su
            self.sumk_to_solver[ish] = su2so
            self.solver_to_sumk_block[ish] = so2su_block
            self.deg_shells[ish] = []

    def create_gf(self, ish=0, gf_function=GfImFreq, space='solver', **kwargs):
        """ Create a zero BlockGf having the correct structure.

        For ``space='solver'``, the structure is according to
        ``gf_struct_solver``, else according to ``gf_struct_sumk``.

        When using GfImFreq as gf_function, typically you have to
        supply beta as keyword argument.

        Parameters
        ----------
        ish : int
            shell index
            If ``space='solver'``, the index of the of the inequivalent correlated shell,
            if ``space='sumk'``, the index of the correlated shell
        gf_function : constructor
            function used to construct the Gf objects constituting the
            individual blocks; default: GfImFreq
        space : 'solver' or 'sumk'
            which space the structure should correspond to
        **kwargs :
            options passed on to the Gf constructor for the individual
            blocks
        """

        return self._create_gf_or_matrix(ish, gf_function, BlockGf, space, **kwargs)

    def create_matrix(self, ish=0, space='solver', dtype=np.complex_):
        """ Create a zero matrix having the correct structure.

        For ``space='solver'``, the structure is according to
        ``gf_struct_solver``, else according to ``gf_struct_sumk``.

        Parameters
        ----------
        ish : int
            shell index
            If ``space='solver'``, the index of the of the inequivalent correlated shell,
            if ``space='sumk'``, the index of the correlated shell
        space : 'solver' or 'sumk'
            which space the structure should correspond to
        """

        def gf_function(indices):
            return np.zeros((len(indices), len(indices)), dtype=dtype)

        def block_function(name_list, block_list):
            d = dict()
            for i in range(len(name_list)):
                d[name_list[i]] = block_list[i]
            return d

        return self._create_gf_or_matrix(ish, gf_function, block_function, space)

    def _create_gf_or_matrix(self, ish=0, gf_function=GfImFreq, block_function=BlockGf, space='solver', **kwargs):
        if space == 'solver':
            gf_struct = self.gf_struct_solver
        elif space == 'sumk':
            gf_struct = self.gf_struct_sumk_dict
        else:
            raise Exception(
                "Argument space has to be either 'solver' or 'sumk'.")

        names = list(gf_struct[ish].keys())
        blocks = []
        for n in names:
            G = gf_function(indices=gf_struct[ish][n], **kwargs)
            blocks.append(G)
        G = block_function(name_list=names, block_list=blocks)
        return G

    def check_gf(self, G, ish=None, space='solver'):
        """ check whether the Green's function G has the right structure

        This throws an error if the structure of G is not the same
        as ``gf_struct_solver`` (for ``space=solver``) or
        ``gf_struct_sumk`` (for ``space=sumk``)..

        Parameters
        ----------
        G : BlockGf or list of BlockGf
            Green's function to check
            if it is a list, there should be as many entries as there
            are shells, and the check is performed for all shells (unless
            ish is given).
        ish : int
            shell index
            default: 0 if G is just one Green's function is given,
            check all if list of Green's functions is given
        space : 'solver' or 'sumk'
            which space the structure should correspond to
        """

        return self._check_gf_or_matrix(G, ish, space)

    def check_matrix(self, G, ish=None, space='solver'):
        """ check whether the matrix G has the right structure

        This throws an error if the structure of G is not the same
        as ``gf_struct_solver`` (for ``space=solver``) or
        ``gf_struct_sumk`` (for ``space=sumk``)..

        Parameters
        ----------
        G : dict of matrices or list of dict of matrices
            matrix to check
            if it is a list, there should be as many entries as there
            are shells, and the check is performed for all shells (unless
            ish is given).
        ish : int
            shell index
            default: 0 if G is just one matrix is given,
            check all if list of dicts is given
        space : 'solver' or 'sumk'
            which space the structure should correspond to
        """

        return self._check_gf_or_matrix(G, ish, space)

    def _check_gf_or_matrix(self, G, ish=None, space='solver'):
        if space == 'solver':
            gf_struct = self.gf_struct_solver
        elif space == 'sumk':
            gf_struct = self.gf_struct_sumk_dict
        else:
            raise Exception(
                "Argument space has to be either 'solver' or 'sumk'.")

        if isinstance(G, list):
            assert len(G) == len(gf_struct),\
                "list of G does not have the correct length"
            if ish is None:
                ishs = list(range(len(gf_struct)))
            else:
                ishs = [ish]
            for ish in ishs:
                self.check_gf(G[ish], ish=ish, space=space)
            return

        if ish is None:
            ish = 0

        if isinstance(G, BlockGf):
            for block in gf_struct[ish]:
                assert block in G.indices,\
                    "block " + block + " not in G (shell {})".format(ish)
            for block, gf in G:
                assert block in gf_struct[ish],\
                    "block " + block + " not in struct (shell {})".format(ish)
                assert list(gf.indices) == 2 * [list(map(str, gf_struct[ish][block]))],\
                    "block " + block + \
                    " has wrong indices (shell {})".format(ish)
        else:
            for block in gf_struct[ish]:
                assert block in G,\
                    "block " + block + " not in G (shell {})".format(ish)
            for block, gf in list(G.items()):
                assert block in gf_struct[ish],\
                    "block " + block + " not in struct (shell {})".format(ish)
                assert list(range(len(gf))) == 2 * [list(map(str, gf_struct[ish][block]))],\
                    "block " + block + \
                    " has wrong indices (shell {})".format(ish)

    def convert_operator(self, O, ish=0):
        """ Converts a second-quantization operator from sumk structure
            to solver structure.

        Parameters
        ----------
        O : triqs.operators.Operator
            Operator in sumk structure

        ish : int
            shell index on which the operator acts
        """

        from triqs.operators import Operator, c, c_dag

        T = self.transformation[ish]
        sk2s = self.sumk_to_solver[ish]

        O_out = Operator(0)

        for monomial in O:
            coefficient = monomial[-1]
            new_monomial = Operator(1)
            #if coefficient > 1e-10:
            for single_operator in monomial[0]:
                new_single_operator = Operator(0)
                daggered = single_operator[0]

                blockname = single_operator[1][0]
                i = single_operator[1][1]
                for j in range(len(T[blockname])):
                    if sk2s[(blockname, j)] != (None, None):
                        if daggered:
                            new_single_operator +=  (T[blockname][j,i] * c_dag(*sk2s[(blockname, j)]))
                        else:
                            new_single_operator +=  (T[blockname][j,i].conjugate() * c(*sk2s[(blockname, j)]))

                new_monomial *= new_single_operator

            O_out += new_monomial * coefficient
        return O_out



    def convert_gf(self, G, G_struct=None, ish_from=0, ish_to=None, show_warnings=True,
                   G_out=None, space_from='solver', space_to='solver', ish=None, **kwargs):
        """ Convert BlockGf from its structure to this structure.

        .. warning::

            Elements that are zero in the new structure due to
            the new block structure will be just ignored, thus
            approximated to zero.

        Parameters
        ----------
        G : BlockGf
            the Gf that should be converted
        G_struct : BlockStructure or str
            the structure of that G or None (then, this structure
            is used)
        ish_from : int
            shell index of the input structure
        ish_to : int
            shell index of the output structure; if None (the default),
            it is the same as ish_from
        show_warnings : bool or float
            whether to show warnings when elements of the Green's
            function get thrown away
            if float, set the threshold for the magnitude of an element
            about to be thrown away to trigger a warning
            (default: 1.e-10)
        G_out : BlockGf
            the output Green's function (if not given, a new one is
            created)
        space_from : 'solver' or 'sumk'
            whether the Green's function ``G`` corresponds to the
            solver or sumk structure of ``G_struct``
        space_to : 'solver' or 'sumk'
            whether the output Green's function should be according to
            the solver of sumk structure of this structure
        **kwargs :
            options passed to the constructor for the new Gf
        """

        if ish is not None:
            warn(
                'The parameter ish in convert_gf is deprecated. Use ish_from and ish_to instead.')
            ish_from = ish
            ish_to = ish
        return self._convert_gf_or_matrix(G, G_struct, ish_from, ish_to,
                                          show_warnings, G_out, space_from, space_to, **kwargs)

    def convert_matrix(self, G, G_struct=None, ish_from=0, ish_to=None, show_warnings=True,
                       G_out=None, space_from='solver', space_to='solver'):
        """ Convert matrix from its structure to this structure.

        .. warning::

            Elements that are zero in the new structure due to
            the new block structure will be just ignored, thus
            approximated to zero.

        Parameters
        ----------
        G : dict of numpy array
            the matrix that should be converted
        G_struct : BlockStructure or str
            the structure of that G or None (then, this structure
            is used)
        ish_from : int
            shell index of the input structure
        ish_to : int
            shell index of the output structure; if None (the default),
            it is the same as ish_from
        show_warnings : bool or float
            whether to show warnings when elements of the Green's
            function get thrown away
            if float, set the threshold for the magnitude of an element
            about to be thrown away to trigger a warning
            (default: 1.e-10)
        G_out : dict of numpy array
            the output numpy array (if not given, a new one is
            created)
        space_from : 'solver' or 'sumk'
            whether the matrix ``G`` corresponds to the
            solver or sumk structure of ``G_struct``
        space_to : 'solver' or 'sumk'
            whether the output matrix should be according to
            the solver of sumk structure of this structure
        **kwargs :
            options passed to the constructor for the new Gf
        """

        return self._convert_gf_or_matrix(G, G_struct, ish_from, ish_to,
                                          show_warnings, G_out, space_from, space_to)

    def _convert_gf_or_matrix(self, G, G_struct=None, ish_from=0, ish_to=None, show_warnings=True,
                              G_out=None, space_from='solver', space_to='solver', **kwargs):
        if ish_to is None:
            ish_to = ish_from

        warning_threshold = 1.e-10
        if isinstance(show_warnings, float):
            warning_threshold = show_warnings
            show_warnings = True

        if G_struct is None:
            G_struct = self

        if space_from == 'solver':
            gf_struct_from = G_struct.gf_struct_solver[ish_from]
            eff_trans_from = G_struct.effective_transformation_solver[ish_from]
            block_mapping_from = G_struct.sumk_to_solver_block[ish_from]
        elif space_from == 'sumk':
            gf_struct_from = G_struct.gf_struct_sumk_dict[ish_from]
            eff_trans_from = {block: np.eye(len(indices))
                              for block, indices in G_struct.gf_struct_sumk[ish_from]}
            block_mapping_from = {b: [b] for b in gf_struct_from}
        else:
            raise Exception(
                "Argument space_from has to be either 'solver' or 'sumk'.")

        if space_to == 'solver':
            gf_struct_to = self.gf_struct_solver[ish_to]
            eff_trans_to = self.effective_transformation_solver[ish_to]
            block_mapping_to = self.solver_to_sumk_block[ish_to]
        elif space_to == 'sumk':
            gf_struct_to = self.gf_struct_sumk_dict[ish_to]
            eff_trans_to = {block: np.eye(len(indices))
                            for block, indices in self.gf_struct_sumk_list[ish_to]}
            block_mapping_to = {b: b for b in gf_struct_to}
        else:
            raise Exception(
                "Argument space_to has to be either 'solver' or 'sumk'.")

        if isinstance(G, BlockGf):
            # create a Green's function to hold the result
            if G_out is None:
                if not 'mesh' in kwargs and not 'beta' in kwargs:
                    kwargs['mesh'] = G.mesh
                G_out = self.create_gf(ish=ish_to, space=space_to, **kwargs)
            else:
                self.check_gf(G_out, ish=ish_to, space=space_to)
        elif isinstance(G, dict):
            if G_out is None:
                G_out = self.create_matrix(ish=ish_to, space=space_to)
            else:
                self.check_matrix(G_out, ish=ish_to, space=space_to)
        else:
            raise Exception('G is neither BlockGf nor dict.')

        for block_to in list(gf_struct_to.keys()):
            if isinstance(G, BlockGf):
                G_out[block_to].zero()
            else:
                G_out[block_to][:] = 0.0
            block_intermediate = block_mapping_to[block_to]
            block_from = block_mapping_from[block_intermediate]
            T_to = eff_trans_to[block_to]
            g_help = G_out[block_to].copy()
            for block in block_from:
                T_from = eff_trans_from[block]
                if isinstance(G, BlockGf):
                    g_help.from_L_G_R(np.dot(T_to, np.conjugate(np.transpose(T_from))),
                                      G[block],
                                      np.dot(T_from, np.conjugate(np.transpose(T_to))))
                    G_out[block_to] << G_out[block_to] + g_help
                else:
                    g_help = np.dot(np.dot(T_to, np.conjugate(np.transpose(T_from))),
                                    np.dot(G[block],
                                           np.dot(T_from, np.conjugate(np.transpose(T_to)))))
                    G_out[block_to] += g_help

        if show_warnings:
            # we back-transform it
            G_back = G_struct._convert_gf_or_matrix(G_out, self, ish_from=ish_to,
                                                    ish_to=ish_from,
                                                    show_warnings=False,  # else we get an endless loop
                                                    space_from=space_to, space_to=space_from, **kwargs)
            for name, gf in (G_back if isinstance(G, BlockGf) else list(G_back.items())):
                if isinstance(G, BlockGf):
                    maxdiff = np.max(np.abs(G_back[name].data - G[name].data),
                                     axis=0)
                else:
                    maxdiff = G_back[name] - G[name]

                if space_to == 'solver' and self == G_struct: # do comparison in solver (ignore diff. in ignored orbitals)
                    tmp = self.create_matrix(space='sumk')
                    tmp[name] = maxdiff
                    maxdiff = G_struct._convert_gf_or_matrix(tmp, self, ish_from=ish_from,
                                                    ish_to=ish_to,
                                                    show_warnings=False,
                                                    space_from=space_from, space_to=space_to, **kwargs)

                    for block in maxdiff:
                        maxdiff_b = maxdiff[block]
                        if np.any(maxdiff_b > warning_threshold):
                            warn('Block {} maximum difference:\n'.format(name) + str(maxdiff))


                elif np.any(maxdiff > warning_threshold):
                    warn('Block {} maximum difference:\n'.format(name)
                         + str(maxdiff))

        return G_out

    def approximate_as_diagonal(self):
        """ Create a structure for a GF with zero off-diagonal elements.

        .. warning::

            In general, this will throw away non-zero elements of the
            Green's function. Be sure to verify whether this approximation
            is justified.
        """

        self.gf_struct_solver=[]
        self.solver_to_sumk=[]
        self.solver_to_sumk_block=[]
        for ish in range(len(self.sumk_to_solver)):
            self.gf_struct_solver.append({})
            self.solver_to_sumk.append({})
            self.solver_to_sumk_block.append({})
            for frm,to in list(self.sumk_to_solver[ish].items()):
                if to[0] is not None:
                    self.gf_struct_solver[ish][frm[0]+'_'+str(frm[1])]=[0]
                    self.sumk_to_solver[ish][frm]=(frm[0]+'_'+str(frm[1]),0)
                    self.solver_to_sumk[ish][(frm[0]+'_'+str(frm[1]),0)]=frm
                    self.solver_to_sumk_block[ish][frm[0]+'_'+str(frm[1])]=frm[0]

    def __eq__(self,other):
        def compare(one,two):
            if type(one)!=type(two):
                if not (isinstance(one, (bool, np.bool_)) and isinstance(two, (bool, np.bool_))):
                    return False
            if one is None and two is None:
                return True
            if isinstance(one,list) or isinstance(one,tuple):
                if len(one) != len(two):
                    return False
                for x,y in zip(one,two):
                    if not compare(x,y):
                        return False
                return True
            elif isinstance(one,(int,bool, str, np.bool_)):
                return one==two
            elif isinstance(one,np.ndarray):
                return np.all(one==two)
            elif isinstance(one,dict):
                if set(one.keys()) != set(two.keys()):
                    return False
                for k in set(one.keys()).intersection(list(two.keys())):
                    if not compare(one[k],two[k]):
                        return False
                return True
            warn('Cannot compare {}'.format(type(one)))
            return False

        for prop in [ "gf_struct_sumk", "gf_struct_solver",
                "solver_to_sumk", "sumk_to_solver", "solver_to_sumk_block",
                "deg_shells","transformation", "corr_to_inequiv"]:
            if not compare(getattr(self,prop),getattr(other,prop)):
                return False
        return True

    def copy(self):
        return copy.deepcopy(self)

    def __reduce_to_dict__(self):
        """ Reduce to dict for HDF5 export."""

        ret = {}
        for element in [ "gf_struct_sumk", "gf_struct_solver",
                         "solver_to_sumk_block","deg_shells",
                         "transformation", "corr_to_inequiv"]:
            ret[element] = getattr(self,element)
            if ret[element] is None:
                ret[element] = 'None'

        if ret["transformation"] is None:
            ret["transformation"] = "None"

        def construct_mapping(mapping):
            d = []
            for ish in range(len(mapping)):
                d.append({})
                for k,v in list(mapping[ish].items()):
                    d[ish][repr(k)] = repr(v)
            return d

        ret['solver_to_sumk']=construct_mapping(self.solver_to_sumk)
        ret['sumk_to_solver']=construct_mapping(self.sumk_to_solver)
        return ret

    @classmethod
    def __factory_from_dict__(cls,name,D) :
        """ Create from dict for HDF5 import."""

        def reconstruct_mapping(mapping):
            d = []
            for ish in range(len(mapping)):
                d.append({})
                for k,v in list(mapping[ish].items()):
                    # literal_eval is a saje alternative to eval
                    d[ish][literal_eval(k)] = literal_eval(v)
            return d

        for elem in D:
            if D[elem]=="None":
                D[elem] = None

        D['solver_to_sumk']=reconstruct_mapping(D['solver_to_sumk'])
        D['sumk_to_solver']=reconstruct_mapping(D['sumk_to_solver'])
        return cls(**D)

    def __str__(self):
        s=''
        s+= "corr_to_inequiv "+str(self.corr_to_inequiv)+'\n'
        s+= "gf_struct_sumk "+str(self.gf_struct_sumk)+'\n'
        s+= "gf_struct_solver "+str(self.gf_struct_solver)+'\n'
        s+= "solver_to_sumk_block "+str(self.solver_to_sumk_block)+'\n'
        for el in ['solver_to_sumk','sumk_to_solver']:
            s+=el+'\n'
            element=getattr(self,el)
            for ish in range(len(element)):
                s+=' shell '+str(ish)+'\n'
                def keyfun(el):
                    return '{}_{:05d}'.format(el[0],el[1])
                keys = sorted(list(element[ish].keys()),key=keyfun)
                for k in keys:
                    s+='  '+str(k)+str(element[ish][k])+'\n'
        s += "deg_shells\n"
        for ish in range(len(self.deg_shells)):
            s+=' shell '+str(ish)+'\n'
            for l in range(len(self.deg_shells[ish])):
                s+='  equivalent group '+str(l)+'\n'
                if isinstance(self.deg_shells[ish][l],dict):
                    for key, val in list(self.deg_shells[ish][l].items()):
                        s+='   '+key+('*' if val[1] else '')+':\n'
                        s+='    '+str(val[0]).replace('\n','\n    ')+'\n'
                else:
                    for key in self.deg_shells[ish][l]:
                        s+='   '+key+'\n'
        s += "transformation\n"
        s += str(self.transformation)
        return s

from h5.formats import register_class
register_class(BlockStructure)
