from ..h5.archive import register_class
from .gf import Gf, LazyExpression
from copy import deepcopy
import numpy as np
import operator

class BlockGf:
    """
    Generic Green's Function by Block.
    """
    def __init__(self, **kwargs):
        """
        * BlockGf(name_list = list of names, block_list = list of blocks, make_copies = False, name = '')

           * ``name_list``: list of the names of the blocks, e.g. ["up","down"].
           * ``block_list``: list of blocks of Green's functions.
           * ``make_copies``: If True, it makes a copy of the blocks and build the Green's function from these copies.
           * ``gf_struct``: 
           * ``mesh``: 
           * ``name_block_generator``: 
         """
        self.name = kwargs.pop('name', 'G')

        if 'name_block_generator' in kwargs:
            self.block_names = []
            self.g_list = []
            for name, block in kwargs['name_block_generator']:
                self.block_names.append(name)
                self.g_list.append(block)
        else:
            if 'name_list' in kwargs:
                self.block_names = kwargs['name_list']
            else:
                self.block_names = [bl[0] for bl in kwargs['gf_struct']]
    
            if 'block_list' in kwargs:
                self.g_list = kwargs['block_list']
            else:
                mesh = kwargs['mesh']
                indices = [bl[1] for bl in kwargs['gf_struct']]
                self.g_list = [
                    Gf(beta=mesh.beta, statistic=mesh.statistic,
                       mesh=mesh, indices=indices_)
                    for indices_ in indices
                ]
        
        if isinstance(self.g_list, Gf):
            self.g_list = [self.g_list]

        make_copies = kwargs.pop('make_copies', False)
        if make_copies:
            self.g_list = [deepcopy(b) for b in self.g_list]

        self.g_dict = {k: v for k, v in zip(self.block_names, self.g_list)}

        if isinstance(self.block_names, tuple):
            self.block_names = list(self.block_names)
        if isinstance(self.g_list, tuple):
            self.g_list = list(self.g_list)
        assert isinstance(self.block_names, list)
        for name in self.block_names:
            assert isinstance(name, str)
        assert isinstance(self.g_list, list)
        for g in self.g_list:
            assert isinstance(g, Gf)
        
        self._sanity_check()
    
    def _sanity_check(self):
        for ib, name in enumerate(self.block_names):
            assert id(self.g_dict[name]) == id(self.g_list[ib])
    
    @property
    def indices(self):
        return self.block_names
    
    def __iter__(self):
        return zip(self.block_names, self.g_list)

    def __getitem__(self, key):
        return self.g_dict[key]

    def __setitem__(self, key, val):
        self.g_dict[key] << val

    def _first(self):
        return self.g_list[0]

    def __len__(self):
        return len(self.g_list)

    def __lshift__(self, other):
        if isinstance(other, BlockGf):
            for name, g in self:
                g << other[name]
        elif isinstance(other, LazyExpression):
            for name, g in self:
                g << other
        else:
            raise RuntimeError("Invalid other!")

    def __iadd__(self, other):
        assert type(other) in [BlockGf, list]
        for bl, bl2 in zip(self.g_list, other):
            bl2_ = bl2
            if isinstance(bl2_, tuple) and len(bl2_) == 2:
                bl2_ = bl2_[1]
            bl += bl2_
        return self

    def __add__(self, other):
        return self.__add_sub__(other, operator.iadd)

    def __sub__(self, other):
        return self.__add_sub__(other, operator.isub)

    def __add_sub__(self, other, op):
        assert type(other) in [BlockGf, list], f"Invalid type{type(other)}"
        res = self.copy()
        for bl, bl2 in zip(res.g_list, other):
            bl2_ = bl2
            if isinstance(bl2_, tuple) and len(bl2_) == 2:
                bl2_ = bl2_[1]
            op(bl, bl2_)
        return res
    
    def __isub__(self, other):
        assert type(other) in [BlockGf, list]
        for bl, bl2 in zip(self.g_list, other):
            bl2_ = bl2
            if isinstance(bl2_, tuple) and len(bl2_) == 2:
                bl2_ = bl2_[1]
            bl -= bl2_
        return self

    def __mul__(self, other):
        if not np.isscalar(other):
            return NotImplemented
        bg = self.copy()
        for name, bl in bg:
            bl.data *= other
        return bg
    
    __rmul__ = __mul__

    def __imul__(self, other):
        if not np.isscalar(other):
            return NotImplemented
        for g in self.g_list:
            g.data *= other
        return self

    @property
    def beta(self):
        return self._first().beta

    @property
    def mesh(self):
        return self._first().mesh

    @property
    def n_blocks(self):
        """ Number of blocks"""
        return len(self.g_list)

    def copy(self):
        """ Return a deep copy of self """
        return deepcopy(self)
    
    def zero(self):
        """ Return fill all blocks with zero """
        for g in self.g_list:
            g.zero()

    def __write_hdf5__(self, group, key):
        """ Write to a HDF5 file"""
        group.create_group(key)
        group[key]['block_names'] = self.block_names
        for name, g in self:
            group[key][name] = g
        group[key].write_attr('Format', 'BlockGf')

    @classmethod
    def __factory_from_dict__(cls, key, dict) :
        block_names = dict['block_names']
        return cls(
            name_list = block_names,
            block_list = [dict[name] for name in block_names]
        )

    def invert(self):
        """ Invert in place """
        for g in self.g_list:
            g.invert()
    
    def inverse(self):
        block_list = [g.inverse() for g in self.g_list]
        return BlockGf(name_list=self.block_names, block_list=block_list)

    def density(self):
        return {name: bl.density() for name, bl in self}

    def total_density(self):
        dense_ = self.density()
        return np.real(np.sum([np.trace(v) for k, v in dense_.items()]))

register_class (BlockGf)