from ..h5.archive import register_class
from .gf import Gf
from copy import deepcopy
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
         """
        self.name = kwargs.pop('name', 'G')
        self.block_names = kwargs['name_list']
        self.g_list = kwargs['block_list']
        self.g_dict = {k: v for k, v in zip(self.block_names, self.g_list)}
        make_copies = kwargs.pop('make_copies', False)

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
        
        if make_copies:
            self.g_list = [deepcopy(b) for b in self.g_list]
    
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
        assert isinstance(other, BlockGf)
        for name, g in self:
            g << other[name]
    
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


register_class (BlockGf)