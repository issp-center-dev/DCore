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
        make_copies = kwargs.pop('make_copies', False)

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
        return self.g_list[key]

    def __setitem__(self, key, val):
        self.g_list[key] << val

    def _first(self):
        return self.g_list[0]

    def __len__(self):
        return len(self.g_list)
    
    @property
    def beta(self):
        return self._first().beta

    @property
    def n_blocks(self):
        """ Number of blocks"""
        return len(self.g_list)
    
    #def __reduce_to_dict__(self):
        #d = dict(self)
        #d['block_names'] = list(self.indices)
        #return d

    def __write_hdf5__(self, group, key):
        """ Write to a HDF5 file"""
        group[key].write_attr('Format', 'BlockGf')
        group[key]['block_names'] = self.block_names
        for name, g in self:
            group[key][name] = g


from dcore.triqs_compat.h5.formats import register_class
register_class (BlockGf)