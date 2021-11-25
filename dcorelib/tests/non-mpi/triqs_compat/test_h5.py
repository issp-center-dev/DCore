import numpy as np
from dcorelib.triqs_compat.h5 import HDFArchive as HDFArchive2
import pytest

triqs_h5_available = False
try:
    from h5 import HDFArchive
    triqs_h5_available = True
except ImportError:
    pass


# https://numpy.org/doc/stable/user/basics.types.html
int_types = [int, np.int_, bool, np.bool_]

@pytest.mark.skipif(not triqs_h5_available, reason="TRIQS h5 is not installed.")
class TestClass:
    @pytest.mark.parametrize("dtype", [np.float64, np.complex128])
    def test_rw_array(self, dtype):
        shape = (3, 2)
        a = np.array(np.random.randn(*shape), dtype=dtype)
    
        with HDFArchive('test.h5', 'w') as f:
            f['data'] = a
    
        with HDFArchive2('test2.h5', 'w') as f:
            f['data'] = a
    
        with HDFArchive2('test.h5', 'r') as f:
            assert np.array_equal(a, f['data'])
    
        with HDFArchive('test2.h5', 'r') as f:
                assert np.array_equal(a, f['data'])
    
    def test_rw_list(self):
        val = [1, 2, 3]
        with HDFArchive('test.h5', 'w') as f:
            f['data'] = val
    
        with HDFArchive2('test2.h5', 'w') as f:
            f['data'] = val
    
        for _Archive in [HDFArchive, HDFArchive2]:
            for _file in ['test.h5', 'test2.h5']:
                with _Archive(_file, 'r') as f:
                    val_reconst = f['data']
                    assert val == val_reconst
    
    
    def test_rw_dict(self):
        val = {'1':-1, 'A': 'b'}
        for t in int_types:
            val[t.__name__] = t(1)
        with HDFArchive('test_dict.h5', 'w') as f:
            f['data'] = val
    
        with HDFArchive2('test_dict2.h5', 'w') as f:
            f['data'] = val
    
        for _Archive in [HDFArchive, HDFArchive2]:
            for _file in ['test_dict.h5', 'test_dict2.h5']:
                with _Archive(_file, 'r') as f:
                    val_reconst = f['data']
                    assert val == val_reconst
    
def test_write_strings():
    with HDFArchive2('test_strings.h5', 'w') as f:
        f['data'] = np.array(['a', 'b']).astype('S')

def test_subgroup():
    with HDFArchive2('test.h5', 'w') as f:
        f.create_group('subgrp')
        print(type(f['subgrp']))

#def test_bool():
    #with h5py.File('test_bool.h5', 'w') as f:
        ##f['bool'] = True