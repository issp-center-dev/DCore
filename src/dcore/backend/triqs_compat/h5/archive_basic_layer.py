# Copyright (c) 2019-2020 Simons Foundation
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http:#www.apache.org/licenses/LICENSE-2.0.txt
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


from .archive import *
import h5py
import numpy as np

def _has_dataset(group, key):
    """ Return if group has a dataset with name 'key' """
    if key not in group:
        return False
    return isinstance(group[key], h5py.Dataset)

def _has_subgroup(group, key):
    """ Return if group has a subgroup with name 'key' """
    if key not in group:
        return False
    return isinstance(group[key], h5py.Group)

def _from_numpy_type(s):
    if isinstance(s, bytes):
        return s.decode('utf-8')
    elif isinstance(s, np.integer):
        return int(s)
    elif isinstance(s, np.floating):
        return float(s)
    elif isinstance(s, np.complexfloating):
        return complex(s)
    else:
        return s
 
class HDFArchiveGroupBasicLayer:
    _class_version = 1

    def __init__(self, parent, subpath ):
        """  """
        self.options = parent.options
        self._group = parent._group[subpath] if subpath else parent._group
        self.ignored_keys = [] 
        self.cached_keys = list(self._group.keys())

    def _init_root(self, LocalFileName, open_flag, libver=('v110')):
        try :
            if open_flag == 'r':
                fich = h5py.File(LocalFileName, open_flag, libver=libver, swmr=True)
            else:
                fich = h5py.File(LocalFileName, open_flag, libver=libver)
                try:
                    fich.swmr_mode = True
                except:
                    pass
        except :
            print("Cannot open the HDF file %s"%LocalFileName)
            raise
        self._group = fich['/']

    def is_group(self,p) :
        """Is p a subgroup ?"""
        assert len(p)>0 and p[0]!='/'
        return p in self.cached_keys and _has_subgroup(self._group, p)

    def is_data(self,p) :
        """Is p a leaf ?"""
        assert len(p)>0 and p[0]!='/'
        return p in self.cached_keys and _has_dataset(self._group, p)

    def write_attr (self, key, val) :
        self._group.attrs[key] = val

    def read_attr (self, key) :
        return _from_numpy_type(self._group.attrs[key])

    def _read (self, key):
        if '__complex__' in self._group[key].attrs and \
            int(self._group[key].attrs['__complex__']) == 1:
            # For compatibility with TRIQS
            val = self._group[key][()]
            assert val.shape[-1] == 2
            return val.view(np.complex128).reshape(val.shape[:-1])
        val = _from_numpy_type(self._group[key][()])
        return val

    def _write(self, key, val) :
        if isinstance(val, np.ndarray) and np.iscomplexobj(val):
            # For compatibility with TRIQS
            self._group[key] = val.view(float).reshape(val.shape +(2,))
            self._group[key].attrs['__complex__'] = 1
        elif isinstance(val, bool):
            self._group[key] = numpy.bool_(val)
        elif isinstance(val, int) or issubclass(type(val), np.integer):
            self._group[key] = np.int_(val)
        else:
            self._group[key] = val
        self.cached_keys.append(key)

    def _raw_write(self, key, val) :
        self._group[key] = val
        self.cached_keys.append(key)
    
    def _flush(self):
        pass

    def create_group (self,key):
        self._group.create_group(key)
        self.cached_keys.append(key)

    def keys(self) :
        return self.cached_keys

    def _clean_key(self,key, report_error=False) :
        if report_error and key not in self.keys() :
             raise KeyError("Key %s is not in archive !!"%key)
        if key in self.cached_keys :
          del self._group[key]
          self.cached_keys.remove(key)
        else: raise KeyError("Key %s is not in archive !!"%key)

    def read_hdf5_format_from_key(self, key):
        if 'Format' in self._group[key].attrs:
            return _from_numpy_type(self._group[key].attrs['Format'])
        elif 'TRIQS_HDF5_data_scheme' in self._group[key].attrs:
            return _from_numpy_type(self._group[key].attrs['TRIQS_HDF5_data_scheme'])
        else:
            return ''

