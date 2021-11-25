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


import sys,numpy
from importlib import import_module
from .archive_basic_layer import HDFArchiveGroupBasicLayer
from .formats import register_class, register_backward_compatibility_method, get_format_info
import h5py

vls_dt = h5py.string_dtype(encoding='utf-8')

# -------------------------------------------
#
#  Various wrappers for basic python types.
#
# --------------------------------------------
class List:
    def __init__(self,ob) :
        self.ob = ob
    def __reduce_to_dict__(self) :
        return {str(n):v for n,v in enumerate(self.ob)}
    @classmethod
    def __factory_from_dict__(cls, name, D) :
        return [x for n,x in sorted([(int(n), x) for n,x in list(D.items())])]

class Tuple:
    def __init__(self,ob) :
        self.ob = ob
    def __reduce_to_dict__(self) :
        return {str(n):v for n,v in enumerate(self.ob)}
    @classmethod
    def __factory_from_dict__(cls, name, D) :
        return tuple(x for n,x in sorted([(int(n), x) for n,x in list(D.items())]))

class Dict:
    def __init__(self,ob) :
        self.ob = ob
    def __reduce_to_dict__(self) :
        return {str(n):v for n,v in list(self.ob.items())}
    @classmethod
    def __factory_from_dict__(cls, name, D) :
        return {n:x for n,x in list(D.items())}

register_class(List)
register_backward_compatibility_method('PythonListWrap', 'List')

register_class(Tuple)
register_backward_compatibility_method('PythonTupleWrap', 'Tuple')

register_class(Dict)
register_backward_compatibility_method('PythonDictWrap', 'Dict')

# -------------------------------------------
#
#  A view of a subgroup of the archive
#
# --------------------------------------------

class HDFArchiveGroup(HDFArchiveGroupBasicLayer):
    """
    """
    _wrappedType = {
        list : List,
        tuple : Tuple,
        dict : Dict
    }
    _MaxLengthKey = 500

    def __init__(self, parent, subpath) :
        # We want to hold a reference to the parent group, if we are not at the root
        # This will prevent a premature destruction of the root HDFArchive object
        if not self is parent: self.parent = parent
        self.options = parent.options
        HDFArchiveGroupBasicLayer.__init__(self, parent, subpath)
        self.options = parent.options
        self.key_as_string_only = self.options['key_as_string_only']
        self._reconstruct_python_objects = self.options['reconstruct_python_object']
        self.is_top_level = False

    #-------------------------------------------------------------------------
    def __contains__(self,key) :
        return key in list(self.keys())

    #-------------------------------------------------------------------------
    def values(self) :
        """
        Generator returning the values in the group
        """
        def res() :
            for name in list(self.keys()) :
                yield self[name]
        return res()

   #-------------------------------------------------------------------------
    def items(self) :
        """
        Generator returning couples (key, values) in the group.
        """
        def res() :
            for name in list(self.keys()):
                yield name, self[name]
        return res()

    #-------------------------------------------------------------------------
    def __iter__(self) :
        """Returns the keys, like a dictionary"""
        def res() :
            for name in list(self.keys()) :
                yield name
        return res()

    #-------------------------------------------------------------------------
    def __len__(self) :
        """Returns the length of the keys list """
        return  len(list(self.keys()))

    #-------------------------------------------------------------------------
    def update(self,object_with_dict_protocol):
        for k,v in list(object_with_dict_protocol.items()) : self[k] = v

    #-------------------------------------------------------------------------
    def __delitem__(self,key) :
        self._clean_key(key,True)

    #-------------------------------------------------------------------------
    def __setitem__(self,key,val) :
        assert '/' not in key, "/ can not be part of a key"

        if key in list(self.keys()) :
            self._clean_key(key) # clean things

        # Transform list, dict, etc... into a wrapped type that will allow HDF reduction
        if type(val) in self._wrappedType: val = self._wrappedType[type(val)](val)

        # write the attributes
        def write_attributes(g) :
           """Use the _hdf5_format_ if it exists otherwise the class name"""
           ds = val._hdf5_format_ if hasattr(val,"_hdf5_format_") else val.__class__.__name__
           try :
             get_format_info(ds)
           except :
             err = """
               You are trying to store an object of type "%s", with the format "%s".
               This format is not registered, so you will not be able to reread the class.
               Didn't you forget to register your class in h5.formats?
               """ %(val.__class__.__name__,ds)
             raise IOError(err)
           g.write_attr("Format", ds)

        if hasattr(val,'__write_hdf5__') : # simplest protocol
            try:
                val.__write_hdf5__(self, key)
            except Exception as e:
                raise RuntimeError(
                    f"Error in writing a Python object of type {type(val)}! " +
                    f"Thrown error was '{e}'"
                )
            self.cached_keys.append(key) # I need to do this here
        elif hasattr(val,'__reduce_to_dict__') : # Is it a HDF_compliant object
            self.create_group(key) # create a new group
            d = val.__reduce_to_dict__()
            if not isinstance(d,dict) : raise ValueError(" __reduce_to_dict__ method does not return a dict. See the doc !")
            SubGroup = HDFArchiveGroup(self,key)
            for k, v in list(d.items()) : SubGroup[k] = v
            write_attributes(SubGroup)
        elif isinstance(val,numpy.ndarray) : # it is a numpy
            try :
               self._write( key, numpy.array(val,copy=1,order='C') )
            except RuntimeError:
               print("HDFArchive is in trouble with the array %s"%val)
               raise
        elif isinstance(val, HDFArchiveGroup) : # will copy the group recursively
            # we could add this for any object that has .items() in fact...
            SubGroup = HDFArchiveGroup(self, key)
            for k,v in list(val.items()) : SubGroup[k]=v
        else : # anything else... expected to be a scalar
            try :
               self._write( key, val)
               self.cached_keys.append(key) # I need to do this here
            except:
               raise #ValueError, "Value %s\n is not of a type suitable to storage in HDF file"%val
        self._flush()

    #-------------------------------------------------------------------------
    #def get_raw (self,key):
        #"""Similar to __getitem__ but it does NOT reconstruct the python object,
        #it presents it as a subgroup"""
        #return self.__getitem1__(key,False)

    #-------------------------------------------------------------------------
    def __getitem__(self,key) :
        """Return the object key, possibly reconstructed as a python object if
        it has been properly set up"""
        # If the key contains /, grabs the subgroups
        if '/' in key:
            a,l =self, key.split('/')
            for s in l[:-1]: a = a.get_raw(s)
            return a[l[-1]]
        return self.__getitem1__(key,self._reconstruct_python_objects)

    #-------------------------------------------------------------------------
    def __getitem1__(self, key, reconstruct_python_object, hdf5_format = None) :

        if key not in self :
            raise KeyError("Key %s does not exist."%key)

        if self.is_group(key) :
            SubGroup = HDFArchiveGroup(self,key) # View of the subgroup
            bare_return = lambda: SubGroup
        elif self.is_data(key) :
            bare_return = lambda: self._read(key)
        else :
            raise KeyError("Key %s is of unknown type !!"%key)

        if not reconstruct_python_object : return bare_return()

        # try to find the format
        if hdf5_format is None:
            hdf5_format = self.read_hdf5_format_from_key(key)
            if hdf5_format == "":
                return bare_return()
        
        try :
            fmt_info = get_format_info(hdf5_format)
        except KeyError:
            print("Warning : The hdf5 format %s is not recognized. Returning as a group. Hint : did you forgot to import this python class ?"%hdf5_format)
            return bare_return()

        r_class_name  = fmt_info.classname
        r_module_name = fmt_info.modulename
        r_readfun = fmt_info.read_fun
        if not (r_class_name and r_module_name) : return bare_return()
        try:
            r_class = getattr(import_module(r_module_name),r_class_name)
        except KeyError:
            raise RuntimeError("I cannot find the class %s to reconstruct the object !"%r_class_name)
        if r_readfun:
            return r_readfun(self._group, key)
        if hasattr(r_class,"__factory_from_dict__"):
            assert self.is_group(key), "__factory_from_dict__ requires a subgroup"
            reconstruct = lambda k: SubGroup.__getitem1__(k, reconstruct_python_object, fmt_info.backward_compat.get(k, None))
            values = {k: reconstruct(k) for k in SubGroup}
            return r_class.__factory_from_dict__(key, values)

        raise ValueError("Impossible to reread the class %s for group %s and key %s"%(r_class_name,self, key))

    #---------------------------------------------------------------------------
    def __str__(self) :
        def pr(name) :
            if self.is_group(name) :
                return "%s : subgroup"%name
            elif self.is_data(name) : # can be an array of a number
                return "%s : data "%name
            else :
                raise ValueError("oopps %s"%name)

        s= "HDFArchive%s with the following content:\n"%(" (partial view)" if self.is_top_level else '')
        s+='\n'.join([ '  '+ pr(n) for n in list(self.keys()) ])
        return s

    #-------------------------------------------------------------------------
    def __repr__(self) :
        return self.__str__()

    #-------------------------------------------------------------------------
    def apply_on_leaves (self,f) :
        """
           For each named leaf (name,value) of the tree, it calls f(name,value)
           f should return :
            - `None`                    : no action is taken
            - an `empty tuple` ()       : the leaf is removed from the tree
            - an hdf-compliant value    : the leaf is replaced by the value
        """
        def visit_tree(n,d):
          for k in d:# Loop over the subgroups in d
              if d.is_group(k) : visit_tree(k,d[k])
              else :
                  r = f(k,d[k])
                  if not r is None : d[k] = r
                  elif r == () : del d[k]
        visit_tree('/',self['/'])

    # These two methods are necessary for "with"
    def __enter__(self): return self
    def __exit__(self, type, value, traceback): pass

# -------------------------------------------
#
#  The main class
#
# --------------------------------------------

class HDFArchive(HDFArchiveGroup):
    """
    """
    _class_version = 1

    def __init__(self, url_name, open_flag = 'a', key_as_string_only = True,
            reconstruct_python_object = True, init = {}):
        r"""
           Parameters
           -----------
           url_name : string
             The url of the hdf5 file.

                  * If url is a simple string, it is interpreted as a local file name

                  * If url is a remote url (e.g. `http://ipht.cea.fr/triqs/data/single_site_bethe.output.h5` )
                    then the h5 file is downloaded in temporary file and opened.
                    In that case, ``open_flag`` must be 'r', read-only mode.
                    The temporary file is deleted at exit.
           open_flag : Legal modes: r, w, a (default)
           key_as_string_only : True (default)
           init : any generator of tuple (key,val), e.g. a dict.items().
             It will fill the archive with these values.

           Attributes
           ----------
           LocalFileName : string
             the name of the file or of the local downloaded copy
           url_name : string
             the name of the Url

           Examples
           --------
           >>> # retrieve a remove archive (in read-only mode) :
           >>> h = HDFArchive( 'http://ipht.cea.fr/triqs/data/single_site_bethe.output.h5')
           >>>
           >>> # full copy of an archive
           >>> HDFArchive( f, 'w', init = HDFArchive(fmp,'r').items())  # full
           >>>
           >>> # partial copy of file of name fmp, with only the key 'G'
           >>> HDFArchive( f, 'w', init = [ (k,v) for (k,v) in HDFArchive(fmp,'r') if k in ['G'] )
           >>>
           >>> # faster version : the object are only retrieved when needed (list comprehension vs iterator comprehension)
           >>> HDFArchive( f, 'w', init = ( (k,v) for (k,v) in HDFArchive(fmp,'r') if k in ['G'] ) )
           >>>
           >>> # partial copy with processing on the fly with the P function
           >>> HDFArchive( f, 'w', init = ( (k,P(v)) for (k,v) in HDFArchive(fmp,'r') if k in ['G'] ) )
           >>>
           >>> # another variant with a filtered dict
           >>> HDFArchive( f, 'w', init = HDFArchive(fmp,'r').items(lambda k :  k in ['G'] ))

        """
        import os,os.path
        assert open_flag in ['r','w','a'], "Invalid mode"
        assert isinstance(url_name,str), "url_name must be a string"

        # If it is a url, retrieve it and check mode is read only
        #import urllib.request
        #try:
            #LocalFileName, http_message = urllib.request.urlretrieve(url_name)
            ## a url must be read only
            #assert open_flag == 'r', "You retrieve a distant Url %s which is not local, so it must be read-only. Use 'r' option"%url_name
        #except ValueError: # Not a valid URL -> Local File
            #LocalFileName, http_message = url_name, None

        LocalFileName, http_message = url_name, None

        if open_flag == 'w':
            # destroys the file, ignoring errors
            try: os.remove(os.path.abspath(LocalFileName))
            except OSError: pass

        self._init_root(LocalFileName, open_flag)
        self.options = {'key_as_string_only' : key_as_string_only,
                        'do_not_overwrite_entries' : False,
                        'reconstruct_python_object': reconstruct_python_object,
                        'UseAlpsNotationForComplex'  : True
                        }
        HDFArchiveGroup.__init__(self,self,"")
        self.is_top_level = True
        for k,v in init : self[k]=v

    def __del__(self):
      # We must ensure the root group is closed before closing the file
      if hasattr(self, '_group'):
          self._flush()
          del self._group

    # These two methods are necessary for "with"
    def __enter__(self): return self

    def __exit__(self, type, value, traceback):
      self._flush()
      del self._group

#--------------------------------------------------------------------------------

class HDFArchiveInert:
    """
    A fake class for the node in MPI. It does nothing, but
    permits to write simply :
       a= mpi.bcast(H['a']) # run on all nodes
    -[] : __getitem__ returns self so that H['a']['b'] is ok...
    - setitem : does nothing.
    """
    def HDFArchive_Inert(self):
        pass
    def __getitem__(self,x)   : return self
    def __setitem__(self,k,v) : pass

#--------------------------------------------------------------------------------


