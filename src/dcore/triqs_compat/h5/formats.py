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

import re
from copy import deepcopy

class FormatInfo:
    """
    This class encapsulates essential information about for a particular h5 format.
    The information includes classname, modulename, documentaiton, the function to read.
    Further it provides information about the hdf5_format strings of subgroup keys,
    which is relevant for providing backward compatible reads.
    """
    def __init__(self, classname, modulename, doc, hdf5_format, read_fun) :
        self.classname, self.modulename, self.doc, self.read_fun = classname, modulename, doc, read_fun
        self.format_name = hdf5_format
        self.backward_compat = {} # key ->  hdf5_format

    def __str__(self) :
        return """
        Name of the class : %s
        Name of the module : %s
        Documentation : %s"""%(self.classname,self.modulename,self.doc)

# Dictionary containing the FormatInfo for all hdf5_format strings
_formats_dict= {}

_formats_backward_compat = [] # List of (regex, format, lambda)

def register_class (cls, doc = None, read_fun = None, hdf5_format = None):
    """
     For each class, register it with::

         from h5.formats import register_class
         register_class (GfImFreq, doc= doc_if_different_from cls._hdf5_format_doc_ )

    """
    hdf5_format = hdf5_format or (cls._hdf5_format_ if hasattr(cls,"_hdf5_format_") else cls.__name__)
    assert hdf5_format not in _formats_dict, "class %s is already registered"%hdf5_format
    doc = doc if doc else (cls._hdf5_format_doc_ if hasattr(cls,"_hdf5_format_doc_") else {})
    _formats_dict[hdf5_format] = FormatInfo(cls.__name__, cls.__module__, doc, hdf5_format, read_fun)


def register_backward_compatibility_method(regex, clsname, fun = lambda s: {}):
    """
    regex : the regular expression to match the hdf5_format (e.g. "Gf" for GfImfreq_x_whatever....)
    clsname : the class name that it corresponds to
    fun : a lambda taking hdf5_format and returning a dict
          field_name -> hdf5_format
          to read old data where not every subobjects have a hdf5_format.
    """
    _formats_backward_compat.append((regex, clsname, fun))


def get_format_info(hdf5_format):
    """
    Given an hdf5_format string, return the associated FormatInfo object.
    """
    # If present exactly, we return it
    if hdf5_format in _formats_dict:
        return _formats_dict[hdf5_format]

    # Enter compatibility mode.
    match_lst = [(regex,clsname,fun) for (regex,clsname,fun) in _formats_backward_compat if re.match(regex,hdf5_format)]
    if len(match_lst) == 0:
        raise KeyError("H5 Format %s is not registered and no backward compatibility found"%hdf5_format)
    if len(match_lst) > 1:
        raise KeyError("H5 Format %s : ambiguous backward compatibility layers : %s"%([regex for (regex,clsname,fun) in match_lst]))
    regex,clsname,fun = match_lst[0]

    # Make a copy of the associated Format object and supplement it with backward compatibility information
    fmt = deepcopy(_formats_dict[clsname])
    fmt.backward_compat = fun(hdf5_format)

    return fmt
