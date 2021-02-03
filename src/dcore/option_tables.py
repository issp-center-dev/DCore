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


import sys
import numpy
from .program_options import *

from io import StringIO

def readable_type_string(t):
    if t == int:
        return "Integer"
    elif t == float:
        return "Float"
    elif t == str:
        return "String"
    elif t == bool:
        return "Bool"
    else:
        raise RuntimeError("Unsupported type " + str(t))


def max_length(str_list):
    """
    Compute the size of the longest string in a given list of strings
    """
    return numpy.amax([len(s) for s in str_list])

def generate_description(p, section):
    """
    Generate descriptions of all options in a given section

    :param p: parser
    :param section: str
        section name
    :return description: str
        Description of all options
    """

    output = StringIO()

    # Check length of strings
    name_strings = p.get_predefined_options(section)
    type_strings = []
    default_value_strings = []
    description_strings = []
    for option in p.get_predefined_options(section):
        type_strings.append(readable_type_string(p.get_type(section, option)))
        default_value_strings.append(str(p.get_default_value(section, option)))
        description_strings.append(p.get_description(section, option))

    width = [max_length(sl)+4 for sl in [name_strings, type_strings, default_value_strings, description_strings]]

    def print_one_line(*str_list):
        for column in range(3):
            print(str_list[column].ljust(width[column], ' '), end=' ', file=output)
        print(str_list[3].ljust(width[3], ' '), file=output)

    print_one_line('='*width[0], '='*width[1], '='*width[2], '='*width[3])
    print_one_line('Name', 'Type', 'Default', 'Description')
    print_one_line('='*width[0], '='*width[1], '='*width[2], '='*width[3])
    for i in range(len(name_strings)):
        print_one_line(name_strings[i], type_strings[i], default_value_strings[i], description_strings[i])
    print_one_line('='*width[0], '='*width[1], '='*width[2], '='*width[3])

    return output.getvalue()

def generate_all_description():
    p = create_parser()
    desc = []
    for section in p.get_predefined_sections():
        desc.append("\n[{}]\n".format(section))
        desc.append(generate_description(p, section))
    return ''.join(desc)

if __name__ == '__main__':
    import sys
    p = create_parser()

    if len(sys.argv) != 2:
        raise RuntimeError("Invalid number of arguments!")
    prefix = sys.argv[1]
    print("Writing tables of runtime options into ", prefix)

    for section in p.get_predefined_sections():
        desc = generate_description(p, section)
        with open(prefix+'/'+section+'_desc.txt', 'w') as f:
            print(desc, file=f)
