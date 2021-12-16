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
import os
import copy
from enum import Enum
from warnings import warn
from collections import OrderedDict

try:
    import configparser
except ImportError:
    import configparser as configparser

class TypedTuple(object):
    def __init__(self, data, elem_type):
        if isinstance(data, str):
            chars_ignore = ['(', ')', '[', ']', ' ']
            for c in chars_ignore:
                data = data.replace(c, '')
            data = data.strip(',')

            self._data = tuple([elem_type(x) for x in data.split(',')])
        elif isinstance(data, tuple):
            for x in data:
                if not isinstance(x, elem_type):
                    raise RuntimeError("Invalid element in tuple.")
            self._data = data
        else:
            raise RuntimeError("Invalid data for Tuple")

    def __repr__(self):
        """
        Return a string representation like (1, 2, 3).

        :return: str
        """
        return '(' + ' , '.join([str(x) for x in self._data]) + ')'

    def to_tuple(self):
        assert isinstance(self._data, tuple)
        return self._data

class IntTuple(TypedTuple):
    def __init__(self, data):
        if isinstance(data, IntTuple):
            self._data = data._data
            assert isinstance(self._data, tuple)
        else:
            super(IntTuple, self).__init__(data, int)

class FloatTuple(TypedTuple):
    def __init__(self, data):
        if isinstance(data, FloatTuple):
            self._data = data._data
            assert isinstance(self._data, tuple)
        else:
            super(FloatTuple, self).__init__(data, float)


class OptionStatus(Enum):
    VALID = 0
    DEPRECATED = 1
    RETIRED = 2

def cast(value_type, string):
    if value_type == bool:
        if string in ['true', 'True']:
            return True
        elif string in ['false', 'False']:
            return False
        else:
            raise ValueError("Cannot cast string "+string+" to bool.")
    else:
        return value_type(string)


class TypedParser(object):
    """
    Parser of an init file. One is able to define options and sections (data types and default values).
    """
    def __init__(self, sections_to_be_used=[]):
        """
        Initializer
        :rtype:
        """
        self.__config_parser = configparser.ConfigParser()
        self.__config_parser.optionxform = str

        self.__definitions = OrderedDict()
        self.__results = OrderedDict()
        self.__section_to_be_used = sections_to_be_used

        self.__read = False

        # Names of sections in which we accept options that are not predefined by add_option().
        self.__allow_undefined_options = []

    def add_option(self, section, option, dtype, default, string, status = OptionStatus.VALID):
        """
        :param section: section name
        :param option: option name
        :param dtype: data type
        :param default: default value
        :param string: short description
        :param status: VALID, DEPRECATED, or RETIRED
        """

        if not section in self.__section_to_be_used:
            return

        if self.__read:
            raise RuntimeError("Do not add option after an input file has been read!")

        if section in self.__definitions and option in self.__definitions[section]:
            raise RuntimeError("Redefinition of an option is not allowed!")


        if section not in self.__definitions:
            self.__definitions[section] = OrderedDict()

        if section not in self.__results:
            self.__results[section] = OrderedDict()

        self.__definitions[section][option] = {'dtype' : dtype,
                                               'description' : string,
                                               'default' : dtype(default),
                                               'status' : status}
        self.__results[section][option] = dtype(default)

    def allow_undefined_options(self, section):
        """
        Allow undefined options for the given section.
        Otherwise, it throws when an undefined option is encountered during reading an ini file.

        :param section:  section name
        :return:
        """
        if type(section) == str and not (section in self.__allow_undefined_options):
            self.__allow_undefined_options.append(section)
        else:
            raise ValueError("section must be a str!")

    def read(self, in_file):
        """
        Read an init file. This function must not be called more than once.

        :param in_file:
        :return:
        """
        if self.__read:
            raise RuntimeError("An input file has been already read!")

        self.__read = True

        if not os.path.exists(in_file):
            raise RuntimeError("Not found "+in_file)
        self.__config_parser.read(in_file)

        for sect in self.__config_parser.sections():
            if sect not in self.__section_to_be_used:
                continue

            if sect not in self.__results:
                self.__results[sect] = {}

            for opt in self.__config_parser.options(sect):
                value = self.__config_parser.get(sect, opt).strip('\'').strip('"')
                if sect in self.__definitions and opt in self.__definitions[sect]:
                    # if an option is pre-defined.
                    if self.__definitions[sect][opt]['status'] == OptionStatus.DEPRECATED:
                        warn("Parameter {0} [{1}] is deprecated.".format(opt, sect))
                    if self.__definitions[sect][opt]['status'] == OptionStatus.RETIRED:
                        warn("Parameter {0} [{1}] is not used anymore.".format(opt, sect))
                    self.__results[sect][opt] = cast(self.__definitions[sect][opt]['dtype'], value)
                else:
                    # if an option is not pre-defined, use the value given in the input file
                    if sect in self.__allow_undefined_options:
                        self.__results[sect][opt] = value
                    else:
                        raise RuntimeError("Undefined option " + opt + " is not allowed in section " + sect + "!")

    def get(self, sect, opt):
        """
        Get the value of the given option in the given section

        :param sect: section name
        :param opt:  option name
        :return: value
        """
        return self.__results[sect][opt]

    def get_type(self, sect, opt):
        """
        Get the type of a given option
        """
        return self.__definitions[sect][opt]['dtype']

    def get_description(self, sect, opt):
        """
        Get the description of a given option
        """
        return self.__definitions[sect][opt]['description']

    def get_default_value(self, sect, opt):
        """
        Get the default value of a given option
        """
        return self.__definitions[sect][opt]['default']

    def get_predefined_sections(self):
        """
        Get a list of sections predefined by add_option()
        """
        return list(self.__definitions.keys())

    def get_predefined_options(self, section):
        """
        Get a list of options predefined by add_option() in a given section
        """
        return list(self.__definitions[section].keys())

    def as_dict(self):
        """
        Convert all options and their values into a dict object

        :return: dict object
        """

        # convert OrderedDict to ordinary dict recursively
        def convert_ordered_dict_to_dict(obj):
            if isinstance(obj, OrderedDict):
                r = {}
                for key, val in list(obj.items()):
                    r[key] = convert_ordered_dict_to_dict(val)
                return r
            else:
                return copy.deepcopy(obj)

        return convert_ordered_dict_to_dict(self.__results)


    def print_options(self):
        """
        Print a list of all options and sections

        :return: None
        """

        print("\n  @ Defined options")
        for section_name, section_data in list(self.__definitions.items()):
            print("")
            print(("   ["+section_name+"] block"))
            for option_name, option_data in list(section_data.items()):
                print(("     option =  {0} : type = {1} : description = \"{2}\" : default value = {3}".format(
                    option_name, option_data[0].__name__, option_data[1], option_data[2])))
