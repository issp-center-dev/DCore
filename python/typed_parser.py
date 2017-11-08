import os, copy

try:
    import configparser
except ImportError:
    import ConfigParser as configparser


class TypedParser(object):
    """
    Parser of an init file. One is able to define options and sections (data types and default values).
    """
    def __init__(self):
        """
        Initializer
        """
        self.__config_parser = configparser.ConfigParser()
        self.__config_parser.optionxform = str

        self.__definitions= {}
        self.__results = {}

        self.__read = False

        # Names of sections in which we accept options that are not predefined by add_option().
        self.__allow_undefined_options = []

    def add_option(self, section, option, dtype, default, string):
        """
        :param section: section name
        :param option: option name
        :param dtype: data type
        :param default: default value
        :param string: short description
        """
        if self.__read:
            raise RuntimeError("Do not add option after an input file has been read!")

        if not section in self.__definitions:
            self.__definitions[section] = {}

        if not section in self.__results:
            self.__results[section] = {}

        self.__definitions[section][option] = (dtype, string, dtype(default))
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
            if not sect in self.__results:
                self.__results[sect] = {}

            for opt in self.__config_parser.options(sect):
                value = self.__config_parser.get(sect, opt)
                if sect in self.__definitions and opt in self.__definitions[sect]:
                    # if an option is pre-defined.
                    self.__results[sect][opt] = self.__definitions[sect][opt][0](value)
                else:
                    # if an option is not pre-defined, use the value in the input file
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

    def as_dict(self):
        """
        Convert all options and their values into a dict object

        :return: dict object
        """
        return copy.deepcopy(self.__results)

    def print_options(self):
        """
        Print a list of all options and sections

        :return:
        """

        # To be implemented
        pass