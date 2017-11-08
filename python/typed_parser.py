import os, copy

try:
    import configparser
except ImportError:
    import ConfigParser as configparser


class TypedParser(object):
    def __init__(self):
        # Make case sensitive parser
        self.__config_parser = configparser.ConfigParser()
        self.__config_parser.optionxform = str

        self.__definitions= {}
        self.__results = {}

        self.__read = False

        # Names of sections in which we accept options that are not predefined by add_option().
        self.__allow_undefined_options = []

    def add_option(self, section, option, dtype, default, string):
        """
        section, option, type, default value, string
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
        if type(section) == str and not (section in self.__allow_undefined_options):
            self.__allow_undefined_options.append(section)
        else:
            raise ValueError("section must be a str!")

    def read(self, in_file):
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
        return self.__results[sect][opt]

    def as_dict(self):
        return copy.deepcopy(self.__results)