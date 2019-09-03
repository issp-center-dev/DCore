from __future__ import print_function

import os
from collections import OrderedDict

try:
    # python 3
    from configparser import ConfigParser
except ImportError:
    # pytnon 2
    from ConfigParser import SafeConfigParser as ConfigParser

file_desc = 'examples_info.md'
info_file = 'info'
info_section = 'info'
keys = ['author', 'date', 'version', 'description']


def get_info(filename):
    cfg = ConfigParser()
    cfg.read(filename)

    def get(key):
        try:
            val = cfg.get(info_section, key)
        except:
            val = ""
        return val

    _info = OrderedDict({key: get(key) for key in keys})
    return _info


def get_all_info():
    _infos = OrderedDict()

    for _dir in sorted(os.listdir('.')):
        file_path = os.path.join(_dir, info_file)
        if not os.path.exists(file_path):
            continue

        print(_dir)

        try:
            info = get_info(file_path)
            # print(info)
            _infos[_dir] = info
        except:
            continue

    return _infos


def gen_description(list_str):
    _desc = "|"
    for entry in list_str:
        _desc += " " + entry + " |"
    return _desc


if __name__ == '__main__':

    infos = get_all_info()
    # print(infos)

    with open(file_desc, 'w') as f:
        desc = gen_description(['directory',] + keys)
        print(desc, file=f)
        desc = gen_description(['-----',] * (len(keys)+1))
        print(desc, file=f)

        for direc, info in infos.items():
            desc = gen_description([direc, ] + [info[key] for key in keys])
            print(desc, file=f)
