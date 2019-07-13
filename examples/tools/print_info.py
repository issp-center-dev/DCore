from __future__ import print_function

import os
import ConfigParser
from collections import OrderedDict


file_desc = 'examples_desc.txt'
info_file = 'info'
info_section = 'info'
keys = ['author', 'date', 'version', 'description']


def get_info(filename):
    cfg = ConfigParser.SafeConfigParser()
    cfg.read(filename)

    def get(key):
        # val = cfg.get('info', key)
        try:
            val = cfg.get(info_section, key)
        except:
            val = ""
        return val

    # info = OrderedDict()
    # for key in keys:
    #     info[key] = get(key)

    info = OrderedDict({key: get(key) for key in keys})

    return info


def get_all_info():
    infos = {}

    for dir in os.listdir('.'):
        file = os.path.join(dir, info_file)
        # print(file)
        if not os.path.exists(file):
            continue

        # if not os.path.isdir(dir):
        #     continue

        # os.chdir(dir)
        print(dir)

        try:
            info = get_info(file)
            print(info)
            infos[dir] = info
        except:
            continue

    return infos


def gen_description(list_str):
    desc = "|"
    for entry in list_str:
        desc += " " + entry + " |"
    return desc


if __name__ == '__main__':

    infos = get_all_info()
    print(infos)

    with open(file_desc, 'w') as f:
        desc = gen_description(['directory',] + keys)
        print(desc, file=f)
        desc = gen_description(['-----',] * (len(keys)+1))
        print(desc, file=f)

        for dir, info in infos.items():
            desc = gen_description([dir,] + [info[key] for key in keys])
            print(desc, file=f)