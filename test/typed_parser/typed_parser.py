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
from __future__ import print_function

from pytriqs.applications.dcore.typed_parser import TypedParser


def read_file():
    p = TypedParser()
    p.add_option("sectionA", "a", int, -1000, "a in sectionA")
    p.allow_undefined_options("sectionB")

    params = p.as_dict()
    assert params["sectionA"]["a"] == -1000

    p.read("parser.in")

    params = p.as_dict()
    assert params["sectionA"]["a"] == 1
    assert params["sectionB"]["b"] == 'B'


# Detect undefined option?
def detect_undefined_option():
    p2 = TypedParser()
    with open('parser_test_2.in', 'w') as f:
        print("[sectionAA]", file=f)
        print("aa = 2", file=f)

    thrown = False
    try:
        p2.read("parser_test_2.in")
    except RuntimeError:
        thrown = True
    assert thrown


read_file()
detect_undefined_option()
