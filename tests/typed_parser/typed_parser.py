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


from dcore.typed_parser import *
import pytest

@pytest.fixture(scope="function")
def test_read_file(request):
    p = TypedParser(['sectionA', 'sectionB'])

    p.add_option("sectionA", "a", int, -1000, "a in sectionA", OptionStatus.RETIRED)
    p.allow_undefined_options("sectionB")

    # SectionC must be ignored.
    p.add_option("sectionC", "c", int, -1000, "c in sectionC")

    params = p.as_dict()
    assert params["sectionA"]["a"] == -1000

    p.read(request.fspath.dirname+"/parser.in")

    params = p.as_dict()
    assert params["sectionA"]["a"] == 1
    assert params["sectionB"]["b"] == 'B'

    assert "sectionC" not in params


# Detect undefined option?
@pytest.fixture(scope="function")
def test_detect_undefined_option(request):
    p2 = TypedParser(['sectionAA'])
    with open(request.fspath.dirname+'/parser_test_2.in', 'w') as f:
        print("[sectionAA]", file=f)
        print("aa = 2", file=f)

    thrown = False
    try:
        p2.read(request.fspath.dirname+"/parser_test_2.in")
    except RuntimeError:
        thrown = True
    assert thrown


def test_float_tuple():
    t = FloatTuple('(1, 2, 3,)')
    assert str(t) == '(1.0 , 2.0 , 3.0)'

    t2 = FloatTuple(t)
    assert str(t2) == str(t)

def test_int_tuple():
    t = IntTuple('(1, 2, 3,)')
    assert str(t) == '(1 , 2 , 3)'

    t2 = IntTuple((10, 10, 10))
