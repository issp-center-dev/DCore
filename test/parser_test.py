from __future__ import print_function

from pytriqs.applications.pydmft.typed_parser import TypedParser

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
    except:
        thrown = True
    assert thrown

read_file()
detect_undefined_option()
