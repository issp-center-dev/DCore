from __future__ import print_function

from pytriqs.applications.pydmft.typed_parser import TypedParser


p = TypedParser()
p.add_option("sectionA", "a", int, -1000, "a in sectionA")

params = p.as_dict()
assert params["sectionA"]["a"] == -1000

p.read("parser.in")

params = p.as_dict()
assert params["sectionA"]["a"] == 1
assert params["sectionB"]["b"] == 'B'
