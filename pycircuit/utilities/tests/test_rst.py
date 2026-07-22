import pycircuit.utilities.rst as rst

def test_table():
    s = rst.table([['v1', 'i'], ['V', 'A'], [1,2e-3], [3,4e-3]], headerrows=2)
    assert s == """==== =======
v1   i      
V    A      
==== =======
   1   0.002
   3   0.004
==== ======="""
