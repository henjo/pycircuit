import pycircuit.utilities.rst as rst
from nose.tools import *

def test_table():
    s = rst.table([['v1', 'i'], ['V', 'A'], [1,2e-3], [3,4e-3]], headerrows=2)
    assert_equal(s, """==== =======
v1   i      
V    A      
==== =======
   1   0.002
   3   0.004
==== =======""")
