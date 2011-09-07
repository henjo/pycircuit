# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import re
import os
import doctest
import pycircuit

from pycircuit.post.cds.psf import *
## Unit tests for psf module
    
class toPSFASCTests(unittest.TestCase):
    def setUp(self):
        self.dir = os.path.dirname(__file__)
    def testString(self):
        s=String(value="test1234")
        self.assertEqual(s.toPSFasc(), '"test1234"')
    def testFloat64(self):
        v=Float64(value=9.731734322e-31)
        self.assertEqual(v.toPSFasc(), '9.73173e-31')
    def testComplexFloat64(self):
        v=ComplexFloat64(value=complex(0.000291959, -2.81639e-09))
        self.assertEqual(v.toPSFasc(), '(0.000291959 -2.81639e-09)')
    def testProperty(self):
        p=PropertyFloat64(name="tolerance.relative", value=1e-3)
        self.assertEqual(p.toPSFasc(prec=9), '"tolerance.relative" 0.00100000')
        p=PropertyString(name="date", value="2:30:26 PM, Thur Sep 27, 2007")
        self.assertEqual(p.toPSFasc(),'"date" "2:30:26 PM, Thur Sep 27, 2007"')

for section in [("HEADER", "header"), ("TYPE", "types"), ("SWEEP", "sweeps"), ("TRACE", "traces"), ("VALUE", "values")]:
    def method(self):
        psfasc = PSFASCSplitter(self.dir + "/psfasc/dc.dc.asc")
        expected=psfasc.sections[section[0]].strip()
        psf = PSFReader(self.dir + "/psf/dc.dc")
        psf.open()
        self.assertEqual(getattr(psf, section[1]).toPSFasc(), expected)
    setattr(toPSFASCTests, "test%sSection" % section[1], method)

class PSFTests(unittest.TestCase):
    def setUp(self):
        self.dir = os.path.dirname(__file__)
    def testPSFasc(self):
        psf=PSFReader(self.dir + "/psf/dc.dc")
        psf.open()
        
        psfascfile=open(self.dir + "/psfasc/dc.dc.asc")
        
        for actual, expected in zip(psf.toPSFasc().split("\n"), psfascfile.readlines()):
            self.assertEqual(actual, expected.strip())
        
class PSFASCSplitter:
    def __init__(self, filename):
        self.sections={}
        f=open(filename)
        buffer=""
        section=None
        for line in f:
            if line.strip() in ("HEADER", "TYPE", "SWEEP", "TRACE", "VALUE", "END"):
                if section:
                    self.sections[section] = buffer
                section=line.strip()
                buffer=line
            else:
                buffer+=line
        if section:
            self.sections[section] = buffer

def psfAscAdjust(str):
    """There are minor differences how the psf util handles floats, -0.0000 is
    shown as 0.0000 in Python. This function corrects for that.
    >>> psfAscAdjust("-0.00000")
    '0.00000'
    >>> psfAscAdjust("bla bla bla -0.00 0.00 -0.001 1.00000e-05")
    'bla bla bla 0.00 0.00 -0.001 1.00000e-05'
    """
    str = re.sub("-(0\.0*$)", '\\1', str)
    return re.sub("-(0\.0*\D)", '\\1', str)
    
# def load_tests(loader, tests, ignore):
#     tests.addTests(doctest.DocTestSuite(pycircuit.post.cds.test.test_psfasc))
#     #tests.addTests(doctest.DocTestSuite(pycircuit.post.cds.psf))
#     return tests


