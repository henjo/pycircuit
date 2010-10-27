# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import re

from pycircuit.post.cds.psf import *
## Unit tests for psf module
    
class toPSFASCTests(unittest.TestCase):
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
        p=PropertyFloat64(None, name="tolerance.relative", value=1e-3)
        self.assertEqual(p.toPSFasc(prec=9), '"tolerance.relative" 0.00100000')
        p=PropertyString(None, name="date", value="2:30:26 PM, Thur Sep 27, 2007")
        self.assertEqual(p.toPSFasc(),'"date" "2:30:26 PM, Thur Sep 27, 2007"')
    def testHeaderSection(self):
        psfasc = PSFASCSplitter("examples/dc.dc.asc")
        expected=psfasc.sections["HEADER"].strip()
        psf = PSFReader("examples/dc.dc")
        psf.open()
        self.assertEqual(psf.header.toPSFasc(), expected)

    def testTypesSection(self):
        psfasc = PSFASCSplitter("examples/dc.dc.asc")
        expected=psfasc.sections["TYPE"].strip()
        f=open("examples/dc.dc")
        psf = PSFReader("examples/dc.dc")
        psf.open()
        self.assertEqual(psf.types.toPSFasc(), expected)

    def testSweepSection(self):
        psfasc = PSFASCSplitter("examples/dc.dc.asc")
        expected=psfasc.sections["SWEEP"].strip()
        f=open("examples/dc.dc")
        psf = PSFReader("examples/dc.dc")
        psf.open()
        self.assertEqual(psf.sweeps.toPSFasc(), expected)

    def testTraceSection(self):
        psfasc = PSFASCSplitter("examples/dc.dc.asc")
        expected=psfasc.sections["TRACE"].strip()
        f=open("examples/dc.dc")
        psf = PSFReader("examples/dc.dc")
        psf.open()
        self.assertEqual(psf.traces.toPSFasc(), expected)

    def testValueSection(self):
        psfasc = PSFASCSplitter("examples/dc.dc.asc")
        expected=psfasc.sections["VALUE"].strip()
        f=open("examples/dc.dc")
        psf = PSFReader("examples/dc.dc")
        psf.open()
        self.assertEqual(psf.values.toPSFasc(), expected)


class PSFTests(unittest.TestCase):
    def testPSFasc(self):
        filename = "examples/dc.dc"

        psf=PSFReader(filename)
        psf.open()
        
        psfascfile=open(filename+".asc")
        
        for actual, expected in zip(psf.toPSFasc().split("\n"), psfascfile.readlines()):
#            print actual, "==", expected.strip()
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
    

if __name__ == '__main__':
    import doctest
    doctest.testmod()

#    unittest.main()

