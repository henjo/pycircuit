# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

# Read different kind of PSF files and compare with output from the psf command

from pycircuit.post.cds.psf import PSFReader
import os, os.path, sys
import unittest
from pycircuit.post.cds.test.test_psfasc import psfAscAdjust

class readAndExportTest():
    psffile = None

    def runTest(self):
        if not self.psffile:
            return
        psf = PSFReader(self.psffile)
        psf.open()

        psfascfile = open(self.psfascfile)

        for actual, expected in zip(psf.toPSFasc().split("\n"), psfascfile.readlines()):
            expected = psfAscAdjust(expected)
            actual = psfAscAdjust(actual)
##            print actual, "==", expected.strip()
            self.assertEqual(actual, expected.strip())

    def __str__(self):
        return "Readin PSF (%s), export to PSFASC, comparing with psf tool output."%self.psffile
        

def makeSuite(suite):
    # Find PSF files with .asc counterparts
    dirname=os.path.dirname(__file__)
    if dirname == "":
        dirname = "."

    psfdir = os.path.join(dirname, "psf")
    psfascdir = os.path.join(dirname, "psfasc")
        
    psffiles = [f for f in os.listdir(os.path.join(dirname, "psf")) if f[0] != "."]

    for filename in psffiles:
        class TestPsf(unittest.TestCase, readAndExportTest):
            pass
        t = TestPsf()
        t.psffile = os.path.join(psfdir, filename)
        t.psfascfile = os.path.join(psfascdir, filename+".asc")
        suite.addTest(t)
    return suite

#def load_tests(loader, tests, ignore):
#    return makeSuite(tests)
