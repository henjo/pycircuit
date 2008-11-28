# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

# Read different kind of PSF files and compare with output from the psf command

from psf import PSFReader
import os, os.path, sys
import unittest
from psfasctests import psfAscAdjust

class readAndExportTest(unittest.TestCase):
    def __init__(self, psffile, psfascfile):
        unittest.TestCase.__init__(self)
        self.psffile = psffile
        self.psfascfile = psfascfile

    def runTest(self):
        f = open(self.psffile)
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
        

def makeSuite():
    return suite

if __name__ == '__main__':
    suite = unittest.TestSuite()

    # Find PSF files with .asc counterparts
    dirname=os.path.dirname(sys.argv[0])
    if dirname == "":
        dirname = "."

    psfdir = os.path.join(dirname, "psf")
    psfascdir = os.path.join(dirname, "psfasc")
        
    psffiles = [f for f in os.listdir(os.path.join(dirname, "psf")) if f[0] != "."]

    for file in psffiles:
        suite.addTest(readAndExportTest(os.path.join(psfdir, file), os.path.join(psfascdir, file+".asc")))

    unittest.main(defaultTest="makeSuite")
#    makeSuite()

            
        
        
