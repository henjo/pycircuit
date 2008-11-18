#!/usr/bin/python
import sys
import psf

if len(sys.argv) < 2:
    print "Usage: psftoasc PSFFILE"
else:
    filename = sys.argv[1]

    psfobj = psf.PSFFile()
    psfobj.deSerializeFile(open(filename))

    print psfobj.toPSFasc()
