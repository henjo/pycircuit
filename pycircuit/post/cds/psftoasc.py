# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

#!/usr/bin/python
import sys
from pycircuit.post.cds import psf

if len(sys.argv) < 2:
    print("Usage: psftoasc PSFFILE")
else:
    filename = sys.argv[1]

    psfobj = psf.PSFData.fromFile(open(filename))

    print(psfobj.toPSFasc())
