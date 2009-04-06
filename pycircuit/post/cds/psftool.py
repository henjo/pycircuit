# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

#!/usr/bin/python
import getopt, sys
import psf

def usage():
    print "Usage: psftoasc [-OPTION] PSFFILE"
    print "Options:"
    print "  -h, --help\tPrint this help"
    print "  -i, --info\tPrint information about the PSF"
    print "  -a, --asc\tExport to PSFasc (default)"

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ihav", ["info" "help", "output="])
    except getopt.GetoptError:
        # print help information and exit:
        sys.exit(2)
    output = None
    verbose = False

    cmd = "psfasc"
    
    for o, a in opts:
        if o == "-v":
            verbose = True
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-i", "--info"):
            cmd = "info"
        if o in ("-a", "--asc"):
            cmd = "info"
        

    if len(args) < 1:
        usage()
        sys.exit(2)
        
    filename = args.pop()

    psfobj = psf.PSFReader(filename)
    psfobj.verbose = verbose
    psfobj.open()

    if cmd == "psfasc":
        print psfobj.toPSFasc()
    elif cmd == "info":
        print psfobj.info()
        if verbose:
            psfobj.printme()
    


if __name__ == "__main__":
    main()

import sys

