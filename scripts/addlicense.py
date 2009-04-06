#!/bin/env python

import os
import fnmatch
import logging

def add_license_header(dir, licenseheader, excludedfiles, oldheaders=[]):
    def excluded(filename):
        return filename in excludedfiles

    def check_license(filename):
        logging.debug('Checking file: ' + filename)

        f = open(filename)
        
        lines = f.readlines()

        ## Remove old headers
        for header in oldheaders:
            if lines[:len(header)] == header:
                lines = lines[len(header):]            
                logging.info('found old header in ' + filename)

        if not lines[:len(licenseheader)] == licenseheader:
            logging.info('adding license header to ' + filename)
            f = open(filename, 'w')
            f.writelines(licenseheader + lines )

    for dirpath, dirnames, filenames in os.walk(dir, topdown=True):
        ## Remove dot files
        for dirname in reversed(dirnames):
            if dirname[0] == '.':
                dirnames.remove(dirname)

        for filename in filenames:
            if filename.endswith('.py') and not filename.startswith('__') \
                    and not excluded(filename):
                check_license(os.path.join(dirpath, filename))

oldheaders = [
    ['# Copyright (c) 2008 Pycircuit Development Team\n',
     '# See LICENSE for details.\n',
     '\n']
    ]

dir = '../pycircuit'

licenseheader = list(open('licenseheader.txt', 'r').readlines())

logging.basicConfig(level=logging.INFO)

add_license_header(dir, licenseheader, ['DESolver.py', 'test_DE'], 
                   oldheaders=oldheaders)
