# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import logging

from pycircuit.post import InternalResultDict, IVResultDict, Waveform

import numpy as np
import re

## Result parser
class NGspiceResult(InternalResultDict, IVResultDict):
    def __init__(self, s, xlabels=None, ylabel=None):
        super(NGspiceResult, self).__init__()

        self.parse(s, xlabels)

    def parse(self, s, xlabels, sweep_dimensions=1):
        lines = s.split('\n')

        header = lines[0][1:]
        names = re.split("\s+", header.rstrip())
        data = {}

        if xlabels == None:
            xlabels = names[:sweep_dimensions]
        else:
            xlabels = list(xlabels)

        logging.debug('Result header: ' + header)

        for name in names:
            data[name] = []

        for line in lines[1:]:
            vals = map(sp2float, line.split())

            for name, val in zip(names,vals):
                data[name].append(val)

        logging.debug('Result data items: ' + str(data))

        ## Store result in ResultDict
        xlist = [v for k, v in data.items() if k in xlabels]
 
        for k, v in data.items():
            if k not in xlabels:
                self[k] = Waveform(xlist, np.array(data[k]),
                                   xlabels = xlabels, ylabel=k)

    def v(self, plus, minus=None):
        """Returns the voltage between the plus and minus node or potential of plus node"""
        vplus = self['v(' + str(plus) + ')'] 
        vminus = 0
        if minus != None:
            vminus = self['v(' + str(minus) + ')']
        return vplus - vminus

    def i(self, terminal):
        """Returns the current flowing into the given terminal or branch object"""
        raise NotImplementedError()


def sp2float(text):    
        """
        converts spice units to float
        """
        text=text.strip()
        tbl={"m":"e-3","M":"e-3","Meg":"e6","K":"e3","G":"e9","u":"e-6",
             "n":"e-9","p":"e-12","f":"e-15"}
        
        keys=tbl.keys()
        for key in keys:
            if text.endswith(key):
                text=text.replace(key,tbl[key])
                
        return float(text)

