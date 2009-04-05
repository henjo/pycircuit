import logging

from pycircuit.post import InternalResultDict, Waveform

import numpy as np

## Result parser
class GnucapResult(InternalResultDict):
    def __init__(self, filename, xlabels=None, ylabel=None):
        super(GnucapResult, self).__init__()

        self.parse(filename, xlabels)

    def parse(self, filename, xlabels):
        if xlabels == None:
            xlabels = []
        else:
            xlabels = list(xlabels)

        ## For debug
        f = open(filename)
        logging.debug('gnucap result file: ' + f.read())

        f = open(filename)
        header=f.readline()[1:]
        names=header.split()
        data={}

        logging.debug('Result header: ' + header)

        for name in names:
            data[name]=[]

        for line in f.readlines():
            vals=map(sp2float, line.split())
            
            if len(vals) > len(names):
                names.insert(0, 'x')
                data['x'] = []
                xlabels.append('x')

            for name, val in zip(names,vals):
                data[name].append(val)

        f.close()                

        logging.debug('Result data items: ' + str(data))

        ## Store result in ResultDict
        xlist = [v for k, v in data.items() if k in xlabels]
 
        for k, v in data.items():
            if k not in xlabels:
                self[k] = Waveform(xlist, np.array(data[k]),
                                   xlabels = xlabels, ylabel=k)


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

