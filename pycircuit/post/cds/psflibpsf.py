# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import unittest

import libpsf

class PSFReader(object):
    def __init__(self, filename=None, asc=False):
        if asc:
            raise NotImplemented("Ascii PSF not implemented in psf library")
        self.filename = filename
        self.ds = None

    def open(self):
        """Open a PSF file and read its headers.

        Example:
        Trying to open a valid psf file
        >>> psf=PSFReader('./test/psf/srcSweep')
        >>> psf.open()
        """

        self.ds = libpsf.PSFDataSet(self.filename)

        self.hprops = self.ds.get_header_properties()

    def getNSweepPoints(self):
        """Returns number of sweeps. 0 if not swept.

        >>> psf=PSFReader('./test/psf/srcSweep')
        >>> psf.open()
        >>> psf.getNSweepPoints()
        4
        """
        if self.ds == None:
            ValueError("Please open the PSF file first")
        return self.ds.get_sweep_npoints()

    def getNSweeps(self):
        """Returns the number of nested sweeps

        >>> psf=PSFReader('./test/psf/srcSweep')
        >>> psf.open()
        >>> psf.getNSweeps()
        1
        """
        return 1

    def __len__(self):
        return len(self.ds.get_signal_names())

    def getValueNames(self):
        """Returns a tuple of the names of the traces

        >>> psf=PSFReader('./test/psf/srcSweep')
        >>> psf.getValueNames()
        >>> psf.open()
        >>> psf.getValueNames()
        ('VOUT', 'VIN', 'R0')

        >>> psf=PSFReader('./test/resultdirs/simple/opBegin')
        >>> psf.open()
        >>> psf.getValueNames()
        ('R0', 'V1', 'V0', 'E0', 'VIN', 'NET9', 'VOUT')

        """
        if self.ds == None:
            raise ValueError("Please open the PSF file first")
        return tuple(self.ds.get_signal_names())
    
    def getSweepParamValues(self, dim=0):
        """Returns a numpy.array of sweep parameter values for sweep dimension dim.

        >>> psf=PSFReader('./test/psf/srcSweep')
        >>> psf.open()
        >>> psf.getSweepParamValues(0)
        array([ 1.,  2.,  3.,  4.])

        windowed result
        >>> psf=PSFReader('./test/psf/timeSweep')
        >>> psf.open()
        >>> psf.getSweepParamValues(0)[:3]
        array([  0.00000000e+00,   2.00000000e-11,   5.33333333e-11])

        """
        if self.ds == None:
            raise ValueError("Please open the PSF file first")
        return self.ds.get_sweep_values()

    def getValuePropertiesByName(self, name):
        """Returns the properties associated with value
        
        >>> psf=PSFReader('./test/psf/opBegin')
        >>> psf.open()
        >>> psf.getValuePropertiesByName("XIRXRFMIXTRIM0.XM1PDAC1.XMN.MAIN")["Region"]
        'subthreshold'

        """
        if self.ds == None:
            raise ValueError("Please open the PSF file first")
        return self.ds.get_signal_properties(name)

    def getValuesByName(self, name):
        """Returns a numpy.array of trace values for swept results and a scalar for non swept.

        Example:
        swept psf file
        >>> psf=PSFReader('./test/psf/srcSweep')
        >>> psf.open()
        >>> psf.getValuesByName("VOUT")
        array([-6., -4., -2.,  0.])
        >>> psf.getValuesByName("VIN")
        array([ 1.,  2.,  3.,  4.])

        swept psf with complex numbers
        >>> psf=PSFReader('./test/psf/frequencySweep')
        >>> psf.open()
        >>> res = psf.getValuesByName("ANT_CM")
        >>> len(res)
        123
        >>> res[:3]
        array([ 0.6+0.j,  0.0+0.j,  0.0+0.j])
        
        swept windowed psf file 
        >>> psf=PSFReader('./test/psf/timeSweep')
        >>> psf.open()
        >>> psf.getValuesByName("INP")[0:3]
        array([ 0.6       ,  0.62486899,  0.66211478])

        non-swept psf file
        >>> psf=PSFReader('./test/psf/dcOpInfo.info')
        >>> psf.open()
        >>> psf.getValuesByName("IREG21U_0.MP5.b1")['betadc']
        4.7957014499434756

        swept psf file withouth groups
        >>> psf=PSFReader('./test/resultdirs/parsweep2/C=1e-12,R=1e-12/psf/ac.ac')
        >>> psf.open()
        >>> psf.getValuesByName("net3")
        array([ 0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,
                0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,
                0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,
                0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,
                0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,
                0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,
                0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,
                0.+0.j,  0.+0.j])

        """
        if self.ds == None:
            raise ValueError("Please open the PSF file first")
        return self.ds.get_signal(name)
        
    def nTraces(self):
        """Returns number of traces

        >>> psf=PSFReader('./test/psf/srcSweep')
        >>> psf.open()
        >>> psf.nTraces()
        3
        """
        if self.ds == None:
            raise ValueError("Please open the PSF file first")
        return self.hprops['PSF traces']

    def info(self):
        s="Number of sweeps: %d\n"%self.getNSweeps()
        if self.getNSweeps() > 0:
            s+="Number of sweep points: %d\n"%self.getNSweepPoints()
        s+="Number of traces: %d"%self.nTraces()
        return s

if __name__ == "__main__":
    import doctest
    doctest.testmod()


