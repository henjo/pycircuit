import os

from nose.tools import *

from pycircuit.post.cds import PSFResultSet
from pycircuit.post import Waveform
from pycircuit.post.testing import *

psfresultset_testvectors = {
    'dcop.raw':
        (('dcOp-dc', 'vout', 2.5), 
         ('dcOp-dc', 'vin' , 5.0)),
    'dcsweep.raw':
        (('dc1-dc', 'out', 
          Waveform(([  1.        ,   3.66666667,   6.33333333,  9.        ],), 
                    [  2.        ,   4.66666667,   7.33333333,  10.        ])),),
}

withlibpsf = 'withlibpsf'
withoutlibpsf = 'withoutlibpsf'

def test_psfresultset():
    for uselibpsf in False, True:
        for rawdir in psfresultset_testvectors.keys():
            if uselibpsf:
                yield psfresultset, withlibpsf, rawdir
            else:
                yield psfresultset, withoutlibpsf, rawdir

def psfresultset(libpsf, rawdir):
    if libpsf == withlibpsf:
        os.environ['USELIBPSF'] = '1'
    else:
        os.environ['USELIBPSF'] = '0'

    rs = PSFResultSet(rawdir)

    for resultname, signal, refvalue in psfresultset_testvectors[rawdir]:
        value = rs[resultname][signal]
        if isinstance(refvalue, Waveform):
            assert_waveform_almost_equal(value, refvalue)
        else:
            assert_almost_equal(value, refvalue)


    if libpsf == withlibpsf:
        ## Verify that libpsf was really used if libpsf was set
        assert bool(int(os.environ['USELIBPSF'])), 'libpsf could not be loaded'

