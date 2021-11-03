import os
from pycircuit.post.cds.psf import PSFReader

def test_read_sweep_psf():
    filename = os.path.join(
        os.path.dirname(__file__),
        "psf/bwswp_acbw.sweep"
    )
    psf = PSFReader(filename) 
    psf.open()
