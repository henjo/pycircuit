import os
from pycircuit.post.cds.psf import PSFReader

def test_read_sweep_psf():
    filename = os.path.join(
        os.path.dirname(__file__),
        "psf/bwswp_acbw.sweep"
    )
    psf = PSFReader(filename) 
    psf.open()
    ## (closed: the reader sits in a reference cycle, and its file
    ## otherwise closes inside whichever test the collector runs in)
    psf.file.close()
