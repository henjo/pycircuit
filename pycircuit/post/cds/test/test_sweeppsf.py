from pycircuit.post.cds.psf import PSFReader

def test_read_sweep_psf():
    filename = "psf/bwswp_acbw.sweep"
    psf = PSFReader(filename)
    psf.open()
