# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Spectre PSF results through polars-psf.

``PSFResultSet(resultdir)`` keeps the ResultDict interface of earlier pycircuit versions::

    rs = PSFResultSet('sim.raw')
    rs['dc1-dc']['out']                        # swept signal: a polars_waveform.Waveform
    rs['dcOp-dc']['vout']                      # operating point: a number
    rs['designParamVals-info']['top-level']    # info: a dict

Result keys are ``<name>-<type>`` as in the logFile; plain names (``rs['dc1']``) work too.
Parametric sweeps come back as families (one curve per parameter combination). For more, such
as lazy tables, filtering and schematic names, use :func:`polars_psf.open` directly
(``rs.dataset``).

Virtuoso/SKILL integration: skillbridge (https://github.com/unihd-cag/skillbridge) for a running
Virtuoso; simdeck.virtuoso (https://github.com/henjo/simdeck) for a headless session, SKILL files
and parsing SKILL values.
"""

import polars_psf

from pycircuit.post.result import ResultDict

__all__ = ['PSFResultSet', 'PSFResult']


class PSFResultSet(ResultDict):
    """The results of a Spectre result directory (or a single PSF file)."""

    def __init__(self, resultdir):
        self.dataset = polars_psf.open(resultdir)
        self._rows = self.dataset.results.select('name', 'type').rows()

    def keys(self):
        return [f'{name}-{type_}' for name, type_ in self._rows]

    def __len__(self):
        return len(self._rows)

    def __getitem__(self, key):
        for name, type_ in self._rows:
            if key in (name, f'{name}-{type_}'):
                return PSFResult(self.dataset.result(name))
        raise KeyError(f'no result {key!r} in {self.keys()}')


class PSFResult(ResultDict):
    """One result: signal name -> waveform (swept) or value (not swept)."""

    def __init__(self, result):
        self.result = result

    def keys(self):
        return self.result.names

    def __len__(self):
        return len(self.keys())

    def __getitem__(self, name):
        if name not in self.keys():
            raise KeyError(f'no signal {name!r} in {self.result.name}')
        if self.result.sweep_name is None:
            return self.result.value(name)
        return self.result.v(name)

    def v(self, plus, minus=None):
        w = self.result.v(plus)
        return w if minus is None else w - self.result.v(minus)

    def i(self, terminal):
        return self.result.i(terminal)
