# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import pycircuit.sim

from pycircuit.sim.gnucap.circuit import Circuit
from pycircuit.sim.gnucap.result import GnucapResult
from pycircuit.sim.gnucap.session import GnucapSessionPexpect, GnucapSessionDirect

import logging
import tempfile
import numpy as np
import types

class Simulation(pycircuit.sim.Simulation):
    def __init__(self, circuit, direct = False, executable = None):
        ## If circuit is a string, assume it's a netlist and create a circuit
        if type(circuit) == types.StringType:
            circuit = Circuit(circuit)
            
        super(Simulation, self).__init__(circuit)

        if direct:
            self.session = GnucapSessionDirect()
        else:
            self.session = GnucapSessionPexpect(executable)

        if circuit != None and not isinstance(circuit, Circuit):
            raise ValueError('Circuit instance must be a %s instance'%
                             str(Circuit))
        
        if circuit != None:
            self.update_netlist()

    def command(self, command, parse_result = False):
        """Send gnucap command to session
        
        Returns the result from gnucap or a GnucapResult object if parse_result
        is True
        """
        result = self.session.command(command)

        if parse_result:
            return GnucapResult(result)
        
        return result

    def run_analysis(self, analysis):
        return analysis.run()

    def update_netlist(self):
        self.send_netlist(str(self.circuit))

    def send_netlist(self, netlist):
        netlistfile = tempfile.NamedTemporaryFile()

        ## In the development version one has to send
        ## a build command before the netlist
        netlistfile.write('The netlist\n' + netlist)

        netlistfile.flush()

        self.session.command('get ' + netlistfile.name)

        netlistfile.close()


## Exeptions
class EngineError(Exception): pass
class GnucapError(Exception): pass

