# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import pycircuit.sim
import pycircuit.utilities.which as which

from pycircuit.sim.gnucap.circuit import Circuit
from pycircuit.sim.gnucap.result import GnucapResult

import os, os.path
import pexpect
import logging
import tempfile
import numpy as np

class Simulation(pycircuit.sim.Simulation):
    prompt = 'gnucap>'

    def __init__(self, circuit, executable = None):
        super(Simulation, self).__init__(circuit)

        if circuit != None and not isinstance(circuit, Circuit):
            raise ValueError('Circuit instance must be a %s instance'%
                             str(Circuit))
        
        if executable == None:
            if 'GNUCAP' in os.environ:
                executable = os.environ['GNUCAP']
            else:
                executable = which.which('gnucap')

        if not os.access(executable, os.X_OK):
            raise EngineError('gnucap executable not found')
            
        self.gnucap = executable
        
        self.setup()

        if circuit != None:
            self.update_netlist()

    def run_analysis(self, analysis):
        return analysis.run()

    def setup(self):
        session = pexpect.spawn(self.gnucap, timeout=2)
        
        session.expect(self.prompt)

        firstline = session.before.split('\n')[0]
        
        self.version = tuple((int(x) for x in 
                              firstline.split(' ')[1].split('.')))
        
        session.setecho(False)
        
        logging.info('Successfully established connection with gnucap')

        self.session = session
        
    def send(self, line):
        logging.debug('Sending: ' + line)
        self.session.sendline(line)
        
        self.session.expect(self.prompt)

        reply = self.session.before.strip()

        ## Strip first line 
        if self.version[0] < 2007:
            reply = '\n'.join(reply.split('\n')[1:]).strip()
        
        logging.debug('Got: ' + reply)

        return reply

    def send_command(self, command, xlabels = None):
        """Send gnucap command and parse results"""
        resultfile = tempfile.NamedTemporaryFile() 
        
        self.send(command + ' > %s'%resultfile.name)

        return GnucapResult(resultfile.name, xlabels)

    def update_netlist(self):
        self.send_netlist(str(self.circuit))

    def send_netlist(self, netlist):
        netlistfile = tempfile.NamedTemporaryFile()

        ## In the development version one has to send
        ## a build command before the netlist
        if self.version[0] > 2007:
            netlistfile.write('clear\nbuild\n' + netlist)
        else:
            netlistfile.write(netlist)

        netlistfile.flush()

        reply = self.send('include ' + netlistfile.name)

        if len(reply) > 0:
            raise GnucapError(reply)
        
        netlistfile.close()

    def close(self):
        self.session.sendline('exit')
        self.session.expect(pexpect.EOF)

    def __del__(self):
        self.close()


## Exeptions
class EngineError(Exception): pass
class GnucapError(Exception): pass

