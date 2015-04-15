import re
import logging
import pexpect
import pycircuit.utilities.which as which
import os, os.path
import sys, StringIO

class NGspiceSession(object):
    def command(self, cmd):
        pass
    def lock(self):
        pass
    def unlock(self):
        pass

class NGspiceSessionPexpect(NGspiceSession):
    prompt = 'ngspice '
    
    def __init__(self, executable = None):
        if executable == None:
            if 'NGSPICE' in os.environ:
                executable = os.environ['NGSPICE']
            else:
                executable = which.which('ngspice')

        if not os.access(executable, os.X_OK):
            raise EngineError('ngspice executable not found')
            
        self.ngspice = executable
        
        self._setup()

    def command(self, line):
        logging.debug('Sending: ' + line)
        self.session.sendline(line)
        
        self.session.expect(self.prompt)

        reply = self.session.before.strip()

        ## remove linefeeds
        reply = re.sub("\r","",reply)
        
        reply = '\n'.join(reply.split('\n')[1:]).strip()
        
        logging.debug('Got: ' + reply)

        return reply

    def _setup(self):
        session = pexpect.spawn(self.ngspice, timeout=2)
        
        session.expect(self.prompt)

        firstlines = session.before.split('\n')
        
        self.version = firstlines[1].strip('** ').split(' :')[0]


        #self.version = tuple((int(x) for x in 
        #                      firstline.split(' ')[1].split('.')))
        
        session.setecho(False)
        
        logging.info('Successfully established connection with ngspice version '+
                     self.version)

        self.session = session
        

    def __del__(self):
        self.session.sendline('exit')
        self.session.expect(pexpect.EOF)
