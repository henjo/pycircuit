import re
import logging
import pexpect
import pycircuit.utilities.which as which
import os, os.path
import sys, StringIO

class GnucapSession(object):
    def command(self, cmd):
        pass
    def lock(self):
        pass
    def unlock(self):
        pass

class GnucapSessionPexpect(GnucapSession):
    prompt = 'gnucap>'
    
    def __init__(self, executable = None):
        if executable is None:
            if 'GNUCAP' in os.environ:
                executable = os.environ['GNUCAP']
            else:
                executable = which.which('gnucap')

        if not os.access(executable, os.X_OK):
            raise EngineError('gnucap executable not found')
            
        self.gnucap = executable
        
        self._setup()

    def command(self, line):
        logging.debug('Sending: ' + line)
        self.session.sendline(line)
        
        self.session.expect(self.prompt)

        reply = self.session.before.strip()

        ## Strip first line 
#        if self.version[0] < 2007:
        
        ## remove linefeeds
        reply = re.sub("\r","",reply)
        
        reply = '\n'.join(reply.split('\n')[1:]).strip()
        
        logging.debug('Got: ' + reply)

        return reply

    def _setup(self):
        session = pexpect.spawn(self.gnucap, timeout=2)
        
        session.expect(self.prompt)

        firstline = session.before.split('\n')[0]
        
        self.version = tuple((int(x) for x in 
                              firstline.split(' ')[1].split('.')))
        
        session.setecho(False)
        
        logging.info('Successfully established connection with gnucap version '+
                     str(self.version))

        self.session = session
        

    def __del__(self):
        self.session.sendline('exit')
        self.session.expect(pexpect.EOF)

class GnucapSessionDirect(GnucapSession):
    def __init__(self):
        self.__class__.locked = False
        
        import gnucap
        self.gnucap = gnucap
        
        gnucap.command("set lang=acs")
        self.runmode = gnucap.SET_RUN_MODE(gnucap.rBATCH)
        
    def lock(self):
        if self.__class__.locked == True:
            raise(Exception("Session is already locked"))
        else:
            self.__class__.locked = True
            
    def unlock(self):
        self.__class__.locked = False
        
    def command(self, command): 
        logging.debug('Sending: ' + command)
        reply = self.gnucap.command(command).strip()
        logging.debug('Got: ' + reply)
        return reply

