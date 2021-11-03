# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import pycircuit.post.cds.skill as skill
import sys
from pycircuit.post.cds.cds import CadenceSession, find_virtuoso
import socket

class CadenceSessionTCP(CadenceSession):
    """Cadence session using a TCP/IP skill command server
    
    >>> s = CadenceSessionTCP()
    >>> s.send("(list 1 2 3)")
    (list 1 2 3)
    >>> s.callfunc('list', 1,2,3)
    (list 1 2 3)
    >>> s.list(1,2,3)
    (list 1 2 3)
    
    """
    def __init__(self, host='localhost', port=50008, verbose=False):
        self.verbose = verbose
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.connect((host, port))

    def send(self, expr):
        if self.verbose:
            print("Sending: "+expr)

        self.sock.send(expr)

        reply = self.sock.recv(16384)

        if reply.startswith('(') and reply.endswith(')'):
            result = skill.parse(reply[1:-1])
        elif not reply.startswith('(') and reply.endswith('()'):
            result = skill.parse(reply[:-2])
        else:
            raise Exception("Invalid reply from server: %s" % reply)

        return result

    def __del__(self):
        self.sock.send('bye')
        self.sock.close()

if __name__ == "__main__":
    import doctest
    doctest.testmod()
