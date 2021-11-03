# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import pycircuit.post.cds.skill as skill
import sys
from pycircuit.post.cds.cds import CadenceSession

class CadenceSessionFile(CadenceSession):
	"""Class to handle that writes to a skill file instead of real cadence session.

	>>> s = CadenceSessionFile(file=sys.stdout)
	>>> s.send("(list 1 2 3)")
	(list 1 2 3)
	>>> s.callfunc('list', 1,2,3)
	(list 1 2 3)
	>>> s.list(1,2,3)
	(list 1 2 3)
	>>> s.testfunc(1,2,3,apa=3,b="test",c=False)
	(testfunc 1 2 3 ?c nil ?apa 3 ?b "test")
	
	"""
	def __init__(self, file=sys.stdout, timeout=30,verbose=False):
		self.file = file
		self.verbose = verbose

	def send(self, expr):
		if self.verbose:
			print "Sending: "+expr
		self.file.write(expr+'\n')
		return None
	def __del__(self):
		pass


if __name__ == "__main__":
    import doctest
    doctest.testmod()

