# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import pycircuit.post.cds.skillparser as skillparser
from types import *

"""
Python module for various Cadence Skill interfacing functions
"""

class SkillObject(object):
    """Class that executes an expression and keeps track of result in the Skill parser"""
    def __init__(self, session):
        self.session = session
        self._varname = 'tmpvar_%s'%id(self)
        self.valid = False
        self.value = None
    	
    @property
    def varname(self):
        if not self.valid:
            raise ValueError('SkillObject referenced before eval was called')
        return self._varname
    
    def eval(self, expr):
        self.session.send('%s=%s'%(self._varname, expr), parse=False)
        self.value = self.session.send('%s' % self._varname, 
    				   parse=True, parseall=True)
        self.valid = True
        return self.value

    def __eq__(self, a): return self.value == a
    def __str__(self): return str(self.value)
    def __nonzero__(self): return bool(self.value)
    def __repr__(self): return repr(self.value)

    def __del__(self):
        ## FIXME
        ## self.varname = 'unbound
        pass



class Symbol(object):
    def __init__(self, name):
        self.name = name
    def __nonzero__(self):
        if self.name=='nil':
            return False
        elif self.name=='t':
            return True
        else:
            raise ValueError('Cannot convert Symbol to boolean')

    def __eq__(self, a):
        return isinstance(a, Symbol) and self.name == a.name

    def __hash__(self, a):
        return hash(self.name)

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

def toSkill(x):
    """Convert python type to skill code

	>>> toSkill([1,2,3,"apa"])
	'(list 1 2 3 "apa")'
	>>> toSkill([1,False,"apa",True])
	'(list 1 nil "apa" t)'
        
	"""
    if type(x) is ListType or type(x) is TupleType:
        return "(list " + " ".join([toSkill(e) for e in x]) +")"
    elif type(x) is StringType:
        return '"' + x + '"'
    elif type(x) is BooleanType and x == True:
        return 't'
    elif type(x) is BooleanType and x == False or x is None:
        return 'nil'
    elif isinstance(x, Symbol):
        return "'%s"%x.name
    elif isinstance(x, SkillObject):
        return x.varname
    else:
        return str(x)

def parse(str):
	"""
	Parse skill expression.
	@return Python representation of str

	>>> parse("1")==1
	True
	>>> parse("1.0")==1.0
	True
	>>> parse('"apa"')=="apa"
	True
	>>> parse('apa')!="apa"
	True
	>>> parse("(1 2 3)")==[1,2,3]
	True
	>>> parse("nil").name == 'nil'
	True
	>>> parse("nil")==Symbol('nil')
	True
	>>> parse("(1 ?a)")==[1,Symbol('?a')]
	True

	"""
	return skillparser.parse('expr', str)

## Define some useful constants
nil = Symbol('nil')

if __name__ == "__main__":
    import doctest
    doctest.testmod()
