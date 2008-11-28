# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import skillparser
from types import *

"""
Python module for various Cadence Skill interfacing functions
"""

class Symbol:
    def __init__(self, name):
        self.name = name
    def __bool__(self):
        if self.name=='nil':
            return False
        elif self.name=='t':
            return True
        else:
            raise ValueError('Cannot convert Symbol to boolean')

    def __eq__(self, a):
        return self.name == a.name

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
    elif type(x) is BooleanType and x == False:
        return 'nil'
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

if __name__ == "__main__":
    import doctest
    doctest.testmod()
