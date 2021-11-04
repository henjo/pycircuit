# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import pycircuit.post.cds.psf as psf

def createValue(psfobj, name, typename=None, value=None):
    if psfobj.getNSweeps() == 0:
        return psf.NonSweepValue(psfobj, name=name, typeid=psfobj.types.getTypeByName(typename).id, value=value)


def updatepsf(psfobj):
        psfobj.updateHeader()
        if psfobj.getNSweeps()==0:
            psfobj.values = psf.ValuesSectionNonSweep(psfobj)

# Begin -- grammar generated by Yapps
import sys, re
from yapps import runtime

class psfascScanner(runtime.Scanner):
    patterns = [
        ('"VALUE"', re.compile('VALUE')),
        ('"SWEEP"', re.compile('SWEEP')),
        ('"GROUP"', re.compile('GROUP')),
        ('"TRACE"', re.compile('TRACE')),
        ('"ARRAY"', re.compile('ARRAY')),
        ('"STRUCT\\("', re.compile('STRUCT\\(')),
        ('"\\*"', re.compile('\\*')),
        ('"STRING"', re.compile('STRING')),
        ('"BYTE"', re.compile('BYTE')),
        ('"INT"', re.compile('INT')),
        ('"DOUBLE"', re.compile('DOUBLE')),
        ('"FLOAT"', re.compile('FLOAT')),
        ('"TYPE"', re.compile('TYPE')),
        ('r"\\)"', re.compile('\\)')),
        ('r"\\("', re.compile('\\(')),
        ('"PROP"', re.compile('PROP')),
        ('"HEADER"', re.compile('HEADER')),
        ('"END"', re.compile('END')),
        ('\\s+', re.compile('\\s+')),
        ('FLOATNUM', re.compile('-?[0-9]+\\.[0-9e+-]*')),
        ('INTNUM', re.compile('-?[0-9]+')),
        ('STR', re.compile('"([^\\\\"]+|\\\\.)*"')),
        ('EOF', re.compile('$')),
        ('ARRAYLEN', re.compile('(\\*|[0-9])')),
    ]
    def __init__(self, str,*args,**kw):
        runtime.Scanner.__init__(self,None,{'\\s+':None,},str,*args,**kw)

class psfasc(runtime.Parser):
    Context = runtime.Context
    def psfasc(self, _parent=None):
        _context = self.Context(_parent, self._scanner, 'psfasc', [])
        psfobj = psf.PSFReader()
        headersec = self.headersec(psfobj, _context)
        psfobj.header = headersec
        if self._peek('"TYPE"', '"SWEEP"', '"TRACE"', '"END"', '"VALUE"', 'STR', context=_context) == '"TYPE"':
            typesec = self.typesec(psfobj, _context)
        if self._peek('"SWEEP"', '"TRACE"', '"END"', '"VALUE"', 'STR', context=_context) == '"SWEEP"':
            sweepsec = self.sweepsec(_context)
        if self._peek('"TRACE"', '"END"', '"VALUE"', 'STR', context=_context) == '"TRACE"':
            tracesec = self.tracesec(_context)
        updatepsf(psfobj)
        if self._peek('"END"', '"VALUE"', 'STR', context=_context) == '"VALUE"':
            valuesec = self.valuesec(psfobj, _context)
        self._scan('"END"', context=_context)
        EOF = self._scan('EOF', context=_context)
        return psfobj

    def name(self, _parent=None):
        _context = self.Context(_parent, self._scanner, 'name', [])
        STR = self._scan('STR', context=_context)
        return eval(STR)

    def headersec(self, psfobj, _parent=None):
        _context = self.Context(_parent, self._scanner, 'headersec', [psfobj])
        self._scan('"HEADER"', context=_context)
        header = psf.HeaderSection(psfobj)
        while self._peek('STR', '"TYPE"', '"SWEEP"', '"TRACE"', '"END"', '"VALUE"', context=_context) == 'STR':
            property = self.property(_context)
            header.addProperty(property)
        return header

    def property(self, _parent=None):
        _context = self.Context(_parent, self._scanner, 'property', [])
        name = self.name(_context)
        _token = self._peek('FLOATNUM', 'INTNUM', 'STR', context=_context)
        if _token == 'FLOATNUM':
            FLOATNUM = self._scan('FLOATNUM', context=_context)
            return psf.PropertyFloat64(name=name, value=FLOATNUM)
        elif _token == 'INTNUM':
            INTNUM = self._scan('INTNUM', context=_context)
            return psf.PropertyUInt(name=name, value=INTNUM)
        else: # == 'STR'
            STR = self._scan('STR', context=_context)
            return psf.PropertyString(name=name, value=eval(STR))

    def proplist(self, _parent=None):
        _context = self.Context(_parent, self._scanner, 'proplist', [])
        self._scan('"PROP"', context=_context)
        self._scan('r"\\("', context=_context)
        props = []
        while self._peek('r"\\)"', 'STR', context=_context) == 'STR':
            property = self.property(_context)
            props.append(property)
        self._scan('r"\\)"', context=_context)
        return props

    def typesec(self, psfobj, _parent=None):
        _context = self.Context(_parent, self._scanner, 'typesec', [psfobj])
        self._scan('"TYPE"', context=_context)
        typesec = psf.TypeSection(psfobj)
        while self._peek('STR', '"SWEEP"', '"TRACE"', '"END"', '"VALUE"', context=_context) == 'STR':
            type = self.type(psfobj, _context)
            psfobj.types.addType(type)
        return typesec

    def type(self, psfobj, _parent=None):
        _context = self.Context(_parent, self._scanner, 'type', [psfobj])
        name = self.name(_context)
        typedef = self.typedef(psfobj, _context)
        datatypedef=psf.DataTypeDef(psfobj, name=name, datatypeid=typedef[0], structdef=typedef[1])
        if self._peek('"PROP"', 'STR', 'r"\\)"', '"SWEEP"', '"TRACE"', '"END"', '"VALUE"', context=_context) == '"PROP"':
            proplist = self.proplist(_context)
            datatypedef.properties = proplist
        return datatypedef

    def typedef(self, psfobj, _parent=None):
        _context = self.Context(_parent, self._scanner, 'typedef', [psfobj])
        _token = self._peek('"STRUCT\\("', '"FLOAT"', '"INT"', '"ARRAY"', '"STRING"', context=_context)
        if _token == '"STRUCT\\("':
            typedefStruct = self.typedefStruct(psfobj, _context)
            return psf.TYPESTRUCT, typedefStruct
        elif _token == '"FLOAT"':
            typedefFloat = self.typedefFloat(_context)
            return psf.TYPEFLOATDOUBLE, None
        elif _token == '"INT"':
            typedefIntByte = self.typedefIntByte(_context)
            return psf.TYPEINTBYTE, None
        elif _token == '"ARRAY"':
            typedefArray = self.typedefArray(_context)
            return psf.TYPEARRAY, typedefArray
        else: # == '"STRING"'
            typedefString = self.typedefString(_context)
            return psf.TYPESTRING, None

    def typedefFloat(self, _parent=None):
        _context = self.Context(_parent, self._scanner, 'typedefFloat', [])
        self._scan('"FLOAT"', context=_context)
        self._scan('"DOUBLE"', context=_context)

    def typedefIntByte(self, _parent=None):
        _context = self.Context(_parent, self._scanner, 'typedefIntByte', [])
        self._scan('"INT"', context=_context)
        self._scan('"BYTE"', context=_context)

    def typedefString(self, _parent=None):
        _context = self.Context(_parent, self._scanner, 'typedefString', [])
        self._scan('"STRING"', context=_context)
        self._scan('"\\*"', context=_context)

    def typedefStruct(self, psfobj, _parent=None):
        _context = self.Context(_parent, self._scanner, 'typedefStruct', [psfobj])
        self._scan('"STRUCT\\("', context=_context)
        structdef = psf.StructDef()
        while self._peek('r"\\)"', 'STR', context=_context) == 'STR':
            type = self.type(psfobj, _context)
            structdef.children.append(type)
        self._scan('r"\\)"', context=_context)
        return structdef

    def typedefArray(self, _parent=None):
        _context = self.Context(_parent, self._scanner, 'typedefArray', [])
        self._scan('"ARRAY"', context=_context)
        self._scan('r"\\("', context=_context)
        ARRAYLEN = self._scan('ARRAYLEN', context=_context)
        self._scan('r"\\)"', context=_context)
        typedef = self.typedef(_context)
        return (typedef[0], ARRAYLEN)

    def tracesec(self, _parent=None):
        _context = self.Context(_parent, self._scanner, 'tracesec', [])
        self._scan('"TRACE"', context=_context)
        while self._peek('STR', '"END"', '"VALUE"', '"TRACE"', 'r"\\)"', '"SWEEP"', context=_context) == 'STR':
            tgroupdef = self.tgroupdef(_context)
        return None
        while self._peek('STR', '"END"', '"VALUE"', '"TRACE"', 'r"\\)"', '"SWEEP"', context=_context) == 'STR':
            tracedef = self.tracedef(_context)
        return None

    def tgroupdef(self, _parent=None):
        _context = self.Context(_parent, self._scanner, 'tgroupdef', [])
        name = self.name(_context)
        self._scan('"GROUP"', context=_context)
        INTNUM = self._scan('INTNUM', context=_context)
        while self._peek('STR', '"END"', '"VALUE"', '"TRACE"', 'r"\\)"', '"SWEEP"', context=_context) == 'STR':
            tracedef = self.tracedef(_context)

    def tracedef(self, _parent=None):
        _context = self.Context(_parent, self._scanner, 'tracedef', [])
        name = self.name(_context)
        STR = self._scan('STR', context=_context)
        if self._peek('"PROP"', 'STR', '"END"', '"VALUE"', '"TRACE"', 'r"\\)"', '"SWEEP"', context=_context) == '"PROP"':
            proplist = self.proplist(_context)

    def sweepdef(self, _parent=None):
        _context = self.Context(_parent, self._scanner, 'sweepdef', [])
        name = self.name(_context)
        STR = self._scan('STR', context=_context)
        if self._peek('"PROP"', 'STR', '"TRACE"', '"END"', '"VALUE"', 'r"\\)"', '"SWEEP"', context=_context) == '"PROP"':
            proplist = self.proplist(_context)

    def sweepsec(self, _parent=None):
        _context = self.Context(_parent, self._scanner, 'sweepsec', [])
        self._scan('"SWEEP"', context=_context)
        while self._peek('STR', '"TRACE"', '"END"', '"VALUE"', 'r"\\)"', '"SWEEP"', context=_context) == 'STR':
            sweepdef = self.sweepdef(_context)
        return None

    def valuesec(self, psfobj, _parent=None):
        _context = self.Context(_parent, self._scanner, 'valuesec', [psfobj])
        self._scan('"VALUE"', context=_context)
        while 1:
            value = self.value(psfobj, _context)
            psfobj.values.addValue(value)
            if self._peek('STR', '"END"', context=_context) != 'STR': break

    def value(self, psfobj, _parent=None):
        _context = self.Context(_parent, self._scanner, 'value', [psfobj])
        opttypename = None
        name = self.name(_context)
        if self._peek('STR', 'FLOATNUM', 'INTNUM', 'r"\\("', context=_context) == 'STR':
            STR = self._scan('STR', context=_context)
            opttypename = eval(STR)
        value = createValue(psfobj, name, opttypename)
        valuedata = self.valuedata(_context)
        value.setValue(valuedata)
        if self._peek('"PROP"', 'STR', '"END"', context=_context) == '"PROP"':
            proplist = self.proplist(_context)
            value.properties = proplist
        return value

    def nonsweptvalue(self, _parent=None):
        _context = self.Context(_parent, self._scanner, 'nonsweptvalue', [])
        name = self.name(_context)
        STR = self._scan('STR', context=_context)
        while 1:
            _token = self._peek('FLOATNUM', 'INTNUM', 'r"\\("', context=_context)
            if _token == 'r"\\("':
                structarrayvalue = self.structarrayvalue(_context)
            elif _token == 'FLOATNUM':
                FLOATNUM = self._scan('FLOATNUM', context=_context)
            else: # == 'INTNUM'
                INTNUM = self._scan('INTNUM', context=_context)
            if self._peek('FLOATNUM', 'INTNUM', 'r"\\("', '"PROP"', 'STR', 'r"\\)"', '"END"', '"VALUE"', '"TRACE"', '"SWEEP"', context=_context) not in ['FLOATNUM', 'INTNUM', 'r"\\("']: break
        if self._peek('"PROP"', 'STR', '"END"', '"VALUE"', '"TRACE"', 'r"\\)"', '"SWEEP"', context=_context) == '"PROP"':
            proplist = self.proplist(_context)
        return name

    def valuedata(self, _parent=None):
        _context = self.Context(_parent, self._scanner, 'valuedata', [])
        _token = self._peek('FLOATNUM', 'INTNUM', 'STR', 'r"\\("', context=_context)
        if _token == 'r"\\("':
            structarrayvalue = self.structarrayvalue(_context)
            return structarrayvalue
        elif _token == 'FLOATNUM':
            FLOATNUM = self._scan('FLOATNUM', context=_context)
            return float(FLOATNUM)
        elif _token == 'INTNUM':
            INTNUM = self._scan('INTNUM', context=_context)
            return int(INTNUM)
        else: # == 'STR'
            STR = self._scan('STR', context=_context)
            return eval(STR)

    def structarrayvalue(self, _parent=None):
        _context = self.Context(_parent, self._scanner, 'structarrayvalue', [])
        self._scan('r"\\("', context=_context)
        structval = []
        while self._peek('r"\\)"', 'FLOATNUM', 'INTNUM', 'STR', 'r"\\("', context=_context) != 'r"\\)"':
            valuedata = self.valuedata(_context)
            structval.append(valuedata)
        self._scan('r"\\)"', context=_context)
        return structval


def parse(rule, text):
    P = psfasc(psfascScanner(text))
    return runtime.wrap_error_reporter(P, rule)

def is_psfasc(filename):
    """Return true if a file is a PSF ascii file"""
    print(filename)
    with open(filename) as psffile:
        try:
            tmp = psffile.read(6)
        except UnicodeDecodeError:
            tmp = None
        if tmp == 'HEADER':
            return True
        else:
            return False

if __name__ == '__main__':
    from sys import argv, stdin
    if len(argv) >= 2:
        if len(argv) >= 3:
            f = open(argv[2],'r')
        else:
            f = stdin
        print(parse(argv[1], f.read()))
    else: sys.stderr.write('Args:  <rule> [<filename>]')
# End -- grammar generated by Yapps
