# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import unittest
import struct, os, re
import operator
import numpy
import psfasc
from copy import copy

from struct import unpack, pack

class PSFInvalid(Exception):
    pass

def warning(str):
    print "Warning: "+str

def indent(str, n=2):
    return "\n".join([' '*n+s for s in str.split("\n")])

class PSFData(object):
    @classmethod
    def fromFile(cls, file):
        obj = cls()
        obj.deSerializeFile(file)
        return obj

    size=None
    
    def __init__(self, value=None, extarg=None):
        self.value = value
        self.extarg = extarg
    def setValue(self, value):
        self.value = value        
    def __eq__(self, a):
        return self.value == a
    def __cmp__(self, a):
        return cmp(self.value, a)
    def __hash__(self):
        return hash(self.value)
    def deSerializeFile(self, file):
        pass
    def getSize(self):
        self.size
    def getValue(self):
        return self.value        
    def __str__(self):
        return str(self.value)
    def toPSFasc(self, prec=None):
        return str(self)
    def __repr__(self):
        return self.value.__repr__()

class PSFNumber(PSFData):
    def __int__(self):
        return self.value
    def __add__(self, a):
        return UInt32(self.value+int(a))
    def __mul__(self, a):
        return UInt32(self.value*int(a))
    def __radd__(self, a):
        return UInt32(self.value+int(a))
    def __sub__(self, a):
        return UInt32(self.value-int(a))
    def __rsub__(self, a):
        return UInt32(int(a)-self.value)
    def __div__(self, a):
        return UInt32(self.value/int(a))
    def __rdiv__(self, a):
        return UInt32(int(a)/self.value)
    def __floordiv__(self, a):
        return UInt32(self.value//int(a))
    def __rfloordiv__(self, a):
        return UInt32(int(a)//self.value)
    def __mod__(self, a):
        return UInt32(self.value%int(a))

class Int8(PSFNumber):
    size=4
    def deSerializeFile(self, file, size=None):
        data=file.read(self.size)
        self.value = unpack("b",data[3])[0]
class UInt8(PSFNumber):
    size=4
    def deSerializeFile(self, file, size=None):
        data=file.read(self.size)
        self.value = unpack("B",data[3])[0]
class Int32(PSFNumber):
    size=4
    def deSerializeFile(self, file, size=None):
        self.value = unpack(">i",file.read(self.size))[0]
class UInt32(PSFNumber):
    size=4
    def deSerializeFile(self, file, size=None):
        self.value = unpack(">I",file.read(self.size))[0]
class Int64(PSFNumber):
    size=8
    def __int__(self):
        return self.value
    def deSerializeFile(self, file, size=None):
        self.value = unpack(">q",file.read(self.size))[0]
class UInt64(PSFNumber):
    size=8
    def __int__(self):
        return self.value
    def deSerializeFile(self, file, size=None):
        self.value = unpack(">Q",file.read(self.size))[0]

class Float64(PSFNumber):
    size=8
    def __float__(self):
        return float(self.value)

    def toPSFasc(self, prec=6):
        if prec:
            fmt=('%%#%dg'%prec)
        else:
            fmt='%#g'
        return fmt%self.value

    def deSerializeFile(self, file, size=None):
        self.value = unpack(">d",file.read(self.size))[0]

class Float32(PSFNumber):
    size=4
    def __float__(self):
        return float(self.value)
    def deSerializeFile(self, file, size=None):
        self.value = unpack(">f",file.read(self.size))[0]

class ComplexFloat64(PSFNumber):
    size=16
    def toPSFasc(self, prec=6):
        if prec:
            fmt=('%%#%dg'%prec)
        else:
            fmt='%#g'
            
        return "(" + fmt%self.value.real + " " + fmt%self.value.imag + ")"
    
    def deSerializeFile(self, file, size=None):
        re,im = unpack(">dd",file.read(self.size))
        self.value = complex(re,im)

class String(PSFData):
    def __str__(self):
        return self.value
    def deSerializeFile(self, file, size=None):
        self.len = unpack(">I",file.read(4))[0]
        if self.len < 0x100:
            self.value = file.read(self.len)
            # Pad to 32-bit boundary
            file.read((4-self.len)%4)
        else:
            raise Exception("String too long %d"%self.len)

    def toPSFasc(self, prec=None):
        return "\""+str(self.value)+"\""

class Struct(PSFData):
    def __init__(self, structdef, value=None):
        self.structdef = structdef
        self.value = {}

        if value:
            self.setValue(value)

    def __getitem__(self, key):
        return self.value[key]

    def getValue(self):
        return dict([(k,v.getValue()) for k,v in self.value.items()])

    def setValue(self, value):
        assert(value != None and len(value) == len(self.structdef.children))
        for element, val in zip(self.structdef.children, value):
            valueobj = element.getDataObj()
            valueobj.setValue(val)
            self.value[element.name] = valueobj            

    def deSerializeFile(self, file):
        for element in self.structdef.children:
            value = element.getDataObj()
            value.deSerializeFile(file)
            self.value[element.name] = value

    def toPSFasc(self, prec=None):
        s="(\n"
        for element in self.structdef.children:
            s+=self.value[element.name].toPSFasc(prec)+"\n"
        s+=")"
        return s

    def __repr__(self):
        return "\n".join([indent(s) for s in map(repr,self.value.items())]) + "\n"

class Array(PSFData):
    def setValue(self, value):
        dataclass, length = self.extarg

        if value != None:
            self.children = [dataclass(value=val) for val in value]
        else:
            self.children = [dataclass(value=None) for val in range(length)]

    def getValue(self):
        return [v.getValue() for v in self.children]

    def __iter__(self):
        return self.children.__iter__()

    def __tuple__(self):
        return tuple(self.children)

    def __repr__(self):
        return "\n".join([indent(s) for s in map(str,self.children)]) + "\n"

class Chunk:
    """Base class for chunk"""
    def __init__(self, psf=None, type=None):
        self.psf = psf
        self.fileoffset=None
        if not hasattr(self.__class__, 'type'):
            self.type = type
        self.verbose = False
        self.name = ""
            
    def deSerializeFile(self, file):
        self.fileoffset = file.tell()

        type = UInt32.fromFile(file)
        if (self.type != None) and self.type != type:
            file.seek(-UInt32.size, 1)
            raise IncorrectChunk(type, self.type)

    def __repr__(self):
        return self.__class__.__name__

class NextSectionType(Chunk):
    type=1
class NextSectionSweep(Chunk):
    type=2
class NextSectionTrace(Chunk):
    type=3
class NextSectionValues(Chunk):
    type=4
class EndOfStructDef(Chunk):
    type=18
    
NextSectionClasses = [NextSectionType, NextSectionSweep, NextSectionTrace, NextSectionValues]

class Property(Chunk):
    type=None
    valueclass=None
    def __init__(self, name=None, value=None):
        Chunk.__init__(self)
        self.name = String(name)
        self.value = self.valueclass(value)
    def deSerializeFile(self, file):
        Chunk.deSerializeFile(self, file)

        self.name = String.fromFile(file)
        self.value = self.valueclass.fromFile(file)
    def toPSFasc(self, prec=9):
        return self.name.toPSFasc() + " " + self.value.toPSFasc(prec=prec)
    def __repr__(self):
        return self.__class__.__name__+"("+str(self.name)+","+str(self.value)+")"
class PropertyString(Property):
    type=33
    valueclass=String
class PropertyUInt(Property):
    type=34
    valueclass=UInt32
class PropertyFloat64(Property):
    type=35
    valueclass=Float64
    
PropertyClasses = [PropertyString, PropertyUInt, PropertyFloat64]

TYPEFLOATDOUBLE = 11
TYPEINTBYTE = 1
TYPECOMPLEXDOUBLE = 12
TYPESTRUCT = 16
TYPESTRING = 2 ## Incorrect number
TYPEARRAY = 3 ## Incorrect number
TYPEINTLONG = 5
class DataTypeDef(Chunk):
    """Class representing data type of waveform data"""
    type=16
    ClassDict = {
        TYPEFLOATDOUBLE: Float64,
        TYPEINTBYTE: Int8,
        TYPECOMPLEXDOUBLE: ComplexFloat64,
        TYPESTRING: String,
        TYPEARRAY: Array,
        TYPEINTLONG: Int32
    }
    PSFASCDict = {
        TYPEFLOATDOUBLE: "FLOAT DOUBLE",
        TYPEINTBYTE: "INT BYTE",
        TYPECOMPLEXDOUBLE: "COMPLEX DOUBLE",
        TYPESTRING: "STRING *",
        TYPEINTLONG: "INT LONG"
    }
    
    def __init__(self, psf, id=0, name=None, datatypeid=0, structdef=None):
        Chunk.__init__(self, psf, type)
        self.id = id
        self.name = name
        self.datatypeid = datatypeid

        self.structdef = structdef
        self.properties = []

    def getDataObj(self):
        """Get a data object described by the DataType"""
        if self.datatypeid == TYPESTRUCT:
            return self.structdef.getDataObj()
        elif self.datatypeid == TYPEARRAY:
            return Array(extarg=(self.ClassDict[self.structdef[0]], self.structdef[1]))
        else:
            return self.ClassDict[self.datatypeid](extarg=self.structdef)

    def toPSFasc(self, prec=None):
        r=self.name.toPSFasc(prec) + " "
        if self.datatypeid == TYPESTRUCT:
            r+=self.structdef.toPSFasc(prec)
        elif self.datatypeid == TYPEARRAY:
            r+="ARRAY ( %s ) "%str(self.structdef[1])+self.PSFASCDict[self.structdef[0]]
        else:
            r+= self.PSFASCDict[self.datatypeid]
        if len(self.properties)>0:
            r+=" PROP(\n"
            r+="\n".join([prop.toPSFasc(prec) for prop in self.properties])
            r+="\n)"
        return r

    def getDataSize(self):
        if self.datatypeid == TYPESTRUCT:
            return self.structdef.getDataSize()
        else:
            return self.ClassDict[self.datatypeid].size
            
    def deSerializeFile(self, file):
        start = file.tell()
        Chunk.deSerializeFile(self, file)
        self.id = UInt32.fromFile(file)
        self.name = String.fromFile(file)

        arraytype = UInt32.fromFile(file)

        self.datatypeid = UInt32.fromFile(file)

        if arraytype != 0:
            self.datatypeid, self.structdef = TYPEARRAY, (UInt32.fromFile(file), self.datatypeid)
            
        if self.datatypeid == 16:
            self.structdef = StructDef.fromFile(file, self.psf)

        # Read possible property objects that belongs to the type by peeking ahead
        while True:
            oldpos = file.tell()
            try:
                prop = readChunk(self.psf, file, expectedclasses=PropertyClasses)
                self.properties.append(prop)
            except ValueError:
                file.seek(oldpos)
                break

                
    def __repr__(self):
        return self.__class__.__name__+"("+str({"name":self.name,"id":"0x%x"%self.id, "datatypeid":self.datatypeid, 
                                                "properties":self.properties})+")"

class DataTypeRef(Chunk):
    type=16
    """Class representing link to data type"""
    def __init__(self, psf, type=None):
        Chunk.__init__(self, psf, type)
        self.id = None
        self.name = None
        self.datatypeid = 0
        self.properties = []

    def getDataObj(self):
        """Get a data object described by the DataType"""
        return self.psf.types.idMap[self.datatypeid].getDataObj()

    def toPSFasc(self, prec=None):
        r=self.name.toPSFasc(prec) + " "
        r+=self.psf.types.idMap[self.datatypeid].name.toPSFasc()
        if len(self.properties)>0:
            r+=" PROP(\n"
            r+="\n".join([prop.toPSFasc(prec) for prop in self.properties])
            r+="\n)"
        return r

    def getDataSize(self):
        return self.psf.types.idMap[self.datatypeid].getDataSize()
        
    def deSerializeFile(self, file):
        start = file.tell()
        Chunk.deSerializeFile(self, file)
        self.id = UInt32.fromFile(file)
        self.name = String.fromFile(file)

        self.datatypeid = UInt32.fromFile(file)

        assert(self.datatypeid != 0)

        # Read possible property objects that belongs to the type by peeking ahead
        while True:
            oldpos = file.tell()
            try:
                prop = readChunk(self.psf, file, expectedclasses=PropertyClasses)
                self.properties.append(prop)
            except ValueError:
                file.seek(oldpos)
                break

    def __repr__(self):
        return self.__class__.__name__+"("+str({"name":self.name,"id":"0x%x"%self.id, "datatypeid":self.datatypeid,
                                                "properties":self.properties})+")"

class StructDef(PSFData):
    """Class representing struct definition"""
    @classmethod
    def fromFile(cls, file, psf):
        obj = cls()
        obj.deSerializeFile(file, psf)
        return obj

    def __init__(self):
        self.children = []

    def getDataObj(self):
        return Struct(self)

    def getDataSize(self):
        return sum([child.getDataSize() for child in self.children])

    def toPSFasc(self, prec=None):
        s="STRUCT(\n"
        for child in self.children:
            s+=child.toPSFasc(prec)+"\n"
        s+=")"
        return s
    
    def deSerializeFile(self, file, psf):
        while True:
            chunk = readChunk(psf, file, expectedclasses=[DataTypeDef, EndOfStructDef])
            if isinstance(chunk, EndOfStructDef):
                break
            else:
                self.children.append(chunk)

    def __repr__(self):
        return self.__class__.__name__ + "(\n"+\
               "\n".join(map(str,self.children))+\
               ")\n"

class SimpleContainer(Chunk):
    type = 21
    def __init__(self, psf, type=None, childrenclslist=None, childrenclsignore=None):
        Chunk.__init__(self, psf, type)
        self.section = None
        self.children = []
        self.childrenclslist = childrenclslist
        self.childrenclsignore = childrenclsignore
        self.endpos = None

    def getChunks(self):
        return  self.children
        
    def deSerializeFile(self, file):
        Chunk.deSerializeFile(self, file)
        self.endpos = UInt32.fromFile(file).value
        self.children = []

        while file.tell() < self.endpos:
            chunk = readChunk(self.psf, file, expectedclasses=self.childrenclslist+self.childrenclsignore)
            if chunk.__class__ in self.childrenclslist:
                self.children.append(chunk)

        # Read trailing bytes
        if self.endpos-file.tell() != 0:
            warning("%d trailing bytes in %s"%(self.endpos-file.tell(), self.__class__.__name__))
            
        self.tail = file.read(self.endpos-file.tell())
        file.seek(self.endpos)

    def __repr__(self):
        s=""
        if self.fileoffset:
            s+= "0x%x"%self.fileoffset+ ":"
        s+= self.__class__.__name__  + "(" + str(self.type) +")"
        if self.endpos and self.fileoffset:
            s+= "size="+str(self.endpos-self.fileoffset)
        s+= "\n" + "\n".join([indent(s) for s in map(str,self.children)]) + "\n"
        return s

class Container22(Chunk):
    type=22
    def __init__(self, psf, type=None, n=None, childrenclslist=None):
        Chunk.__init__(self, psf, 22)
        self.section = None
        self.children = []
        self.childrenclslist = childrenclslist
        self.endpos = None

    def getChunks(self):
        return  self.children

    def deSerializeFile(self, file):
        Chunk.deSerializeFile(self, file)
        self.endpos = UInt32.fromFile(file).value # Save end position of Container

        self.children = []
        while file.tell() < self.endpos:
            chunk = readChunk(self.psf, file,
                              expectedclasses=self.childrenclslist)
            self.children.append(chunk)

        # Read trailing bytes
        if self.endpos-file.tell() != 0:
            warning("%d trailing bytes in %s"%(self.endpos-file.tell(), self.__class__.__name__))
            
        self.tail = file.read(self.endpos-file.tell())
        file.seek(self.endpos)

    def __repr__(self):
        return "0x%x"%self.fileoffset +":" + self.__class__.__name__  +\
               "(" + str(self.type) +")" + "\n" + "\n".join([indent(s) for s in map(str,self.children)]) + "\n"


class ZeroPad(Chunk):
    type = 20
    def deSerializeFile(self, file):
        Chunk.deSerializeFile(self, file)
        size = UInt32.fromFile(file).value
        self.endpos = file.tell() + size
        file.seek(self.endpos)

class HashTable(Chunk):
    type = 19
    """Class representing offset of trace data"""
    def __init__(self, psf, n=None):
        Chunk.__init__(self, psf, type)
        self.children = []
        self.extra=[]
    def deSerializeFile(self, file):
        Chunk.deSerializeFile(self, file)

        startpos = file.tell()
        size = UInt32.fromFile(file)

        for i in range(0, size/8):
            id = UInt32.fromFile(file)
            offset = UInt32.fromFile(file)
            self.children.append((id, offset))

    def __repr__(self):
        return self.__class__.__name__+"\n"+ "\n".join(["  0x%x: 0x%x"%(k,v.value) for k,v in self.children])+")"

class HashTableTrace(Chunk):
    type = 19
    """Class representing offset of trace data"""
    def __init__(self, psf):
        Chunk.__init__(self, psf, type)
        self.children = []
    def deSerializeFile(self, file):
        Chunk.deSerializeFile(self, file)
        
        self.size = UInt32.fromFile(file)

        for i in range(0, self.size.value/16):
            id = UInt32.fromFile(file)
            offset = UInt32.fromFile(file)
            data1 = UInt32.fromFile(file).value
            data2 = UInt32.fromFile(file).value
            self.children.append((id,offset,data1,data2))

    def __repr__(self):

        return self.__class__.__name__+"\n"+ "\n".join(["  %s: 0x%x 0x%x 0x%x"%(pack(">I",k.value),v.value,d1,d2) for k,v,d1,d2 in self.children])+")"

class HashContainer(Chunk):
    type=21
    hashclass = HashTable
    def __init__(self, psf, childrenclslist=None, childrenclsignore=None):
        Chunk.__init__(self, psf, type)
        self.section = None
        self.children = []
        self.childrenclslist = childrenclslist
        self.childrenclsignore = childrenclsignore
        self.endpos = None
        self.hashtable = None

    def __len__(self):
        return len(self.children)

    def getChunks(self):
        return  self.children

    def deSerializeFile(self, file):
        Chunk.deSerializeFile(self, file)
        self.endpos = UInt32.fromFile(file).value
        self.children = []

        self.data = Container22(self.psf,
                                childrenclslist=self.childrenclslist)
        self.data.deSerializeFile(file)

        self.hashtable = self.hashclass(self.psf)
        self.hashtable.deSerializeFile(file)

        # Copy children reference from data
        self.children = self.data.children

        self.section = UInt32.fromFile(file)
        
        # Read trailing bytes
        if self.endpos-file.tell() != 0:
            warning("%d trailing bytes in %s"%(self.endpos-file.tell(), self.__class__.__name__))
            
        self.tail = file.read(self.endpos-file.tell())
        file.seek(self.endpos)

    def __repr__(self):
        s=""
        if self.fileoffset:
            s += "0x%x"%self.fileoffset +":" 
        s += self.__class__.__name__  + "(" + str(self.type) +")"
        if self.endpos:
            s+=" size="+str(self.endpos-self.fileoffset) + "\n"
        s += "\n".join([indent(s) for s in map(str,(self.children, self.hashtable))]) + "\n"
        return s
    
class HeaderSection(SimpleContainer):
    type=21
    def __init__(self, psf, n=None):
        SimpleContainer.__init__(self,psf, childrenclslist=PropertyClasses,
                                 childrenclsignore=NextSectionClasses)
        self.properties = {}

    def addProperty(self, prop):
        """Add property to header"""
        self.children.append(prop)
        self.properties[prop.name] = prop.value

    def deSerializeFile(self, file):
        SimpleContainer.deSerializeFile(self, file)

        # Read header properties
        self.properties = {}
        for prop in self.children:
            self.properties[prop.name] = prop.value

    def toPSFasc(self, prec=None):
        r="HEADER\n"
        r+='"PSFversion" "1.00"\n'
        r+="\n".join([child.toPSFasc(prec) for child in self.children \
                      if not child.name.value[0:3].upper() == 'PSF'])
        return r    

class SweepSection(SimpleContainer):
    type=21
    def __init__(self, psf):
        SimpleContainer.__init__(self, psf, childrenclslist=[DataTypeRef],
                                 childrenclsignore=NextSectionClasses)

    def deSerializeFile(self, file):
        SimpleContainer.deSerializeFile(self, file)
        # Read header properties
        self.idMap = {}

        for chunk in self.children:
            self.idMap[chunk.id] = chunk

    def getSweep(self, id):
        return self.idMap[id]
    
    def getNames(self):
        return tuple([str(child.name) for child in self.children])

    def toPSFasc(self, prec=None):
        r="SWEEP\n"
        r+="\n".join([child.toPSFasc(prec) for child in self.children])
        return r


class TypeSection(HashContainer):
    def __init__(self, psf):
        HashContainer.__init__(self, psf, childrenclslist=[DataTypeDef],
                               childrenclsignore=NextSectionClasses)
        
        self.idMap = {}
        self.nameMap = {}
        
    def addType(self, type):
        type.id = self.psf.allocId()
        self.children.append(type)
        self.idMap[type.id] = type
        self.nameMap[type.name] = type
        
    def getType(self, id):
        return self.idMap[id]

    def getTypeByName(self, name):
        return self.nameMap[name]
                
    def deSerializeFile(self, file):
        HashContainer.deSerializeFile(self, file)

        # Read header properties
        self.idMap = {}

        for chunk in self.children:
            self.idMap[chunk.id] = chunk
            self.nameMap[chunk.name] = type
            
    def toPSFasc(self, prec=None):
        r="TYPE\n"
        r+="\n".join([child.toPSFasc(prec) for child in self.children])
        return r


class TraceSection(HashContainer):
    hashclass = HashTableTrace
    def __init__(self, psf):
        HashContainer.__init__(self, psf, childrenclslist=[GroupDef, DataTypeRef])
        self.idMap = {}
        self.nameIndex = {}
        
    def deSerializeFile(self, file):
        HashContainer.deSerializeFile(self, file)

        self.idMap = {}
        
        for index, chunk in enumerate(self.children):
            self.idMap[chunk.id] = chunk
            if isinstance(chunk, GroupDef):
                self.nameIndex.update(dict([(par, (index,)+value) for par,value in chunk.getNameIndex().items()]))
            else:
                self.nameIndex[chunk.name] = (index,)

    def getNameIndex(self):
        return self.nameIndex
            
    def toPSFasc(self, prec=None):
        r="TRACE\n"
        r+="\n".join([child.toPSFasc(prec) for child in self.children])
        return r
    def getTraceNames(self):
        result = []
        for trace in self.children:
            if isinstance(trace,GroupDef):
                result += trace.getNames()
            else:
                result.append(trace.name)
        return tuple(map(str, result))
    def getTraceIndexByName(self, name):
        """Returns an index to the given trace name
        
        The index is hierarchical so if if the traces are divided into 2 groups the index (0,1) means
        child 1 of group 0
        
        >>> psf=PSFReader('./test/psf/srcSweep')
        >>> psf.open()
        >>> psf.traces.getTraceIndexByName("VIN")
        (0, 1)

        >>> psf=PSFReader('./test/resultdirs/parsweep2/C=1e-12,R=1e-12/psf/ac.ac')
        >>> psf.open()
        >>> psf.traces.getTraceIndexByName("net3")
        (0,)

        """
        return self.nameIndex[name]

class ValuesSectionNonSweep(HashContainer):
    type=21
    def __init__(self, psf):
        HashContainer.__init__(self, psf, childrenclslist=[NonSweepValue])
        self.idMap={}
        self.nameMap={}

    def addValue(self, value):
        value.id = self.psf.allocId()
        if not isinstance(value, NonSweepValue):
            raise ValueError("Value should be a NonSweepValue")
        self.idMap[value.id] = value
        self.nameMap[value.name] = value
        self.children.append(value)

    def deSerializeFile(self, file):
        HashContainer.deSerializeFile(self, file)

        for child in self.children:
            self.nameMap[child.name] = child

    def getValuePropertiesByName(self, name):
        return dict([(prop.name, prop.value) for prop in self.nameMap[name].properties])

    def getValueByName(self, name):
        return self.nameMap[name].getValue()

    def getValueNames(self):
        return tuple([child.name for child in self.children])

    def toPSFasc(self, prec=None):
        r="VALUE\n"
        r+="\n".join([child.toPSFasc(prec) for child in self.children])
        return r

class ValuesSectionSweep(SimpleContainer):
    type=21
    def deSerializeFile(self, file):
        Chunk.deSerializeFile(self, file)
        self.endpos = UInt32.fromFile(file).value

        windowedsweep = self.psf.header.properties.has_key('PSF window size')

        if windowedsweep:
            el = ZeroPad(self.psf)
            el.deSerializeFile(file)

        isweep=0
        while isweep < self.psf.header.properties['PSF sweep points']:
            if windowedsweep:
                value = SweepValueWindowed(self.psf)
            else:
                value = SweepValueSimple(self.psf)

            isweep += value.deSerializeFile(file, n=self.psf.header.properties['PSF sweep points']-isweep)

            self.children.append(value)

        self.section = UInt32.fromFile(file)

        # Read trailing bytes
        if self.endpos-file.tell() != 0:
            warning("%d trailing bytes in %s"%(self.endpos-file.tell(), self.__class__.__name__))
            self.tail = file.read(self.endpos-file.tell())

        file.seek(self.endpos)

    def getSweepParamValues(self):
        return reduce(operator.__add__, [child.getSweepParamValues() for child in self.children])

    def getValueNames(self):
        return self.psf.traces.getTraceNames()

    def __len__(self):
        return len(self.psf.traces)

    def getValueByName(self, name):
        windowedsweep = self.psf.header.properties.has_key('PSF window size')

        index = self.psf.traces.getTraceIndexByName(name)

        result = []
        for child in self.children:
            obj=child
            for i in index:
                obj = obj.children[i]

            # If windowed sweep, each child will be a list of values in the window
            if windowedsweep:
                result += [v.getValue() for v in obj]
            else:
                result.append(obj.getValue())
                
        return numpy.array(result)

    def toPSFasc(self, prec=None):
        r="VALUE\n"
        r+="\n".join([child.toPSFasc(prec) for child in self.children])
        return r

class NonSweepValue(Chunk):
    type=16
    def __init__(self, psf, id=None, typeid=None, name=None, value=None):
        Chunk.__init__(self, psf, type)
        self.id = id
        self.name = name
        self.typeid = typeid

        if typeid:
            self.valuetype = self.psf.types.idMap[self.typeid]
        else:
            self.valuetype = None

        if value:
            self.value = value
        elif self.valuetype:
            self.value = self.valuetype.getDataObj()
        else:
            self.value = None
        self.properties = []

    def getValue(self):
        return self.value.getValue()

    def setValue(self, value):
        self.value.setValue(value)

    def deSerializeFile(self, file):
        startpos = file.tell()
        Chunk.deSerializeFile(self, file)
        self.id = UInt32.fromFile(file)
        self.name = String.fromFile(file)

        self.typeid = UInt32.fromFile(file)

        assert(self.typeid != 0)

        self.valuetype = self.psf.types.idMap[self.typeid]

        self.value = self.valuetype.getDataObj()
        
        self.value.deSerializeFile(file)

        print [ddef.datatypeid for ddef in self.valuetype.structdef.children]
        
        # Read possible property objects that belongs to the type by peeking ahead
        while True:
            oldpos = file.tell()
            try:
                prop = readChunk(self.psf, file, expectedclasses=PropertyClasses)
                self.properties.append(prop)
            except ValueError:
                file.seek(oldpos)
                break
        
    def toPSFasc(self, prec=None):
        r = self.name.toPSFasc(prec) + " " + self.valuetype.name.toPSFasc(prec) + " " + self.value.toPSFasc(prec)
        if len(self.properties)>0:
            r+=" PROP(\n"
            r+="\n".join([prop.toPSFasc(prec) for prop in self.properties])
            r+="\n)"
        return r

    def __repr__(self):
        return self.__class__.__name__+"("+str({"name":self.name, "id":"0x%x"%self.id, "typeid":"0x%x"%self.typeid,
                                                "properties":self.properties,"value":self.value})+")"



class SweepValue(Chunk):
    """Class representing waveform data"""
    type = 16
    def __init__(self, psf, type=None):
        Chunk.__init__(self, psf, type)
        self.id = None
        self.linktypeid = UInt32()
        self.datatypeid = UInt32()
        self.paramtype = None
        self.paramvalue = None
        self.children = []
        self.properties = []

    def deSerializeFile(self, file, n=None):
        pass

    def getSweepParamValues(self):
        pass
    

    def __len__(self):
        return len(self.children)
    
    def __repr__(self):
        return self.__class__.__name__ + "(" + str(self.paramtype.name) + "=" + str(self.paramvalue) +","+ \
               "children="+str(self.children) +")\n"

class SweepValueSimple(SweepValue):
    def deSerializeFile(self, file, n=None):
        Chunk.deSerializeFile(self, file)

        self.paramtypeid = UInt32.fromFile(file)

        self.paramtype = self.psf.sweeps.getSweep(self.paramtypeid)
        self.paramvalue = self.paramtype.getDataObj()
        self.paramvalue.deSerializeFile(file)

        for datatype in self.psf.traces.children:
            datatypeid = UInt32.fromFile(file)

            if datatypeid in (17,16):
                valuetypeid = UInt32.fromFile(file)

                if valuetypeid != datatype.id:
                    raise Exception("Unexpected trace trace type id %d should be %d"%(valuetypeid, datatype.id))

                value = datatype.getDataObj()
                value.deSerializeFile(file)
                self.children.append(value)
            else:
                raise Exception("Datatypeid unknown 0x%x"%self.datatypeid)

        return 1

    def getSweepParamValues(self):
        return [self.paramvalue.getValue()]

    def toPSFasc(self, prec=None):
        r=self.paramtype.name.toPSFasc(prec) + " " +self.paramvalue.toPSFasc(prec)+"\n"
        r+="\n".join([valuetype.name.toPSFasc(prec) + " "  + value.toPSFasc(prec) \
                      for valuetype, value in zip(self.psf.traces.children, self.children)])
        return r
            

class SweepValueWindowed(SweepValue):
    def deSerializeFile(self, file, n=None):
        bufferstart = file.tell()

        Chunk.deSerializeFile(self, file)

        self.paramtypeid = UInt32.fromFile(file)

        assert(len(self.psf.sweeps.children) == 1)
        self.paramtype=self.psf.sweeps.children[0]
        self.paramvalue = []

        # Get sweep parameter values
        paramvaluesize = self.paramtype.getDataSize()
        windowsize = self.psf.header.properties['PSF window size'].value
        leftinwindow = (file.tell()//windowsize + 1)*windowsize - file.tell()

        windowlen = leftinwindow//paramvaluesize; 

        if n > windowlen:
            n = windowlen

        for j in xrange(n):
            paramvalue = self.paramtype.getDataObj()
            paramvalue.deSerializeFile(file)
            if j < n:
                self.paramvalue.append(paramvalue)

        # Get trace values
        for trace in self.psf.traces.children:
            value = trace.getDataObj()
            value.deSerializeFile(file, count=n,
                                  windowsize=self.psf.header.properties['PSF window size'].value)
            self.children.append(value)

        # Skip trailing padding bytes
        padsize = int((self.psf.header.properties['PSF buffer size'] - (file.tell()-bufferstart))% \
                  self.psf.header.properties['PSF buffer size'])
        file.seek(padsize, 1)

        return n

    def getSweepParamValues(self):
        return [v.getValue() for v in self.paramvalue]

    def toPSFasc(self, prec=None):
        r=''
        for i, paramvalue in enumerate(self.paramvalue):
            r+=self.paramtype.name.toPSFasc(prec) + " " + paramvalue.toPSFasc(prec) + "\n"
            r+="\n".join([trace.name.toPSFasc(prec) + " " + value.toPSFasc(prec=prec, index=i) \
                          for trace,value in zip(self.psf.traces.children, self.children)])
            if i < len(self.paramvalue)-1:
                r+="\n"
        return r

class GroupData(PSFData):
    def __init__(self, groupdef):
        PSFData.__init__(self)
        self.groupdef = groupdef
        self.children = []
    def deSerializeFile(self, file, count=None, windowsize=None):
        for element in self.groupdef.children:
            if count==None:
                value = element.getDataObj()
                value.deSerializeFile(file)
                self.children.append(value)
            else:
                valuearray=[]

                # If a window is used in the PSF file, the entire window is stored
                # and the data is aligned to the end of the window. So we need
                # to skip window size - data size
                file.seek(int(windowsize - count*element.getDataSize()), 1)

                for i in xrange(0,count):
                    value = element.getDataObj()
                    value.deSerializeFile(file)
                    valuearray.append(value)

                self.children.append(valuearray)

    def toPSFasc(self, prec=None, index=None):
        if index != None:
            return "\n".join([v[index].toPSFasc(prec) for v in self.children])
        else:
            return "\n".join([v.toPSFasc(prec) for v in self.children])

    def getSize(self):
        return self.groupdef.getDataSize()

    def __repr__(self):
        return "GroupData" + "\n" + "\n".join([indent(s) for s in map(repr,self.children)]) + "\n"
            
class GroupDef(Chunk):
    type=17
    """Class representing group of traces"""
    def __init__(self, psf):
        Chunk.__init__(self, psf)
        self.children=[]
        self.datasize=None

    def getDataObj(self):
        return GroupData(self)
    def deSerializeFile(self, file):
        Chunk.deSerializeFile(self, file)
        self.id = UInt32.fromFile(file)
        self.name = String.fromFile(file)
        self.nchildren = UInt32.fromFile(file)

        # Read children
        self.children = []
        self.datasize = 0
        for i in range(0, self.nchildren):
            child = DataTypeRef(self.psf)
            child.deSerializeFile(file)
            self.children.append(child)
            self.datasize += child.getDataSize()

    def getNameIndex(self):
        return dict([(v.name, (i,)) for i,v in enumerate(self.children)])

    def toPSFasc(self, prec=None):
        s=self.name.toPSFasc(prec) + " GROUP %d\n"%len(self.children)
        s+="\n".join([child.toPSFasc(prec) for child in self.children])
        return s

    def getDataSize(self):
        return self.datasize

    def getNames(self):
        return [str(child.name) for child in self.children]

    def __repr__(self):
        return "0x%x"%self.fileoffset +":" + self.__class__.__name__+ "(id=0x%x"%self.id+", nchildren=%d"%self.nchildren+")\n" + "\n".join([indent(s) for s in map(str,self.children)]) + "\n"


class UnknownChunk(Exception):
    def __init__(self, chunktype):
        self.type = chunktype
    def __str__(self):
        return "Unknown chunk of type: %d"%self.type

class InvalidChunk(Exception):
    def __init__(self, chunk):
        self.chunk = chunk
    def __str__(self):
        return "Invalid %s"%(self.chunk.__class__.__name__)

class IncorrectChunk(Exception):
    def __init__(self, type, expectedtype):
        self.type = type
        self.expectedtype = expectedtype
    def __str__(self):
        return "Incorrect chunk type %d (should be %d)"%(self.type, self.expectedtype)

class LastValue(Exception):
    pass
    
def readChunk(psf, file, expectedclasses=None):
    type = UInt32.fromFile(file)
    file.seek(-4, 1) # Rewind one word since the type will be read again by the deSerializeFile function

    if expectedclasses:
        if not type in [cls.type for cls in expectedclasses]:
            raise ValueError("Unexpected type %d, not in "%type + str([cls.type for cls in expectedclasses]))
        for cls in expectedclasses:
            if type == cls.type:
                chunk = cls(psf)
    else:
        raise Exception("Use expectedclasses!")
        if type == 21:
            chunk = Section(psf)
        elif type == 20:
            chunk = ZeroPad(psf)
        elif type == 22:
            chunk = Container22(psf, type, n=n)
        elif type == 33:
            chunk = PropertyString(psf)
        elif type == 34:
            chunk = PropertyUInt(psf)
        elif type == 35:
            chunk = PropertyFloat64(psf)
        elif type == 16:
            chunk = DataTypeDef(psf,type)
        elif type == 17:
            chunk = GroupDef(psf)
        elif type == 19:
            chunk = HashTable(psf, n=n)
        elif type in (1,2,3,4):
            file.seek(4,1)
            return None
        else:
            warning("Unknown chunk %d"%type)
            raise UnknownChunk(type)

    chunk.deSerializeFile(file)
    
    return chunk

class PSFReader(object):
    def __init__(self, filename=None, asc=None):
        self.header = None
        self.types = TypeSection(self)
        self.sweeps = None
        self.traces = None
        self.lastid = 0x10000000;
        self.verbose = False
        self.filename = filename
        self.file = None
        self.values = None
        self.asc = asc
        
    def open(self):
        """Open a PSF file and read its headers.

        Example:
        Trying to open a valid psf file
        >>> psf=PSFReader('./test/psf/srcSweep')
        >>> psf.open()
        """
        
        if self.asc == None:
            self.asc = psfasc.is_psfasc(self.filename)

        if not self.asc:
            self.file = open(self.filename, "rb")
            
            if self.validate():
                self.deSerializeFile(self.file)
            else:
                raise PSFInvalid("Invalid PSF file")
        else:
            newpsfobj = psfasc.parse("psfasc", open(self.filename).read())
            self.header = newpsfobj.header                
            self.types = newpsfobj.types
            self.sweeps = newpsfobj.sweeps
            self.traces = newpsfobj.traces
            self.values = newpsfobj.values
            self.lastid = newpsfobj.lastid
            self.verbose = newpsfobj.verbose
            
    def validate(self):
        """Check if the PSF file is valid.
        Returns True if valid, False otherwise

        >>> psf=PSFReader('./test/psf/srcSweep')
        >>> psf.validate()
        True
        >>> psf=PSFReader('./test/psfasc/srcSweep.asc')
        >>> psf.validate()
        False
        """
        if self.file == None:
            file = open(self.filename, "rb")
        else:
            file = self.file
            
        # Read Clarissa signature
        file.seek(-4-8,2)
        clarissa = file.read(8)
        return clarissa == "Clarissa"

    def getNSweepPoints(self):
        """Returns number of sweeps. 0 if not swept.

        >>> psf=PSFReader('./test/psf/srcSweep')
        >>> psf.open()
        >>> psf.getNSweepPoints()
        4
        """
        if self.file == None:
            ValueError("Please open the PSF file first")
        return self.header.properties['PSF sweep points']

    def getNSweeps(self):
        """Returns the number of nested sweeps

        >>> psf=PSFReader('./test/psf/srcSweep')
        >>> psf.open()
        >>> psf.getNSweeps()
        1
        """
        if self.file == None:
            ValueError("Please open the PSF file first")
        return self.header.properties['PSF sweeps']

    def __len__(self):
        return len(self.values)

    def getValueNames(self):
        """Returns a tuple of the names of the traces

        >>> psf=PSFReader('./test/psf/srcSweep')
        >>> psf.getValueNames()
        >>> psf.open()
        >>> psf.getValueNames()
        ('VOUT', 'VIN', 'R0')

        >>> psf=PSFReader('./test/resultdirs/simple/opBegin')
        >>> psf.open()
        >>> psf.getValueNames()
        ('R0', 'V1', 'V0', 'E0', 'VIN', 'NET9', 'VOUT')

        """
        if self.values:
            return self.values.getValueNames()
    
    def getSweepParamValues(self, dim=0):
        """Returns a numpy.array of sweep parameter values for sweep dimension dim.

        >>> psf=PSFReader('./test/psf/srcSweep')
        >>> psf.open()
        >>> psf.getSweepParamValues(0)
        array([ 1.,  2.,  3.,  4.])

        windowed result
        >>> psf=PSFReader('./test/psf/timeSweep')
        >>> psf.open()
        >>> psf.getSweepParamValues(0)[:3]
        array([  0.00000000e+00,   2.00000000e-11,   5.33333333e-11])

        """
        return numpy.array(self.values.getSweepParamValues())

    def getValuePropertiesByName(self, name):
        """Returns the properties associated with value
        
        >>> psf=PSFReader('./test/psf/opBegin')
        >>> psf.open()
        >>> psf.getValuePropertiesByName("XIRXRFMIXTRIM0.XM1PDAC1.XMN.MAIN")["Region"]
        'subthreshold'

        """
        return self.values.getValuePropertiesByName(name)

    def getValuesByName(self, name):
        """Returns a numpy.array of trace values for swept results and a scalar for non swept.

        Example:
        swept psf file
        >>> psf=PSFReader('./test/psf/srcSweep')
        >>> psf.open()
        >>> psf.getValuesByName("VOUT")
        array([-6., -4., -2.,  0.])
        >>> psf.getValuesByName("VIN")
        array([ 1.,  2.,  3.,  4.])

        swept psf with complex numbers
        >>> psf=PSFReader('./test/psf/frequencySweep')
        >>> psf.open()
        >>> res = psf.getValuesByName("ANT_CM")
        >>> len(res)
        123
        >>> res[:3]
        array([ 0.6+0.j,  0. +0.j,  0. +0.j])
        
        swept windowed psf file 
        >>> psf=PSFReader('./test/psf/timeSweep')
        >>> psf.open()
        >>> psf.getValuesByName("INP")[0:3]
        array([ 0.6       ,  0.62486899,  0.66211478])

        non-swept psf file
        >>> psf=PSFReader('./test/psf/dcOpInfo.info')
        >>> psf.open()
        >>> psf.getValuesByName("IREG21U_0.MP5.b1")['betadc']
        4.7957014499434756

        swept psf file withouth groups
        >>> psf=PSFReader('./test/resultdirs/parsweep2/C=1e-12,R=1e-12/psf/ac.ac')
        >>> psf.open()
        >>> psf.getValuesByName("net3")
        array([ 0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,
                0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,
                0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,
                0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,
                0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,
                0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,
                0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,
                0.+0.j,  0.+0.j])

        """
        return self.values.getValueByName(name)
        
    def nTraces(self):
        """Returns number of traces

        >>> psf=PSFReader('./test/psf/srcSweep')
        >>> psf.open()
        >>> psf.nTraces()
        3
        """
        if self.file == None:
            ValueError("Please open the PSF file first")
        return self.header.properties['PSF traces']

    def allocId(self):
        self.lastid+=1
        return self.lastid-1

    def info(self):
        s="Number of sweeps: %d\n"%self.getNSweeps()
        if self.getNSweeps() > 0:
            s+="Number of sweep points: %d\n"%self.getNSweepPoints()
        s+="Number of traces: %d"%self.nTraces()
        return s

    def updateHeader(self):
        if self.sweeps:
            sweeps = len(self.sweeps.children)
        else:
            sweeps=0
        self.header.addProperty(PropertyUInt("PSF sweeps", sweeps))

    def deSerializeFile(self, file):
        # Find filesize
        file.seek(0,2)
        filesize = file.tell()
        
        # Last word contains the size of the data
        file.seek(-4,2)
        datasize = UInt32.fromFile(file).value
        if self.verbose:
            print "Total data size: ",datasize

        # Read Clarissa signature
        file.seek(-4-8,2)
        clarissa = file.read(8)
        if not clarissa == "Clarissa":
            raise ValueError("Clarissa signature not found")

        # Read section index table
        sectionoffsets = {}
        file.seek(-4-8-8,2)
        pos = file.tell()

        sectionnums = []
        while file.tell() >= datasize:
            sectionnum = UInt32.fromFile(file)
            sectionnums.insert(0,sectionnum.value)
            offset = UInt32.fromFile(file)
            sectionoffsets[sectionnum] = offset
            pos -= 8
            file.seek(pos)

        offsets = [sectionoffsets[secnum] for secnum in sectionnums]
        sizes = map(operator.sub, offsets[1:]+[datasize], offsets)
        sectionsizes = dict(zip(sectionnums, sizes))

        if self.verbose:
            print sectionoffsets, sectionsizes
        
        file.seek(0)

        self.unk1 = UInt32.fromFile(file)
        if self.verbose:
            print "First word: 0x%x"%self.unk1

        # Load headers
        file.seek(int(sectionoffsets[0]))
        self.header = HeaderSection(self)
        self.header.deSerializeFile(file)
        if self.verbose:
            print "HEADER"
            print self.header
        

        if sectionoffsets.has_key(1):
            file.seek(int(sectionoffsets[1]))
            self.types.deSerializeFile(file)

            if self.verbose:
                print "TYPE"
                print self.types

        if sectionoffsets.has_key(2):
            file.seek(int(sectionoffsets[2]))
            self.sweeps = SweepSection(self)
            self.sweeps.deSerializeFile(file)

            if self.verbose:
                print "SWEEPS"
                print self.sweeps

        if sectionoffsets.has_key(3):
            file.seek(int(sectionoffsets[3]))
            self.traces = TraceSection(self)
            self.traces.deSerializeFile(file)

        if sectionoffsets.has_key(4):
            file.seek(int(sectionoffsets[4]))
            # Load data
            if self.sweeps:
                self.values = ValuesSectionSweep(self)
            else:
                self.values = ValuesSectionNonSweep(self)
            self.values.deSerializeFile(file)

    def printme(self):
        print "HEADER"
        print self.header
        print "TYPES"
        print self.types
        if self.sweeps:
            print "SWEEP"
            print self.sweeps
        if self.traces:
            print "TRACE"
            print self.traces
        print "VALUES"
        print self.values

    def toPSFasc(self, prec=None):
        """Export to PSF ascii"""
        sections = [self.header.toPSFasc(prec), self.types.toPSFasc(prec)]
        if self.sweeps:
            sections.append(self.sweeps.toPSFasc(prec))
        if self.traces:
            sections.append(self.traces.toPSFasc(prec))
        if self.values:
            sections.append(self.values.toPSFasc(prec))
        r="\n".join(sections) + "\n"
        r+="END\n"
        return r

    def __repr__(self):
        return "\n".join(map(str, (self.header, self.types, self.sweeps, self.traces, self.values)))


if __name__ == "__main__":
    import doctest
    doctest.testmod()


