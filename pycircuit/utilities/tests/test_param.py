# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import unittest
from nose.tools import *
from numpy.testing import assert_array_almost_equal, assert_array_equal

from pycircuit.utilities.param import *

def test_parameter():
    p = Parameter('name', desc='Description', unit='V', default=1)

    assert_equal(str(p), 'name')
    assert_equal(repr(p),
                 "Parameter('name', desc='Description', unit='V', default=1)")

    p2 = p.copy()
    
    assert_equal(p, p2)

def test_subclass_parameter():
    p = Parameter('name', desc='Description', unit='V', default=1)

    class MyParameter(Parameter):
        pass

    p3 = MyParameter('name', desc='Description', unit='V', default=1)
    
    assert_not_equal(p,p3)

    assert_equal(repr(p3),
                 "MyParameter('name', desc='Description', unit='V', default=1)")

class ParameterTest(unittest.TestCase):
    """Test Parameter class"""

    def testInstance(self):
        """Testing instantiating Parameter objects"""
        param1 = Parameter("gm", "Transconductance", "A/V", default=1e-3)
        param2 = Parameter(name="gm", desc="Transconductance", unit="A/V", default=1e-3)
        for param in param1,param2:
            self.assertEqual(param.name, 'gm')
            self.assertAlmostEqual(param.default, 1e-3)
            self.assertEqual(param.desc, 'Transconductance')
            self.assertEqual(param.unit, 'A/V')

    def testIncorrectInstantiation(self):
        def testIncorrect():
            param = Parameter(apa=1)

        self.failUnlessRaises(TypeError, testIncorrect)
        
class ParameterDictTest(unittest.TestCase):
    def testInstantiate(self):
        paramdict = ParameterDict()
        gm = Parameter(name="gm", desc="Transconductance", unit="A/V", 
                       default=1e-3)
        gds = Parameter(name="gds", desc="Output conductance", unit="A/V", 
                        default=1e-6)
        paramdict = ParameterDict(gm,gds,
                                  gm=2e-3
                                  )
        self.assertEqual(paramdict.parameters, [gm,gds])
        self.assertEqual(paramdict.gm, 2e-3)
        self.assertEqual(paramdict.gds, 1e-6)

    def testAppendParameter(self):
        """Test appending a parameter"""
        paramdict = ParameterDict()
        self.assertEqual(len(paramdict), 0)
        gmparam = Parameter(name="gm", desc="Transconductance", unit="A/V")
        paramdict.append(gmparam)
        self.assertEqual(len(paramdict), 1)
        self.assertEqual(paramdict['gm'], gmparam)

    def testGetItem(self):
        """Test getting parameter"""
        paramdict = ParameterDict()
        gmparam = Parameter(name="gm", desc="Transconductance", unit="A/V", 
                            default=2.0e-6)
        paramdict.append(gmparam)
        self.assertEqual(paramdict['gm'], gmparam)

    def testGetSetValue(self):
        paramdict = ParameterDict()
        gmparam = Parameter(name="gm", desc="Transconductance", unit="A/V", 
                            default=2.0e-6)
        paramdict.append(gmparam)

        self.assertEqual(paramdict.get('gm'), 2e-6)

        self.assertEqual(paramdict.get(gmparam), 2e-6)

        self.assertEqual(paramdict.gm, 2e-6)

        paramdict.gm = 3e-6
        self.assertEqual(paramdict.gm, 3e-6)

    def testSet(self):
        paramdict = ParameterDict()
        paramdict.append(Parameter(name="gm", desc="Transconductance", 
                                   unit="A/V", default=2.0e-6))
        paramdict.set(gm=4e-6)
        self.assertEqual(paramdict.gm, 4e-6)
        def test():
            paramdict.set(gds=3e-6)
        self.failUnlessRaises(KeyError, test)
        
    def testContains(self):
        paramdict = ParameterDict()
        gmparam = Parameter(name="gm", desc="Transconductance", unit="A/V", 
                            default=2.0e-6)
        paramdict.append(gmparam)
        self.assertTrue('gm' in paramdict)
        self.assertTrue(gmparam in paramdict)
        self.assertFalse('gds' in paramdict)
        self.assertFalse(Parameter(name='gds') in paramdict)

    def testGetParameters(self):
        paramdict = ParameterDict()
        gmparam = Parameter(name="gm", desc="Transconductance", unit="A/V", 
                            default=2.0e-6)
        gdsparam = Parameter(name="gds")
        paramdict.append(gmparam)
        self.assertEqual(paramdict.parameters, [gmparam])
        paramdict.append(gdsparam)
        self.assertEqual(paramdict.parameters, [gmparam, gdsparam])

    def testCopy(self):
        param1 = Parameter(name="vth0", desc="Threshold voltage", unit="V", 
                           default=0.3)
        paramdict = ParameterDict(param1)

        paramdict2 = paramdict.copy()
        paramdict2.vth0 = 0.4
        paramdict3 = paramdict2.copy(vth0=0.5)
        self.assertEqual(paramdict2.vth0, 0.4)
        self.assertEqual(paramdict.vth0, 0.3)
        self.assertEqual(paramdict3.vth0, 0.5)

def test_parameter_dict_symbolic():
    
    a,b = Parameter('a'), Parameter('b')

    gmparam = Parameter(name="gm", desc="Transconductance", unit="A/V", 
                        default=2.0e-6)
    gdsparam = Parameter(name="gds")

    pdict_expr = ParameterDict(gmparam,gdsparam)
    pdict_ab = ParameterDict(a,b)

    pdict_expr.gm = "2*a + b"
    pdict_expr.gds = "-3*a + 2*b"

    pdict_ab.a = 3
    pdict_ab.b = 4

    pdict_values = pdict_expr.eval_expressions((pdict_ab,))

    assert_equal(pdict_values.gm, 10)
    assert_equal(pdict_values.gds, -1)
    
    class Variable(Parameter):
        pass
    
    vara = Variable('a')

    pdict_vars = ParameterDict(vara)

    pdict_vars.a = 10

    pdict_expr.gds = "-3*a + 2*b"

    pdict_values = pdict_expr.eval_expressions((pdict_ab, pdict_vars))

    assert_equal(pdict_values.gm, 24)
    assert_equal(pdict_values.gds, -22)

    ## Try to use parameter not defined in pdict_ab
    pdict_expr.gm = "2*a + b + c"

    assert_raises(EvalError, lambda: pdict_expr.eval_expressions((pdict_ab,)))
    
def test_update_parameterdict():
    paramdict1 = ParameterDict()
    paramdict2 = ParameterDict()

    gmparam = Parameter(name="gm", desc="Transconductance", unit="A/V")
    paramdict1.append(gmparam)
    paramdict2.append(gmparam)

    paramdict1.gm = 10
    paramdict2.gm = 30

    paramdict1.update_values(paramdict2)
    
    assert_equal(paramdict1.gm, 30)
    assert_equal(paramdict2.gm, 30)
    
