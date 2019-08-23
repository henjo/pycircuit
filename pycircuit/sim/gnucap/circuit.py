# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import pycircuit.sim

class Circuit(pycircuit.sim.Circuit):
    def __init__(self, netlist=None):
        self.netlist = netlist

        self.instances = {}

    def __eq__(self, a): return str(self) == str(a)

    def __setitem__(self, key, value):
        if not isinstance(value, Instance):
            raise ValueError('instance must be Instance object')

        if value.name_prefix is not None and not key.startswith(value.name_prefix):
            raise ValueError('instance name must start with %s'%
                             value.name_prefix)

        self.instances[key] = value

    def __str__(self):
        if self.netlist is None:
            res = ''
            for k, v in self.instances.items():
                res += k + ' ' + str(v) + '\n'

            return res
        else:
            return self.netlist


class Instance(object):
    def __init__(self, *args, **kvargs):
        self.args = args
        self.kvargs = kvargs

    def __str__(self):
        return ' '.join([str(arg) for arg in self.args] +
                        [str(k) + ' ' + str(v) for k,v in self.kvargs.items()])

    name_prefix = None

class VS(Instance):
    name_prefix = 'V'
class R(Instance):
    name_prefix = 'R'
class C(Instance):
    name_prefix = 'C'

