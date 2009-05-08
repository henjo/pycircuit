# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from pycircuit.sim import Variable
from pycircuit.utilities.param import Parameter, ParameterDict
from pycircuit.utilities.misc import indent, inplace_add_selected, \
    inplace_add_selected_2d, create_index_vectors
from constants import *
from copy import copy
import numpy as np
import types
import numeric

default_toolkit = numeric

class Node(object):
    """A Node object represents a point in an electric circuit"""
    def __init__(self, name=None, isglobal = False):
        if name.endswith('!'):
            name = name[:-1]
            isglobal = True
            
        self.name = name
        self.isglobal = isglobal

    def __hash__(self): return hash(self.name)

    def __eq__(self, a): 
        try:
            return self.name == a.name
        except:
            return False        

    @property
    def V(self):
        return Quantity('V', self)

    def __str__(self):
        name = self.name
        if self.isglobal:
            name = name + '!'
        return name 

    def __repr__(self):
        if self.isglobal:
            return self.__class__.__name__ + '(' + repr(self.name) + ', ' \
                'isglobal=True)'
        else:
            return self.__class__.__name__ + '(' + repr(self.name) + ')'            

class Branch(object):
    """A branch connects two nodes.
    
    A branch is used in modified nodal analysis to describe components that 
    defines a voltage between two nodes as a function of current flowing 
    between these nodes. Examples are voltage sources and inductors.
    Positive current through a branch is defined as a current flowing from 
    plus to minus.a
    
    """
    def __init__(self, plus, minus, name=None):
        """Initiate a branch

        Arguments:
        plus -- Node object connected to the postive terminal of the branch
        minus -- Node object connected to the negative terminal of the branch

        Keyword arguments:
        name -- branch name

        """

        self.plus = plus
        self.minus = minus
        self.name = name

    def __hash__(self): return hash(self.plus) ^ hash(self.minus)

    def __eq__(self, a): 
        try:
            return self.plus == a.plus and self.minus == a.minus
        except:
            return False        

    @property
    def V(self):
        return Quantity('V', self)

    @property
    def I(self):
        return Quantity('I', self)

    def __repr__(self):
        return 'Branch('+repr(self.plus)+','+repr(self.minus)+')'

### Default reference node
gnd = Node("gnd", isglobal=True)

defaultepar = ParameterDict(
    Parameter("T", "Temperature", unit="K", default = 300))

class Circuit(object):
    """Basic circuit class 

    The circuit class models electric circuits but could be used for
    any conservative system. A circuit object contains a list of nodes and
    branches which are associated with node voltages and branch currents in the
    modelled system. Note that the values of these quantities are
    stored in separate analysis classes, never inside a circuit object.

    The nodes are connected to the outside through terminals. When the circuit
    is instanciated, the outside nodes are passed to the object via the 
    terminals.

    **Attributes**
        *nodes*
          A list that contains Node objects. The first k nodes are terminal 
          nodes where k is the number of terminals.

        *branches*
          list of Branch objects. The solver will solve for the 
          currents through the branches.

        *terminals*
          list of terminal names

        *instparams*
          A list of valid instance parameters (Parameter objects)

        *mpar*
          A class variable with a ParameterDict containing model specific 
          parameters

        *ipar*
          A ParameterDict containing instance specific parameters

        *nodenames*
          A dictionary that maps a local node name to the node object in
          nodes. If the node is connnected to superior hierarchy levels
          through a terminal the terminal name must be the same as the
          local node name

        *terminalhook*
          Temporary storage of information about what nodes in the superior
          hierarchy level the terminals should be connected to during 
          instantiation. It is a dictionary where the keys are terminal names
          and the values are node objects. The value is None when it is not 
          used. The only reason for this attribute is to allow for bottom-up 
          instantiations like: cir1['I1'] = R('n1', gnd)

        *linear* 
          A boolean value that is true if i(x) and q(x) are linear 
          functions

    """

    
    nodes = []
    branches = []
    terminals = []
    mpar = ParameterDict()
    instparams = []
    linear = True
    
    def __init__(self, *args, **kvargs):
        if 'toolkit' in kvargs:
            self.toolkit = kvargs['toolkit']
            del kvargs['toolkit']
        else:
            self.toolkit = default_toolkit

        self.nodenames = {}
        self.ipar = ParameterDict(*self.instparams, **kvargs)
        
        ## Set up ipar expressions
        self.ipar_expressions = ParameterDict(*self.instparams, **kvargs)
        for par in self.instparams:
            if par.name not in kvargs:
                setattr(self.ipar_expressions, par.name, None)

        ## Add terminal nodes
        for terminal in self.terminals:
            self.append_node(Node(terminal))

        ## Set temporary terminal mapping information for use by instantiation
        ## method in higher hierarchy
        self.terminalhook = dict(zip(self.terminals, args))

        ## Subscribe to updates of instance parameters
        if hasattr(self, 'update'):
            self.ipar.attach(self)
            self.update(self.ipar)

    def __eq__(self, a):
        return self.__class__ == a.__class__ and \
            self.nodes == a.nodes and \
            self.nodenames == a.nodenames and self.branches == a.branches and \
            self.ipar == a.ipar
        
    def __copy__(self):
        newc = self.__class__()
        newc.nodes = copy(self.nodes)    
        newc.nodenames = copy(self.nodenames)    
        newc.branches = copy(self.branches)    
        newc.ipar = copy(self.ipar)
        return newc

    def add_nodes(self, *names):
        """Create internal nodes in the circuit and return the new nodes

        >>> c = Circuit()
        >>> n1, n2 = c.add_nodes("n1", "n2")
        >>> c.nodes
        [Node('n1'), Node('n2')]
        >>> 'n1' in c.nodenames and 'n2' in c.nodenames
        True
        
        """
        newnodes = []
        for name in names:
            newnodes.append(Node(name))
            self.append_node(newnodes[-1])

        return tuple(newnodes)

    def add_node(self, name):
        """Create and internal node in the circuit and return the new node"""
        return self.add_nodes(name)[0]

    def append_node(self, node):
        """Append node object to circuit"""
        ## Make a copy of node list so the class is unchanged
        if self.__class__.nodes is self.nodes:
            self.nodes = list(self.nodes)

        if node not in self.nodes:
            self.nodes.append(node)
        self.nodenames[node.name] = node

    def append_branches(self, *branches):
        """Append node object to circuit"""
        ## Make a copy of node list so the class is unchanged
        if self.__class__.branches is self.branches:
            self.branches = list(self.branches)

        self.branches.extend(branches)

    def get_terminal_branch(self, terminalname):
        """Find the branch that is connected to the given terminal

        If no branch is found or if there are more branches than one, None is 
        returned

        Returns
        -------
        Tuple of Branch object and an integer indicating if terminal is 
        connected to the positive or negative side of the branch. 
        1 == positive side, -1 == negative side
        
        >>> from elements import *
        >>> net1 = Node('net1')
        >>> net2 = Node('net2')
        >>> VS = VS(net1, net2)
        >>> VS.get_terminal_branch("minus")
        (Branch(Node('plus'),Node('minus')), -1)
        
        """
        plusbranches = [] ## Branches with positive side connected to the
                          ## terminal
        minusbranches = [] ## Branches with negative side connected to the
                           ## terminal
        for branch in self.branches:
            if branch.plus == self.nodenames[terminalname]:
                plusbranches.append(branch)
            elif branch.minus == self.nodenames[terminalname]:
                minusbranches.append(branch)

        if len(plusbranches + minusbranches) != 1:
            return None
        elif len(plusbranches) == 1:
            return plusbranches[0], 1
        elif len(minusbranches) == 1:
            return minusbranches[0], -1            

    def get_node_index(self, node):
        """Get row in the x vector of a node instance"""

        if not isinstance(node, Node):
            node = Node(str(node))

        return self.nodes.index(node)

    def get_branch_index(self, branch):
        """Get row in the x vector of a branch instance"""
        if branch in self.branches:
            return len(self.nodes) + self.branches.index(branch)
        else:
            raise ValueError('Branch %s is not present in circuit (%s)'%
                             (str(branch), str(self.branches)))

    def get_node(self, name):
        """Find a node by name.
        
        >>> c = Circuit()
        >>> n1 = c.add_node("n1")
        >>> c.get_node('n1')
        Node('n1')
        
        """
        return self.nodenames[name]

    def get_node_name(self, node):
        """Find the name of a node
        
        >>> c = Circuit()
        >>> n1 = c.add_node("n1")
        >>> c.get_node_name(n1)
        'n1'

        """
        for k, v in self.nodenames.items():
            if v == node:
                return k
        
    def connect_terminals(self, **kvargs):
        """Connect nodes to terminals by using keyword arguments

        """
        for terminal, node in kvargs.items():
            ## Sanity check
            if type(terminal) is not types.StringType:
                raise Exception("%s should be string"%str(terminal))
            if terminal not in self.terminals:
                raise ValueError('terminal '+str(terminal)+' is not defined')
            
            if not isinstance(node, Node):
                node = Node(str(node))
            
            if terminal in self.nodenames:
                oldterminalnode = self.nodenames[terminal]
                if oldterminalnode != node:
                    self.nodes.remove(self.nodenames[terminal])
            
            if node not in self.nodes:
                self.nodes.insert(self._nterminalnodes, node)
            
            self.nodenames[terminal] = node            
            
    def save_current(self, terminal):
        """Returns a circuit where a current probe is added at a terminal
        
        >>> from elements import *
        >>> cir = R(Node('n1'), gnd, r=1e3)
        >>> newcir = cir.save_current('plus')
        >>> newcir.G(np.zeros(4))
        array([[ 0.   ,  0.   ,  0.   ,  1.   ],
               [ 0.   ,  0.001, -0.001,  0.   ],
               [ 0.   , -0.001,  0.001, -1.   ],
               [ 1.   ,  0.   , -1.   ,  0.   ]])
        """
        
        if self.get_terminal_branch(terminal) == None:
            return ProbeWrapper(self, terminals = (terminal,))
        else:
            return self            

    @property
    def n(self):
        """Return size of x vector"""
        return len(self.nodes) + len(self.branches)

    def terminal_nodes(self):
        """Return a list of all terminal nodes"""
        return self.nodes[0:self._nterminalnodes]

    def non_terminal_nodes(self, instancename = None):
        """Return a list of all non-terminal nodes. 

        If the instancename is set, the local nodes
        will have a instancename<dot> prefix added to the node name

        """
        if instancename == None:
            return self.nodes[self._nterminalnodes:]
        else:
            result = []
            for node in self.nodes[len(self.terminals):]:
                if node.isglobal:
                    result.append(node)
                else:
                    result.append(Node(instancename + '.' + node.name))
            return result

    def G(self, x, epar=defaultepar):
        """Calculate the G (trans)conductance matrix given the x-vector"""
        return self.toolkit.zeros((self.n, self.n))

    def C(self, x, epar=defaultepar):
        """Calculate the C (transcapacitance) matrix given the x-vector"""
        return self.toolkit.zeros((self.n, self.n))

    def u(self, t=0.0, epar=defaultepar, analysis=None):
        """Calculate the u column-vector of the circuit at time t

        Arguments
        ---------

        epar -- ParameterDict with environment parameters such as temperature
        analysis -- This argument gives the possibility to have analysis 
                    dependent sources.
                    for normal time dependent and dc sources this argument 
                    should be None
        
        """
        return self.toolkit.zeros(self.n)

    def i(self, x, epar=defaultepar):
        """Calculate the i vector as a function of the x-vector

        For linear circuits i(x(t)) = G*x
        """
        return self.toolkit.dot(self.G(x), x)

    def q(self, x, epar=defaultepar):
        """Calculate the q vector as a function of the x-vector

        For linear circuits q(x(t)) = C*x
        """
        return self.toolkit.dot(self.C(x), x)

    def CY(self, x, w, epar=defaultepar):
        """Calculate the noise sources correlation matrix

        Arguments
        ---------
        x -- (numpy array) the state vector
        w -- Angular frequency
        epar -- (ParameterDict) Environment parameters

        """
        return self.toolkit.zeros((self.n, self.n))

    def next_event(self, t):
        """Returns the time of the next event given the current time t"""
        return inf
    
    def name_state_vector(self, x, analysis=''):
        """Return a dictionary of the x-vector keyed by node and branch names

        >>> from elements import *
        >>> c = SubCircuit()
        >>> n1 = c.add_node('net1')
        >>> c['is'] = IS(gnd, n1, i=1e-3)
        >>> c['R'] = R(n1, gnd, r=1e3)
        >>> c.name_state_vector(np.array([[1.0]]))
        {'net1': 1.0}

        >>> 

        """
        result = {}
        for xvalue, node in zip(x[:len(self.nodes)][0], self.nodes):
            result[self.get_node_name(node)] = xvalue

        nnodes = len(self.nodes)
        for i, xvalue, branch in enumerate(zip(x[nnodes:], self.branches)):
            result['i' + analysis + str(i) + ')'] = xvalue            

        return result

    def extract_v(self, x, nodep, noden=None, refnode=gnd, 
                  refnode_removed=False):
        """Extract voltage between nodep and noden from the given x-vector.

        If noden is not given the voltage is taken between nodep and refnode. 
        x-vectors with the reference node removed can be handled by setting 
        the refnode_removed to True.

        *x*
          x-vector

        *nodep*
          Node object or node reference in text format of positive node

        *noden*
          Node object or node reference in text format of negative node

        *refnode*
          reference node

        *refnode_removed*
          If set the refernce node is expected to be removed from the x-vector

        >>> from elements import *
        >>> c = SubCircuit()
        >>> n1, n2 = c.add_nodes('n1','n2')
        >>> c['R1'] = R(n1, n2, r=1e3)
        >>> c['R2'] = R(n2, gnd, r=1e3)
        >>> c.extract_v(np.array([1.0, 0.5, 0.0]), 'n1', 'n2')
        0.5
        >>> c.extract_v(np.array([1.0, 0.5, 0.0]), c.nodes[0])
        1.0
        >>> c.extract_v(np.array([1.0, 0.5]), c.nodes[0], refnode_removed = True)
        1.0
        
        """
        v = []
        for node in nodep, noden:
            if type(node) is types.StringType:
                node = self.get_node(node)
            elif node == None:
                node = refnode

            nodeindex = self.get_node_index(node)

            refnodeindex = self.get_node_index(refnode)

            if refnode_removed:
                if nodeindex > refnodeindex:
                    nodeindex -= 1

            if nodeindex != refnodeindex:
                v.append(x[nodeindex])
            else:
                v.append(0)

        return v[0] - v[1]

        

    def extract_i(self, x, branch_or_term, xdot = None,
                  refnode = gnd, refnode_removed = False,
                  t = 0,
                  linearized = False, xdcop = None):
        """Extract branch or terminal current from the given x-vector.

        *x* 
           x-vector

        *branch_or_term*
           Branch object or terminal name

        *xdot*
           dx/dt vector. this is needed if there is no branch defined at the
           terminal

        *refnode*
           reference node

        *refnode_removed*
           If set the refernce node is expected to be removed from the x-vector

        *t*
           Time when the sources are to be evaluated
        
        *linearized*
           Set to True if the AC current is wanted

        *xdcop*
           *xcdop* is the DC operation point x-vector if linearized == True

        >>> from elements import *
        >>> c = SubCircuit()
        >>> net1 = c.add_node('net1')
        >>> c['vs'] = VS(net1, gnd)
        >>> c.extract_i(np.array([1.0, 0, -1e-3]), 'vs.minus')
        0.001
        >>> c.extract_i(np.array([1.0, -1e-3]), 'vs.minus', refnode_removed = True)
        0.001
        
        """
        dot = self.toolkit.dot
        
        if type(branch_or_term) is types.StringType:
            ## Calculate current going in to the terminal as
            ## self.i(x)[terminal_node] + u(t) + dq(x)/dt. 
            ## This will work since i(x) returns
            ## the sum of all the currents going out from the
            ## terminal node that originates from devices within
            ## the circuit. According to Kirchoff's current law
            ## of the terminal node
            ## -I_external + sum(I_internal_k) = 0
            ## Where I_external represents the current coming from
            ## outside the circuit going *in* to the terminal node,
            ## I_internal_k represents one of the currents that flows
            ## from the terminal node to a device within the circuit.
            ## So now we can calculate I_external as 
            ## I_external = sum(I_internal_k) = 
            ## self.I(x)[terminal_node] + u(t) + dq(x)/dt =
            ## self.I(x)[terminal_node] + u(t) + sum(dq(x)/dx_k * dx_k/dt) =
            ## self.I(x)[terminal_node] + u(t) + C(x) * dx/dt

            branch_sign = self.get_terminal_branch(branch_or_term)

            if branch_sign != None:
                branch, sign = branch_sign
            else:
                terminal_node = self.nodenames[branch_or_term]
                
                if xdot == None:
                    raise ValueError('xdot argument must not be None if no' 
                                     'branch is connected to the terminal')

                terminal_node_index = self.get_node_index(terminal_node)

                if linearized:
                    return dot(self.G(xdcop)[terminal_node_index], x) + \
                        dot(self.C(xdcop)[terminal_node_index], xdot) + \
                        self.u(t, analysis = 'ac')[terminal_node_index]
                        
                else:
                    return self.i(x)[terminal_node_index] + \
                        dot(self.C(x)[terminal_node_index], xdot) + \
                        self.u(t)[terminal_node_index]
        else:
            branch = branch_or_term
            sign = 1

        branchindex = self.get_branch_index(branch)

        if refnode_removed:
            branchindex -= 1

        return sign * x[branchindex]      

    def update_ipar(self, parent_ipar, variables=None):
        """Calculate numeric values of instance parameters"""
        substvalues = ((Variable, variables),
                       (Parameter, parent_ipar))
        newipar = self.ipar_expressions.eval_expressions(substvalues)
        self.ipar.update(newipar)

    def __repr__(self):
        return self.__class__.__name__ + \
               '(' + \
               ','.join([repr(self.nodenames[term].name) for term in self.terminals] +
                        ['%s=%s'%(par.name, self.ipar.get(par)) 
                         for par in self.ipar.parameters]) + ')'

    def _instance_nodes(self, instancenodes, instance, instancename):
        """Return circuit nodes from instance nodes
        """
        for instancenode in instancenodes:
            if instancenode.isglobal:
                yield instancenode
            elif instancenode.name in instance.terminals:
                terminal = instancenode.name
                yield self.term_node_map[instancename][terminal]
            else:
                yield Node(instancename + '.' + instancenode.name)

    def _instance_branches(self, instance, instancename, 
                           instancebranches = None):
        """Return circuit branches from instance branches
        """
        if instancebranches == None:
            instancebranches = instance.branches

        for instancebranch in instancebranches:
            plus, minus = self._instance_nodes([instancebranch.plus, 
                                                instancebranch.minus],
                                               instance, instancename)
            yield Branch(plus,minus)

    @property
    def _nterminalnodes(self):
        """Return number of terminal nodes"""
        return len(self.terminals)

class SubCircuit(Circuit):
    """
    SubCircuit is container for circuit instances

    **Attributes**
        *elements* 
          dictionary of Circuit objects keyed by its instance name

        *elementnodemap*
          list of translation lists which translate between node indices of the
          elements and node indices of the parent

        *term_node_map* 
          dictionary of instance terminal to node object maps keyed by the 
          instance name. The maps are themselves dictionaries keyed by 
          the terminal names.

    """
    def __init__(self, *args, **kvargs):
        self.elements = {}
        self.elementnodemap = {}
        self.term_node_map = {}
        self._mapmatrix = {}
        Circuit.__init__(self, *args, **kvargs)

    def __eq__(self, a):
        return super(SubCircuit, self).__eq__(a) and \
            self.elements == a.elements and \
            self.elementnodemap == a.elementnodemap and \
            self.term_node_map == a.term_node_map

    def __copy__(self):
        newc = Circuit.__copy__(self)
        
        newc.elements = copy(self.elements)
        newc.elementnodemap = copy(self.elementnodemap)
        newc.term_node_map = copy(self.term_node_map)
        newc._rep_nodemap_list = copy(self._rep_nodemap_list)
        
        return newc

    def netlist(self, top = True):
        """
        >>> from elements import *
        >>> a = SubCircuit()
        >>> a['R1'] = R(1,2)
        >>> print a.netlist()
        R1 1 2 R r=1000.0 noisy=True
    
        """
        out = []

        if top:
            subcircuits = set([instance.__class__
                           for instance in self.elements.values() 
                           if isinstance(instance, SubCircuit)])

            for subcircuit_class in subcircuits:
                out.append(subcircuit_class().netlist(top = False))
        else:
            subcktdef = '.subckt ' + self.__class__.__name__ + ' ' + \
                ' '.join(self.terminals)

            out.append(subcktdef)

        for instancename, instance in self.elements.items():
            termnodes = self._instance_nodes(instance.nodes, instance, 
                                             instancename)
            nodes = ' '.join([str(self.term_node_map[instancename][terminal])
                              for terminal in instance.terminals])
            
            parameters = ' '.join(['%s=%s'%(par.name, instance.ipar.get(par)) 
                                   for par in instance.ipar.parameters])

            if top:
                n_indent = 0
            else:
                n_indent = 2
            
            out.append(indent(instancename + ' ' + nodes + ' ' + 
                              instance.__class__.__name__ + ' ' +
                              parameters, n = n_indent))

        if not top:
            out.append('.ends')
            
        return '\n'.join(out)


    def __str__(self):
        return self.netlist(top=False)

    def add_instance(self, instancename, instance, **connection):
        """Add instance to the circuit.
        
        optional named arguments are used to map terminals of the 
        instance to nodes of the circuit
        """

        if instance.toolkit is not self.toolkit:
            raise ValueError('Instance must use the same toolkit as parent')

        if instancename in self.elements:
            del self[instancename]

        self.elements[instancename] = instance

        ## Add local nodes and branches from new instance
        for node in instance.non_terminal_nodes(instancename):
            self.append_node(node)

        ## Create term_node_map entry for the new instance
        term_node_map = self.term_node_map[instancename] = {}

        ## Connect terminal to node
        for terminal, node in connection.items():
            ## Create a node object if it is not already
            if not isinstance(node, Node):
                node = Node(str(node))
            
            ## Add node
            self.append_node(node)

            ## Update terminal-node map
            term_node_map[terminal] = node            

        ## Add branches
        newbranches = self._instance_branches(instance, instancename)
        self.append_branches(*newbranches)

        ## Update circuit node - instance map
        self.update_node_map()

    def __setitem__(self, instancename, element):
        """Adds an instance to the circuit"""

        self.add_instance(instancename, element, **element.terminalhook)

        element.terminalhook = None


    def __delitem__(self, instancename):
        """Removes instance from circuit
        
        >>> from elements import *
        >>> c = SubCircuit()
        >>> c['V'] = VS(gnd, gnd)
        >>> del c['V']
        >>> c.branches
        []
        
        """
        element = self.elements.pop(instancename)

        ## Remove floating terminal nodes and internal nodes
        othernodes = set(self.terminal_nodes())
        for instance_name, e in self.elements.items():
            terminal_nodes = set([self.term_node_map[instance_name][term]
                                  for term in e.terminals])
            othernodes.update(terminal_nodes)
        internal_nodes = set(element.non_terminal_nodes(instancename))
        terminal_nodes = set([self.term_node_map[instancename][term]
                              for term in element.terminals])
        floating_terminal_nodes =  terminal_nodes - othernodes
        removed_nodes = internal_nodes | floating_terminal_nodes

        for node in removed_nodes:
            self.nodes.remove(node)
            del self.nodenames[node.name] 

        for branch in self._instance_branches(element, instancename):
            self.branches.remove(branch)

        del self.term_node_map[instancename]

        self.update_node_map()

    def __getitem__(self, instancename):
        """Get local or hierarchical instance by name"""
        
        instname_parts = [part for part in instancename.split('.')]

        if instancename == '':
            return self
        if len(instname_parts) == 1:
            return self.elements[instname_parts[0]]
        else:
            top_instancename = instname_parts[0]
            sub_instancename = '.'.join(instname_parts[1:])
            return self.elements[top_instancename][sub_instancename]

    def get_terminal_branch(self, terminalname):
        """Find the branch that is connected to the given terminal

        If no branch is found or if there are more branches than one, None is 
        returned

        The name can be a hierachical name and the notation is 'I1.I2.plus' for
        the terminal 'plus' of instance I2 of instance I1

        Returns
        -------
        Tuple of Branch object and an integer indicating if terminal is 
        connected to the positive or negative side of the branch. 
        1 == positive side, -1 == negative side

        >>> from elements import *
        >>> c = SubCircuit()
        >>> net1, net2 = c.add_nodes('net1', 'net2')
        >>> c['vs'] = VS(net1, net2)
        >>> c.get_terminal_branch("vs.minus")
        (Branch(Node('net1'),Node('net2')), -1)
        
        """
        hierlevels = [part for part in terminalname.split('.')]

        if len(hierlevels)==1:
            return Circuit.get_terminal_branch(self, hierlevels[0])
        else:
            instancename = hierlevels[0]
            topelement = self.elements[instancename]
            branch_sign = \
                topelement.get_terminal_branch('.'.join(hierlevels[1:]))

            if branch_sign:
               ## Add prefix to node names in branch
               branch_gen = self._instance_branches(topelement, instancename, 
                                                    (branch_sign[0],))

               return branch_gen.next(), branch_sign[1]
        
    def get_node(self, name):
        """Find a node by name.
        
        >>> from elements import *
        >>> c = SubCircuit()
        >>> c['V1'] = VS(1,0)
        >>> c.get_node('V1.plus')
        Node('1')
        
        """
        if name in self.nodenames:
            return self.nodenames[name]
        else:
            path = name.split('.')

            top = path[0]

            if len(path) < 2:
                return ValueError('Node name %s not found'%name)
            elif len(path) > 2:
                return top + '.' + self[top].get_node('.'.join(path[1:]))
            else:
                element_node_index = self[top].get_node_index('.'.join(path[1:]))
                node_index = self.elementnodemap[top][element_node_index]
                return self.nodes[node_index]
            
    def get_node_name(self, node):
        """Find the name of a node
        
        >>> c1 = SubCircuit()
        >>> c2 = SubCircuit()
        >>> c1['I1'] = c2
        >>> n1 = c2.add_node("net1")
        >>> c1.get_node_name(n1)
        'net1'

        """

        ## Use name of node object if present
        if node.name != None:
            return node.name
        
        ## First search among the local nodes
        name = Circuit.get_node_name(self, node)
        if name != None:
            return name
        
        ## Then search in the circuit elements
        for instname, element in self.elements.items():
            name =  element.get_node_name(node)
            if name != None:
                return instname + '.' + name
        
    def update_node_map(self):
        """Update the elementnodemap attribute"""

        self.elementnodemap = {}
        self._rep_nodemap_list = {}
        
        for instance_name, element in self.elements.items():
            nodemap = self.term_node_map[instance_name]
            element_nodes = [nodemap[terminal] for terminal in element.terminals]

            for node in element.non_terminal_nodes(instance_name):
                element_nodes.append(node)
        
            element_branches = self._instance_branches(element, instance_name)

            nodemap = \
                [self.nodes.index(node) for node in element_nodes] + \
                [self.branches.index(branch) + len(self.nodes) 
                 for branch in element_branches]

            self.elementnodemap[instance_name] = nodemap

            ## Create mapping matrix
            if len(nodemap) > 0:
                mapmatrix = self.toolkit.zeros((self.n, len(nodemap)),
                                               dtype = np.integer)
            
                for inst_node_index, node_index in enumerate(nodemap):
                    mapmatrix[node_index, inst_node_index] = 1
                self._mapmatrix[instance_name] = mapmatrix
            else:
                self._nodemap = None

    def update_ipar(self, parent_ipar, variables=None):
        """Calculate numeric values of instance parameters"""
        super(SubCircuit, self).update_ipar(parent_ipar, variables)
        
        ## Update ipar in elements
        for element in self.elements.values():
            element.update_ipar(self.ipar, variables)
        
    def G(self, x, epar=defaultepar):
        return self._add_element_submatrices('G', x, (epar,))

    def C(self, x, epar=defaultepar):
        return self._add_element_submatrices('C', x, (epar,))

    def u(self, t=0.0, epar=defaultepar, analysis=None):
        dtype = None
        if analysis == 'ac':
            dtype = self.toolkit.ac_u_dtype

        return self._add_element_subvectors('u', None, (t,epar,analysis), 
                                            dtype=dtype)

    def i(self, x, epar=defaultepar):
        return self._add_element_subvectors('i', x, (epar,))

    def CY(self, x, w, epar=defaultepar):
        """Calculate composite noise source correlation matrix

        The noise sources in one element are assumed to be uncorrelated 
        with the noise sources in the other elements.

        """
        return self._add_element_submatrices('CY', x, (w, epar,))

    def save_current(self, terminal):
        """Returns a circuit where the given terminal current is saved
        
        The terminal can be a hierarchical name and the notation is I1.term
        for terminal 'term' of instance 'I1'
        """
        
        hierterm = [part for part in terminal.split('.')]

        if len(hierterm) == 1:
            return Circuit.save_current(self, terminal)
        elif len(hierterm) >= 2:
            base = self
            for instance in hierterm[:-2]:
                base = base[instance]

            base.add_instance(hierterm[0], base[hierterm[0]].save_current(hierterm[1]),
                              **base.term_node_map[hierterm[0]])
        else:
            raise Exception('Invalid terminal name: %s'%terminal)
        
        return self

    def extract_i(self, x, branch_or_term, xdot = None,
                  refnode = gnd, refnode_removed = False, 
                  linearized = False, xdcop = None):
        if type(branch_or_term) is types.StringType:
            if self.get_terminal_branch(branch_or_term) == None:

                hierlevels = [part for part in branch_or_term.split('.')]

                if len(hierlevels) > 1:
                    instance_name = hierlevels[0]
                    instance = self.elements[instance_name]
                    terminal_name = '.'.join(hierlevels[1:])

                    ## Get slice of x-vector for the instance
                    nodemap = self.elementnodemap[instance_name]

                    subx = x[nodemap]

                    if xdot != None:
                        subxdot = xdot[nodemap]
                    else:
                        subxdot = None

                    return instance.extract_i(subx, terminal_name, 
                                              xdot = subxdot, 
                                              refnode=refnode,
                                              refnode_removed=refnode_removed,
                                              linearized = linearized, 
                                              xdcop = xdcop)


        return Circuit.extract_i(self, x, branch_or_term, xdot = xdot,
                                 refnode=refnode, 
                                 refnode_removed=refnode_removed,
                                 linearized = linearized, xdcop = xdcop)
        
    def update(self, subject):
        """This is called when an instance parameter is updated"""
        for element in self.elements.values():
            element.update_ipar(self.ipar)
        
    def _add_element_submatrices(self, methodname, x, args):
        dot = self.toolkit.dot
        n = self.n
        lhs = self.toolkit.zeros((n,n))

        for instance, element in self.elements.items():
            nodemap = self.elementnodemap[instance]

            if x != None:
                subx = x[nodemap]
                try:
                    rhs = getattr(element, methodname)(subx, *args)
                except Exception, e:
                    raise e.__class__(str(e) + ' at element ' + str(element) 
                                      + ', args='+str(args))
            else:
                rhs = getattr(element, methodname)(*args)
                
            T = self._mapmatrix[instance]
        
            lhs += dot(dot(T, rhs), T.T)

        return lhs

    def _add_element_subvectors(self, methodname, x, args, dtype=None):
        n = self.n
        lhs = self.toolkit.zeros(n, dtype=dtype)

        for instance,element in self.elements.items():
            if x != None:
                subx = x[self.elementnodemap[instance]]
                rhs = getattr(element, methodname)(subx, *args)
            else:
                rhs = getattr(element, methodname)(*args)

            T = self._mapmatrix[instance]
            
            lhs += self.toolkit.dot(T, rhs)

        return lhs

    @property
    def xflatelements(self):
        """Iterator over all elements and subelements"""
        for e in self.elements.values():
            if not isinstance(e, SubCircuit):
                yield e
            else:
                for sube in e.xflatelements:
                    yield sube

    def translate_branch(self, branch, instance):
        """Return branch from a local branch in an instance"""
        return Branch(self.get_node(instance + '.' + branch.plus.name),
                      self.get_node(instance + '.' + branch.minus.name))
                    
class ProbeWrapper(SubCircuit):
    """Circuit wrapper that adds voltage sources for current probing"""
    def __init__(self, circuit, terminals = ()):

        ## Copy nodes, branches, terminals and parameters
        self.terminals = circuit.terminals

        terminalnodes = [circuit.get_node(terminal)
                         for terminal in circuit.terminals]
        super(ProbeWrapper, self).__init__(*terminalnodes)
        
        self.add_instance('wrapped', circuit, 
                          **dict(zip(circuit.terminals, terminalnodes)))
        
        for terminal in terminals:
            self.save_current(terminal)

    def save_current(self, terminal):
        """Returns a circuit where the given terminal current is saved"""
        
        ## Add probe to terminal if it does not already exists
        if self.get_terminal_branch(terminal) == None:
            node = self.get_node(terminal)
            ## Add internal node
            internal_node = self.add_node(terminal + '_internal')

            ## Create zero-voltage voltage source between terminal and
            ## internal node
            self[terminal + '_probe'] = IProbe(node, internal_node)

            ## Re-connect wrapped circuit to internal node
            term_node_map = self.term_node_map['wrapped']
            circuit = self['wrapped']
            del self['wrapped']
            term_node_map[terminal] = internal_node
            self.add_instance('wrapped', circuit, **term_node_map)

        return self
    

class CircuitProxy(Circuit):
    def __init__(self, circuit, parent=None, instance_name=None):
        super(CircuitProxy, self).__init__(self)
        self.device = circuit
        self.terminals = circuit.terminals
        self.nodes = circuit.nodes
        self.nodenames = circuit.nodenames
        self.branches = circuit.branches
        self.ipar = circuit.ipar
        
        ## Find out how this instance was connected to its parent
        ## and set terminalhook accordingly
        if isinstance(parent, SubCircuit) and instance_name != None:
            self.terminalhook = parent.term_node_map[instance_name]

    def G(self, x, epar=defaultepar): return self.device.G(x,epar)
    def C(self, x, epar=defaultepar): return self.device.C(x,epar)
    def u(self, t=0.0, epar=defaultepar, analysis=None): 
        return self.device.u(x,epar)
    def i(self, x, epar=defaultepar): return self.device.i(x,epar)
    def q(self, x, epar=defaultepar): return self.device.q(x,epar)
    def CY(self, x, w, epar=defaultepar): return self.device.CY(x,epar)

def instjoin(*instnames):
    """Return hierarchical instance names from instance name components
    
    >>> instjoin('I1','VS')
    'I1.VS'
    >>> instjoin('','VS')
    'VS'
    
    """
    return '.'.join([part for part in instnames if len(part) > 0])
    
class IProbe(Circuit):
    """Zero voltage independent voltage source used for current probing"""
    terminals = ['plus', 'minus']

    def __init__(self, plus, minus, **kvargs):
        Circuit.__init__(self, plus, minus, **kvargs)
        self.append_branches(Branch(self.nodenames['plus'], 
                                    self.nodenames['minus']))

    def G(self, x, epar=defaultepar):
        return self.toolkit.array([[0 , 0, 1],
                                   [0 , 0, -1],
                                   [1 , -1, 0]])

    @property
    def branch(self):
        """Return the branch (plus, minus)"""
        return self.branches[0]

class Quantity(object):
    """Reference to voltage or current of a branch or node

    The quantity can be used in behavioural modelling or post processing
    when one want to refer to a voltage or current of a branch or node
    
    >>> a, b = Node('a'), Node('b')
    >>> Quantity('I', Branch(a,b))
    I(a,b)
    >>> Quantity('V', a)
    V(a)

    """
    
    def __init__(self, quantity, branch_or_node):
        """The quantity can be 'V' or 'I' which corresponds to voltage or
        current of the Branch or Node object branch_or_node"""

        if quantity not in ('V', 'I'):
            raise ValueError("quantity must be either 'V' or 'I'")
        if not isinstance(branch_or_node, (Node, Branch)):
            raise ValueError('branch_or_node must be a Branch or Node object')
        
        if quantity == 'I' and isinstance(branch_or_node, Node):
            raise ValueError('Current can only be taken on branches')

        self.quantity = quantity
        self.branch_or_node = branch_or_node

    @property
    def isnode(self): return isinstance(self.branch_or_node, Node)

    @property
    def isbranch(self): return isinstance(self.branch_or_node, Branch)
        
    def __repr__(self):
        if isinstance(self.branch_or_node, Branch):
            return self.quantity + '(' + str(self.branch_or_node.plus.name) + \
                ',' + str(self.branch_or_node.minus.name) + ')'
        else:
            return self.quantity + '(' + str(self.branch_or_node.name) + ')'

if __name__ == "__main__":
    import doctest
    doctest.testmod()
