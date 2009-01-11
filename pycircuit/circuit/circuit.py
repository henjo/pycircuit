# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from numpy import array, zeros, concatenate, dot, exp, inf
from pycircuit.utilities.param import Parameter, ParameterDict
from pycircuit.utilities.misc import indent, inplace_add_selected, \
    inplace_add_selected_2d, create_index_vectors
from constants import *
from copy import copy
import types


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

    Attributes
    ----------
    *nodes*
      A list that contains Node objects. The first k nodes are terminal nodes
      where k is the number of terminals.

    *branches*
      list of Branch objects. The solver will solve for the 
      currents through the branches.

    *terminals*
      list that contains terminal names

    *instparams*
      A list of valid instance parameters (Parameter objects)

    *mpar*
      A class variable with a ParameterDict containing model specific parameters

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
      and the values are node objects. The value is None when it is not used.
      The only reason for this attribute is to allow for bottom-up 
      instantiations like: cir1['I1'] = R('n1', gnd)

    *linear* 
      A boolean value that is true if i(x) and q(x) are linear 
      functions

    """
    terminals = []
    mpar = ParameterDict()
    instparams = []
    linear = True
    
    def __init__(self, *args, **kvargs):
        self.nodes = []
        self.internalnodes = []
        self.nodenames = {}
        self.branches = []
        self.ipar = ParameterDict(*self.instparams, **kvargs)

        ## Add terminal nodes
        for terminal in self.terminals:
            self.append_node(Node(terminal))

        ## Set temporary terminal mapping information for use by instantiation
        ## method in higher hierarchy
        self.terminalhook = dict(zip(self.terminals, args))
        
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
        if node not in self.nodes:
            self.nodes.append(node)
        self.nodenames[node.name] = node

    def get_terminal_branch(self, terminalname):
        """Find the branch that is connected to the given terminal

        If no branch is found or if there are more branches than one, None is 
        returned

        Returns
        -------
        Tuple of Branch object and an integer indicating if terminal is 
        connected to the positive or negative side of the branch. 
        1 == positive side, -1 == negative side
        
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
        return self.nodes.index(node)

    def get_branch_index(self, branch):
        """Get row in the x vector of a branch instance"""        
        return len(self.nodes) + self.branches.index(branch)

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
        
        >>> cir = R(Node('n1'), gnd, r=1e3)
        >>> newcir = cir.save_current('plus')
        >>> newcir.G(zeros(4))
        array([[0, 0, 0, 1],
               [0, 0.001, -0.001, 0],
               [0, -0.001, 0.001, -1],
               [1, 0, -1, 0]], dtype=object)
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
        return zeros((self.n, self.n), dtype=object)

    def C(self, x, epar=defaultepar):
        """Calculate the C (transcapacitance) matrix given the x-vector"""
        return zeros((self.n, self.n), dtype=object)

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
        return zeros(self.n, dtype=object)

    def i(self, x, epar=defaultepar):
        """Calculate the i vector as a function of the x-vector

        For linear circuits i(x(t)) = G*x
        """
        return dot(self.G(x), x)

    def q(self, x, epar=defaultepar):
        """Calculate the q vector as a function of the x-vector

        For linear circuits q(x(t)) = C*x
        """
        return dot(self.C(x), x)

    def CY(self, x, w, epar=defaultepar):
        """Calculate the noise sources correlation matrix

        Arguments
        ---------
        x -- (numpy array) the state vector
        w -- Angular frequency
        epar -- (ParameterDict) Environment parameters

        """
        return zeros((self.n, self.n), dtype=object)

    def next_event(self, t):
        """Returns the time of the next event given the current time t"""
        return inf
    
    def name_state_vector(self, x, analysis=''):
        """Return a dictionary of the x-vector keyed by node and branch names

        >>> c = SubCircuit()
        >>> n1 = c.add_node('net1')
        >>> c['is'] = IS(gnd, n1, i=1e-3)
        >>> c['R'] = R(n1, gnd, r=1e3)
        >>> c.name_state_vector(array([[1.0]]))
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
        
        >>> c = SubCircuit()
        >>> n1, n2 = c.add_nodes('n1','n2')
        >>> c['R1'] = R(n1, n2, r=1e3)
        >>> c['R2'] = R(n2, gnd, r=1e3)
        >>> c.extract_v(array([1.0, 0.5, 0.0]), 'n1', 'n2')
        0.5
        >>> c.extract_v(array([1.0, 0.5, 0.0]), c.nodes[0])
        1.0
        >>> c.extract_v(array([1.0, 0.5]), c.nodes[0], refnode_removed = True)
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

        >>> c = SubCircuit()
        >>> net1 = c.add_node('net1')
        >>> c['vs'] = VS(net1, gnd)
        >>> c.extract_i(array([1.0, 0, -1e-3]), 'vs.minus')
        0.001
        >>> c.extract_i(array([1.0, -1e-3]), 'vs.minus', refnode_removed = True)
        0.001
        
        """
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

    def _instance_branches(self, instance, instancename, instancebranches = None):
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

    Attributes
    ----------
    *elements* 
      dictionary of Circuit objects keyed by its instance name
    
    *elementnodemap*
      list of translation lists that translate between node indices of the
      elements to the node index in the SubCircuit object for each element

    *term_node_map* 
      dictionary of instance terminal to node object maps keyed by the instance
      name. the map is itself a dictionary keyed by the terminal names

    """
    def __init__(self, *args, **kvargs):
        Circuit.__init__(self, *args, **kvargs)
        self.elements = {}
        self.elementnodemap = {}
        self._rep_nodemap_list = {}
        self.term_node_map = {}
        
    def netlist(self, top = True):
        """
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
        self.branches.extend(newbranches)

        ## Update circuit node - instance node map
        self.update_node_map()

    def __setitem__(self, instancename, element):
        """Adds an instance to the circuit"""

        self.add_instance(instancename, element, **element.terminalhook)

        element.terminalhook = None


    def __delitem__(self, instancename):
        """Removes instance from circuit
        
        >>> c = SubCircuit()
        >>> c['V'] = VS(gnd, gnd)
        >>> del c['V']
        >>> c.branches
        []
        
        """
        element = self.elements.pop(instancename)

        othernodes = set([])
        for instance_name, e in self.elements.items():
            othernodes.update(set(e.non_terminal_nodes(instance_name)))

        for node in set(element.non_terminal_nodes(instancename)) - othernodes:
            self.nodes.remove(node)
            del self.nodenames[node.name] 

        for branch in self._instance_branches(element, instancename):
            self.branches.remove(branch)

        del self.term_node_map[instancename]

        self.update_node_map()

    def __getitem__(self, instancename):
        """Get instance"""
        return self.elements[instancename]

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
            self._rep_nodemap_list[instance_name] = create_index_vectors(nodemap)
                
    def G(self, x, epar=defaultepar):
        return self._add_element_submatrices('G', x, (epar,))

    def C(self, x, epar=defaultepar):
        return self._add_element_submatrices('C', x, (epar,))

    def u(self, t=0.0, epar=defaultepar, analysis=None):
        return self._add_element_subvectors('u', None, (t,epar,analysis))

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
        
    def _add_element_submatrices(self, methodname, x, args):
        n = self.n
        lhs = zeros((n,n), dtype=object)

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
                
            inplace_add_selected_2d(lhs, self._rep_nodemap_list[instance], rhs)

        return lhs

    def _add_element_subvectors(self, methodname, x, args):
        n = self.n
        lhs = zeros(n, dtype=object)

        for instance,element in self.elements.items():
            if x != None:
                subx = x[self.elementnodemap[instance]]
                rhs = getattr(element, methodname)(subx, *args)
            else:
                rhs = getattr(element, methodname)(*args)

            inplace_add_selected(lhs, self._rep_nodemap_list[instance], rhs)

        return lhs

    @property
    def xflatelements(self):
        """Iterator over all elements and subelements"""
        for e in self.elements.values():
            if not isinstance(e, SubCircuit):
                yield e
            else:
                yield e.xflatelements

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
    
class R(Circuit):
    """Resistor element

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> c['R']
    R('plus','minus',r=1000.0,noisy=True)
    >>> c.G(zeros(2))
    array([[0.001, -0.001],
           [-0.001, 0.001]], dtype=object)
    >>> c = SubCircuit()
    >>> n2=c.add_node('2')
    >>> c['R'] = R(n1, n2, r=1e3)
    >>> c.G(zeros(2))
    array([[0.001, -0.001],
           [-0.001, 0.001]], dtype=object)

    """
    terminals = ['plus', 'minus']
    instparams = [Parameter(name='r', desc='Resistance', unit='ohm', 
                            default=1e3),
                  Parameter(name='noisy', desc='No noise', unit='', 
                            default=True),
                  ]

    def G(self, x, epar=defaultepar):
        g = 1/self.ipar.r
        return  array([[g, -g],
                        [-g, g]], dtype=object)

    def CY(self, x, w, epar=defaultepar):
        if self.ipar.noisy:
            if 'kT' in epar:
                iPSD = 4*epar.kT / self.ipar.r
            else:
                iPSD = 4*kboltzmann * epar.T / self.ipar.r
        else:
            iPSD = 0

        return  array([[iPSD, -iPSD],
                       [-iPSD, iPSD]], dtype=object)

class G(Circuit):
    """Conductor element

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['G'] = G(n1, gnd, g=1e-3)
    >>> c['G']
    G('plus','minus',g=0.001,nonoise=False)
    >>> c.G(zeros(2))
    array([[0.001, -0.001],
           [-0.001, 0.001]], dtype=object)

    """
    terminals = ['plus', 'minus']
    instparams = [Parameter(name='g', desc='Conductance', unit='S', 
                            default=1e-3),
                  Parameter(name='nonoise', 
                            desc='If true the conductance is noiseless', 
                            unit='', default=False),
                  ]

    def G(self, x, epar=defaultepar):
        g = self.ipar.g
        return  array([[g, -g],
                        [-g, g]], dtype=object)

    def CY(self, x, w, epar=defaultepar):
        if not self.ipar.nonoise:
            if 'kT' in epar:
                iPSD = 4*epar.kT*self.ipar.g
            else:
                iPSD = 4*kboltzmann*epar.T*self.ipar.g
            return  array([[iPSD, -iPSD],
                           [-iPSD, iPSD]], dtype=object)
        else:
            return super(G, self).CY(x, w, epar=epar)
        

class C(Circuit):
    """Capacitor

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['C'] = C(n1, gnd, c=1e-12)
    >>> c.G(zeros(2))
    array([[0, 0],
           [0, 0]], dtype=object)
    >>> c.C(zeros(2))
    array([[1e-12, -1e-12],
           [-1e-12, 1e-12]], dtype=object)

    """

    terminals = ['plus', 'minus']    
    instparams = [Parameter(name='c', desc='Capacitance', 
                            unit='F', default=1e-12)]

    def C(self, x, epar=defaultepar):
        return array([[self.ipar.c, -self.ipar.c],
                      [-self.ipar.c, self.ipar.c]])

class L(Circuit):
    """Inductor

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['L'] = L(n1, gnd, L=1e-9)
    >>> c.G(zeros(3))
    array([[0.0, 0.0, 1.0],
           [0.0, 0.0, -1.0],
           [1.0, -1.0, 0.0]], dtype=object)
    >>> c.C(zeros(3))
    array([[0, 0, 0],
           [0, 0, 0],
           [0, 0, 1e-09]], dtype=object)
    """
    terminals = ['plus', 'minus']    
    instparams = [Parameter(name='L', desc='Inductance', 
                            unit='H', default=1e-9)]

    _G = array([[0.0 , 0.0, 1.0],
                [0.0 , 0.0, -1.0],
                [1.0 , -1.0, 0.0]])
    def __init__(self, plus, minus, L=0.0):
        Circuit.__init__(self, plus, minus, L=L)
        self.branches.append(Branch(plus, minus))
    def G(self, x, epar=defaultepar):
        return self._G
    def C(self, x, epar=defaultepar):
        n = self.n
        C = zeros((n,n), dtype=object)
        C[-1,-1] = self.ipar.L
        return C

class IProbe(Circuit):
    """Zero voltage independent voltage source used for current probing"""
    terminals = ['plus', 'minus']

    def __init__(self, plus, minus, **kvargs):
        Circuit.__init__(self, plus, minus, **kvargs)
        self.branches.append(Branch(self.nodenames['plus'], self.nodenames['minus']))

    def G(self, x, epar=defaultepar):
        return array([[0 , 0, 1],
                      [0 , 0, -1],
                      [1 , -1, 0]], dtype=object)

    @property
    def branch(self):
        """Return the branch (plus, minus)"""
        return self.branches[0]


class VS(Circuit):
    """Independent voltage source

    >>> from analysis import DC
    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> DC(c).solve(refnode=gnd).x
    array([ 1.5   ,  0.    , -0.0015])
    
    """
    terminals = ['plus', 'minus']
    instparams = [Parameter(name='v', desc='Source voltage', 
                            unit='V', default=1),
                  Parameter(name='vac', desc='AC analysis amplitude', 
                            unit='V', default=1),
                  Parameter(name='noisePSD', 
                            desc='Voltage noise power spectral density', 
                            unit='V^2/Hz', default=0)]

    def __init__(self, plus, minus, **kvargs):
        Circuit.__init__(self, plus, minus, **kvargs)
        self.branches.append(Branch(self.nodenames['plus'], self.nodenames['minus']))

    def G(self, x, epar=defaultepar):
        return array([[0 , 0, 1],
                      [0 , 0, -1],
                      [1 , -1, 0]], dtype=object)

    def u(self, t=0.0, epar=defaultepar, analysis=None):
        if analysis == 'ac':
            return array([0, 0, -self.ipar.vac], dtype=object)
        elif analysis == None:
            return array([0, 0, -self.ipar.v], dtype=object)
        else:
            return super(VS, self).u(t,epar,analysis)

    def CY(self, x, w, epar=defaultepar):
        CY = super(VS, self).CY(x, w)
        CY[2, 2] = self.ipar.noisePSD
        return CY

    @property
    def branch(self):
        """Return the branch (plus, minus)"""
        return self.branches[0]

class IS(Circuit):
    """Independent DC current source

    >>> from analysis import DC, gnd as gnd2
    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['is'] = IS(gnd, n1, i=1e-3)
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> DC(c).solve(refnode=gnd).x
    array([ 1.,  0.])

    """
    instparams = [Parameter(name='i', desc='DC Current', 
                            unit='A', default=1e-3),
                  Parameter(name='iac', desc='Small signal current amplitude', 
                            unit='A', default=0),
                  Parameter(name='noisePSD', 
                            desc='Current noise power spectral density', 
                            unit='A^2/Hz', default=0.0)]
    terminals = ['plus', 'minus']

    def u(self, t=0.0, epar=defaultepar, analysis=None):
        if analysis == None:
            return array([self.ipar.i, -self.ipar.i])
        elif analysis == 'ac':
            return array([self.ipar.iac, -self.ipar.iac])
        else:
            return super(IS, self).u(t,epar,analysis)

    def CY(self, x, w, epar=defaultepar):
        return  array([[self.ipar.noisePSD, -self.ipar.noisePSD],
                       [-self.ipar.noisePSD, self.ipar.noisePSD]], dtype=object)

class VCVS(Circuit):
    """Voltage controlled voltage source

    >>> from analysis import DC
    >>> c = SubCircuit()
    >>> n1, n2 =c.add_nodes('1', '2')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['vcvs'] = VCVS(n1, gnd, n2, gnd, g=2)
    >>> c.nodes
    [Node('1'), Node('2'), Node('gnd', isglobal=True)]
    >>> c.branches
    [Branch(Node('1'),Node('gnd', isglobal=True)), Branch(Node('2'),Node('gnd', isglobal=True))]
    >>> c['vcvs'].G(zeros(4))
    array([[0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0],
           [0, 0, 0, 0, 1],
           [0, 0, 0, 0, -1],
           [2, -2, -1, 1, 0]], dtype=object)

    """
    instparams = [Parameter(name='g', desc='Voltage gain', 
                            unit='V/V', default=1)]
    terminals = ['inp', 'inn', 'outp', 'outn']
    def __init__(self, *args, **kvargs):
        Circuit.__init__(self, *args, **kvargs)
        self.branches.append(Branch(self.nodenames['outp'], 
                                    self.nodenames['outn']))
        
    def G(self, x, epar=defaultepar):
        G = super(VCVS, self).G(x)
        branchindex = -1
        inpindex, innindex, outpindex, outnindex = \
           (self.nodes.index(self.nodenames[name]) 
            for name in ('inp', 'inn', 'outp', 'outn'))
        G[outpindex, branchindex] += 1
        G[outnindex, branchindex] += -1
        G[branchindex, outpindex] += -1
        G[branchindex, outnindex] += 1
        G[branchindex, inpindex] += self.ipar.g
        G[branchindex, innindex] += -self.ipar.g
        return G

class VCCS(Circuit):
    """Voltage controlled current source

    >>> from analysis import DC
    >>> c = SubCircuit()
    >>> n1,n2 = c.add_nodes('1', '2')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['vccs'] = VCCS(n1, gnd, n2, gnd, gm=1e-3)
    >>> c['rl'] = R(n2, gnd, r=1e3)
    >>> DC(c).solve(refnode=gnd).x
    array([ 1.5, -1.5,  0. ,  0. ])

    """
    terminals = ['inp', 'inn', 'outp', 'outn']
    instparams = [Parameter(name='gm', desc='Transconductance', 
                            unit='A/V', default=1e-3)]
    
    def G(self, x, epar=defaultepar):
        G = super(VCCS, self).G(x)
        gm=self.ipar.gm
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name]) 
             for name in ('inp', 'inn', 'outp', 'outn'))
        G[outpindex, inpindex] += gm
        G[outpindex, innindex] += -gm
        G[outnindex, inpindex] += -gm
        G[outnindex, innindex] += gm
        return G

class Nullor(Circuit):
    """Nullor

    From Wikipedia, the free encyclopedia

     A nullor is a theoretical two-port network comprised of a nullator at 
     its input and a norator at its output.[1]
     Nullors represent an ideal amplifier, having infinite current, 
     voltage, transconductance and transimpedance gain.[2] 
     Its transmission parameters are all zero.

     1. The name "nullor" was introduced by H.J. Carlin, 
      Singular network elements, 
      IEEE Trans. Circuit Theory, March 1965, vol. CT-11, pp. 67-72.
 
     2. Verhoeven C J M van Staveren A Monna G L E Kouwenhoven, 
       M H L & Yildiz E (2003). 
       Structured electronic design: negative feedback amplifiers.
       Boston/Dordrecht/London: Kluwer Academic, §2.2.2 pp. 32-34. 
       ISBN 1402075901.

    """
    terminals = ('inp', 'inn', 'outp', 'outn')

    def __init__(self, *args, **kvargs):
        Circuit.__init__(self, *args, **kvargs)
        self.branches.append(Branch(self.nodenames['outp'], 
                                    self.nodenames['outn']))

    def G(self, x, epar=defaultepar):
        G = super(Nullor, self).G(x)
        branchindex = -1
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name]) 
             for name in ('inp', 'inn', 'outp', 'outn'))

        G[outpindex, branchindex] += 1
        G[outnindex, branchindex] += -1
        G[branchindex, inpindex] += 1
        G[branchindex, innindex] += -1
        return G


class Transformer(Circuit):
    """Ideal transformer

    >>> from analysis import DC
    >>> c = SubCircuit()
    >>> n1, n2 = c.add_nodes('1', '2')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['vcvs'] = Transformer(n1, gnd, n2, gnd, n=2)
    >>> c['vcvs'].nodes
    [Node('inp'), Node('inn'), Node('outp'), Node('outn')]
    >>> c['vcvs'].branches
    [Branch(Node('outp'),Node('outn'))]
    >>> c['vcvs'].G(zeros(4))
    array([[0, 0, 0, 0, 2],
           [0, 0, 0, 0, -2],
           [0, 0, 0, 0, 1],
           [0, 0, 0, 0, -1],
           [-1, 1, 2, -2, 0]], dtype=object)

    """
    instparams = [Parameter(name='n', desc='Winding ratio', unit='', default=1)]
    terminals = ['inp', 'inn', 'outp', 'outn']
    def __init__(self, *args, **kvargs):
        Circuit.__init__(self, *args, **kvargs)
        self.branches.append(Branch(self.nodenames['outp'], 
                                    self.nodenames['outn']))

    def G(self, x, epar=defaultepar):
        G = super(Transformer, self).G(x)
        branchindex = -1
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name]) 
             for name in ('inp', 'inn', 'outp', 'outn'))
        G[inpindex, branchindex] += self.ipar.n
        G[innindex, branchindex] += -self.ipar.n
        G[outpindex, branchindex] += 1
        G[outnindex, branchindex] += -1
        G[branchindex, outpindex] += self.ipar.n
        G[branchindex, outnindex] += -self.ipar.n
        G[branchindex, inpindex] += -1
        G[branchindex, innindex] += 1
        return G

class Diode(Circuit):
    terminals = ['plus', 'minus']
    mpar = Circuit.mpar.copy( 
        Parameter(name='IS', desc='Saturation current', 
                  unit='A', default=1e-13))
    linear = False
    def G(self, x, epar=defaultepar):
        VD = x[0]-x[1]
        VT = kboltzmann*epar.T / qelectron
        g = self.mpar.IS*exp(VD/VT)/VT
        return array([[g, -g],
                      [-g, g]], dtype=object)

    def i(self, x, epar=defaultepar):
        """
        
        """
        VD = x[0]-x[1]
        VT = kboltzmann*epar.T / qelectron
        I = self.mpar.IS*(exp(VD/VT)-1.0)
        return array([I, -I])


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
