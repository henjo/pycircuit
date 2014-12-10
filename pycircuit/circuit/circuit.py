# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from pycircuit.utilities.param import Parameter, ParameterDict#, EvalError
from pycircuit.utilities.misc import indent, inplace_add_selected, \
    inplace_add_selected_2d, create_index_vectors
from copy import copy
import types
import numeric
import numpy as np

default_toolkit = numeric

timedomain_analyses = ('dc', 'tran')

def makenode(node_or_name):
    if node_or_name is None:
        return None
    if isinstance(node_or_name, Node):
        return node_or_name
    else:
        return Node(str(node_or_name))

def basename(hiername):
    sepchar = '.'
    return hiername.split(sepchar)[-1]

def instname(hiername):
    sepchar = '.'
    return sepchar.join(hiername.split(sepchar)[:-1])

class Node(object):
    """A Node object represents a point in an electric circuit"""
    def __init__(self, name=None, isglobal = False):
        name = str(name)

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

    @property
    def I(self):
        return Quantity('I', self)

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
    
    A branch is used to tell the nodal analyzer to enforce a potential or define a flow
    between its plus and minus nodes. It may also just be used by other branches to 
    measure potential difference or current through a potential branch. A branch may also
    be 
    

    **Attributes**
        *plus*
          Positive node
        *minus*
          Negative node
        *i*
          Flow between its plus and minus node
        *q*
          Time integral
        *v*
          Potential between plus and 
        *G*
          i jacobian, di/d(v or i of input branches). The derivatives are stored as a list
        *C*
          q jacobian, dq/d(v or i of input branches). The derivatives are stored as a list
        *potential*
          If set the branch the solver will enforce a potential between it's plus and minus nodes. Otherwise
          it is a flow branch.
        *output* 
          out is a string that indicates if the circuit may set its q, u and/or i value. For example, if the
          string contains the letter 'q' it may set the q value. The string 'qui' indicates that the circuit
          may set all three values.
        *input* 
          If set, the circuit may use the voltage or current value of the branch
        *indirect*
          If set, the output quantity is not set directly but rather it takes the value that is required
          to solve an equation. 
        *linear*
          If true, the voltage or current of output branches is a linear function
        *noisy*
          If true, the voltage or current of output branches is noisy    
        *noise_correlated*
          If true, the branch noise is correlated with noise from other branches
        *noise_constant*
          If true, the branch noise does not vary with operating point, the noise may 
          with the absolute temperature and still be considered constant
    
    """

    ## Possibly add inormation regarding memory(q), linearity etc.
    
    potential = False
    input = False
    linear = True
    output = ''
    
    def __init__(self, plus, minus, potential=None, 
                 output='', input=False, linear=True, 
                 indirect=False,
                 noisy=False, noise_correlated = False, noise_constant = False
                 ): # default is 'flow' branch, not 'potential' branch
        if potential is not None:
            self.potential = potential
        self.plus             = makenode(plus)
        self.minus            = makenode(minus)
        self.output           = output
        self.input            = input
        self.indirect         = indirect
        self.linear           = linear
        self.noisy            = noisy
        self.noise_correlated = noise_correlated
        self.noise_constant   = noise_constant

    def __repr__(self):
        return '%s(%s, %s)' % (self.__class__.__name__, self.plus, self.minus)

    def G(self,wrt):
        '''Method for calculating conductance
        '''
        pass

    def C(self,wrt):
        '''Method for calculating derivative of q wrt to branch potential or flow
        '''
        pass
    
class BranchI(Branch):
    '''Flow branch

    i = dq/dt
    '''
    potential = False

class BranchV(Branch):
    '''Potential type branch

    v = dq/dt
    '''
    potential = True

class BranchRef(Branch):
    """Reference to a branch in an instance in a SubCircuit"""
    def __init__(self, instname, inst, branch):
        self.instname    = instname
        self.inst        = inst
        self.branch      = branch

        for attr in ('plus', 'minus', 'potential', 'input',  'output', 'linear', 
                     'indirect', 'noisy', 'noise_constant', 'noise_correlated'):
            setattr(self, attr, getattr(self.branch, attr))

    def plusnode(self, cir):
        """Return circuit node name of branch plus node"""
        return cir.get_node(instjoin(self.instname, self.plus.name))

    def minusnode(self, cir):
        """Return circuit node name of branch plus node"""
        return cir.get_node(instjoin(self.instname, self.minus.name))

    def __repr__(self):
        return 'BranchRef(%s, %s, %s)' % (self.instname, self.inst, repr(self.branch))

    def __hash__(self):
        return hash(self.branch) ^ hash(self.instname)

    def __eq__(self, a):
        if a.__class__ is BranchRef:
            return (self.branch == a.branch) and (self.instname == a.instname)
        else:
            return False

### Default reference node
gnd = Node("gnd", isglobal=True)

defaultepar = ParameterDict(
    Parameter("T", "Temperature", unit="K", default = 300,),
    Parameter("t", "Time", unit="s", default = 0,),
    Parameter("w", "Angular frequency used by frequency domain analyses", 
              unit="rad/s", default = 0,),
    Parameter("analysis", "Analysis", default = 'ac',),
    )

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
          list of Branch objects.

        *terminals*
          list of terminal names

        *instparams*
          A list of valid instance parameters (Parameter objects)

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
    instparams = []
    linear = True
    
    def __init__(self, *args, **kvargs):
        if 'toolkit' in kvargs:
            self.toolkit = kvargs['toolkit']
            del kvargs['toolkit']
        else:
            self.toolkit = default_toolkit

        self.nodenames = {}

        ## Add terminal nodes
        self.add_terminals(self.terminals)

        ## Set temporary terminal mapping information for use by instantiation
        ## method in higher hierarchy
        self.terminalhook = dict(zip(self.terminals, args))

        ## Create instance parameters
        self.iparv = ParameterDict(*self.instparams)
        self.ipar = ParameterDict(*self.instparams)

        ## Subscribe to changes on ipar 
        self.ipar.attach(self, updatemethod='_ipar_changed')

        ## set instance parameters from arguments
        self.ipar.set(**kvargs)
        
        ## Subscribe to updates of instance parameters
        if hasattr(self, 'update'):
            self.iparv.attach(self)
            self.update(self.ipar)

        ## Copy branches from class variables
        ## Is this needed or not!!
#        self.branches = map(copy, self.branches)

    def __eq__(self, a):
        return self.__class__ == a.__class__ and \
            self.nodes == a.nodes and \
            self.nodenames == a.nodenames and self.branches == a.branches and \
            self.iparv == a.iparv

    ## Change/remove this ??
    def __copy__(self):
        newc = self.__class__()
        newc.toolkit = self.toolkit
        newc.nodes = copy(self.nodes)    
        newc.nodenames = copy(self.nodenames)    
        newc.branches = copy(self.branches)    
        newc.instparams = copy(self.instparams)
        newc.ipar = copy(self.ipar)        
        newc.ipar.detach(self)
        newc.ipar.attach(newc, updatemethod='_ipar_changed')
        newc.iparv = copy(self.iparv)
        if hasattr(newc, 'update'):
            newc.iparv.detach(self)
            newc.iparv.attach(newc)
            newc.update(newc.ipar)
        newc.linear = copy(self.linear)        
        newc.terminals = copy(self.terminals)
        return newc

    def _ipar_changed(self, subject):
        self.update_iparv(ignore_errors=True)

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

        if node is None:
            node_name = None
        else:
            node_name = node.name

        self.nodenames[node_name] = node

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
        
        terminalnode      = self.get_node(terminalname)
        terminal_instname = instname(terminalname)

        result = None
        ## Keep track of branches connected to terminal and 
        ## wether they are inside or outside instanace of terminal
        inside_branches  = [] 
        outside_branches = []
        for branch in self.xflatbranches():
            ## Is branch connected to terminal?
            connected_branch = (branch.plusnode(self) == terminalnode) or \
                               (branch.minusnode(self) == terminalnode)
            if connected_branch:
                ## Is branch in same element as terminal?
                in_same_element = branch.instname.startswith(terminal_instname)

                if in_same_element:
                    inside_branches.append(branch)
                else:
                    outside_branches.append(branch)

        if len(outside_branches) == 1 and outside_branches[0].potential:
            sign = -1
            branch = outside_branches[0]
        elif len(inside_branches ) == 1 and inside_branches[0].potential:
            sign = 1
            branch = inside_branches[0]
        else:
            return None

        if branch.plusnode(self) == terminalnode:
            result = (branch, sign)
        elif branch.minusnode(self) == terminalnode:
            result = (branch, -sign)

        return result

    ## Need to be moved, since circuit doesn't know about x-vector!!
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

    def add_terminals(self, terminals):
        """Add terminals to circuit 

        >>> c = Circuit()
        >>> c.add_terminals(["n1"])
        >>> c.terminals
        ['n1']

        """

        if self.__class__.terminals is self.terminals:
            self.terminals = list(self.terminals)

        for terminal in terminals:
            # add terminal to terminal list if it is not included
            if terminal not in self.terminals:
                self.terminals.append(terminal)

            ## If no node with terminal name exists create node
            if not self.nodenames.has_key(terminal):                
                self.add_node(terminal) 

            node = self.nodenames[terminal]

            ## move node to position k in nodes as it is
            ## now a terminal node
            if not self.nodes.index(node) < self._nterminalnodes:
                self.nodes.remove(node)
                self.nodes.insert(self._nterminalnodes-1, node)
  
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
        >>> import numpy as np
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
    ## Should be moved to NA!!
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

    def next_event(self, t):
        """Returns the time of the next event given the current time t"""
        return inf
    
    def update_iparv(self, parent_ipar=None, globalparams=None,
                     ignore_errors=False):
        """Calculate numeric values of instance parameters"""
        
        substvalues = tuple(p for p in (globalparams, parent_ipar) if p)
            
        newipar = self.ipar.eval_expressions(substvalues, 
                                             ignore_errors=ignore_errors)

        self.iparv.update_values(newipar)

    def eval_iqu(self, inputvalues, epar):
        return ()

    def eval_iqu_and_der(self, inputvalues, epar):
        return self.eval_iqu_and_der_func(self, inputvalues, epar)

    def eval_noise(self, epar):
        return ()

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

    ## Branch related functions
    def xflatbranches(self, instancefilter=None, branchfilter=None, **matchall):
        """Returns branch reference object to all branches in circuit

        The instancefilter is a function that takes an instance name as argument
        and returns True if the instance is to be included.
        The branchfilter argument is a function that takes a branch object
        as argument and returns True if the branch is to be included.
        Finally for simpler branch filtering tasks the function also takes
        keyword arguments with branch properties that all have to match certain 
        values. For example c.xflatbranches(potential=True, input=True) will
        only yield input potential branches. Note, keyword arguments only work 
        if no branchfilter is supplied.
        """
        if branchfilter is None:
            branchfilter = matchbranches(**matchall)

        for branch in self.branches:
            if branchfilter(branch):
                yield BranchRef('', self, branch)

    @property
    def inputbranches(self):
        return [branch for branch in self.branches if branch.input]

    @property
    def _nterminalnodes(self):
        """Return number of terminal nodes"""
        return len(self.terminals)

    @property
    def I(self):
        if len(self.branches) == 1:
            return Quantity('I', self.branches[0])
        else:
            raise Exception('The I property only works for circuits with one branch')

    @property
    def V(self):
        if len(self.branches) == 1:
            return Quantity('V', self.branches[0])
        else:
            raise Exception('The V property only works for circuits with one branch')

class SubCircuit(Circuit):
    """
    SubCircuit is container for circuit instances

    **Attributes**
        *elements* 
          dictionary of Circuit objects keyed by its instance name

        Maybe we can remove this!!
        *elementnodemap* 
          list of translation lists which translate between node indices of the
          elements and node indices of the parent

        *term_node_map* 
          dictionary of instance terminal to node object maps keyed by the 
          instance name. The maps are themselves dictionaries keyed by 
          the terminal names.

    """
    elements = {}
    elementnodemap = {}
    term_node_map = {}

    def __init__(self, *args, **kvargs):
        super(SubCircuit, self).__init__(*args, **kvargs)
        self.elements = {}
        self.elementnodemap = {}
        self.term_node_map = {}
        self._mapmatrix = {}

    def __eq__(self, a):
        return super(SubCircuit, self).__eq__(a) and \
            self.elements == a.elements and \
            self.elementnodemap == a.elementnodemap and \
            self.term_node_map == a.term_node_map

    def __copy__(self):
        newc = super(SubCircuit, self).__copy__()        
        newc.elements = {}
        for instance_name, element in self.elements.items():
            newc.elements[instance_name] = copy(self.elements[instance_name])
        newc.elementnodemap = copy(self.elementnodemap)
        newc.term_node_map = copy(self.term_node_map)
        
        newc.update_node_map()

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

        ## Propagate parent toolkit to instance
        instance.toolkit = self.toolkit

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
                if node is None:
                    node = None
                else:
                    node = Node(str(node))
            
            ## Add node
            self.append_node(node)

            ## Update terminal-node map
            term_node_map[terminal] = node            

        ## Update circuit node - instance map
        self.update_node_map()

        ## update iparv
        instance.update_iparv(self.iparv, ignore_errors=True)

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
                element_node = self[top].get_node('.'.join(path[1:]))
                element_node_index = self[top].nodes.index(element_node)
                node_index = self.elementnodemap[top][element_node_index]
                return self.nodes[node_index]
            else:
                element_node = Node('.'.join(path[1:]))
                element_node_index = self[top].nodes.index(element_node)
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
        
            nodemap = [self.nodes.index(node) for node in element_nodes]

            self.elementnodemap[instance_name] = nodemap

            ## Create mapping matrix
            if len(nodemap) > 0:
                mapmatrix = np.zeros((self.n, len(nodemap)),
                                     dtype = np.integer)
                for inst_node_index, node_index in enumerate(nodemap):
                    mapmatrix[node_index, inst_node_index] = 1
                self._mapmatrix[instance_name] = mapmatrix
            else:
                self._nodemap = None

    def update_iparv(self, parent_ipar=None, globalparams=None, 
                     ignore_errors = False):
        """Calculate numeric values of instance parameters"""
        super(SubCircuit, self).update_iparv(parent_ipar, globalparams,
                                             ignore_errors=ignore_errors)

        ## Update ipar in elements
        for element in self.elements.values():
            element.update_iparv(self.iparv, globalparams,
                                 ignore_errors=ignore_errors)
        
    def save_current(self, terminal):
        """Returns a circuit where the given terminal current is saved
        
        The terminal can be a hierarchical name and the notation is I1.term
        for terminal 'term' of instance 'I1'
        """
        
        if self.get_terminal_branch(terminal) is not None:
            return self

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

    def update(self, subject):
        """This is called when an instance parameter is updated"""
        for element in self.elements.values():
            element.update_iparv(self.iparv, ignore_errors=True)

    def find_class_instances(self, instance_class):
        '''Find all elements of this type

        All resistors in a circuit for example.
        '''
        instances = []        
        for instanceName, element in self.elements.items():
            if isinstance(element, instance_class):
                instances.append(instanceName)
            elif isinstance(element, SubCircuit):
                for instance in element.find_class_instances(instance_class):
                    instances.append(instanceName + '.' + instance)
        return instances
    
    def xflatinstances(self, instancefilter=None, instnameprefix=''):
        """Iterator over all elements and subelements
        
        If instancefilter is given is used to filter out instances by name.
        The parameter should contain a sequence of instance names.
        """
        for instname, e in self.elements.items():
            full_instname = instjoin(instnameprefix, instname)

            ## Filter out instances by name
            if instancefilter is not None and \
               full_instname not in instancefilter:
                continue
            
            if not isinstance(e, SubCircuit):
                yield full_instname, e
            else:
                for subinst in e.xflatinstances(instnameprefix=full_instname):
                    yield subinst

    def get_branch_global_nodes(self, branch, instname=None):
        """Return node indices of plus and minus node of branch
        
        instname is used to get node map of branches of instances
        """
        if instname is not None:
            return (self.get_node(instname + '.' + branch.plus.name),
                    self.get_node(instname + '.' + branch.minus.name))
        else:
            return (branch.plus.name, branch.minus.name)

    def xflatbranches(self, instancefilter=None, branchfilter=None, **matchall):
        """Returns branch reference object to all branches in circuit

        The instancefilter is a function that takes an instance name as argument
        and returns True if the instance is to be included.
        The branchfilter argument is a function that takes a branch object
        as argument and returns True if the branch is to be included.
        Finally for simpler branch filtering tasks the function also takes
        keyword arguments with branch properties that all have to match certain 
        values. For example c.xflatbranches(potential=True, input=True) will
        only yield input potential branches. Note, keyword arguments only work 
        if no branchfilter is supplied.
        """
        if branchfilter is None:
            branchfilter = matchbranches(**matchall)

        for instname, inst in self.xflatinstances(instancefilter=instancefilter):
            for branch in inst.branches:
                if branchfilter(branch):
                    yield BranchRef(instname, inst, branch)

    def get_connections(self, node):
        node = self.get_node(node)
        
        connections = []

        for instname, inst in self.xflatinstances():
            for terminal in inst.terminals:
                fullterminalname = instjoin(instname, terminal)
                if node == self.get_node(fullterminalname):
                    connections.append(fullterminalname)

        return connections
        
    ## Remove or move or change !!
    # def translate_branch(self, branch, instance):
    #     """Return branch from a local branch in an instance"""
    #     return Branch(self.get_node(instance + '.' + branch.plus.name),
    #                   self.get_node(instance + '.' + branch.minus.name))
                    
class ProbeWrapper(SubCircuit): #!!
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
    

class CircuitProxy(Circuit): #!!
    def __init__(self, circuit, parent=None, instance_name=None):
        super(CircuitProxy, self).__init__(self)
        self.device = circuit
        self.terminals = circuit.terminals
        self.nodes = circuit.nodes
        self.nodenames = circuit.nodenames
        self.branches = circuit.branches
        self.iparv = circuit.iparv
        
        ## Find out how this instance was connected to its parent
        ## and set terminalhook accordingly
        if isinstance(parent, SubCircuit) and instance_name != None:
            self.terminalhook = parent.term_node_map[instance_name]

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
    branches = (BranchV('plus', 'minus'),)
    terminals = ['plus', 'minus']

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

        self.quantity = quantity

        self.branch_or_node = branch_or_node

    def __eq__(self, a):
        return self.quantity == a.quantity and \
               self.branch_or_node == a.branch_or_node

    @property
    def is_i(self): return 'I' == self.quantity

    @property
    def is_v(self): return 'V' == self.quantity

    @property
    def isnode(self): return isinstance(self.branch_or_node, Node)

    @property
    def isbranch(self): return isinstance(self.branch_or_node, Branch)
        
    def __repr__(self):
        if isinstance(self.branch_or_node, Branch):
            branch = self.branch_or_node
            if isinstance(branch, BranchRef) and self.quantity == 'I':
                return self.quantity + '(' + branch.instname + ')'
            else:
                return self.quantity + '(' + str(branch.plus.name) + \
                    ',' + str(branch.minus.name) + ')'
        else:
            node = self.branch_or_node
            return self.quantity + '(' + str(node.name) + ')'

def matchbranches(**attr_and_values):
    """Create filter where given values of branch properties all have to match
    """
    def alltrue(branch):
        return all([getattr(branch,k) == v 
                    for k,v in attr_and_values.items()])
    return alltrue
        
if __name__ == "__main__":
    import doctest
    doctest.testmod()
