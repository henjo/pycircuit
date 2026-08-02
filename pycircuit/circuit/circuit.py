# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from pycircuit.utilities.param import Parameter, ParameterDict#, EvalError
from pycircuit.utilities.misc import indent, inplace_add_selected, \
    inplace_add_selected_2d, create_index_vectors
from copy import copy
import contextlib
from .toolkit import numeric
from .toolkit import symbolic
import numpy as np

## The process-wide fallback toolkit used when a circuit is built without an
## explicit toolkit= argument.
##
## Assigning this global directly is deprecated: it is mutable shared state that
## leaks between analyses and tests (a construction-time choice that is easy to
## get wrong).  Prefer passing toolkit= explicitly, or scope it with the
## use_toolkit() context manager below.
default_toolkit = numeric


def set_default_toolkit(toolkit):
    """Set the process-wide default toolkit.

    Deprecated in favour of passing ``toolkit=`` explicitly or using
    :func:`use_toolkit`; provided so the global is only mutated through a
    documented entry point.
    """
    global default_toolkit
    default_toolkit = toolkit


@contextlib.contextmanager
def use_toolkit(toolkit):
    """Temporarily select the default toolkit within a ``with`` block.

    Restores the previous value on exit, so the choice does not leak to later
    code.  This is the recommended replacement for assigning
    ``circuit.default_toolkit`` directly::

        with use_toolkit(symbolic):
            cir = SubCircuit()
            ...
    """
    global default_toolkit
    previous = default_toolkit
    default_toolkit = toolkit
    try:
        yield toolkit
    finally:
        default_toolkit = previous


timedomain_analyses = ('dc', 'tran')

class Node():
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

class Branch():
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
    Parameter("T", "Temperature", unit="K", default = 300),
    ## Present so devices have a defined value even when evaluated outside an
    ## analysis.  `Analysis.__init__` overwrites it on its own epar copy; devices
    ## previously fell through to a hard-coded 1e-12 in each of three places when
    ## the attribute was missing, which is how `bypass=True` came to do nothing at
    ## all in a transient (the transient never passed its epar, so every device saw
    ## this dict and took the missing-attribute branch).
    ##
    ## -1.0 means "never bypass", matching the `bypass=False` default: the device
    ## test is `abs(vnew - vold) < bypasstol`, which no non-negative difference can
    ## satisfy against a negative tolerance.  The old implicit 1e-12 meant a small
    ## amount of bypassing always happened, unasked.
    Parameter("bypasstol", "Device bypass tolerance", unit="V", default = -1.0))

class Circuit():
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

    
    nodes = ()
    branches = ()
    terminals = ()
    instparams = ()
    linear = True
    
    def __init__(self, *args, **kvargs):
        if 'toolkit' in kvargs:
            self.toolkit = kvargs['toolkit']
            del kvargs['toolkit']
        else:
            self.toolkit = default_toolkit

        if self.nodes is self.__class__.nodes:
            self.nodes = list(self.__class__.nodes)
        if self.branches is self.__class__.branches:
            self.branches = list(self.__class__.branches)
        if self.terminals is self.__class__.terminals:
            self.terminals = list(self.__class__.terminals)
        if self.instparams is self.__class__.instparams:
            self.instparams = list(self.__class__.instparams)

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
            
    def __eq__(self, a):
        return self.__class__ == a.__class__ and \
            self.nodes == a.nodes and \
            self.nodenames == a.nodenames and self.branches == a.branches and \
            self.iparv == a.iparv
        
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
        if node not in self.nodes:
            self.nodes.append(node)
        self.nodenames[node.name] = node

    def append_branches(self, *branches):
        """Append node object to circuit"""
        ## Make a copy of branch list so the class is unchanged
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
        
        >>> from pycircuit.circuit.elements import *
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

    def get_node_index(self, node, refnode=None):
        """Get row in the x vector of a node instance

           If the refnode argument is given the reference node
           is assumed to be removed
        """

        if not isinstance(node, Node):
            node = Node(str(node))
        if refnode and not isinstance(refnode, Node):
            refnode = Node(str(refnode))

        if node in self.nodes:
            index = self.nodes.index(node)
            if refnode is not None:
                irefnode = self.nodes.index(refnode)
                if index == irefnode:
                    return None
                if index > irefnode:
                    return index - 1
            return index
        else:
            raise ValueError('Node %s is not in circuit node list (%s)'%
                             (str(node), str(self.nodes)))

    def instance_branch_indices(self, instancename):
        """Rows of the solution vector holding ``instancename``'s branch currents.

        STAGE 10.3.  Returns a list, empty for an element that declares no
        branches.  Raises ``KeyError`` for an unknown instance rather than
        returning an empty list, because "this element has no branches" and "there
        is no such element" are different answers and only one of them is a typo.

        Resolved through the recorded span rather than by searching
        ``self.branches``: that search is ambiguous for parallel elements, whose
        branches compare equal.
        """
        start, stop = self._instance_branch_span[instancename]
        offset = len(self.nodes)
        return [offset + i for i in range(start, stop)]

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

        for terminal in terminals:
            # add terminal to terminal list if it is not included
            if terminal not in self.terminals:
                self.terminals.append(terminal)

            ## If no node with terminal name exists create node
            if terminal not in self.nodenames:                
                self.add_node(terminal) 

            node = self.nodenames[terminal]

            ## move node to position k in nodes as it is
            ## now a terminal node
            if not self.nodes.index(node) < self._nterminalnodes:
                self.nodes.remove(node)
                self.nodes.insert(self._nterminalnodes-1, node)
                
    def accept_step(self, t, x, epar):
        """Called by the transient solver when a time step is accepted.
        This allows elements with state history (e.g. T-Lines) to update their internal history buffers.
        """
        pass

    def reset_state(self, epar=None):
        """Discard any state carried over from a previous analysis.

        STAGE 8(d).  Elements that accumulate history across a run -- `TLine` is
        the only one today -- kept it forever, so a second analysis on the same
        circuit object started from the first one's leftovers.  Measured: a `TLine`
        DC solve gave **v(b) = 0.500000 before** a transient and **0.000000 after**,
        and a second transient began with 12 stale history points and ended with
        73.

        Called at the start of every analysis, so which analysis ran previously
        cannot change the answer.  A no-op for stateless elements, which is all of
        them but one.
        """
        pass

    def max_timestep(self):
        """The largest timestep this element can be integrated with, or ``None``.

        STAGE 8(d).  A delay element is only meaningful if the solver samples
        finely enough to resolve its delay: `TLine` interpolates its history at
        ``t - TD``, and with ``dt`` comparable to ``TD`` there is nothing to
        interpolate between.  Measured under `fixed_timestep`, TD = 1e-9: the
        observed propagation delay came out **2.00x TD at dt = 1e-9, 4.00x at
        2e-9 and 8.00x at 5e-9** -- silently.

        ``None`` means "no opinion", which is every element except `TLine`.
        """
        return None
  
    def connect_terminals(self, **kvargs):
        """Connect nodes to terminals by using keyword arguments

        """
        for terminal, node in kvargs.items():
            ## Sanity check
            if not isinstance(terminal, str):
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
        
        >>> from pycircuit.circuit.elements import *
        >>> import numpy as np
        >>> cir = R(Node('n1'), gnd, r=1e3)
        >>> newcir = cir.save_current('plus')
        >>> newcir.G(np.zeros(4))
        array([[ 0.   ,  0.   ,  0.   ,  1.   ],
               [ 0.   ,  0.001, -0.001,  0.   ],
               [ 0.   , -0.001,  0.001, -1.   ],
               [ 1.   ,  0.   , -1.   ,  0.   ]])
        """
        
        if self.get_terminal_branch(terminal) is None:
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
        if instancename is None:
            return self.nodes[self._nterminalnodes:]
        else:
            result = []
            for node in self.nodes[len(self.terminals):]:
                if node.isglobal:
                    result.append(node)
                else:
                    result.append(Node(instancename + '.' + node.name))
            return result

    def G(self, x, epar=defaultepar, params_tree=None):
        """Calculate the G (trans)conductance matrix given the x-vector"""
        return self.toolkit.zeros((self.n, self.n))

    def C(self, x, epar=defaultepar, params_tree=None):
        """Calculate the C (transcapacitance) matrix given the x-vector"""
        return self.toolkit.zeros((self.n, self.n))

    def u(self, t=0.0, epar=defaultepar, analysis=None, params_tree=None):
        return self.toolkit.zeros(self.n)

    def dudt(self, t=0.0, epar=defaultepar, analysis=None, params_tree=None):
        """Time derivative of the source vector, for Fang's ``p = df_ckt/dh``.

        STAGE 12B.  With the step size an unknown (DAC 2013 eq 11) the residual
        is differentiated with respect to it, and the sources are evaluated at
        ``t_{m-1} + h`` -- so ``du/dt`` enters ``p`` directly, and on a driven
        circuit it is usually the larger of the two contributions.

        Zero here is not a placeholder: it is correct for every element whose
        contribution does not depend on time, which is most of them. Elements
        carrying a `TimeFunction` override it; an element whose time dependence
        cannot be differentiated should RAISE rather than inherit this, so a
        coupled run fails loudly instead of solving a subtly different problem.
        """
        return self.toolkit.zeros(self.n)

    def i(self, x, epar=defaultepar, params_tree=None):
        """Calculate the i vector as a function of the x-vector"""
        return self.toolkit.dot(self.G(x), x)

    def q(self, x, epar=defaultepar, params_tree=None):
        """Calculate the q vector as a function of the x-vector"""
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
        return self.toolkit.inf
    
    def name_state_vector(self, x, analysis=''):
        """Return a dictionary of the x-vector keyed by node and branch names

        >>> from pycircuit.circuit.elements import *
        >>> import numpy as np
        >>> c = SubCircuit()
        >>> n1 = c.add_node('net1')
        >>> c['vs'] = VS(n1, gnd, v=1.0)
        >>> c['R'] = R(n1, gnd, r=1e3)
        >>> c.name_state_vector(np.array([1.0, 0.0, -1e-3]))
        {'net1': np.float64(1.0), 'gnd': np.float64(0.0), 'i0)': np.float64(-0.001)}

        """
        result = {}
        for xvalue, node in zip(x[:len(self.nodes)], self.nodes):
            result[self.get_node_name(node)] = xvalue

        nnodes = len(self.nodes)
        for i, (xvalue, branch) in enumerate(zip(x[nnodes:], self.branches)):
            result['i' + analysis + str(i) + ')'] = xvalue

        return result

    def stamp_v(self, x, value, nodep, noden=None, refnode=None):
        """Stamp value in vector such that x[nodep] += value, x[noden] -= value
        
        If refnode is not None the reference node is assumed to be removed
        from vector
        """
        x[self.get_node_index(nodep, refnode)] += value
        x[self.get_node_index(noden, refnode)] -= value

    def remove_refnode(self, matrices, refnode):
        """Remove refnode from vectors or matrices"""
        n = self.get_node_index(refnode)
        result = []
        
        for A in matrices:
            for axis in range(len(A.shape)):
                A=self.toolkit.delete(A, [n], axis=axis)
            result.append(A)
        return tuple(result)

    def check_kcl(self, x, t=0.0, iq=None, abstol=1e-12):
        """Manually evaluate KCL at each node by summing the currents and ensuring 
        the residual is within < abstol.
        """
        f = self.i(x) + self.u(t)
        if iq is not None:
            f += iq
        # Check only the node rows (len(self.nodes))
        nnodes = len(self.nodes)
        node_residuals = f[:nnodes]
        import numpy as np
        max_residual = np.max(np.abs(node_residuals))
        if max_residual > abstol:
            raise ValueError("KCL check failed! Max residual {} > {}".format(max_residual, abstol))
        return max_residual

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

        >>> from pycircuit.circuit.elements import *
        >>> import numpy as np
        >>> c = SubCircuit()
        >>> n1, n2 = c.add_nodes('n1','n2')
        >>> c['R1'] = R(n1, n2, r=1e3)
        >>> c['R2'] = R(n2, gnd, r=1e3)
        >>> c.extract_v(np.array([1.0, 0.5, 0.0]), 'n1', 'n2')
        np.float64(0.5)
        >>> c.extract_v(np.array([1.0, 0.5, 0.0]), c.nodes[0])
        np.float64(1.0)
        >>> c.extract_v(np.array([1.0, 0.5]), c.nodes[0], refnode_removed = True)
        np.float64(1.0)
        
        """
        v = []
        for node in nodep, noden:
            if type(node) is str:
                node = self.get_node(node)
            elif node is None:
                node = refnode

            if refnode_removed:
                nodeindex = self.get_node_index(node, refnode)
            else:
                nodeindex = self.get_node_index(node, None)

            if nodeindex is None: ## When node == refnode
                v.append(0)
                continue
                    
            v.append(x[nodeindex])

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
           dx/dt vector. this is needed if dx/dt is non-zero and there is no branch defined at the
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

        >>> from pycircuit.circuit.elements import *
        >>> import numpy as np
        >>> c = SubCircuit()
        >>> net1 = c.add_node('net1')
        >>> c['vs'] = VS(net1, gnd)
        >>> c.extract_i(np.array([1.0, 0, -1e-3]), 'vs.minus')
        np.float64(0.001)
        >>> c.extract_i(np.array([1.0, -1e-3]), 'vs.minus', refnode_removed = True)
        np.float64(0.001)
        
        """
        dot = self.toolkit.dot
        
        if type(branch_or_term) is str:
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

            if branch_sign is not None:
                branch, sign = branch_sign
            else:
                terminal_node = self.nodenames[branch_or_term]
                
                terminal_node_index = self.get_node_index(terminal_node)

                if xdot is not None:
                    if linearized:
                        return dot(self.G(xdcop)[terminal_node_index], x) + \
                            dot(self.C(xdcop)[terminal_node_index], xdot) + \
                            self.u(t, analysis = 'ac')[terminal_node_index]

                    else:
                        return self.i(x)[terminal_node_index] + \
                            dot(self.C(x)[terminal_node_index], xdot) + \
                            self.u(t)[terminal_node_index]
                else:
                    if linearized:
                        return dot(self.G(xdcop)[terminal_node_index], x) + \
                            self.u(t, analysis = 'ac')[terminal_node_index]

                    else:
                        return self.i(x)[terminal_node_index] + \
                            self.u(t)[terminal_node_index]

        else:
            branch = branch_or_term
            sign = 1

        branchindex = self.get_branch_index(branch)

        if refnode_removed:
            branchindex -= 1

        return sign * x[branchindex]      

    def update_iparv(self, parent_ipar=None, globalparams=None,
                     ignore_errors=False):
        """Calculate numeric values of instance parameters"""
        
        substvalues = tuple(p for p in (globalparams, parent_ipar) if p)
            
        newipar = self.ipar.eval_expressions(substvalues, 
                                             ignore_errors=ignore_errors)

        self.iparv.update_values(newipar)

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
        if instancebranches is None:
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
    elements = {}
    term_node_map = {}

    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        self.elements = {}
        ## STAGE 10.3 -- which slice of `self.branches` each instance owns.
        ##
        ## `add_instance` appends an element's branches contiguously, so the span
        ## is exact and O(1) to look up.  Recorded rather than reconstructed
        ## because reconstruction is silently WRONG for parallel elements:
        ## `Branch.__eq__` compares node pairs, so two inductors between the same
        ## two nodes produce EQUAL branches and `branches.index()` returns the
        ## first for both.  Measured -- an initial current given to the second of
        ## two parallel inductors landed on the first one's unknown, with nothing
        ## to indicate it.
        ##
        ## The span is in BRANCH-LIST coordinates, not solution-vector ones: a
        ## branch's row is `len(self.nodes) + offset`, and nodes are still being
        ## added while elements are, so storing the row directly would go stale.
        self._instance_branch_span = {}
        self._elementnodemap = {}
        self._rep_nodemap_list = {}
        self._map_indices_1d = {}
        self._map_indices_2d = {}
        ## STAGE 2+.5 -- rebuild the node map LAZILY.
        ##
        ## `add_instance` used to call `update_node_map()`, which rebuilds the map
        ## for EVERY element.  N insertions therefore did O(N^2) element-work
        ## before a single analysis ran: building an 800-element ladder took 24.8 s,
        ## and a 1600-element one could not be built inside a ten-minute budget --
        ## which is what stopped stage 7b measuring its n=2000 case at all.
        ##
        ## Deferring is safe in a way that incremental appending is NOT: a branch
        ## index is `branch_offset + len(self.nodes)`, so adding any node later
        ## shifts every branch index already computed.  Appending one element's
        ## entry would silently leave the earlier ones stale; rebuilding once, when
        ## the map is first read, cannot.
        self._nodemap_dirty = True
        self.term_node_map = {}
        self._mapmatrix = {}

    def _ensure_node_map(self):
        if self._nodemap_dirty:
            self.update_node_map()

    @property
    def _map_indices_1d(self):
        self._ensure_node_map()
        return self.__map_indices_1d

    @_map_indices_1d.setter
    def _map_indices_1d(self, value):
        self.__map_indices_1d = value

    @property
    def _map_indices_2d(self):
        self._ensure_node_map()
        return self.__map_indices_2d

    @_map_indices_2d.setter
    def _map_indices_2d(self, value):
        self.__map_indices_2d = value

    @property
    def elementnodemap(self):
        self._ensure_node_map()
        return self._elementnodemap

    @elementnodemap.setter
    def elementnodemap(self, value):
        ## An explicit assignment is taken as authoritative -- `__copy__` does this.
        self._elementnodemap = value
        self._nodemap_dirty = False

    def __eq__(self, a):
        ## elementnodemap values are numpy arrays, so a plain dict ``==``
        ## raises "truth value of an array is ambiguous"; compare element-wise.
        def _nodemap_eq(m1, m2):
            if m1.keys() != m2.keys():
                return False
            return all(np.array_equal(m1[k], m2[k]) for k in m1)
        return super().__eq__(a) and \
            self.elements == a.elements and \
            _nodemap_eq(self.elementnodemap, a.elementnodemap) and \
            self.term_node_map == a.term_node_map

    def __copy__(self):
        newc = super().__copy__()        
        newc.elements = {}
        for instance_name, element in self.elements.items():
            newc.elements[instance_name] = copy(self.elements[instance_name])
        newc.elementnodemap = copy(self.elementnodemap)
        newc.term_node_map = copy(self.term_node_map)
        
        newc.update_node_map()

        return newc

    def netlist(self, top = True):
        """
        >>> from pycircuit.circuit.elements import *
        >>> a = SubCircuit()
        >>> a['R1'] = R(1,2)
        >>> print(a.netlist())
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
                node = Node(str(node))
            
            ## Add node
            self.append_node(node)

            ## Update terminal-node map
            term_node_map[terminal] = node            

        ## Add branches
        branch_start = len(self.branches)
        newbranches = self._instance_branches(instance, instancename)
        self.append_branches(*newbranches)
        self._instance_branch_span[instancename] = (branch_start,
                                                    len(self.branches))

        ## Update circuit node - instance map.  STAGE 2+.5: mark, do not rebuild --
        ## the rebuild happens once, when the map is first read.
        self._nodemap_dirty = True

        ## update iparv
        instance.update_iparv(self.iparv, ignore_errors=True)

    def __setitem__(self, instancename, element):
        """Adds an instance to the circuit"""

        self.add_instance(instancename, element, **element.terminalhook)

        element.terminalhook = None


    def __delitem__(self, instancename):
        """Removes instance from circuit
        
        >>> from pycircuit.circuit.elements import *
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

        >>> from pycircuit.circuit.elements import *
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

               return next(branch_gen), branch_sign[1]
        
    def get_node(self, name):
        """Find a node by name.
        
        >>> from pycircuit.circuit.elements import *
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
                node = self[top].get_node('.'.join(path[1:]))
                index = self[top].get_node_index(node)
                index2 = self.elementnodemap[top][index]
                node2 = self.nodes[index2]
                return node2.name
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
        if node.name is not None:
            return node.name
        
        ## First search among the local nodes
        name = Circuit.get_node_name(self, node)
        if name is not None:
            return name
        
        ## Then search in the circuit elements
        for instname, element in self.elements.items():
            name =  element.get_node_name(node)
            if name is not None:
                return instname + '.' + name
        
    def update_node_map(self):
        """Update the elementnodemap attribute"""

        ## STAGE 2+.5 -- ONE DICT INSTEAD OF A LINEAR SCAN PER NODE PER ELEMENT.
        ##
        ## This method looked every node up with `self.nodes.index(node)`, an O(n)
        ## scan whose every step calls `Node.__eq__`.  Building an 800-element
        ## ladder ran `list.index` 642,400 times and `Node.__eq__` **129,042,200**
        ## times, for 27.6 s and 17.1 s respectively.  `Node` already defines
        ## `__hash__` as `hash(self.name)`, consistent with its name-based
        ## `__eq__`, so the nodes were hashable the whole time.
        ##
        ## FIRST OCCURRENCE WINS, deliberately: `list.index` returns the first
        ## match, and a dict comprehension would keep the LAST.  If two `Node`
        ## objects in `self.nodes` ever compare equal, the naive comprehension
        ## would silently renumber matrix rows -- which is exactly what gate
        ## 2+.5-1 exists to catch.
        node_index = {}
        for _i, _node in enumerate(self.nodes):
            if _node not in node_index:
                node_index[_node] = _i

        self._elementnodemap = {}
        self._rep_nodemap_list = {}
        self._map_indices_1d = {}
        self._map_indices_2d = {}
        ## Cleared FIRST: the rebuild below reads `self.nodes`, not the map, so
        ## nothing here can re-enter `_ensure_node_map`.
        self._nodemap_dirty = False
        
        branch_offset = 0
        
        for instance_name, element in self.elements.items():
            nodemap = self.term_node_map[instance_name]
            element_nodes = [nodemap[terminal] for terminal in element.terminals]

            for node in element.non_terminal_nodes(instance_name):
                element_nodes.append(node)
        
            element_branches = list(self._instance_branches(element, instance_name))

            branch_indices = []
            for _ in element_branches:
                branch_indices.append(branch_offset + len(self.nodes))
                branch_offset += 1

            nodemap = \
                [node_index[node] for node in element_nodes] + branch_indices

            import numpy as np
            nodemap = np.array(nodemap, dtype=int)
            self._elementnodemap[instance_name] = nodemap

            ## Create mapping coordinates instead of dense matrices
            if len(nodemap) > 0:
                # For 1D vector stamping (e.g. u, i, q)
                self._map_indices_1d[instance_name] = nodemap
                
                # For 2D matrix stamping (e.g. G, C)
                # We need all (row, col) combinations
                # meshgrid(nodemap, nodemap, indexing='ij') creates the grid
                rows, cols = np.meshgrid(nodemap, nodemap, indexing='ij')
                self._map_indices_2d[instance_name] = (rows.flatten(), cols.flatten())
            else:
                self._nodemap = None

    def update_iparv(self, parent_ipar=None, globalparams=None, 
                     ignore_errors = False):
        """Calculate numeric values of instance parameters"""
        super().update_iparv(parent_ipar, globalparams,
                                             ignore_errors=ignore_errors)

        ## Update ipar in elements
        for element in self.elements.values():
            element.update_iparv(self.iparv, globalparams,
                                 ignore_errors=ignore_errors)
        
        self._build_evaluation_groups()

    def _build_evaluation_groups(self):
        """Ask the toolkit which element classes it will batch.

        The grouping is an optimisation only a differentiating toolkit can
        use, so the toolkit decides and the base returns nothing -- this
        module needs no knowledge of JAX to offer the opportunity.
        """
        self._eval_groups = self.toolkit.evaluation_groups(self)
        
    def G(self, x, epar=defaultepar, params_tree=None):
        return self._add_element_submatrices('G', x, (epar,), params_tree=params_tree)

    def C(self, x, epar=defaultepar, params_tree=None):
        return self._add_element_submatrices('C', x, (epar,), params_tree=params_tree)

    def u(self, t=0.0, epar=defaultepar, analysis=None, params_tree=None):
        dtype = None
        if analysis == 'ac':
            dtype = self.toolkit.ac_u_dtype

        return self._add_element_subvectors('u', None, (t,epar,analysis), params_tree=params_tree, 
                                            dtype=dtype)

    def dudt(self, t=0.0, epar=defaultepar, analysis=None, params_tree=None):
        ## No `ac` dtype branch: this exists for the transient coupled solve, and
        ## a complex `du/dt` would be meaningless there.
        return self._add_element_subvectors('dudt', None, (t, epar, analysis),
                                            params_tree=params_tree)

    def i(self, x, epar=defaultepar, params_tree=None):
        return self._add_element_subvectors('i', x, (epar,), params_tree=params_tree)

    #This seemed to be missing
    def q(self, x, epar=defaultepar, params_tree=None):
        return self._add_element_subvectors('q', x, (epar,), params_tree=params_tree)

    def limit(self, x, x0, epar=defaultepar):
        """Apply each element's Newton limiter and WRITE THE RESULT BACK.

        STAGE 5(c).  This used to read

            subx = x[self.elementnodemap[instance]]
            element.limit(subx, subx0, epar)

        and return `x` unmodified.  `subx` is fancy indexing, therefore a **copy**,
        so a limiter that writes its argument wrote into a temporary that was then
        discarded -- limiting was a no-op for any device that did not keep private
        state.  `Diode` kept `_vlim` and linearised `G` around it, so the one
        limiter in the tree was the one that could not expose the bug.

        Both conventions are accepted while the tree has both: a limiter may return
        the limited sub-vector (the state-free form, which is what stage 5 adds to
        `Semiconductor` and what a traced backend could ever support), or return
        `None` and keep its own state (`Diode`).
        """
        ## HOISTED -- `elementnodemap` is a property since 2+.5, so reading it
        ## inside the loop is a Python call per element per Newton iteration.  The
        ## stamping paths below already hoist; these two did not, and the suite
        ## went from ~8 to ~16 minutes until they did.
        elementnodemap = self.elementnodemap
        for instance, element in self.elements.items():
            if hasattr(element, 'limit'):
                nodemap = elementnodemap[instance]
                subx = x[nodemap]
                subx0 = x0[nodemap]
                limited = element.limit(subx, subx0, epar)
                if limited is not None:
                    ## The write-back the old code was missing.  Fancy-index
                    ## assignment does reach `x`, unlike the read above.
                    x[nodemap] = limited
        return x

    def CY(self, x, w, epar=defaultepar):
        """Calculate composite noise source correlation matrix

        The noise sources in one element are assumed to be uncorrelated 
        with the noise sources in the other elements.

        """
        return self._add_element_submatrices('CY', x, (w, epar,))
    def next_event(self, t):
        """Returns the time of the next event given the current time t
        by polling all elements in the subcircuit."""
        events = [element.next_event(t) for element in self.elements.values() 
                  if hasattr(element, 'next_event')]
        if events:
            return self.toolkit.maximum(t, min(events)) # Ensure we don't go backwards
        return self.toolkit.inf

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
        if type(branch_or_term) is str:
            if self.get_terminal_branch(branch_or_term) is None:

                hierlevels = [part for part in branch_or_term.split('.')]

                if len(hierlevels) > 1:
                    instance_name = hierlevels[0]
                    instance = self.elements[instance_name]
                    terminal_name = '.'.join(hierlevels[1:])

                    ## Get slice of x-vector for the instance
                    nodemap = self.elementnodemap[instance_name]

                    subx = x[nodemap]

                    if xdot is not None:
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
        
    def accept_step(self, t, x, epar):
        """Propagate accept_step to all child elements"""
        ## Hoisted for the same reason as in `limit` -- see the note there.
        elementnodemap = self.elementnodemap
        for instance, element in self.elements.items():
            if hasattr(element, 'accept_step'):
                subx = x[elementnodemap[instance]]
                element.accept_step(t, subx, epar)

    def reset_state(self, epar=None):
        """Propagate reset_state to all child elements -- see Circuit.reset_state."""
        for element in self.elements.values():
            if hasattr(element, 'reset_state'):
                element.reset_state(epar)

    def max_timestep(self):
        """The tightest cap any child element asks for, or ``None``.

        The MINIMUM, not the first found: two delay lines with different ``TD``
        both have to be resolved, and the shorter one governs.
        """
        caps = [c for c in (element.max_timestep()
                            for element in self.elements.values()
                            if hasattr(element, 'max_timestep'))
                if c is not None]
        return min(caps) if caps else None
            
    def update(self, subject):
        """This is called when an instance parameter is updated"""
        for element in self.elements.values():
            element.update_iparv(self.iparv, ignore_errors=True)
        
    def _add_element_submatrices(self, methodname, x, args, params_tree=None):
        ## STAGE 2b -- see the note on `_add_element_subvectors`.  Same two changes:
        ## the per-element `hasattr`/`getattr` probes are hoisted, and the scatter is
        ## done once at the end instead of once per element.
        import numpy as np

        n = self.n
        toolkit = self.toolkit

        # Check if toolkit prefers sparse assembly directly
        build_sparse = hasattr(toolkit, 'build_sparse')
        has_add_at = hasattr(toolkit, 'add_at')
        groups = getattr(self, '_eval_groups', None)
        idxmap = self._map_indices_2d
        elementnodemap = self.elementnodemap

        if build_sparse:
            all_data, all_rows, all_cols = [], [], []
        else:
            lhs = toolkit.zeros((n,n))

        ## Elements the toolkit can evaluate in bulk are stamped in one go; the
        ## loop below then skips them.  A toolkit that cannot batch returns
        ## None and everything goes through the loop.
        batched = toolkit.batched_contributions(
            self, groups, methodname, x, args, params_tree, ndim=2)
        if batched is not None:
            values, (rows, cols) = batched
            if build_sparse:
                all_data.append(values)
                all_rows.append(rows)
                all_cols.append(cols)
            else:
                lhs = toolkit.add_at(lhs, (rows, cols), values)

        pending_rc = []
        pending_val = []

        for instance, element in self.elements.items():
            if groups and element.__class__ in groups:
                continue # Handled by vectorization above

            nodemap = elementnodemap[instance]

            if x is not None:
                subx = x[nodemap]
                try:
                    rhs = getattr(element, methodname)(subx, *args)
                except Exception as e:
                    raise e.__class__(str(e) + ' at element ' + str(element)
                                      + ', args='+str(args))
            else:
                rhs = getattr(element, methodname)(*((None,) + tuple(args)))

            rc = idxmap.get(instance)
            if rc is None:
                continue
            rows, cols = rc

            if build_sparse:
                all_data.append(np.asarray(rhs).flatten())
                all_rows.append(rows)
                all_cols.append(cols)
            elif has_add_at:
                rhs_flat = toolkit.reshape(rhs, (-1,)).flatten()
                lhs = toolkit.add_at(lhs, (rows, cols), rhs_flat)
            else:
                pending_rc.append((rows, cols))
                pending_val.append(np.asarray(rhs).flatten())

        if pending_rc:
            lhs = self._scatter_2d(lhs, pending_rc, pending_val, n)

        if build_sparse:
            import numpy as np
            if not all_data:
                return self.toolkit.build_sparse([], [], [], shape=(n,n))
            return self.toolkit.build_sparse(
                np.concatenate(all_data),
                np.concatenate(all_rows),
                np.concatenate(all_cols),
                shape=(n,n)
            )
        return lhs

    def _add_element_subvectors(self, methodname, x, args, dtype=None, params_tree=None):
        ## STAGE 2b.  Two changes, both mechanical, neither touching the arithmetic:
        ##
        ## 1. HOIST THE PER-ELEMENT PROBES.  `hasattr(self.toolkit, 'add_at')` and
        ##    `getattr(self, '_eval_groups', None)` were evaluated once per element
        ##    per stamp.  `NumericToolkit` has no `add_at`, so every one of those
        ##    `hasattr` calls went through `Toolkit.__getattr__`, which formats an
        ##    error message and raises -- hundreds of thousands of times per run,
        ##    purely to answer a question whose answer cannot change during a solve.
        ##
        ## 2. SCATTER ONCE.  `np.add.at` is a slow unbuffered ufunc path; calling it
        ##    per element repeats that cost per element.  Collecting the indices and
        ##    values and accumulating once with `np.bincount` is the same sum in the
        ##    same order (bincount accumulates in index order, and duplicate indices
        ##    are summed identically), and it is a single pass in C.
        ##
        ## `bincount` only handles real floating point, so complex and object dtypes
        ## keep the `np.add.at` path.  That is not a fallback for correctness -- both
        ## paths are exact -- it is a dtype restriction of `bincount`.
        import numpy as np

        n = self.n
        toolkit = self.toolkit
        lhs = toolkit.zeros(n, dtype=dtype)
        has_add_at = hasattr(toolkit, 'add_at')
        groups = getattr(self, '_eval_groups', None)
        idxmap = self._map_indices_1d
        elementnodemap = self.elementnodemap

        batched = toolkit.batched_contributions(
            self, groups, methodname, x, args, params_tree, ndim=1)
        if batched is not None:
            values, indices = batched
            lhs = toolkit.add_at(lhs, indices, values)

        pending_idx = []
        pending_val = []

        for instance, element in self.elements.items():
            if groups and element.__class__ in groups and x is not None:
                continue # Handled by vectorization above

            if x is not None:
                subx = x[elementnodemap[instance]]
                rhs = getattr(element, methodname)(subx, *args)
            else:
                rhs = getattr(element, methodname)(*args)

            indices = idxmap.get(instance)
            if indices is None:
                continue

            if has_add_at:
                rhs_flat = toolkit.reshape(rhs, (-1,)).flatten()
                lhs = toolkit.add_at(lhs, indices, rhs_flat)
            else:
                pending_idx.append(indices)
                pending_val.append(np.asarray(rhs).flatten())

        if pending_idx:
            lhs = self._scatter_1d(lhs, pending_idx, pending_val, n)

        return lhs

    @staticmethod
    def _scatter_2d(lhs, pending_rc, pending_val, n):
        """Accumulate collected ((rows, cols), value) pairs into ``lhs`` in one pass.

        ``np.add.at(lhs, (rows, cols), v)`` on a 2-D array is exactly
        ``bincount(rows*n + cols, weights=v)`` reshaped: both walk the index array
        in order and sum duplicates as they are met, so the floating-point result is
        identical, not merely equal to within rounding.

        This is only reached when the toolkit has no ``add_at``, and such a toolkit
        also returns ``None`` from ``batched_contributions`` -- so ``lhs`` is still
        the zero matrix here.  That is what makes the single ``lhs +=`` at the end
        exact rather than an accumulation in a different order.
        """
        import numpy as np

        rows = np.concatenate([rc[0] for rc in pending_rc])
        cols = np.concatenate([rc[1] for rc in pending_rc])
        val = np.concatenate(pending_val)

        if val.dtype == object or lhs.dtype == object:
            if lhs.dtype != object:
                lhs = lhs.astype(object)
            np.add.at(lhs, (rows, cols), val)
            return lhs
        if np.iscomplexobj(val) or np.iscomplexobj(lhs):
            if not np.iscomplexobj(lhs):
                lhs = lhs.astype(complex)
            np.add.at(lhs, (rows, cols), val)
            return lhs

        flat = np.asarray(rows, dtype=np.intp) * n + np.asarray(cols, dtype=np.intp)
        lhs += np.bincount(flat, weights=val, minlength=n * n)[:n * n].reshape(n, n)
        return lhs

    @staticmethod
    def _scatter_1d(lhs, pending_idx, pending_val, n):
        """Accumulate collected (index, value) pairs into ``lhs`` in one pass."""
        import numpy as np

        idx = np.concatenate(pending_idx)
        val = np.concatenate(pending_val)

        ## `bincount` is real-only.  Object dtype additionally has to widen `lhs`,
        ## exactly as the per-element path did.
        if val.dtype == object or lhs.dtype == object:
            if lhs.dtype != object:
                lhs = lhs.astype(object)
            np.add.at(lhs, idx, val)
            return lhs
        if np.iscomplexobj(val) or np.iscomplexobj(lhs):
            if not np.iscomplexobj(lhs):
                lhs = lhs.astype(complex)
            np.add.at(lhs, idx, val)
            return lhs

        lhs += np.bincount(idx, weights=val, minlength=n)[:n]
        return lhs

    def find_class_instances(self, instance_class):
        instances = []        
        for instanceName, element in self.elements.items():
            if isinstance(element, instance_class):
                instances.append(instanceName)
            elif isinstance(element, SubCircuit):
                for instance in element.find_class_instances(instance_class):
                    instances.append(instanceName + '.' + instance)
        return instances

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
        super().__init__(*terminalnodes)
        
        self.add_instance('wrapped', circuit, 
                          **dict(zip(circuit.terminals, terminalnodes)))
        
        for terminal in terminals:
            self.save_current(terminal)

    def save_current(self, terminal):
        """Returns a circuit where the given terminal current is saved"""
        
        ## Add probe to terminal if it does not already exists
        if self.get_terminal_branch(terminal) is None:
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
        super().__init__(self)
        self.device = circuit
        self.terminals = circuit.terminals
        self.nodes = circuit.nodes
        self.nodenames = circuit.nodenames
        self.branches = circuit.branches
        self.iparv = circuit.iparv
        
        ## Find out how this instance was connected to its parent
        ## and set terminalhook accordingly
        if isinstance(parent, SubCircuit) and instance_name is not None:
            self.terminalhook = parent.term_node_map[instance_name]

    def G(self, x, epar=defaultepar, params_tree=None): return self.device.G(x,epar)
    def C(self, x, epar=defaultepar, params_tree=None): return self.device.C(x,epar)
    def u(self, t=0.0, epar=defaultepar, analysis=None, params_tree=None): 
        return self.device.u(x,epar)
    def i(self, x, epar=defaultepar, params_tree=None): return self.device.i(x,epar)
    def q(self, x, epar=defaultepar, params_tree=None): return self.device.q(x,epar)
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

    def G(self, x, epar=defaultepar, params_tree=None):
        return self.toolkit.array([[0 , 0, 1],
                                   [0 , 0, -1],
                                   [1 , -1, 0]])

    @property
    def branch(self):
        """Return the branch (plus, minus)"""
        return self.branches[0]

class Quantity():
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
