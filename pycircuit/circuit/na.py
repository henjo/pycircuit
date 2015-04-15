import scipy.sparse as sp
import numpy as np
from itertools import groupby
import types
from copy import copy 

import circuit
from circuit import Quantity, defaultepar
import numeric
from ..utilities import report

from collections import defaultdict

## Constants
I, Q, U, n_IQU = range(4)

class NA(object):
    """Nodal analysis base class.

    **Attributes**
       *linear*
          If True, only the linear branches will be evaluated
          If False, only the non-linear branches will be evaluated
    """
    ## Output vectors and jacobians
    x = None
    i = None
    q = None
    G = None
    C = None
    CY = None
    u = None
    
    ## Length of x-vector
    n = None 

    def __init__(self, cir, toolkit=numeric, refnode=circuit.gnd, linear=None):
        self.cir = cir
        self.refnode = refnode
        self.branches = []
        self.toolkit = toolkit
        self.linear = linear
        self.dirty = True
        self.circuitdirty = True

    def create_matrices(self):
        toolkit = self.toolkit

        n = self.n
        self.x          = toolkit.zeros(n)
        self.i          = toolkit.zeros(n)
        self.q          = toolkit.zeros(n)
        self.u          = toolkit.zeros(n)
        self.outvectors = [self.i, self.q, self.u]

        self.G          = toolkit.zeros((n,n))
        self.Gconst     = toolkit.zeros((n,n))
        self.C          = toolkit.zeros((n,n))
        self.Cconst     = toolkit.zeros((n,n))
        self.jacobians  = [self.G, self.C]
        self.CY         = toolkit.zeros((n,n))

    def setup_GCconst(self):
        """Evaluate linear part"""
        self.Gconst     = self.toolkit.zeros((self.n,self.n))
        self.Cconst     = self.toolkit.zeros((self.n,self.n))

        ## Set up branch mapping for linear and non-source branches only
        ## and evaluate G
        branchfilter = lambda branch: branch.linear and 'u' not in branch.output
        self.setup_branch_mapping(branchfilter)
        self.dirty = False # unset dirty flag to avoid infinite recursion
        self.update(self.x, defaultepar)
        self.dirty = True

        ## Update Gconst with it
        self.Gconst = self.G.copy()
        self.Cconst = self.C.copy()

    def update(self, x, epar, static=True, dynamic=True, jacobians=True, noise=False):
        raise NotImplementedError()
    
    def extract_v(self, x, nodep, noden=None):
        """Extract voltage between nodep and noden from the given x-vector.

        If noden is not given the voltage is taken between nodep and refnode. 

        *x*
          Solution vector

        *nodep*
          Node object or node reference in text format of positive node

        *noden*
          Node object or node reference in text format of negative node
        """
        raise NotImplementedError()
    
class MNA(NA):
    """Modified nodal analysis"""
    def setup(self):
        if self.circuitdirty:
            self.setup_circuit()
            self.create_matrices()
            self.generate_eval_iqu_and_der()
            self.circuitdirty = False
 
        ## Set up Gconstant
        self.setup_GCconst()

        ## Only evaluate non-linear branches and sources
        ## as all linear non-source branches are already evaluated in
        ## Gconst and Cconst            
        branchfilter = lambda branch: not branch.linear or 'u' in branch.output
        self.setup_branch_mapping(branchfilter)

        ## Now everything is set up
        self.dirty = False

    def setup_circuit(self):
        """Set up circuit related attributes"""
        ## Find mapping between circuit node indices and x index
        if self.refnode is not None:
            self.nodes = [node for node in self.cir.nodes 
                          if node is not self.refnode and node is not None]
        else:
            self.nodes = self.cir.nodes

        ## Find all potential branches and mapping to the x-vector
        self.potentialbranch_x_map = {}
        self.potentialbranches = []
        n_nodes = len(self.nodes)
        branches = self.cir.xflatbranches(potential=True)
        for i, branch in enumerate(branches):
            self.potentialbranch_x_map[branch] = i + n_nodes
            self.potentialbranches.append(branch)

        self.n = n_nodes + len(self.potentialbranch_x_map)

    def describe_x_vector(self, ):
        """Describe what quantities the elements of the x-vector represents

        >>> from elements import *
        >>> import numpy as np
        >>> c = SubCircuit()
        >>> c['VS'] = VS(1, gnd, v=1e-3)
        >>> c['R'] = R(1, gnd, r=1e3)
        >>> na = MNA(c, refnode=gnd, toolkit=symbolic)
        >>> na.describe_x_vector()
        [V(1), I(VS)]

        """
        return [Quantity('V', node)   for node in self.nodes] + \
               [Quantity('I', branch) for branch in self.potentialbranches]

    def describe_i_vector(self):
        """Describe what quantities the elements of the i-vector represents"""

        return [Quantity('I', node)   for node in self.nodes] + \
               [Quantity('V', branch) for branch in self.potentialbranches]

    def get_node_index(self, node, instname=None):
        """Return index to node voltage in solution vector of node"""

        # A node that is None means global ground
        if node is None:
            return None

        node = circuit.makenode(node)

        if instname is not None:
            node = self.cir.get_node(circuit.instjoin(instname, node.name))

        if node in (self.refnode, None):
            return None
        elif node in self.nodes:
            return self.nodes.index(node)
        else:
            from IPython import embed
            embed()
            raise ValueError('Node %s is not in MNA node list (%s)'%
                             (str(node), str(self.nodes)))

    def set_out_vector_branch(self, outtype, instname, branch, value):
        """Add to i, q or u vector elements that corresponds to given branch
        
        This function is used by some analyses to set stimulus without changing the circuit

        if instance name is given and branch is set to None the first branch of that
        instance will be used.
        """
        if instname is None:
            inst = self.cir
        else:
            inst = self.cir[instname]
        

        if branch is None:
            branch = inst.branches[0]

        out = self.outvectors[outtype]

        if branch.potential:
            branchref = circuit.BranchRef(instname, inst, branch)
            index = self.potentialbranch_x_map[branchref]
            out[index] += value
        else:
            indexp, indexn = (self.get_node_index(branch.plus, instname),
                              self.get_node_index(branch.minus, instname)) 

            if indexp is not None:
                out[indexp] += value
            if indexn is not None:
                out[indexn] += value

    def extract_v(self, x, nodep, noden=None):
        """Extract voltage between nodep and noden from the given x-vector.

        If noden is not given the voltage is taken between nodep and refnode. 

        *x*
          Solution vector

        *nodep*
          Node object or node reference in text format of positive node

        *noden*
          Node object or node reference in text format of negative node

        >>> from elements import *
        >>> import numpy as np
        >>> c = SubCircuit()
        >>> n1, n2 = c.add_nodes('n1','n2')
        >>> c['R1'] = R(n1, n2, r=1e3)
        >>> c['R2'] = R(n2, gnd, r=1e3)
        >>> mna = MNA(cir)
        >>> mna.extract_v(np.array([1.0, 0.5, 0.0]), 'n1', 'n2')
        0.5
        >>> mna.extract_v(np.array([1.0, 0.5, 0.0]), c.nodes[0])
        1.0
        
        """
        v = []
        for node in nodep, noden:
            if node is None:
                node = self.refnode

            index = self.get_node_index(node)

            if index is None: ## When node == refnode
                v.append(0)
            else:
                v.append(x[index])

        return v[0] - v[1]

    def extract_i(self, x, branch_or_term, xdot = None):
        """Extract branch or terminal current from the given x-vector.

        *branch_or_term*
           Branch object or terminal name

        *xdot*
           dx/dt vector. this is needed if dx/dt is non-zero and there is no branch defined 
           at the terminal

        >>> from elements import *
        >>> import numpy as np
        >>> c = SubCircuit()
        >>> net1 = c.add_node('net1')
        >>> c['vs'] = VS(net1, gnd)
        >>> c.extract_i(np.array([1.0, -1e-3]), 'vs.minus')
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

            branch_sign = self.cir.get_terminal_branch(branch_or_term)
            
            if branch_sign is not None:
                branch, sign = branch_sign
            else:
                raise ValueError('Current not saved. Use circuit.save_current')

                terminal_node = self.cir.get_node(branch_or_term)
                
                terminal_node_index = self.get_node_index(terminal_node)

                if xdot is not None:
                    if linearized:
                        return dot(self.G[terminal_node_index], x) + \
                            dot(self.C[terminal_node_index], xdot) + \
                            self.u[terminal_node_index]

                    else:
                        return self.i[terminal_node_index] + \
                            dot(self.C[terminal_node_index], xdot) + \
                            self.u[terminal_node_index]
                else:
                    if linearized:
                        return dot(self.G[terminal_node_index], x) + \
                            self.u[terminal_node_index]

                    else:
                        return self.i[terminal_node_index] + \
                            self.u[terminal_node_index]

        else:
            branch = branch_or_term
            sign = 1

        branchindex = self.potentialbranch_x_map[branch] 

        return sign * x[branchindex]      

    def map_branch_input(self, branchref):
        if not branchref.potential:
            return (self.get_node_index(branchref.plus,  branchref.instname),
                    self.get_node_index(branchref.minus, branchref.instname))
        else:
            return self.potentialbranch_x_map[branchref], None

    def map_branch_output(self, branchref):
        if not branchref.potential:
            return (self.get_node_index(branchref.plus,  branchref.instname),
                    self.get_node_index(branchref.minus, branchref.instname))
        else:
            return self.potentialbranch_x_map[branchref], None

    def setup_GCconst(self):
        super(MNA, self).setup_GCconst()
        ## Add gyrators for all direct potential branches
        ## Indirect potential nodes only have a KCL that ensures that the flow
        ## of the branch goes out and in of the plus and minus nodes
        ##
        ## KCL:
        ## i(branch.plus)  = i(branch)
        ## i(branch.minus) = i(branch)
        ## KVL:
        ## v(branch.plus) - v(branch.minus) = 0
        ##
        for i, branch in enumerate(self.cir.xflatbranches()):
            if branch.potential:
                ## Find indices to plus and minus nodes and branch current
                i_x_plus  = self.get_node_index(branch.plus,  branch.instname)
                i_x_minus = self.get_node_index(branch.minus, branch.instname)
                i_x_branch = self.potentialbranch_x_map[branch]

                ## KCL
                if i_x_plus is not None:
                    self.Gconst[i_x_plus,  i_x_branch] += 1
                if i_x_minus is not None:
                    self.Gconst[i_x_minus, i_x_branch] -= 1

                if not branch.indirect:
                    ## KVL
                    ## -V(branch) + i(x) + dq/dt(i(x0) + u(t) = 0
                    if i_x_plus is not None:
                        self.Gconst[i_x_branch, i_x_plus]  -= 1
                    if i_x_minus is not None:
                        self.Gconst[i_x_branch, i_x_minus] += 1

    def setup_branch_mapping(self, branchfilter=None):
        ## Set up mapping between global x-vector and inputs of all input branches.
        ##
        ## Branch inputs are potentials for flow branches and
        ## flows for potential branches
        ## The mapping is calculated as:
        ## x_branch[ix_branch_p] <- x[ix_p]   
        ## x_branch[ix_branch_n] <- x_branch[ix_branch_n] - x[ix_n]
        ## where x_branch is the concatenated inputs of all instances
        ## 
        ## This allows for either direct mapping between the x-vector and branch input
        ## or as the difference between two x-values
        ##
        ## Output mapping is set up between the branch outputs
        ## branch_iqu = [inst0_ioutbranch0, inst0_ioutbranch1, ...
        ##               inst0_qoutbranch0, inst0_qoutbranch1, ...
        ##               inst0_uoutbranch0, inst0_uoutbranch1, ...
        ##               inst1_ioutbranch0, inst1_ioutbranch1, ...
        ##               inst1_qoutbranch0, inst1_qoutbranch1, ...
        ##               inst1_uoutbranch0, inst1_uoutbranch1, ...
        ##               ...
        ## and global i, q, u vectors
        ## as
        ## i[i_index_p] <- i[i_index_p] + branch_i[iqu_i_index_p]
        ## i[i_index_n] <- i[i_index_n] - branch_i[iqu_i_index_n]

        self.out_branch_index_p       = list_of_lists(n_IQU)
        self.out_branch_index_n       = list_of_lists(n_IQU)
        self.out_index_p              = list_of_lists(n_IQU)
        self.out_index_n              = list_of_lists(n_IQU)
        self.jacobian_index_p1        = list_of_lists(n_IQU)
        self.jacobian_index_n1        = list_of_lists(n_IQU)
        self.jacobian_index_p2        = list_of_lists(n_IQU)
        self.jacobian_index_n2        = list_of_lists(n_IQU)
        self.jacobian_branch_index_p1 = list_of_lists(n_IQU)
        self.jacobian_branch_index_n1 = list_of_lists(n_IQU)
        self.jacobian_branch_index_p2 = list_of_lists(n_IQU)
        self.jacobian_branch_index_n2 = list_of_lists(n_IQU)

        inputbranch_counter     = 0
        out_branch_counter      = 0
        jacobian_branch_counter = 0
        self.ix_branch_p        = []
        self.ix_branch_n        = []
        self.ix_p               = []
        self.ix_n               = []
        instances               = []
        branch_input_sizes      = []
        branch_output_sizes     = []
        jacobian_sizes          = []

        ## Iterate over all branches
        branches = self.cir.xflatbranches(branchfilter=branchfilter)
        for (instname, inst), branchrefs in groupby(branches, 
                                                    lambda b: (b.instname, b.inst)):
            ## Generate list from iterator
            branchrefs = list(branchrefs)

            ## Save instances
            instances.append(inst)

            ## Categorize branches in input and output branches
            inbranches =  [branch 
                           for branch in branchrefs if branch.input]
            outbranches = [[branch 
                            for branch in branchrefs if 'i' in branch.output],
                           [branch 
                            for branch in branchrefs if 'q' in branch.output],
                           [branch 
                            for branch in branchrefs if 'u' in branch.output]]

            ## Calculate sizes of instance slices of branch input and output vectors
            branch_output_size = sum(map(len, outbranches))
            branch_input_sizes.append(len(inbranches))
            branch_output_sizes.append(branch_output_size)
            jacobian_sizes.append(branch_output_size * len(inbranches))

            ## Calculate mapping for branch inputs
            for branch in inbranches:
                ixp, ixn = self.map_branch_input(branch)
                if ixp is not None:
                    self.ix_branch_p.append(inputbranch_counter)
                    self.ix_p.append(ixp)

                if ixn is not None:
                    self.ix_branch_n.append(inputbranch_counter)
                    self.ix_n.append(ixn)
                
                inputbranch_counter += 1
            
            ## Calculate mapping between branch outputs and jacobians
            ## output i, q, u vectors and jacobians G, C matrices
            for out_type in I, Q, U:
                for branch in outbranches[out_type]:
                    ioutp, ioutn = self.map_branch_output(branch)
                    ## Set up mapping index vectors between branches and i,q,u
                    if ioutp is not None:
                        self.out_branch_index_p[out_type].append(out_branch_counter)
                        self.out_index_p[out_type].append(ioutp)
                    if ioutn is not None:
                        self.out_branch_index_n[out_type].append(out_branch_counter)
                        self.out_index_n[out_type].append(ioutn)

                    ## Set up mapping index vectors between branches and G, C
                    if out_type in (I,Q):
                        for inbranch in inbranches:
                            ixp, ixn = self.map_branch_input(inbranch)

                            if ioutp is not None and ixp is not None:
                                self.jacobian_index_p1[out_type].append((ioutp, ixp))
                                self.jacobian_branch_index_p1[out_type].append(jacobian_branch_counter)

                            if ioutp is not None and ixn is not None:
                                self.jacobian_index_n1[out_type].append((ioutp, ixn))
                                self.jacobian_branch_index_n1[out_type].append(jacobian_branch_counter)

                            if ioutn is not None and ixp is not None:
                                self.jacobian_index_n2[out_type].append((ioutn, ixp))
                                self.jacobian_branch_index_n2[out_type].append(jacobian_branch_counter)

                            if ioutn is not None and ixn is not None:
                                self.jacobian_index_p2[out_type].append((ioutn, ixn))
                                self.jacobian_branch_index_p2[out_type].append(jacobian_branch_counter)

                            jacobian_branch_counter += 1

                    out_branch_counter += 1

        ## Create branch input, output and jacobian vectors
        self.x_branch   = self.toolkit.zeros(inputbranch_counter)
        self.iqu        = self.toolkit.zeros(out_branch_counter)
        self.jacobian   = self.toolkit.zeros(jacobian_branch_counter)

        ## Slice it over instances
        x_branch_slices = slice_array(self.x_branch, branch_input_sizes)
        iqu_slices      = slice_array(self.iqu,      branch_output_sizes)
        jacobian_slices = slice_array(self.jacobian, jacobian_sizes)

        self.inst_and_slices = zip(instances, 
                                   x_branch_slices,
                                   iqu_slices, 
                                   jacobian_slices)

        ## Set up branch noise to CY matrix mapping
        self.setup_branch_noise_mapping()

    def setup_branch_noise_mapping(self, branchfilter=None):
        # Set up mapping between instance noise correlations and global CY matrix
        
        # Output mapping is set up between the branch noise outputs
        # branch_noise = [inst0_noise_corr_branch_0_0, inst0_noise_corr_branch_0_1, 
        #                 inst0_noise_corr_branch_1_0, inst0_noise_corr_branch_1_1, 
        #                 inst1_noise_corr_branch_0_0, inst1_noise_corr_branch_1_1, 
        #                 ...
        #                ]
        
        # and global CY matrix
        # as
        # CY[cy_indices_add] <- CY[cy_indices_add] + branch_noise[branch_noise_index_add]
        # CY[cy_indices_sub] <- CY[cy_indices_sub] - branch_noise[branch_noise_index_sub]

        cy_branch_i_p, cy_branch_i_n = None, None
        cy_branch_j_p, cy_branch_j_n = None, None
        def append_quad(cy_i_p, cy_i_n, cy_j_p, cy_j_n, branch_i):
            self.cy_indices_add.append((cy_i_p, cy_j_p))
            self.branch_noise_index_add.append(branch_i)

            if (cy_i_n is not None) and (cy_j_n is not None):
                self.cy_indices_add.append((cy_i_n, cy_j_n))
                self.branch_noise_index_add.append(branch_i)

            if cy_j_n is not None:
                self.cy_indices_sub.append((cy_i_p, cy_j_n))
                self.branch_noise_index_sub.append(branch_i)
            
            if cy_i_n is not None:
                self.cy_indices_sub.append((cy_i_n, cy_j_p))
                self.branch_noise_index_sub.append(branch_i)
        
        ## Init index vectors
        self.cy_indices_add = []
        self.cy_indices_sub = []
        self.branch_noise_index_add = []
        self.branch_noise_index_sub = []

        ## Iterate over noisy branches
        branches = self.cir.xflatbranches(branchfilter=lambda branch: branch.noisy)
        noise_sizes = [] # number of noise correlation parameters returned by inst.eval_noise()
        instances = []
        branch_index = 0
        for (instname, inst), branchrefs in groupby(branches, 
                                                    lambda b: (b.instname, b.inst)):
            start_branch_index = branch_index
            instances.append(inst)

            ## Generate list from iterator
            branchrefs = list(branchrefs)

            n_branches = len(branchrefs)
            i_noise = 0 # index in eval_noise method return value
            for i, branchref_i in enumerate(branchrefs):
                cy_branch_i_p, cy_branch_i_n = self.map_branch_output(branchref_i)
                
                if branchref_i.branch.noise_correlated:
                    for j, branchref_j in enumerate(branchrefs):
                        if branchref_j.branch.noise_correlated:
                            cy_branch_j_p, cy_branch_j_n = self.map_branch_output(branch_j)
                            append_quad(cy_branch_i_p, cy_branch_i_n,
                                        cy_branch_j_p, cy_branch_j_n, branch_index)
                            branch_index += 1
                else:
                    append_quad(cy_branch_i_p, cy_branch_i_n,
                                cy_branch_i_p, cy_branch_i_n, branch_index)
                    branch_index += 1

            noise_sizes.append(branch_index - start_branch_index)

        self.branch_noise = self.toolkit.zeros(branch_index)

        noise_slices = slice_array(self.branch_noise, noise_sizes)
        
        self.noise_inst_and_slices = zip(instances, noise_slices)


    def update(self, x, epar, static=True, dynamic=True, jacobians=True, noise=False):
        if self.dirty:
            self.setup()

        ## Init output vectors and matrices
        if x is not None and x is not self.x:
            self.x[:] = x
        self.i[:] = self.toolkit.dot(self.Gconst, self.x)
        self.q[:] = self.toolkit.dot(self.Cconst, self.x)
        self.u[:] = 0
        
        if jacobians:
            self.G = self.Gconst.copy()
            ## the jacobians vector needs to be updated beause G is a new object
            self.jacobians[I] = self.G
            self.C = self.Cconst.copy()
            ## the jacobians vector needs to be updated beause C is a new object
            self.jacobians[Q] = self.C

        ## Update input potentials of flow branches
        x_branch                    = self.x_branch
        x_branch[self.ix_branch_n]  =  0
        x_branch[self.ix_branch_p]  =  self.x[self.ix_p]
        x_branch[self.ix_branch_n]  -= self.x[self.ix_n]
        
        ## Evaluate circuit
        if dynamic:
            for inst_and_slices in self.inst_and_slices:
                inst, x_branch_slice, iqu_slice, der_slice = inst_and_slices

                iqu_slice[:], der_slice[:] = inst.eval_iqu_and_der(x_branch_slice, epar)
        else:
            for inst_and_slices in self.inst_and_slices:
                inst, x_branch_slice, iqu_slice, der_slice = inst_and_slices

                iqu_slice[:] = inst.eval_iqu(x_branch_slice, epar)

        ## Add outputs to i,q,u-vectors
        ## It would have been nice if this could have been done as:
        ## self.i[i_index_p] += iqu[iqu_i_index_p]
        ## self.i[i_index_n] -= iqu[iqu_i_index_n]
        ## Unfortunately i_index_p may have the same index twice which screws upp the addition
        branch_output = self.iqu

        for out_type in I,Q,U:
            inplace_add_from(self.outvectors[out_type],
                             branch_output, 
                             self.out_index_p[out_type], 
                             self.out_branch_index_p[out_type])
            inplace_sub_from(self.outvectors[out_type],
                             branch_output, 
                             self.out_index_n[out_type], 
                             self.out_branch_index_n[out_type])
            
        if dynamic:
            branch_jacobian = self.jacobian

            for out_type in I, Q:
                inplace_add_from(self.jacobians[out_type], branch_jacobian, 
                                 self.jacobian_index_p1[out_type], 
                                 self.jacobian_branch_index_p1[out_type])
                inplace_sub_from(self.jacobians[out_type], branch_jacobian,
                                 self.jacobian_index_n1[out_type], 
                                 self.jacobian_branch_index_n1[out_type])
                inplace_add_from(self.jacobians[out_type], branch_jacobian, 
                                 self.jacobian_index_p2[out_type], 
                                 self.jacobian_branch_index_p2[out_type])
                inplace_sub_from(self.jacobians[out_type], branch_jacobian,
                                 self.jacobian_index_n2[out_type], 
                                 self.jacobian_branch_index_n2[out_type])
            
        if noise:
            for inst, branch_noise_slice in self.noise_inst_and_slices:
                branch_noise_slice[:] = inst.eval_noise(epar)

            inplace_add_from(self.CY, self.branch_noise, 
                             self.cy_indices_add, self.branch_noise_index_add)
            inplace_sub_from(self.CY, self.branch_noise, 
                             self.cy_indices_sub, self.branch_noise_index_sub)

    def generate_eval_iqu_and_der(self):
        for instname, inst in self.cir.xflatinstances():
            if not hasattr(inst, 'eval_iqu_and_der_func'):
                self.toolkit.generate_eval_iqu_and_der(inst)

def list_of_lists(n):    
    return [[] for i in range(n)]

def inplace_add_from(dest, src, idest, isrc):
    """dest[idest] += src[isrc]"""
    for i_dst, i_src in zip(idest, isrc):
        dest[i_dst] += src[i_src]

def inplace_sub_from(dest, src, idest, isrc):
    """dest[idest] += src[isrc]"""
    for i_dst, i_src in zip(idest, isrc):
        dest[i_dst] -= src[i_src]

def slice_array(x, slicesizes):
    """Create a number of consecutive slices from array x with given sizes

    Example:

    array_slices(x, (3, 4, 2))

    would return

    x[0:4], x[4:8], x[8:10]

    """
    indices = np.cumsum([0] + list(slicesizes))
    return [x[indices[i]:indices[i+1]] for i in range(len(slicesizes))]

