import numpy as np
from itertools import groupby

import circuit
import numeric

from collections import defaultdict

## Constants
I, Q, U, n_IQU = range(4)

class NA(object):
    '''Nodal analysis base class.
    '''
    ## Output vectors and jacobians
    i = None
    q = None
    G = None
    C = None
    CY = None
    u = None

    ## Length of x-vector
    n = None 
    
    def __init__(self, cir, toolkit=numeric, refnode=circuit.gnd):
        self.cir = cir
        self.refnode = refnode
        self.branches = []
        self.toolkit = toolkit

        if refnode is not None:
            self.refnode_index = self.cir.get_node_index(refnode)

        self.generate_eval_iqu_and_der()

    def create_matrices(self):
        toolkit = self.toolkit

        n = self.n
        self.i          = toolkit.zeros(n)
        self.q          = toolkit.zeros(n)
        self.u          = toolkit.zeros(n)
        self.outvectors = [self.i, self.q, self.u]

        self.G          = toolkit.zeros((n,n))
        self.Gconst     = toolkit.zeros((n,n))
        self.C          = toolkit.zeros((n,n))
        self.jacobians  = [self.G, self.C]
        self.CY         = toolkit.zeros((n,n))

    def update(x,t):
        raise NotImplementedError()
    
    ## Moved from circuit, Needs to be rewritten!!
    # def extract_v(self, x, nodep, noden=None, refnode=gnd, 
    #               refnode_removed=False):
    #     """Extract voltage between nodep and noden from the given x-vector.

    #     If noden is not given the voltage is taken between nodep and refnode. 
    #     x-vectors with the reference node removed can be handled by setting 
    #     the refnode_removed to True.

    #     *x*
    #       x-vector

    #     *nodep*
    #       Node object or node reference in text format of positive node

    #     *noden*
    #       Node object or node reference in text format of negative node

    #     *refnode*
    #       reference node

    #     *refnode_removed*
    #       If set the refernce node is expected to be removed from the x-vector

    #     >>> from elements import *
    #     >>> import numpy as np
    #     >>> c = SubCircuit()
    #     >>> n1, n2 = c.add_nodes('n1','n2')
    #     >>> c['R1'] = R(n1, n2, r=1e3)
    #     >>> c['R2'] = R(n2, gnd, r=1e3)
    #     >>> c.extract_v(np.array([1.0, 0.5, 0.0]), 'n1', 'n2')
    #     0.5
    #     >>> c.extract_v(np.array([1.0, 0.5, 0.0]), c.nodes[0])
    #     1.0
    #     >>> c.extract_v(np.array([1.0, 0.5]), c.nodes[0], refnode_removed = True)
    #     1.0
        
    #     """
    #     v = []
    #     for node in nodep, noden:
    #         if type(node) is types.StringType:
    #             node = self.get_node(node)
    #         elif node == None:
    #             node = refnode

    #         if refnode_removed:
    #             nodeindex = self.get_node_index(node, refnode)
    #         else:
    #             nodeindex = self.get_node_index(node, None)

    #         if nodeindex == None: ## When node == refnode
    #             v.append(0)
    #             continue
                    
    #         v.append(x[nodeindex])

    #     return v[0] - v[1]

        
    # ## Moved from NA, Needs to be rewritten!!
    # def extract_i(self, x, branch_or_term, xdot = None,
    #               refnode = gnd, refnode_removed = False,
    #               t = 0,
    #               linearized = False, xdcop = None):
    #     """Extract branch or terminal current from the given x-vector.

    #     *x* 
    #        x-vector

    #     *branch_or_term*
    #        Branch object or terminal name

    #     *xdot*
    #        dx/dt vector. this is needed if dx/dt is non-zero and there is no branch defined at the
    #        terminal

    #     *refnode*
    #        reference node

    #     *refnode_removed*
    #        If set the refernce node is expected to be removed from the x-vector

    #     *t*
    #        Time when the sources are to be evaluated
        
    #     *linearized*
    #        Set to True if the AC current is wanted

    #     *xdcop*
    #        *xcdop* is the DC operation point x-vector if linearized == True

    #     >>> from elements import *
    #     >>> import numpy as np
    #     >>> c = SubCircuit()
    #     >>> net1 = c.add_node('net1')
    #     >>> c['vs'] = VS(net1, gnd)
    #     >>> c.extract_i(np.array([1.0, 0, -1e-3]), 'vs.minus')
    #     0.001
    #     >>> c.extract_i(np.array([1.0, -1e-3]), 'vs.minus', refnode_removed = True)
    #     0.001
        
    #     """
    #     dot = self.toolkit.dot
        
    #     if type(branch_or_term) is types.StringType:
    #         ## Calculate current going in to the terminal as
    #         ## self.i(x)[terminal_node] + u(t) + dq(x)/dt. 
    #         ## This will work since i(x) returns
    #         ## the sum of all the currents going out from the
    #         ## terminal node that originates from devices within
    #         ## the circuit. According to Kirchoff's current law
    #         ## of the terminal node
    #         ## -I_external + sum(I_internal_k) = 0
    #         ## Where I_external represents the current coming from
    #         ## outside the circuit going *in* to the terminal node,
    #         ## I_internal_k represents one of the currents that flows
    #         ## from the terminal node to a device within the circuit.
    #         ## So now we can calculate I_external as 
    #         ## I_external = sum(I_internal_k) = 
    #         ## self.I(x)[terminal_node] + u(t) + dq(x)/dt =
    #         ## self.I(x)[terminal_node] + u(t) + sum(dq(x)/dx_k * dx_k/dt) =
    #         ## self.I(x)[terminal_node] + u(t) + C(x) * dx/dt

    #         branch_sign = self.get_terminal_branch(branch_or_term)

    #         if branch_sign != None:
    #             branch, sign = branch_sign
    #         else:
    #             terminal_node = self.nodenames[branch_or_term]
                
    #             terminal_node_index = self.get_node_index(terminal_node)

    #             if xdot != None:
    #                 if linearized:
    #                     return dot(self.G(xdcop)[terminal_node_index], x) + \
    #                         dot(self.C(xdcop)[terminal_node_index], xdot) + \
    #                         self.u(t, analysis = 'ac')[terminal_node_index]

    #                 else:
    #                     return self.i(x)[terminal_node_index] + \
    #                         dot(self.C(x)[terminal_node_index], xdot) + \
    #                         self.u(t)[terminal_node_index]
    #             else:
    #                 if linearized:
    #                     return dot(self.G(xdcop)[terminal_node_index], x) + \
    #                         self.u(t, analysis = 'ac')[terminal_node_index]

    #                 else:
    #                     return self.i(x)[terminal_node_index] + \
    #                         self.u(t)[terminal_node_index]

    #     else:
    #         branch = branch_or_term
    #         sign = 1

    #     branchindex = self.get_branch_index(branch)

    #     if refnode_removed:
    #         branchindex -= 1

    #     return sign * x[branchindex]      


class MNA(NA):
    """Modified nodal analysis"""

    def __init__(self, cir, toolkit=numeric, refnode=circuit.gnd):
        super(MNA, self).__init__(cir, refnode=refnode, toolkit=toolkit)

        ## Find mapping between circuit node indices and x index
        if self.refnode is not None:
            self.node_index_to_x_index = range(self.refnode_index) + [None] + \
                range(self.refnode_index, len(self.cir.nodes)-1)
            self.n  = len(self.cir.nodes) - 1
        else:
            self.node_index_to_x_index = range(self.cir.nodes)        
            self.n  = len(self.cir.nodes)

        ## Find all potential branches and mapping to the x-vector
        self.potentialbranch_x_map = {}
        for i, branchinfo in enumerate(self.cir.xflatbranchmap(potential=True)):
            instname, inst, branch, nodemap = branchinfo

            key = (instname, branch)
            self.potentialbranch_x_map[key] = i + self.n

        self.n += len(self.potentialbranch_x_map)

        ## Create matrices
        self.create_matrices()
 
        ## Set up Gconstant
        self.setup_Gconst()
        
        self.setup_branch_mapping()

    def map_branch_input(self, instname, inst, branch, nodemap):
        if not branch.potential:
            return (self.node_index_to_x_index[nodemap[0]],
                    self.node_index_to_x_index[nodemap[1]])
        else:
            return self.potentialbranch_x_map[(instname, branch)], None

    def map_branch_output(self, instname, inst, branch, nodemap):
        if not branch.potential:
            return (self.node_index_to_x_index[nodemap[0]],
                    self.node_index_to_x_index[nodemap[1]])
        else:
            return self.potentialbranch_x_map[(instname, branch)], None

    def setup_Gconst(self):
        ## Add gyrators for all potential branches
        ## KCL:
        ## i(branch.plus)  = i(branch)
        ## i(branch.minus) = i(branch)
        ## KVL:
        ## v(branch.plus) - v(branch.minus) = 0
        for i, branchinfo in enumerate(self.cir.xflatbranchmap()):
            instname, inst, branch, nodemap = branchinfo

            if branch.potential:
                ## Find indices to plus and minus nodes and branch current
                i_x_plus  = self.node_index_to_x_index[nodemap[0]] 
                i_x_minus = self.node_index_to_x_index[nodemap[1]] 
                i_x_branch = self.potentialbranch_x_map[(instname, branch)]

                ## KCL
                if i_x_plus is not None:
                    self.Gconst[i_x_plus,  i_x_branch] += 1
                if i_x_minus is not None:
                    self.Gconst[i_x_minus, i_x_branch] -= 1
                ## KVL
                if i_x_plus is not None:
                    self.Gconst[i_x_branch, i_x_plus]  += 1
                if i_x_minus is not None:
                    self.Gconst[i_x_branch, i_x_minus] -= 1
            

    def setup_branch_mapping(self):
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
        for (instname, inst), instbranchmaps in groupby(self.cir.xflatbranchmap(), 
                                                        lambda x: x[0:2]):
            ## Generate list from iterator
            instbranchmaps = list(instbranchmaps)

            ## Save instances
            instances.append(inst)

            ## Extract branches and mapping for each instance
            branches   = [branch          
                          for iname, inst, branch, mapping in instbranchmaps]
            branchmaps = dict([(branch, mapping)
                               for iname, inst, branch, mapping in instbranchmaps])

            ## Categorize branches in input and output branches
            inbranches =  [branch 
                           for branch in branches if branch.input]
            outbranches = [[branch 
                            for branch in branches if 'i' in branch.output],
                           [branch 
                            for branch in branches if 'q' in branch.output],
                           [branch 
                            for branch in branches if 'u' in branch.output]]

            ## Calculate sizes of instance slices of branch input and output vectors
            branch_output_size = sum(map(len, outbranches))
            branch_input_sizes.append(len(inbranches))
            branch_output_sizes.append(branch_output_size)
            jacobian_sizes.append(branch_output_size * len(inbranches))

            ## Calculate mapping for branch inputs
            for branch in inbranches:
                ixp, ixn = self.map_branch_input(instname, inst, branch, 
                                                 branchmaps[branch])
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
                    ioutp, ioutn = self.map_branch_output(instname, inst, 
                                                          branch, 
                                                          branchmaps[branch])
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
                            nodemap = branchmaps[inbranch]
                            ixp, ixn = self.map_branch_input(instname,inst,
                                                             inbranch,
                                                             nodemap)
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

    def update(self, x, epar, static=True, dynamic=True, jacobians=True, noise=False):
        ## Init output vectors and matrices
        self.i[:] = self.toolkit.dot(self.Gconst, x)
        self.q[:] = 0
        
        if jacobians:
            self.G[:,:] = self.Gconst[:,:]
            self.C[:,:] = 0

        ## Update input potentials of flow branches
        x_branch                    = self.x_branch
        x_branch[self.ix_branch_n]  =  0
        x_branch[self.ix_branch_p]  =  x[self.ix_p]
        x_branch[self.ix_branch_n]  -= x[self.ix_n]
        
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

    def generate_eval_iqu_and_der(self):
        for inst in self.cir.xflatelements:
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

