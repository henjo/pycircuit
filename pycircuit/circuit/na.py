import numpy as np

import circuit
import numeric

from collections import defaultdict

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
        self.i  = toolkit.zeros(n)
        self.q  = toolkit.zeros(n)
        self.u  = toolkit.zeros(n)

        self.G  = toolkit.zeros((n,n))
        self.C  = toolkit.zeros((n,n))
        self.CY = toolkit.zeros((n,n))
        
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
        for i, branchinfo in enumerate(self.cir.xflatbranchmap()):
            instname, inst, branch, nodemap = branchinfo

            if branch.potential:
                key = (instname, branch)
                self.potentialbranch_x_map[key] = i + self.n

        self.n += len(self.potentialbranch_x_map)

        ## Create matrices
        self.create_matrices()
 
        self.setup_branch_mapping()

    def map_branch_input(self, instname, inst, branch, nodemap):
        if not branch.potential:
            return (self.node_index_to_x_index[nodemap[0]],
                    self.node_index_to_x_index[nodemap[1]])
        else:
            return self.potentialbranch_x_map[(instname, branch)], None

    map_branch_output = MNA.map_branch_input

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
        inputbranch_counter = 0
        self.ix_branch_p    = []
        self.ix_branch_n    = []
        self.ix_p           = []
        self.ix_n           = []
        instances           = []
        n_inputs_inst       = defaultdict(int)
        n_outputs_inst      = defaultdict(int)
        n_iq_outputs_inst   = defaultdict(int)
        for instname, inst, branch, nodemap in self.cir.xflatbranchmap():
            if branch.input:
                ixp, ixn = self.map_branch_input_to_x(instname, inst, branch, 
                                                      nodemap)
                if ixp is not None:
                    self.ix_branch_p.append(inputbranch_counter)
                    self.ix_p.append(ixp)

                if ixn is not None:
                    self.ix_branch_n.append(inputbranch_counter)
                    self.ix_n.append(ixn)
                
                n_inputs_inst[inst] += 1

                inputbranch_counter += 1

            if set(branch.output) & set('iqu'):
                ioutp, ioutn = self.map_branch_output(instname, inst, branch, 
                                                      nodemap)


                if 'i' in branch.output:
                    self.i_index_p.append(ioutp)
                
                n_outputs_inst[inst] += 1

                if 'u' not in branch.output:
                    n_iq_outputs_inst[inst] += 1

        ## List of instances
        instances           = n_inputs_inst.keys()

        ## Create x_branch vector
        self.x_branch       = np.zeros(inputbranch_counter)
        xbranch_slize_sizes = [n_inputs_inst[inst] for inst in instances]

        ## Create iqu and derivative vectors
        iqu_slize_sizes = [n_outputs_inst[inst] for inst in instances]
        der_slize_sizes = [n_iq_outputs_inst[inst] * n_inputs_inst[inst]
                           for inst in instances]
        self.iqu        = np.zeros(np.sum(iqu_slize_sizes))
        self.der        = np.zeros(np.sum(der_slize_sizes))

        ## Slice it over instances
        x_branch_slices = slice_array(self.x_branch, xbranch_slize_sizes)
        iqu_slices      = slice_array(self.iqu,      iqu_slize_sizes)
        der_slices      = slice_array(self.der,      der_slize_sizes)
        
        self.inst_and_slices = zip(instances, x_branch_slices, iqu_slices, der_slices)

        ## Create mapping from iqu branch output vectors to 
        ## circuit i, q and u global vectors
        

    def update(self, x, static=True, dynamic=True, jacobians=True, noise=False):
        ## Update input potentials of flow branches
        x_branch                    = self.x_branch
        x_branch[self.ix_branch_n]  =  0
        x_branch[self.ix_branch_p]  =  x[self.ix_p]
        x_branch[self.ix_branch_n]  -= x[self.ix_n]
        
        ## Evaluate circuit
        if dynamic:
            for inst_and_slices in self.inst_and_slices:
                inst, x_branch_slice, iqu_slice, der_slice = inst_and_slices

                iqu_slice[:], der_slice[:] = inst.eval_iqu_and_der(x_branch_slice)
        else:
            for inst, slices in self.inst_and_slices:
                inst, x_branch_slice, iqu_slice, der_slice = inst_and_slices

                iqu_slice[:] = inst.eval_iqu(x_branch_slice)

        ## Add outputs to i,q,u-vectors
        ## It would have been nice if this could have been done as:
        ## self.i[i_index_p] += iqu[iqu_i_index_p]
        ## self.i[i_index_n] -= iqu[iqu_i_index_n]
        ## Unfortunately i_index_p may have the same index twice which screws upp the addition
        iqu = self.iqu
        inplace_add_from(self.i, iqu, self.i_index_p, self.iqu_i_index_p)
        inplace_sub_from(self.i, iqu, self.i_index_n, self.iqu_i_index_n)
        inplace_add_from(self.q, iqu, self.q_index_p, self.iqu_q_index_p)
        inplace_sub_from(self.q, iqu, self.q_index_n, self.iqu_q_index_n)
        inplace_add_from(self.u, iqu, self.u_index_p, self.iqu_u_index_p)
        inplace_sub_from(self.u, iqu, self.u_index_n, self.iqu_u_index_n)

        if dynamic:
            inplace_add_from(self.G, der, G_index_p1, der_G_index)
            inplace_sub_from(self.G, der, G_index_n1, der_G_index)
            inplace_add_from(self.G, der, G_index_p2, der_G_index)
            inplace_sub_from(self.G, der, G_index_n2, der_G_index)

            inplace_add_from(self.C, der, C_index_p1, der_C_index)
            inplace_sub_from(self.C, der, C_index_n1, der_C_index)
            inplace_add_from(self.C, der, C_index_p2, der_C_index)
            inplace_sub_from(self.C, der, C_index_n2, der_C_index)


    def generate_eval_iqu_and_der(self):
        for inst in self.cir.xflatelements:
            if not hasattr(inst, 'eval_iqu_and_der_func'):
                self.toolkit.generate_eval_iqu_and_der(inst)

    
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

