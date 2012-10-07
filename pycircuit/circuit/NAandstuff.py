## Describe content of file
## Explain update of branches in place etc.
## Use docstring

## TODO:

## Separate linear and non-linear

## Sparse matrix interface (to NA)
#Several classes?
#Switching between different sparse matrices

## Topology analysis
#Done in circuit

## Noise
#Separate update process?
#Shaped noise

## Decorated functions
# Pass decorated function (function with some parameters preset) to
# instances/branches to let them update in correct place.
# Will probably give some overhead due to function calling.
# Will look cleaner and give less for-loops in xNA

class NA(object):
    '''Nodal analysis base class.
    '''
    G = None
    C = None
    i = None
    q = None
    CY = None
    u = None
    
    def __init__(self,cir):
        self.branches = []
        pass
    
    def update(x,t):
        pass
    

class MNA(NA):
    '''Modified Nodal Analysis class.

    Only for flat circuit for now !!
    '''
    
    ## all our branches (as from the circuit)

    def __init__(self):
        '''This will be a massive function
        '''
        # init G, C, i, q, CY and u here
        #self.i = 

    def update(self,x,t):
        '''Concept/example implementation of update !!
        Will look quite different when done. 
        '''
        in_i_branches = [[ntioinst0_branchR1],        # in branches in inst0
                         [inst1_branchR1], ....]  # in branches in inst1
        in_i_branch_indices = [ [[0,1]],          # plus and minus node indices for branches in inst0
                                [[1,0]] ]         # plus and minus node indices for branches in inst1

        i_outbranches = [[inst0_branchR1, inst0_branchR2],        # out branches in inst0
                       [inst1_branchR1, inst1_branchR2], ....]  # out branches in inst1
        i_outbranch_indices = [ [[0,1],[2,3]],          # plus and minus node indices for branches in inst0
                             [[1,0], [3,2]] ]         # plus and minus node indices for branches in inst1

        q_outbranches = [[inst0_branchR1, inst0_branchR2],        # out branches in inst0
                       [inst1_branchR1, inst1_branchR2], ....]  # out branches in inst1
        q_outbranch_indices = [ [[0,1],[2,3]],          # plus and minus node indices for branches in inst0
                             [[1,0], [3,2]] ]         # plus and minus node indices for branches in inst1
        self.i[:] = 0
        for inst, inst_inbranches, inst_inbranch_indices, inst_i_outbranches, \
                inst_q_outbranches, inst_i_outbranch_indices, inst_q_outbranch_indices \
                in zip(instances, inbranches, inbranch_indices, i_outbranches, q_outbranches, \
                           i_outbranch_indices, q_outbranch_indices):
            ## update all branches
            for inbranch, (index_p, index_n) in zip(inst_inbranches, inst_inbranch_indices):
                inbranch.v = x(index_p) - x(index_n)
            ## update_qiu on all instances
            inst.update_qiu(t)
            
            ## for all (out) i_branches in all instances 
            for outbranch, (out_index_p, out_index_n) in zip(inst_i_outbranches, inst_i_outbranch_indices):
                self.i[out_index_p] += outbranch.i
                self.i[out_index_n] -= outbranch.i
                self.u[out_index_p] += outbranch.u
                self.u[out_index_n] -= outbranch.u
                
                for g, (in_index_p, in_index_n) in zip(outbranch.G, inst_inbranch_indices):
                    self.G[out_index_p, in_index_p) += g;
                    self.G[out_index_p, in_index_n) -= g;
                    self.G[out_index_n, in_index_p) -= g;
                    self.G[out_index_n, in_index_n) += g;
                
            ## for all (out) q_branches in all instances 
            for outbranch, (index_p, index_n) in zip(inst_q_outbranches, inst_q_outbranch_indices):
                self.q[index_p] += outbranch.q
                self.q[index_n] -= outbranch.q
            
                for c, (in_index_p, in_index_n) in zip(outbranch.C, inst_inbranch_indices):
                    self.C[out_index_p, in_index_p) += c;
                    self.C[out_index_p, in_index_n) -= c;
                    self.C[out_index_n, in_index_p) -= c;
                    self.C[out_index_n, in_index_n) += c;
                    
            for something in someother:
                self.CY

class Circuit(object):
    '''New circuit class
    
    Copy doc from old class.
    '''

    branches = [] # class variable. will be initialised but never touched

    def __init__(self):
        pass

    def i(self):
        pass

    def u(self,t):
        pass

class Branch(object):
    '''Copy most from existing branch doc.
    '''
    ## Possibly add inormation regarding memory(q), linearity etc.
    
    i = None #flow
    q = None #time integral
    v = None #potential

    def __init__(self, plus, minus, potential=False): # default is 'flow' branch, not 'potential' branch
        self.potential = potential
        self.plus = plus
        self.minus = minus

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

    def __init__(self, plus, minus):
        self.plus = plus
        self.minus = minus        

class BranchV(Branch):
    '''Potential type branch

    v = dq/dt
    '''
    potential = True

    def __init__(self, plus, minus):
        self.plus = plus
        self.minus = minus        

class Pi(Circuit):
    '''Pi network circuit
    '''

    terminals = ('in', 'out', 'gnd')
    branches = [BranchI(1,'in'), BranchI(1,'gnd'), BranchI('out',1)]

    #iout_branches = [ i for i in branches if i.out == True]
    
    def __init__():
        pass

    def update_qiu(self,t):
        '''Update of q, i and u with the help of t
        '''
        
        #Branch aliases
        branch_r1, branch_c1, branch_r2 = self.branches        

        #Branch values
        branch_r1.i = branch_r1.v / 1e3 
        branch_r2.i = branch_r2.v / 2e3
        branch_c1.q = branch_c1.v * 1e-12

        #Branch derivatives
        #These could be lifted out to separate function for use only on nonlinear branches
        branch_r1.G(branch_r1) = 1 / 1e3
        branch_r2.G(branch_r2) = 1 / 2e3
        branch_c1.C(branch_c1) = 1e-12

