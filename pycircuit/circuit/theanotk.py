import theano
import itertools

from .numeric import *

## Import functions from theano
from theano.tensor import cos, sin, tan, cosh, sinh, tanh, log, exp, ceil, \
    floor, sqrt, real, imag

def generate_eval_iqu_and_der(cir):
    """Generate update_qiu_and_der method and attach to circuit instance"""

    n_inbranches    = len(cir.inputbranches)

    n_iqoutbranches = len([None for branch in cir.branches 
                           if 'i' in branch.output or
                              'q' in branch.output])

    def eval_iqu_and_der_func(cir, x, epar):
        if not hasattr(cir, '_theano_eval_iqu_and_der'):
            time = theano.tensor.dscalar('time')
            xsym = theano.tensor.vector('x').reshape((n_inbranches,))

            ## Copy epar and set time as a theano tensor
            myepar = epar.copy()
            myepar.t = time
            
            iqu = cir.eval_iqu(xsym, myepar)
            
            iq = iqu[:n_iqoutbranches]
            
            jacobian = [theano.tensor.grad(output, xsym) for output in iq]
            CG = list(itertools.chain(*jacobian))
            
            y = theano.tensor.concatenate([iqu, CG])
            fg = theano.function([xsym, time], y, on_unused_input='ignore')

            ## If CG is empty, theano.function only returns the function
            ## not the derivative
            if len(CG) == 0:
                cir._theano_eval_iqu_and_der = lambda x, t: (fg(x, t), ())
            else:
                cir._theano_eval_iqu_and_der = fg
            
        return cir._theano_eval_iqu_and_der(x, epar.t)
    
    cir.eval_iqu_and_der_func = eval_iqu_and_der_func
