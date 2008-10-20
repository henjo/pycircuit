import re
from sympy import Symbol, simplify, symbols, series

def approx(expr, patterns, n=2):
    """Approximate an expression using taylor series expansion

    The expression is approximated by doing a substitution of its variables as:

    .. math::

    v_k' = v_k \dot t^k

    where the magnitude of the `v_0, v_1, ..., v_k` variables are assumed to be in falling
    order. After the subsitution the taylor series around `t=0` of the expression is calculated and
    truncated after the n'th term. Finally t is substituted with 1.

    Arguments
    ---------

    expr -- Sympy expression
    patterns -- List of regular expression patterns that match the name of the variables
        in the order ([pattern_v0, pattern_v1, ...])

    Examples
    --------

    >>> a,b = symbols('ab')
    >>> approx(1+a+(a+b)*(2*a+3*b), ['a','b'], 3)
    3*a+3*b**2
    
    """
    t = Symbol('t')
    
    variables = expr.atoms(Symbol)

    sublist = []
    for i, pattern in enumerate(patterns):
        sublist.extend([(var.name, var * t**(i+1))
                        for var in variables if re.match(pattern, var.name)])

    return expr.subs(dict(sublist)).series(t, point=0, n=n, with_order=False)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
