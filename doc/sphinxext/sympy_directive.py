"""
A special directive for rendering results from sympy expressions using the 
math directive from either the pngmath or jsmath extension

Usage
-----

Can be used like this::

    .. sympy::

       x,y = sympify('xy')
       1/sqrt(1+x)

    Keep values of x,y using the persistent option, and use the docstring
    syntax

    .. sympy::
       :persistent:

       >>> 1/sqrt(1+x)

The content is interpreted as doctest formatted if it has a line starting
with ``>>>``.

The ``sympy`` directive supports the options

    include-source : bool
        Whether to display the source code. Default can be changed in conf.py

    persistent : bool
        The python session starts with the namespace from previous sympy 
        section

Configuration options
---------------------

The plot directive has the following configuration options:

    sympy_pre_code
        Code that should be executed before each sympy section

    sympy_include_source
        Default value for the include-source option

TODO
----


"""

import re
import sympy
from docutils.parsers.rst import directives

def setup(app):
    sympy_directive_options = {'include-source': _option_boolean,
                               'persistent': directives.flag}

    app.add_directive('sympy', sympy_directive, True, (0, 1, False),
                      **sympy_directive_options)

    app.add_config_value('sympy_pre_code', default_pre_code, True)
    app.add_config_value('sympy_include_source', True, True)

    app.connect('source-read', purge_sympy_namespace)

#------------------------------------------------------------------------------
# sympy:: directive registration etc.
#------------------------------------------------------------------------------

saved_namespace = None
def sympy_directive(name, arguments, options, content, lineno,
                   content_offset, block_text, state, state_machine):
    global saved_namespace

    document = state_machine.document
    env = document.settings.env
    config = env.config
    
    options.setdefault('include-source', config.sympy_include_source)
    if options['include-source'] is None:
        options['include-source'] = config.plot_include_source

    if 'persistent' in options and saved_namespace is not None:
        ns = saved_namespace
    else:
        ns = {}
        ## Evaluate pre-code
        exec config.sympy_pre_code in ns
 
    rst = ""
    codeblock = ''
    for line in content:
        ## Add line to source code block
        codeblock += line + '\n'

        if only_whitespace(line):
            continue

        ## Evaluate statement
        result = eval_line(unescape_doctest(line), ns)

        if result is not None:
            if options['include-source']:
                rst += rst_codeblock(codeblock)
                codeblock = ''

            ## Add result as math
            if is_sympy_object(result):
                latex_expr = sympy.latex(result, mode='plain')

                rst += '.. math::\n\n' + indent(latex_expr) + '\n'
            else:
                rst += str(result) + '\n'

    ## Flush remaining code block
    if options['include-source'] and codeblock:
        rst += rst_codeblock(codeblock)

    ## Save name space
    saved_namespace = ns        
        
    lines = rst.split("\n")
    if len(lines):
        state_machine.insert_input(
            lines, state_machine.input_lines.source(0))

    return []

def purge_sympy_namespace(app, docname, source):
    saved_namespace = None

def is_sympy_object(o):
    return isinstance(o, (sympy.Basic, sympy.Matrix))

def contains_doctest(text):
    r = re.compile(r'^\s*>>>', re.M)
    m = r.match(text)
    return bool(m)

def only_whitespace(text):
    r = re.compile(r'^\s*$')
    m = r.match(text)
    return bool(m)

def unescape_doctest(text):
    """
    Extract code from a piece of text, which contains either Python code
    or doctests.

    """
    if not contains_doctest(text):
        return text

    code = ""
    for line in text.split("\n"):
        m = re.match(r'^\s*(>>>|...) (.*)$', line)
        if m:
            code += m.group(2) + "\n"
        elif line.strip():
            code += "# " + line.strip() + "\n"
        else:
            code += "\n"
    return code

def indent(s, n=4, notfirstline = False):
    """Indent string

    >>> indent("apa\\nrapa\\nbapa", 4)
    '    apa\\n    rapa\\n    bapa'
    >>> indent("apa\\nrapa\\nbapa", 4, notfirstline=True)
    'apa\\n    rapa\\n    bapa'

    """
    if notfirstline:
        return ('\n' + n*' ').join(s.split('\n'))
    else:
        return '\n'.join([n*' ' + line for line in s.split('\n')])

def rst_codeblock(code):
    return ".. code-block:: python\n\n" + indent(code) + '\n'    

def eval_line(line, ns = {}):
    ## HACK to avoid problems with compile(.., ..., 'single') which
    ## when evaluated pollutes the global variable '_' which later
    ## interferes with Sphinx code
    
    try:
        code = compile(line, '<string>', 'eval')
    except SyntaxError:
        try:
            code = compile(line, '<string>', 'single')
        except SyntaxError:
            code = compile(line, '<string>', 'exec')
    try:
        return eval(code, ns)
    except Exception, e:
        raise Exception('%s raised when evaluating the line: %s'%
                        (repr(e), line))


def _option_boolean(arg):
    if not arg or not arg.strip():
        return None
    elif arg.strip().lower() in ('no', '0', 'false'):
        return False
    elif arg.strip().lower() in ('yes', '1', 'true'):
        return True
    else:
        raise ValueError('"%s" unknown boolean' % arg)

default_pre_code = """from __future__ import division
from sympy import *
x, y, z = symbols('x,y,z')
k, m, n = symbols('k,m,n', integer=True)
f, g, h = map(Function, 'fgh')
"""
