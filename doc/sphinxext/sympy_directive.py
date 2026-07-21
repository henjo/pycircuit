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
import signal
import contextlib
import sympy
from docutils.parsers.rst import directives, Directive
from sphinx.util import logging

logger = logging.getLogger(__name__)


@contextlib.contextmanager
def _time_limit(seconds):
    """Abort the wrapped block with TimeoutError after `seconds` (POSIX only)."""
    def _handler(signum, frame):
        raise TimeoutError('exceeded %d s' % seconds)

    if not hasattr(signal, 'SIGALRM'):
        yield                       # no alarm on this platform; run unbounded
        return
    previous = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)


def _option_boolean(arg):
    if not arg or not arg.strip():
        return None
    elif arg.strip().lower() in ('no', '0', 'false'):
        return False
    elif arg.strip().lower() in ('yes', '1', 'true'):
        return True
    else:
        raise ValueError('"%s" unknown boolean' % arg)


def setup(app):
    app.add_directive('sympy', SympyDirective)

    app.add_config_value('sympy_pre_code', default_pre_code, True)
    app.add_config_value('sympy_include_source', True, True)
    ## Execute sympy:: blocks and render their results as math.  Execution is
    ## robust: a block that raises or exceeds sympy_timeout falls back to
    ## rendering its source (with a warning), so a single stale/slow legacy
    ## example cannot break or hang the whole build.  Set sympy_execute = False
    ## to render every block as source without running anything.
    app.add_config_value('sympy_execute', True, True)
    app.add_config_value('sympy_timeout', 10, True)

    app.connect('source-read', purge_sympy_namespace)

    return {'parallel_read_safe': False, 'parallel_write_safe': True}

#------------------------------------------------------------------------------
# sympy:: directive registration etc.
#------------------------------------------------------------------------------

saved_namespace = None

class SympyDirective(Directive):
    has_content = True
    required_arguments = 0
    optional_arguments = 1
    final_argument_whitespace = False
    option_spec = {'include-source': _option_boolean,
                   'persistent': directives.flag}

    def run(self):
        config = self.state.document.settings.env.config
        state_machine = self.state_machine

        if getattr(config, 'sympy_execute', True):
            try:
                with _time_limit(getattr(config, 'sympy_timeout', 20)):
                    rst = self._execute(config)
            except Exception as exc:
                ## A stale/slow legacy example: fall back to showing its source
                ## rather than break or hang the build.  This is expected,
                ## designed behaviour, so it is logged at info level (the source
                ## still renders) and does not count as a build warning.
                logger.info('sympy:: block not rendered live (%s); '
                            'showing source', exc)
                rst = rst_codeblock('\n'.join(self.content) + '\n')
        else:
            rst = rst_codeblock('\n'.join(self.content) + '\n')

        lines = rst.split("\n")
        if len(lines):
            state_machine.insert_input(
                lines, state_machine.input_lines.source(0))
        return []

    def _execute(self, config):
        """Run the block, returning rST with source blocks and result math.

        Raises on any evaluation error or timeout; run() catches that and falls
        back to source-only rendering for the whole block.
        """
        global saved_namespace

        options = self.options
        options.setdefault('include-source', config.sympy_include_source)
        if options['include-source'] is None:
            options['include-source'] = config.sympy_include_source

        if 'persistent' in options and saved_namespace is not None:
            ns = saved_namespace
        else:
            ns = {}
            exec(config.sympy_pre_code, ns)

        rst = ""
        codeblock = ''
        for line in self.content:
            codeblock += line + '\n'

            if only_whitespace(line):
                continue

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

        saved_namespace = ns
        return rst

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
    except Exception as e:
        raise Exception('%s raised when evaluating the line: %s'%
                        (repr(e), line))


default_pre_code = """from __future__ import division
from sympy import *
x, y, z = symbols('x,y,z')
k, m, n = symbols('k,m,n', integer=True)
f, g, h = map(Function, 'fgh')

# Compatibility shim: fraction_expand was removed from sympy but is used by the
# legacy example pages; it expands the numerator and denominator of a fraction.
def fraction_expand(expr, **hints):
    return expand(expr, frac=True, **hints)
"""
