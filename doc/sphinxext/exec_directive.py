"""A directive that executes Python at build time and inserts its printed reST.

Usage
-----

::

    .. exec-rst::

       print("Some **generated** reStructuredText")
       print("")
       print(".. list-table:: built live")
       ...

The block's ``stdout`` is captured and re-parsed as reStructuredText, so the
content (tables, numbers, figures) is regenerated on every documentation build
rather than being pasted in and going stale.

Robustness mirrors the ``sympy::`` directive: the block runs under a timeout and
any failure falls back to rendering the source (logged at info level), so a
single slow or broken block cannot break or hang the whole build.

Configuration
-------------

    exec_rst_execute : bool
        Run blocks (default ``True``).  Set ``False`` to render source only.
    exec_rst_timeout : int
        Per-block wall-clock limit in seconds (default 60).
"""

import contextlib
import io
import signal

from docutils.parsers.rst import Directive
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


def _indent(text, n=4):
    return '\n'.join((n * ' ' + line if line else line)
                     for line in text.split('\n'))


def setup(app):
    app.add_directive('exec-rst', ExecRstDirective)
    app.add_config_value('exec_rst_execute', True, True)
    app.add_config_value('exec_rst_timeout', 60, True)
    return {'parallel_read_safe': False, 'parallel_write_safe': True}


class ExecRstDirective(Directive):
    has_content = True

    def run(self):
        config = self.state.document.settings.env.config
        code = '\n'.join(self.content)

        if getattr(config, 'exec_rst_execute', True):
            try:
                with _time_limit(getattr(config, 'exec_rst_timeout', 60)):
                    buf = io.StringIO()
                    with contextlib.redirect_stdout(buf):
                        exec(code, {'__name__': '__exec_rst__'})
                    rst = buf.getvalue()
            except Exception as exc:
                ## Fall back to showing the source rather than break/hang the
                ## build; expected, designed behaviour -> info level, no warning.
                logger.info('exec-rst block not rendered live (%s); showing source', exc)
                rst = ".. code-block:: python\n\n" + _indent(code) + "\n"
        else:
            rst = ".. code-block:: python\n\n" + _indent(code) + "\n"

        lines = rst.split("\n")
        self.state_machine.insert_input(
            lines, self.state_machine.input_lines.source(0))
        return []
