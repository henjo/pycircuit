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

A block that fails falls back to rendering its source, so one slow or broken
block cannot break or hang the whole build.  **That fallback is loud**, and
deliberately so: it is logged as a *warning* and an admonition naming the
exception is rendered above the source in the built page.

The reason is worth stating, because the quiet version of this caused a real
problem.  When the fallback logged at info level, a block whose imports were
unresolvable rendered its own source and the build still reported "succeeded" --
so a page that was supposed to prove its numbers by regenerating them was
instead displaying the code that would have produced them, indefinitely and
invisibly.  A live-documentation directive that fails silently is worse than no
directive, because it looks like evidence.

Build with ``sphinx-build -W`` to make any such fallback fatal.  Whether the
docs are allowed to ship with a dead block is then a decision someone makes,
rather than one nobody sees.

Configuration
-------------

    exec_rst_execute : bool
        Run blocks (default ``True``).  Set ``False`` to render source only --
        the one case where falling back is not a failure, so it stays quiet.
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
                ## Fall back to the source so one bad block cannot break the
                ## build -- but say so, in the log *and* on the page.  A live
                ## block that quietly renders source still looks like generated
                ## output to a reader, which is the failure mode this avoids.
                logger.warning(
                    'exec-rst block did not run (%s: %s); rendering its source '
                    'instead -- the figures on this page are NOT live',
                    type(exc).__name__, exc,
                    location=(self.state.document['source'], self.lineno))
                rst = (".. warning::\n\n"
                       "   This block did not run at build time (``%s: %s``), so\n"
                       "   the source is shown instead of its output. Any numbers\n"
                       "   it was meant to generate are missing, not stale.\n\n"
                       % (type(exc).__name__, exc)
                       + ".. code-block:: python\n\n" + _indent(code) + "\n")
        else:
            ## Execution switched off wholesale -- a deliberate choice, not a
            ## failure, so this one stays quiet.
            rst = ".. code-block:: python\n\n" + _indent(code) + "\n"

        lines = rst.split("\n")
        self.state_machine.insert_input(
            lines, self.state_machine.input_lines.source(0))
        return []
