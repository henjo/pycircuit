# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""``InternalResultDict.__len__`` read ``self.results``, but ``__init__``
sets ``self.items`` -- every other method on the class already uses
``self.items`` consistently, so ``len()`` was the one that never worked.
Found while fixing pycircuit/circuit/volterra.py (architecture.md P14),
which calls ``len()`` on the empty result its (deliberately unfinished)
``run()`` returns -- but the bug is in this widely-used class itself
(``feedback.py``, ``nportanalysis.py``, ``shooting.py``, ``analysis_ss.py``
all construct one), not specific to volterra.py.
"""
from pycircuit.post.internalresult import InternalResultDict


def test_len_matches_the_number_of_items():
    empty = InternalResultDict()
    assert len(empty) == 0

    populated = InternalResultDict({'a': 1, 'b': 2})
    assert len(populated) == 2

    populated['c'] = 3
    assert len(populated) == 3
