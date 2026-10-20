## THE ALLOCATOR'S MODE (2026-10-21): with env `PYCIRCUIT_MALLOC_PIN=1` glibc's
## malloc thresholds are pinned before anything allocates an array, so a
## process's speed does not depend on its allocation history (`_malloc`; the
## timing harnesses set it).
from . import _malloc  # noqa: F401 (imported for its effect)
from .circuit import *
from .elements import *
from .analysis import *
from .dcanalysis import *
from .symbolicdc import *
from .analysis_ss import *
from .nportanalysis import *
## STAGE 10.  `Transient` and `PSS` were the only analyses not reachable from the
## package root: `DC`, `DCSweep`, `AC` and `Noise` all arrive through the
## star-imports above, so `from pycircuit.circuit import Transient` -- the import
## every transient example and every doc page uses -- raised ImportError while its
## neighbours worked.
##
## Imported by name rather than with a star, because `transient.py` and
## `shooting.py` both do `from pycircuit.circuit.analysis import *` themselves,
## so a star here would re-export their transitive imports and make the package
## namespace depend on what those modules happen to pull in.
from .transient import Transient
from .shooting import PSS

from .toolkit import numeric, sparse_numeric, symbolic, symbolic_poly