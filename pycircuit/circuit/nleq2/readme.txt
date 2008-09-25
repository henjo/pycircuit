This is an additional package that is distributed unchanged, with permission
from ZIB, under the terms of their own licence and warranty as described in nleq2.f
and *not* the PySCeS GPL licence.

Please take careful notice of the usage conditions, licence  and warranty as
described in nleq2_readme.txt. If you do not agree with these usage conditions
you must either disable the installation of nleq2 or obtain a licence from ZIB.

To disable the co-installation of the nleq2 non-linear solver in
the user configuration section:

-- setup.py --
 # set
 nleq2 = 0

------------------------------------------------------------
Note that a few adaptations may be necessary to utilize the code for
your computer at hand. Please examine the two subroutines D1MACH and
SECOND, which you can locate at the end of the file nleq2.f. Probably,
you need to adapt D1MACH to correspond to the arithmetic of your
computer.
------------------------------------------------------------
====================================================================
For IEEE big (Motorola RiSC) and little (Intel or AMD i386) endian
processors this is done automatically - Brett Olivier 20040421
====================================================================
