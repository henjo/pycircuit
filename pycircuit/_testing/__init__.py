"""Test-suite infrastructure that ships with the package (2026-10-03): the
hdl backend state reader (`state`) and the leak detector (`leaks`), loaded by
the root conftest.  Kept inside the package so every test tree, a pytester
run and a subprocess can import it.  Nothing here is imported by the library
itself."""
