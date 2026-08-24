Example 11: characterising a foundry MOSFET
-------------------------------------------

A production compact model is not a formula you type in — it is a
**model card**, several hundred parameters a foundry extracted from
silicon, plus the geometry scaling that turns those into the values one
transistor of one size actually uses. This example reads IHP's PSP103
card for their open 130 nm process, builds the device, and measures it
the way a bench would: a transfer curve, an output curve, and the
small-signal operating point.

Nothing is fitted. Every number comes off the card.

.. literalinclude:: example11.py

Three things in it are worth pulling out.

**The temperature has to be said twice, and it has to agree.** The
scaling layer scales the card's parameters *to* a temperature; the
element evaluates *at* whatever ``defaultepar.T`` says. They are not
linked. Leaving the two 0.15 K apart is worth 0.04 % of the noise
density — small, invisible in the drain current, and it was the largest
remaining error in this model's comparison for a while.

**Two ways to drive the device.** ``transfer()`` puts it in a circuit
and runs a real DC solve, which is what you want when anything else is
connected to it. ``output()`` evaluates the element directly at a bias
vector, which is faster and exact when the device is the only thing in
the problem.

**Look internal nodes up by name.** The element owns more unknowns than
it has terminals — a gate node behind the gate resistance, and the
auxiliary node the induced gate noise network lives on. So

.. code-block:: python

    G[0, 1]      # 0.0, always

returns exactly zero: the drain current does not depend on the *external*
gate at all, because ``rg`` separates them. The transconductance is the
derivative with respect to the internal gate, and the robust way to find
it is by name:

.. code-block:: python

    n = {str(node): k for k, node in enumerate(e.nodes)}
    gm = G[n['d'], n['gi']]

Live results
````````````

Run when this page is built. The last row is PSP103's own operating
point, recorded from the vendor's compiled binary and committed with the
tests — so the example checks itself against the model it is imitating.

.. exec-rst::

    import os, sys, warnings
    warnings.simplefilter('ignore')
    import numpy as np

    ## The build's working directory is not fixed, so find the example
    ## by walking up from here rather than assuming one.
    d = os.getcwd()
    for _ in range(6):
        cand = os.path.join(d, 'doc', 'src', 'circuit', 'examples')
        for c in (cand, d):
            if os.path.isfile(os.path.join(c, 'example11.py')):
                sys.path.insert(0, c)
                d = None
                break
        if d is None:
            break
        d = os.path.dirname(d)

    try:
        import example11 as ex
        vg = np.linspace(0.4, 1.2, 5)
        idg = ex.transfer(vg)
        vd = np.linspace(0.2, 1.2, 5)
        idd = ex.output(vd)
        op = ex.operating_point()
        ok = True
    except Exception as exc:
        ok = False
        why = str(exc)

    if ok:
        print(".. list-table:: Transfer characteristic, "
              "W/L = 1 / 0.13 um, Vd = 1.2 V")
        print("   :header-rows: 1")
        print("")
        print("   * - Vg")
        print("     - Id")
        for a, b in zip(vg, idg):
            print("   * - %.2f V" % a)
            print("     - %.4e A" % b)
        print("")
        print(".. list-table:: Output characteristic, Vg = 1.2 V")
        print("   :header-rows: 1")
        print("")
        print("   * - Vd")
        print("     - Id")
        for a, b in zip(vd, idd):
            print("   * - %.2f V" % a)
            print("     - %.4e A" % b)
        print("")
        import json
        import pycircuit.circuit
        ref = json.load(open(os.path.join(
            os.path.dirname(pycircuit.circuit.__file__), 'tests', 'data',
            'psp103_ihp_iv.json')))['op']['short']
        pt = next(p for p in ref['points']
                  if (p['vg'], p['vd'], p['vb']) == (1.2, 1.2, 0.0))
        print(".. list-table:: Operating point at Vg = Vd = 1.2 V")
        print("   :header-rows: 1")
        print("")
        print("   * - ")
        print("     - Ids")
        print("     - gm")
        print("     - gds")
        print("     - gmb")
        print("   * - this model")
        print("     - %.5e" % op['ids'])
        print("     - %.5e" % op['gm'])
        print("     - %.5e" % op['gds'])
        print("     - %.5e" % op['gmb'])
        print("   * - PSP103")
        print("     - %.5e" % pt['ids'])
        print("     - %.5e" % pt['gm'])
        print("     - %.5e" % pt['gds'])
        print("     - %.5e" % pt['gmb'])
    else:
        print("*The IHP PDK is not installed here, so this example "
              "cannot be run at build time.* The recorded result is:")
        print("")
        print(".. code-block:: text")
        print("")
        print("    transfer  Vd=1.2: 1.3176e-06 2.4022e-05 1.0103e-04 "
              "2.3480e-04 4.0053e-04")
        print("    output    Vg=1.2: 2.1929e-04 3.2545e-04 3.6395e-04 "
              "3.8525e-04 4.0053e-04")
        print("    op        Ids=4.00526e-04  gm=8.68976e-04  "
              "gds=5.40648e-05  gmb=1.25324e-04")
        print("    psp103    Ids=4.00526e-04  gm=8.68977e-04  "
              "gds=5.40648e-05  gmb=1.25324e-04")

Where to go next
````````````````

:doc:`../compact_models` is the reference page for these devices — what
is in them, what is not, how the card is read, and how far the whole
thing has been measured against the vendor. :doc:`../hdl` is the
language they are written in.
