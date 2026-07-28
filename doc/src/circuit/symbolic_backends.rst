Choosing a symbolic backend
============================

:mod:`pycircuit.circuit` can solve the same circuit through any of six
symbolic toolkits, each answering a different question well. This page is
the map between them — the story that previously lived only in commit
messages and each backend's own doc page. It doesn't replace those pages;
it tells you which one to open.

Two things vary independently, and conflating them is the easiest way to
get lost here:

* **What a number is** — every toolkit on this page uses ``sympy``
  expressions. (The separate ``numeric``/``jax`` toolkits use floats/JAX
  tracers instead; they're not part of this comparison.)
* **How the answer is represented** — a fully expanded ``N(s)/D(s))``, a
  shared graph that is never expanded, or a chain of small assignments.
  This is the axis that actually distinguishes the six backends below.

The single fact that explains most of this page
------------------------------------------------

A **fully symbolic** circuit (every component a distinct symbol) has a
network determinant with exponentially many terms — no representation
yields a genuinely compact closed form, because the information content is
exponential. A circuit with **numeric components and only the frequency**
``s`` **symbolic** has few symbolic generators, and that's where the
fraction-free backends win big. Picking a backend is mostly about which of
these two situations you're in.

Comparison
----------

.. list-table::
   :header-rows: 1
   :widths: 16 20 34 20 10

   * - Toolkit
     - Represents the answer as
     - Best when
     - Status
     - Doc
   * - ``symbolic``
     - Expanded ``N(s)/D(s)`` (sympy ``LUsolve``)
     - Baseline; small circuits; anything not covered below
     - Stable, default
     - —
   * - ``symbolic_poly``
     - Expanded ``N(s)/D(s))``, fraction-free (``DomainMatrix.solve_den``)
     - Numeric components, symbolic ``s`` — the common case
     - Stable, recommended default for that case
     - :doc:`symbolic_poly`
   * - ``ginac_toolkit``
     - Expanded ``N(s)/D(s))``, compiled (GiNaC)
     - Same sweet spot as ``symbolic_poly``, past the size where sympy
       gets slow, *and* the system can be frequency-scaled to O(1)
       coefficients
     - Experimental, conditional win — needs the compiled extension
     - :doc:`ginac_native`
   * - ``symengine_toolkit``
     - Expanded ``N(s)/D(s))``, compiled (symengine)
     - Nothing today — kept as a plug-in point
     - Experimental, **not a speed win** (no fraction-free cancellation)
     - —
   * - ``ddd_toolkit``
     - Shared graph, never expanded
     - Fully symbolic circuits with a lot of internal repetition
       (hierarchy, many nullors) — where "expanded" isn't an option at all
     - Stable, has its own analysis family (two-port, loop gain, noise,
       sensitivity)
     - :doc:`ddd`
   * - ``soe`` (``solve_soe``)
     - Ordered chain of small assignments, never inlined
     - Fully symbolic circuits without DDD's kind of repetition (e.g. a
       plain resistor ladder) — compact through *sharing*, not cancellation
     - Research prototype — **not reachable through any toolkit or
       analysis**; call ``solve_soe`` directly
     - :doc:`soe_symbolic`

Reading the "best when" column as a decision tree: if the circuit is fully
symbolic, ask whether it's *structurally repetitive* (hierarchy, many
identical stages, nullor-heavy) — if so, DDD; if it's more of a flat ladder
or one-off topology, SoE. If most of the circuit is numeric and only ``s``
is symbolic, use ``symbolic_poly`` (or GiNaC once past the point where that
gets slow). ``symbolic`` is always a safe fallback and the right choice
when a circuit doesn't fit either sweet spot. ``symengine`` isn't a live
option today.

Why ``soe`` is different from the other five
----------------------------------------------

Every other row in the table above is a :class:`~pycircuit.circuit.toolkit.Toolkit`
you can pass as ``AC(cir, toolkit=...)`` — DDD in particular has a full
analysis family built on it (``DDDTwoPort``, ``DDDLoopGain``,
``DDDTransimpedance``). ``soe`` is a plain module-level function,
:func:`pycircuit.circuit.soe.solve_soe`, called directly on an assembled
``(A, b)`` system — it has no toolkit wrapper and no ``ac_solution`` path
into ``AC``/``DC``/``Noise``.

That's not an oversight so much as a real architectural mismatch: the
:class:`~pycircuit.circuit.acsolution.ACSolution` interface every other
representation implements assumes a *shared denominator* — ``numerators()``
over one ``denominator()``, ``poles()`` as its roots — because that's true
of an expanded ``N(s)/D(s)`` **and** of a DDD graph (a diagram is still
mathematically ``N(s)/D(s))``, just never written out). SoE's whole point is
that there often *isn't* one shared denominator worth naming — inlining its
assignments back into that shape is exactly the exponential blow-up it
exists to avoid (see :doc:`soe_symbolic`). Giving it a real ``AC``
integration would mean either extending ``ACSolution`` to a representation
that doesn't fit its current contract, or building SoE a separate analysis
path the way DDD has its own family — a design question, not something this
page resolves. Until that's decided, treat ``soe`` as what its own page
already says: a research prototype, useful directly, not part of the
toolkit hierarchy.
