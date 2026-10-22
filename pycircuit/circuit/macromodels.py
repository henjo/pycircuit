# -*- coding: utf-8 -*-
from pycircuit.circuit import SubCircuit, Parameter
from pycircuit.circuit.elements import R, C, VS, BSource

class OpAmp(SubCircuit):
    """
    Practical Operational Amplifier Macro Model.
    
    Includes:
    - Finite DC open-loop gain (Aol)
    - Input resistance (Rin)
    - Output resistance (Rout)
    - Dominant pole (GBW - Gain Bandwidth Product)
    - Output voltage limiting (Vdd / Vss) via non-linear tanh function.
    """
    terminals = ('inp', 'inn', 'outp', 'outn')
    instparams = [
        Parameter('Aol', 'DC open-loop gain', default=1e5),
        Parameter('GBW', 'Gain-Bandwidth product', default=1e6, unit='Hz'),
        Parameter('Rin', 'Input resistance', default=1e6, unit='Ohm'),
        Parameter('Rout', 'Output resistance', default=50.0, unit='Ohm'),
        Parameter('Vdd', 'Positive supply rail limit', default=15.0, unit='V'),
        Parameter('Vss', 'Negative supply rail limit', default=-15.0, unit='V'),
    ]

    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        
        inp = self.nodenames['inp']
        inn = self.nodenames['inn']
        outp = self.nodenames['outp']
        outn = self.nodenames['outn']
        
        ## ⚠ The children take EXPRESSIONS of this element's parameters, and
        ## the two sources' functions are rebuilt by `update` when theirs
        ## change, so a write after construction reaches them.  Until
        ## 2026-10-22 each got the value at construction, and a later write
        ## of `Aol`, `GBW`, `Rin`, `Rout`, `Vdd` or `Vss` changed nothing.

        # 1. Input resistance
        self['Rin'] = R(inp, inn, r='Rin')
        
        # 2. Gain and Dominant Pole Stage
        # The dominant pole frequency fd = GBW / Aol
        # fd = 1 / (2 * pi * R_pole * C_pole)
        # Let's set R_pole = 1.0 ohm for simplicity
        R_pole = 1.0
        
        n_mid = self.add_node('mid')
        
        # The VCCS drives the mid node. Gm = Aol / R_pole = Aol
        # Wait, if VCCS transconductance is Gm, Vmid = Gm * R_pole * (Vinp - Vinn)
        # Since R_pole = 1.0, Gm = Aol.
        gain, limit = self._functions()
        self['Gm1'] = BSource(inp, inn, n_mid, outn, i_func=gain)
        self['R_pole'] = R(n_mid, outn, r=R_pole)
        ## (C_pole = 1 / (2 pi R_pole fd), fd = GBW / Aol: the expression in
        ## the order the arithmetic always ran, numpy's pi spelled out)
        self['C_pole'] = C(n_mid, outn,
                           c='1.0 / (2.0 * 3.141592653589793 * 1.0 * (GBW / Aol))')
        
        # 3. Output limiting and Output Resistance
        # We use a Voltage-Controlled Voltage Source (but implemented as a BSource with a small series resistor, or directly a BSource VCCS pushing into Rout)
        # Actually, V_mid is the unlimited voltage.
        # We need to limit V_mid to between Vss and Vdd.
        # We can use a BSource that evaluates tanh to limit the voltage.
        
        self['B_out'] = BSource(n_mid, outn, outp, outn, i_func=limit)
        self['R_out'] = R(outp, outn, r='Rout')
        self._built_from = self._function_args()

    def _function_args(self):
        p = self.iparv
        return (p.Aol, p.Vdd, p.Vss, p.Rout)

    def _functions(self):
        """The gain stage's and the output limiter's functions, from the
        parameters as they stand."""
        Aol = self.iparv.Aol

        def gain(v):
            return Aol * v

        Vdd, Vss = self.iparv.Vdd, self.iparv.Vss
        Vcenter = (Vdd + Vss) / 2.0
        Vspan = (Vdd - Vss) / 2.0

        Rout = self.iparv.Rout

        # BSource is I = f(Vctrl). Here Vctrl is Vmid.
        def output_limit(v_mid):
            # Try to use JAX if v_mid is a JAX array, else fallback
            try:
                import jax.numpy as jnp
                if isinstance(v_mid, jnp.ndarray) or hasattr(v_mid, 'aval'):
                    tanh = jnp.tanh
                else:
                    import math
                    tanh = math.tanh
            except ImportError:
                import math
                tanh = math.tanh

            v_ideal = Vcenter + Vspan * tanh((v_mid - Vcenter) / Vspan)
            # The current needed to produce v_ideal across Rout is v_ideal / Rout
            # But wait! If we just supply current v_ideal/Rout, and there is a load on the output,
            # the voltage will drop!
            # The voltage source in series with Rout is a Thevenin equivalent.
            # Norton equivalent: current source of V_ideal / Rout, in parallel with Rout!
            return v_ideal / Rout

        return gain, output_limit

    def update(self, subject):
        ## (the sources' functions rebuilt when their parameters changed;
        ## `Circuit.__init__` calls this before the children exist)
        els = self.__dict__.get('elements', {})
        if 'B_out' in els and 'Gm1' in els:
            args = self._function_args()
            if any(a is not b and a != b for a, b in zip(args, self._built_from)):
                gain, limit = self._functions()
                els['Gm1'].ipar.i_func = gain
                els['B_out'].ipar.i_func = limit
                self._built_from = args
        super().update(subject)
