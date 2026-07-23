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
        
        # 1. Input resistance
        self['Rin'] = R(inp, inn, r=self.iparv.Rin)
        
        # 2. Gain and Dominant Pole Stage
        # The dominant pole frequency fd = GBW / Aol
        # fd = 1 / (2 * pi * R_pole * C_pole)
        # Let's set R_pole = 1.0 ohm for simplicity
        Aol = self.iparv.Aol
        GBW = self.iparv.GBW
        R_pole = 1.0
        
        import numpy as np
        fd = GBW / Aol
        C_pole = 1.0 / (2.0 * np.pi * R_pole * fd)
        
        n_mid = self.add_node('mid')
        
        # The VCCS drives the mid node. Gm = Aol / R_pole = Aol
        # Wait, if VCCS transconductance is Gm, Vmid = Gm * R_pole * (Vinp - Vinn)
        # Since R_pole = 1.0, Gm = Aol.
        Gm = Aol
        self['Gm1'] = BSource(inp, inn, n_mid, outn, i_func=lambda v: Gm * v)
        self['R_pole'] = R(n_mid, outn, r=R_pole)
        self['C_pole'] = C(n_mid, outn, c=C_pole)
        
        # 3. Output limiting and Output Resistance
        # We use a Voltage-Controlled Voltage Source (but implemented as a BSource with a small series resistor, or directly a BSource VCCS pushing into Rout)
        # Actually, V_mid is the unlimited voltage.
        # We need to limit V_mid to between Vss and Vdd.
        # We can use a BSource that evaluates tanh to limit the voltage.
        
        Vdd = self.iparv.Vdd
        Vss = self.iparv.Vss
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
            
        self['B_out'] = BSource(n_mid, outn, outp, outn, i_func=output_limit)
        self['R_out'] = R(outp, outn, r=Rout)
