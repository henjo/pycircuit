# Plan for Extending PyCircuit Elements

To support a wider variety of real-world circuits and rival commercial SPICE simulators, `pycircuit`'s element library (`elements.py` and modular files like `mos.py` or `bjt.py`) needs to be expanded. 

Thanks to our new **JAX Auto-Diff** and **Vectorized Engine**, adding non-linear elements is now drastically simpler: developers only need to write the fundamental physics equations (`eval_i` and `eval_q`), and the engine automatically derives Jacobians and vector stamps them into sparse memory matrices!

Here is the structured roadmap for adding new elements.

## Phase 1: Essential Semiconductor Devices
These are the foundational blocks missing for analog design.

1. **Bipolar Junction Transistors (BJT)**
   - **Models**: Implement the standard **Ebers-Moll** (Level 1) model for both NPN and PNP.
   - **Implementation**: Define `eval_i` with forward and reverse active region equations (exponential terms). JAX handles the highly non-linear $g_m$ and $g_\pi$ derivatives seamlessly.
   - **Expansion**: A future upgrade can implement the Gummel-Poon model (including Early effect and high-level injection).

2. **Junction Field Effect Transistor (JFET)**
   - **Models**: Standard Shichman-Hodges JFET model.
   - **Implementation**: `eval_i` utilizing square-law equations for Triode and Saturation regions.

3. **Advanced Diodes**
   - **Zener Diode**: Extend the current `Diode` model to include avalanche breakdown voltages.
   - **Varactor Diode**: A voltage-dependent capacitor. Requires writing an `eval_q(x)` equation where $C(V)$ scales inversely with the depletion width.

## Phase 2: Advanced Sources
Time-domain simulations (Transient) heavily rely on complex input stimulus signals.

1. **Piecewise Linear Sources (PWL)**
   - **Feature**: `VPWL` and `IPWL`.
   - **Implementation**: Interpolation between a provided array of $(Time, Value)$ tuples. Will hook into the `next_event()` timing system so the Transient step-controller explicitly lands on the waveform corners.

2. **Exponential & Modulated Sources**
   - **Feature**: `VExp` (Exponential pulse), `VAM` (Amplitude Modulated), and `VSFFM` (Single Frequency FM).
   - **Implementation**: Pure math functions inherited from `TimeFunction` in `func.py`.

3. **Arbitrary Behavioral Sources (B-Sources)**
   - **Feature**: Non-linear dependent sources (e.g., $V = I^2 \cdot R$ or $V = \sin(V_{in})$).
   - **Implementation**: A general `NonLinearVCVS` using `eval_i(x)` where the user passes a lambda expression. JAX can natively compile user-provided lambda equations!

## Phase 3: Passives & Transmission Lines
For RF and high-speed digital designs.

1. **Mutual Inductance (Coupled Inductors)**
   - **Feature**: Magnetic coupling factor ($K$) between two existing `L` elements.
   - **Implementation**: Introduces off-diagonal elements into the $G$ matrix representing $M = K \sqrt{L_1 L_2}$.

2. **Transmission Lines (T-Lines)**
   - **Feature**: Ideal Lossless Transmission Line.
   - **Implementation**: Requires delay states. Usually implemented using the method of characteristics (Branin's method), effectively creating a time-delayed dependent source.

3. **Voltage/Current Controlled Switches**
   - **Feature**: Ideal switches that trigger based on threshold voltages with hysteresis.
   - **Implementation**: Smooth hyperbolic tangent (`tanh`) transitions to keep the matrix mathematically continuous for Newton-Raphson solvers.

## Phase 4: Macro Models
High-level building blocks to abstract away transistor-level complexity.

1. **Practical Operational Amplifier**
   - **Feature**: A leap beyond the ideal `Nullor`. Includes finite DC gain, dominant-pole bandwidth (GBW), input/output resistance, and voltage limiting.
   - **Implementation**: Composed structurally as a `SubCircuit` containing dependent sources, $R$, $C$, and `VCVS_limited`.

## Engineering Execution Strategy
Because we have `JAXToolkit`:
* We do **not** need to manually compute analytic derivatives for the BJTs, JFETs, or Varactors. 
* We simply translate the textbook current ($I$) and charge ($Q$) equations into `eval_i` and `eval_q` hooks.
* We must ensure parameters (e.g., thermal voltage $V_T$) reference `epar.T` so that temperature sweeps work dynamically.
