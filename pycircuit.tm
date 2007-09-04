<TeXmacs|1.0.6>

<style|article>

<\body>
  <section|The MNA method>

  The MNA method is the most common way of solving the voltages and currents
  of an electric circuit.

  For a circuit with N nodes and M independent voltage sources or inductors,
  the unknown quantities in the MNA method can be written as:

  <\eqnarray*>
    <tformat|<table|<row|<cell|x=<matrix|<tformat|<cwith|1|1|1|1|cell-valign|b>|<table|<row|<cell|v>>|<row|<cell|j>>>>>>|<cell|>|<cell|>>>>
  </eqnarray*>

  where <with|mode|math|v> is a vector of length N that contains the node
  voltages and <with|mode|math|j> is a vector of branch currents through
  independent voltage sources and inductors.

  <subsection|Solving>

  The behaviour of a non-linear circuit with energy storage elements such as
  inductors and capacitors can be found from the solution of a system of
  differential equations:

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|<frac|d|dt><with|prog-scripts|maxima|q(v(t))+i(v(t))+(t)=0><eq-number><label|diffeqnonlin>>>>>
  </eqnarray*>

  Where <with|mode|math|t> denotes time, v<with|mode|math|(t)> is a vector of
  node voltages and terminal currents. <with|mode|math|q(\<cdot\>)> is a
  function that maps <with|mode|math|v(t)> to sums of capacitive charge or
  inductive fluxes at a node, <with|mode|math|i(\<cdot\>)> is a function that
  maps the <with|mode|math|x(t)> vector to a sum of currents at a node and
  finally <with|mode|math|><with|mode|math|u(t)> is the input to the system.

  If the circuit is linear or linearized equation <reference|diffeqnonlin>
  can be rewritten in matrix form as:

  <\eqnarray*>
    <tformat|<table|<row|<cell|C*<frac|d|dt>*x(t)+G*x(t)+u(t)=0>|<cell|>|<cell|>>>>
  </eqnarray*>

  where the matrix elements are given by:

  <\eqnarray*>
    <tformat|<table|<row|<cell|C<rsub|j,k>=<frac|\<partial\>q<rsub|j>|\<partial\>x<rsub|k>>>|<cell|G<rsub|j,k>=<frac|\<partial\>i<rsub|j>|\<partial\>x<rsub|k>>>|<cell|>>>>
  </eqnarray*>

  The <with|mode|math|C> and <with|mode|math|G> matrices can be written as a
  matrix of 4 submatrices as:

  <\eqnarray*>
    <tformat|<table|<row|<cell|C=<matrix|<tformat|<cwith|1|1|1|1|cell-halign|r>|<cwith|1|1|1|1|cell-valign|b>|<table|<row|<cell|C<rsub|G>>|<cell|C<rsub|B>>>|<row|<cell|C<rsub|C>>|<cell|C<rsub|D>>>>>>>|<cell|G=<matrix|<tformat|<cwith|1|1|1|1|cell-halign|r>|<cwith|1|1|1|1|cell-valign|b>|<table|<row|<cell|G<rsub|G>>|<cell|G<rsub|B>>>|<row|<cell|G<rsub|C>>|<cell|G<rsub|D>>>>>>>|<cell|>>>>
  </eqnarray*>

  where:

  <\itemize>
    <item><with|mode|math|G<rsub|G>\<in\>\<cal-R\><rsup|N*\<times\>N>> is
    determined by the (trans)conductances of the circuit

    <item><with|mode|math|G<rsub|B>\<in\>\<cal-R\><rsup|N*\<times\>N>> is
    determined by the voltage sources

    <item><with|mode|math|G<rsub|C>\<in\>\<cal-R\><rsup|N*\<times\>N>> is
    determined by the voltage sources

    <item><with|mode|math|G<rsub|D>\<in\>\<cal-R\><rsup|N*\<times\>N>> can be
    non-zero when the circuit contains dependent sources
  </itemize>

  \;

  \;
</body>

<\initial>
  <\collection>
    <associate|preamble|false>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-2|<tuple|1.1|?>>
    <associate|auto-3|<tuple|1.2|?>>
    <associate|diffeqnonlin|<tuple|1|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>The
      MNA method> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1.5fn>|1.1<space|2spc>Solving
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>
    </associate>
  </collection>
</auxiliary>