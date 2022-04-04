<h1 id="generalizedalpha">GeneralizedAlpha</h1>
<p><span class="math inline"><em>α</em><sub><em>F</em></sub></span> and
<span class="math inline"><em>α</em><sub><em>M</em></sub></span> are
defined differently than in the paper, we use <span
class="math inline"><em>α</em><sub><em>F</em></sub> = (1−<em>α</em><sub><em>f</em></sub>)</span>
and <span
class="math inline"><em>α</em><sub><em>M</em></sub> = (1−<em>γ</em><sub><em>m</em></sub>)</span>
where <span class="math inline"><em>α</em><sub><em>f</em></sub></span>
and <span class="math inline"><em>α</em><sub><em>m</em></sub></span> are
those used in the paper. Like Newmark and all the implicit schemes, the
unconditional stability of this method applies to linear problems. There
are no results showing stability of this method over the wide range of
nonlinear problems that potentially exist. Experience indicates that the
time step for implicit schemes in nonlinear situations can be much
greater than those for explicit schemes.</p>
<ul>
<li><span
class="math inline"><em>α</em><sub><em>M</em></sub> = 1.0, <em>α</em><sub><em>F</em></sub> = 1.0</span>
produces the Newmark Method.</li>
<li><span
class="math inline"><em>α</em><sub><em>M</em></sub> = 1.0</span>
corresponds to the HHT method.</li>
</ul>
<p>The method is second-order accurate provided <span
class="math inline">$\gamma = \tfrac{1}{2} + \alpha_M - \alpha_F$</span>
- The method is unconditionally stable provided <span
class="math inline">$\alpha_M &gt;= \alpha_F &gt;= \tfrac{1}{2},
\beta&gt;=\tfrac{1}{4} +\tfrac{1}{2}(\gamma_M - \gamma_F)$</span></p>
<ul>
<li><p><span class="math inline"><em>γ</em></span> and <span
class="math inline"><em>β</em></span> are optional. The default values
ensure the method is unconditionally stable, second order accurate and
high frequency dissipation is maximized. The defaults are:</p>
<p><span class="math display">$$\gamma = \tfrac{1}{2} + \gamma_M -
\gamma_F$$</span></p>
<p>and</p></li>
</ul>
<p><span class="math display">$$\beta = \tfrac{1}{4}(1 + \gamma_M -
\gamma_F)^2$$</span></p>
<h2 id="theory">Theory</h2>
<p>The Generalized <span class="math inline"><em>α</em></span> method
(sometimes called the <span class="math inline"><em>α</em></span>
method) is a one step implicit method for solving the transient problem
which attempts to increase the amount of numerical damping present
without degrading the order of accuracy. As with the HHT method, the
following Newmark approximations are used:</p>
<p><span
class="math display"><em>U</em><sub><em>t</em> + <em>Δ</em><em>t</em></sub> = <em>U</em><sub><em>t</em></sub> + <em>Δ</em><em>t</em><em>u̇</em><sub><em>t</em></sub> + [(0.5−<em>β</em>)<em>Δ</em><em>t</em><sup>2</sup>]<em>ü</em><sub><em>t</em></sub> + [<em>β</em><em>Δ</em><em>t</em><sup>2</sup>]<em>ü</em><sub><em>t</em> + <em>Δ</em><em>t</em></sub></span></p>
<p><span
class="math display"><em>u̇</em><sub><em>t</em> + <em>Δ</em><em>t</em></sub> = <em>u̇</em><sub><em>t</em></sub> + [(1−<em>γ</em>)<em>Δ</em><em>t</em>]<em>ü</em><sub><em>t</em></sub> + [<em>γ</em><em>Δ</em><em>t</em>]<em>ü</em><sub><em>t</em> + <em>Δ</em><em>t</em></sub></span></p>
<p>but the time-discrete momentum equation is modified:</p>
<p><span
class="math display"><em>R</em><sub><em>t</em> + <em>α</em><sub><em>M</em></sub><em>Δ</em><em>t</em></sub> = <em>F</em><sub><em>t</em> + <em>Δ</em><em>t</em></sub><sup>ext</sup> − <em>M</em><em>ü</em><sub><em>t</em> + <em>α</em><sub><em>M</em></sub><em>Δ</em><em>t</em></sub> − <em>C</em><em>u̇</em><sub><em>t</em> + <em>α</em><sub><em>F</em></sub><em>Δ</em><em>t</em></sub> − <em>F</em><sup>int</sup>(<em>U</em><sub><em>t</em> + <em>α</em><sub><em>F</em></sub><em>Δ</em><em>t</em></sub>)</span></p>
<p>where the displacements and velocities at the intermediate point are
given by:</p>
<p><span
class="math display"><strong>d</strong><sub><em>t</em> + <em>α</em><sub><em>F</em></sub><em>Δ</em><em>t</em></sub> = (1−<em>α</em><sub><em>F</em></sub>)<em>U</em><sub><em>t</em></sub> + <em>α</em><sub><em>F</em></sub><em>U</em><sub><em>t</em> + <em>Δ</em><em>t</em></sub></span></p>
<p><span
class="math display"><strong>v</strong><sub><em>t</em> + <em>α</em><sub><em>F</em></sub><em>Δ</em><em>t</em></sub> = (1−<em>α</em><sub><em>F</em></sub>)<em>u̇</em><sub><em>t</em></sub> + <em>α</em><sub><em>F</em></sub><em>u̇</em><sub><em>t</em> + <em>Δ</em><em>t</em></sub></span></p>
<p><span
class="math display"><strong>a</strong><sub><em>t</em> + <em>α</em><sub><em>M</em></sub><em>Δ</em><em>t</em></sub> = (1−<em>α</em><sub><em>M</em></sub>)<em>ü</em><sub><em>t</em></sub> + <em>α</em><sub><em>M</em></sub><em>ü</em><sub><em>t</em> + <em>Δ</em><em>t</em></sub></span></p>
<p>Following the methods outlined for Newmarks method, linearization of
the nonlinear momentum equation results in the following linear
equations:</p>
<p><span
class="math display"><em>K</em><sub><em>t</em> + <em>Δ</em><em>t</em></sub><sup>*<em>i</em></sup><em>d</em><em>U</em><sub><em>t</em> + <em>Δ</em><em>t</em></sub><sup><em>i</em> + 1</sup> = <em>R</em><sub><em>t</em> + <em>Δ</em><em>t</em></sub><sup><em>i</em></sup></span></p>
<p>where</p>
<p><span class="math display">$$
K_{t+\Delta t}^{*i} = \alpha_F K_t + \frac{\alpha_F \gamma}{\beta \Delta
t} C_t + \frac{\alpha_M}{\beta \Delta t^2} M
$$</span></p>
<p>and</p>
<p><span
class="math display"><em>R</em><sub><em>t</em> + <em>Δ</em><em>t</em></sub><sup><em>i</em></sup> = <em>F</em><sub><em>t</em> + <em>Δ</em><em>t</em></sub><sup>ext</sup> − <em>F</em>(<em>U</em><sub><em>t</em> + <em>α</em><em>F</em><em>Δ</em><em>t</em></sub><sup><em>i</em> − 1</sup>)<sup>int</sup> − <em>C</em><em>u̇</em><sub><em>t</em> + <em>α</em><em>F</em><em>Δ</em><em>t</em></sub><sup><em>i</em> − 1</sup> − <em>M</em><em>ü</em><sub><em>t</em> + <em>α</em><em>M</em><em>Δ</em><em>t</em></sub><sup><em>i</em> − 1</sup></span></p>
<p>The linear equations are used to solve for <span
class="math inline"><em>U</em><sub><em>t</em> + <em>α</em><em>F</em><em>Δ</em><em>t</em></sub>, <em>u̇</em><sub><em>t</em> + <em>α</em><em>F</em><em>Δ</em><em>t</em></sub><em>ü</em><sub><em>t</em> + <em>α</em><em>M</em><em>Δ</em><em>t</em></sub></span>.
Once convergence has been achieved the displacements, velocities and
accelerations at time <span
class="math inline"><em>t</em> + <em>Δ</em><em>t</em></span> can be
computed.</p>
<h2 id="references">REFERENCES</h2>
<p>J. Chung, G.M.Hubert. “A Time Integration Algorithm for Structural
Dynamics with Improved Numerical Dissipation: The Generalized-<span
class="math inline"><em>α</em></span> Method” ASME Journal of Applied
Mechanics, 60, 371:375, 1993.</p>
