
## Ordinary Differential Equations

### Bouc-Wen

Models associated with the Bouc-Wen family are generally given by a relation with the form:

> $$\sigma(\varepsilon) = b E \varepsilon+(1-b) E z(\varepsilon)$$
> $$\dot{z}=\frac{\partial z}{\partial \varepsilon} \frac{\partial \varepsilon}{\partial t} = f(z,\varepsilon,\dot \varepsilon) \dot{\varepsilon}$$

Most variations of the Bouc-Wen model differ in the function $f$.

| Source                     | ^[sym] | ^[slp] | ^[dmg] | ^[]  | Notes  |   
|--------------------------|--|--|--|--|-------------------------------------------------------|
| @bouc1967forced            | - | - | - | - | |
| @wen1976                   | - | - | - | - | |
| @baber1981random           | - | a | - | - | Introduction of degradation to @wen1976 |
| @baber1985random           | - | a | a | - | Addition of pinching to @baber1981random using series composition with a Gaussian.
| @baber1986                 | - | a | b | - |
| @sues1988                  | - | b | - | - | New degradation function proposed for RC, similar to Clough's
| @foliente1995hysteresis    | - | a | c | - | Proposed modification to pinching function of @baber1986 
| @reinhorn1995modeling      |   |   |   |   | 
| @sivaselvan2000hysteretic  |   |   |   |   |
| @charalampakis2009bouc     |   |
| @kottari2014       |
| @pelliciari2019    |
| @liberatore2019    |
| @amir2022statespace        |

<!--
| @song2006
-->

\newpage

Bouc's Hysteresis (1967)
:   ~

  > $$
    f(z,\dot{\varepsilon})=A-\alpha \operatorname{sgn}(\dot{\varepsilon}) z
    %\dot{z}=(A-\alpha \operatorname{sgn}(\dot{\varepsilon}) z) \dot{\varepsilon} \\
  > $$

  Bouc's model yields a separable ODE for $\alpha \ne 0$ with solutions:
  $$
  z(\varepsilon)= \begin{cases}-\frac{1}{\alpha}\left(C_1 e^{-\alpha \varepsilon}-A\right), & \text { if } \dot{\varepsilon}>0 \\ \frac{1}{\alpha}\left(C_2 e^{\alpha \varepsilon}-A\right), & \text { if } \dot{\varepsilon}<0\end{cases}
  $$
  Where $C_1$ and $C_2$ must be calculated piecewise by using
    the final state of the previous branch as the initial condition.


Wen's Hysteresis (1976)

:   Wen extended the original model by Bouc with the following expression:

  > $$
    %\dot{z}=\left\{A-[\alpha \operatorname{sgn}(z \dot{\varepsilon})-\beta]|z|^{n}\right\} \dot{\varepsilon} \\
    f(z,\dot{\varepsilon}) = \frac{d z}{d \varepsilon}=A-[\alpha \operatorname{sgn}(z \dot{\varepsilon})-\beta]|z|^{n}
  > $$


Baber-Wen Degradation (1981)

:  Baber and Wen (1981) present an extension to the Bouc-Wen model
  that includes strength and stiffness deterioration.

  > $$
    f(z,\dot{\varepsilon})=\frac{A-|z|^{n}\{\alpha \operatorname{sgn}(\dot{\varepsilon} z)+\beta\} \nu}{\eta}
  > $$
  > $$
      \begin{aligned}
    A &=A_{o}-\delta_{A} \psi \\
    \nu &=1+\delta_{\nu} \psi \\
    \eta &=1+\delta_{\eta} \psi \\
    \dot{\psi} & =(1-b) E \dot{\varepsilon} z
    \end{aligned}
    > $$
    > $A_{o}$, $\delta_{A}$, $\delta_{\nu}$  and $\delta_{\eta}$ are degradation parameters. 

### Other Models

S&R's Arrangement and Asymmetry
: 


Perfect Plasticity

:  One dimensional perfect plasticity can formulated as the following ODE

  >$$
  \begin{aligned}
  \dot{\sigma} &=E\left[1-\frac{1+\operatorname{sgn}(\sigma \dot{\varepsilon})}{2}\left\{U\left(\sigma_{y}^{+}-\sigma\right) U(\sigma)
  +U\left(-\sigma_{y}^{-}+\sigma\right) U(-\sigma)\right\}\right] \dot{\varepsilon}
  \end{aligned}
  >$$


[sym]: a
[slp]: Slipping
[dmg]: Degradation

### Implementations


OpenSees `BoucWen` and FEDEAS `BoucWen1dMat`

:   Bouc-Wen-Baber (1981) with backward Euler and Newton-Raphson.

  > Application of backward Euler scheme to $\dot e$ and $\dot z$:
  > $$
    \begin{aligned}
    z_{(n+1)}&=z_{(n)}+\Delta t \frac{A_{(n+1)}-\left|z_{(n+1)}\right|^{n}\left\{\gamma+\beta \operatorname{sgn}\left(\frac{\left(\epsilon_{(n+1)}-\epsilon_{(n)}\right)}{\Delta t} z_{(n+1)}\right)\right\} \nu_{(n+1)}}{\eta_{(n+1)}} \frac{\left(\epsilon_{(n+1)}-\epsilon_{(n)}\right)}{\Delta t} \\
    e_{(n+1)}
    &=e_{(n)}+\Delta t(1-\alpha) k_{o} \frac{\left(\epsilon_{(n+1)}-\epsilon_{(n)}\right)}{\Delta t} z_{(n+1)} \\
    \end{aligned}
  > $$
  > 
  > Note that the evolution equations are applied discretely with forms:
  > $$
    \begin{aligned}
  A_{(n+1)} &=A_{o}-\delta_{A} e_{(n+1)} \\
  \nu_{(n+1)} &=1+\delta_{\nu} e_{(n+1)} \\
  \eta_{(n+1)} &=1+\delta_{\eta} e_{(n+1)}
  \end{aligned}
  > $$

### Observations

- Asymmetric load cycles never form a stable closed loop (); ie, asymmetric load cycles always cause drift.
- Stiffness at unloading is always greater than stiffness at loading.
