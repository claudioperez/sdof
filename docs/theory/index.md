# Theory

We are looking to approximate second order differential equations of the form 

$$\mathbf{M} \ddot{\boldsymbol{u}} + \mathbf{C}\dot{\boldsymbol{u}} + \boldsymbol{p}(\boldsymbol{u}, \dot{\boldsymbol{u}}) = \boldsymbol{f}$$

## Newmark Equations

The Newmark-$\beta$ scheme is defined by the following equations:

<!--
$$
\begin{aligned}
d_{n+1} &= d_n + hv_n + h^2 \left(\left(\frac{1}{2}-\beta \right)a_n + \beta a_{n+1}  \right)
\end{aligned}
$$
-->

$$
\begin{aligned}
&\boldsymbol{v}_{n+1}= \boldsymbol{v}_n + (1 - \gamma) \Delta t~ \boldsymbol{a}_n + \gamma \Delta t\,  \boldsymbol{a}_{n+1}\\
&\boldsymbol{d}_{n+1}= \boldsymbol{d}_n + \Delta t~\boldsymbol{v}_n + \tfrac{1}{2}\Delta t^2 \left(\left(1 - 2\beta\right) a_n + 2\beta \boldsymbol{a}_{n+1}\right)\\
%&m  \boldsymbol{a}_{n+1} + c \boldsymbol{v}_{n+1} + p(\boldsymbol{d}_{n+1}) = \boldsymbol{f}_{n+1}^{\textrm{ext}}
\end{aligned}
$$

- Using the extended mean value
  theorem, the Newmark-$\beta$
  method states that the velocity (i.e., first time derivative)
  can be decomposed as,

  $$
  \begin{aligned}
  \dot{\boldsymbol{u}}_{n+1} &= \dot{\boldsymbol{u}}_n + \Delta t~ \ddot{\boldsymbol{u}}_\gamma \\
  \boldsymbol{d}_{n+1}       &= \boldsymbol{u}_n       + \Delta t~ \dot{\boldsymbol{u}}_n + \begin{matrix} \frac 1 2 \end{matrix} \Delta t^2~ \ddot{\boldsymbol{u}}_\beta
  \end{aligned}
  $$

  where

  $$
  \begin{aligned}
   \ddot{\boldsymbol{u}}_\gamma &= (1 - \gamma) \ddot{\boldsymbol{u}}_n + \gamma  \ddot{\boldsymbol{u}}_{n+1}~~~~~~~~0\leq \gamma \leq 1\\
   \ddot{\boldsymbol{u}}_\beta  &= (1 - 2\beta) \ddot{\boldsymbol{u}}_n + 2\beta  \ddot{\boldsymbol{u}}_{n+1}~~~~0\leq 2\beta\leq 1
  \end{aligned}
  $$



- **Explicit central difference scheme** is obtained by setting
  $\gamma=0.5$ and $\beta=0$
  
- **Average constant acceleration (Middle point rule)** is obtained by
  setting $\gamma=0.5$ and $\beta=0.25$


$$
\frac{1}{\Delta t}\left(v_{n+1} - v_n\right) = (1-\gamma) a_n + \gamma a_{n+1}
$$



## Generalized $\alpha$ Scheme

The generalized $\alpha$ scheme imposes equilibrium with the following form:

$$\begin{aligned}
\mathbf{M}\, a _{\alpha_{m}} + \boldsymbol{p}(d_{\alpha_{f}}, v_{\alpha_{f}}) &= \boldsymbol{f}(t_{\alpha_{f}})
\end{aligned}$$

where

$$\begin{aligned}
t_{\alpha_{f}}              &= (1-\alpha_{f}) \, t_{n+1}              + \alpha_{f} \, t_n \\
\boldsymbol{d}_{\alpha_{f}} &= (1-\alpha_{f}) \, \boldsymbol{d}_{n+1} + \alpha_{f} \, \boldsymbol{d}_n \\
\boldsymbol{v}_{\alpha_{f}} &= (1-\alpha_{f}) \, \boldsymbol{v}_{n+1} + \alpha_{f} \, \boldsymbol{v}_n \\
\boldsymbol{a}_{\alpha_{m}} &= (1-\alpha_{m}) \, \boldsymbol{a}_{n+1} + \alpha_{m} \, \boldsymbol{a}_n
\end{aligned}$$


## Implementation


The Newmark equations allow one to express all three unknowns (i.e., $a_{n+1}, v_{n+1}$ or $d_{n+1}$)
in terms of one of these. 

### Linear Newmark Schemes

The most straight-forward way to implement the Newmark scheme, is to form
a single equation in terms of unknown accelerations, $\boldsymbol{a}_{n+1}$.

$$
\mathbf{A} \boldsymbol{a}_{n+1} = \boldsymbol{b}
$$

This simply follows from plugging the Newmark equations into the discrete
equilibrium statement at $t_{n+1}$. However, before this, it is useful rewrite
equations [@eqn:newmark] as follows:

$$
\begin{aligned}
\boldsymbol{d}_{n+1} &= \tilde{\boldsymbol{d}} + c_{da} \, \boldsymbol{a}_{n+1} \\
\boldsymbol{v}_{n+1} &= \tilde{\boldsymbol{v}} + c_{va} \, \boldsymbol{a}_{n+1}
\end{aligned}
$$

where the tilde variables allow us to collect everything that is known at the start of a time step.
Expanding and plugging into the discretized equilibrium equation furnishes
the Newmark scheme in "a-form":

$$\begin{aligned}
\tilde{\boldsymbol{d}} &= \boldsymbol{d}_n + \Delta t \boldsymbol{v}_n + \tfrac{1}{2} \Delta t^2 (1-2\beta) \boldsymbol{a}_n \\ 
\tilde{\boldsymbol{v}} &= \boldsymbol{v}_n + \Delta t (1-\gamma) \boldsymbol{a}_n \\ 
(\mathbf{M} + \gamma\Delta t\, \mathbf{C}+\beta \Delta t^2\, \mathbf{K})\boldsymbol{a}_{n+1} &= \boldsymbol{f}_{n+1} - \mathbf{C}\tilde{\boldsymbol{v}} - \mathbf{K}\tilde{\boldsymbol{d}} \\
\end{aligned}$$

$$\begin{aligned}
\boldsymbol{d}_{n+1} &= \tilde{\boldsymbol{d}} + \beta \Delta t^2 \, \boldsymbol{a}_{n+1} \\
\boldsymbol{v}_{n+1} &= \tilde{\boldsymbol{v}} + \gamma \Delta t \, \boldsymbol{a}_{n+1}
\end{aligned}$$

### Generalized Unknowns

Alternatively, the Newmark equations can be manipulated to produce
a problem in terms of velocity or displacement.
This generalizes as follows:

$$\begin{aligned}
\tilde{\boldsymbol{d}} &= \sum_i b_{di}\boldsymbol{y}_i  \\ 
\tilde{\boldsymbol{v}} &= \sum_i b_{vi}\boldsymbol{y}_i  \\ 
\tilde{\boldsymbol{a}} &= \sum_i b_{ai}\boldsymbol{y}_i  \\ 
(c_a\,\mathbf{M} + c_v\, \mathbf{C} + c_d\, \mathbf{K})\boldsymbol{x}_{n+1} &= \boldsymbol{f}_{n+1} - \left(\mathbf{M}\tilde{\boldsymbol{a}} + \mathbf{C}\tilde{\boldsymbol{v}} + \mathbf{K}\tilde{\boldsymbol{d}} \right) \\
\end{aligned}$$

$$\begin{aligned}
\boldsymbol{d}_{n+1} &= \tilde{\boldsymbol{d}} + c_d  \, \boldsymbol{x}_{n+1} \\
\boldsymbol{v}_{n+1} &= \tilde{\boldsymbol{v}} + c_v  \, \boldsymbol{x}_{n+1} \\
\boldsymbol{a}_{n+1} &= \tilde{\boldsymbol{a}} + c_a  \, \boldsymbol{x}_{n+1}
\end{aligned}$$


### Nonlinear Generalized - $\alpha$

If $\boldsymbol{x}_{n+1}$ denotes our chosen unknown from these
three, we can write the equations in the following form:

$$
\begin{aligned}
\boldsymbol{d}_{\alpha} &= \tilde{\boldsymbol{d}}_x + c_{dx} \, \boldsymbol{x}_{n+1} \\
\boldsymbol{v}_{\alpha} &= \tilde{\boldsymbol{v}}_x + c_{vx} \, \boldsymbol{x}_{n+1} \\
\boldsymbol{a}_{\alpha} &= \tilde{\boldsymbol{a}}_x + c_{ax} \, \boldsymbol{x}_{n+1} \\
\end{aligned}
$$

for some coefficients $c_{yx}$, and 
where $\tilde{\left(\cdot\right)}_x$ variables encapsulate the information that is known at the start of the time step.

Applying this to a nonlinear problem yields:

$$
\begin{aligned}
r(\boldsymbol{x}_{n+1})&= \mathbf{M}\,  \left(\tilde{\boldsymbol{a}} + c_{ax}\,\boldsymbol{x}_{n+1}\right)
                        + \mathbf{C}\,  \left(\tilde{\boldsymbol{v}} + c_{vx}\,\boldsymbol{x}_{n+1}\right)
                        + \boldsymbol{p}\left(\tilde{\boldsymbol{d}} + c_{dx}\,\boldsymbol{x}_{n+1},
                                              \tilde{\boldsymbol{v}} + c_{vx}\,\boldsymbol{x}_{n+1}\right)
                        - \boldsymbol{f}(t_{n+1})
\end{aligned}
$$

Linearizing for $r_\eta = r(\boldsymbol{x}_{n+1}^i + \eta \, d\boldsymbol{x})$ yields

$$
\begin{aligned}
\mathcal{L} \, r_{\eta} &= r(\boldsymbol{x}^i_{n+1}) + \mathbf{A} d \boldsymbol{x} \\
\end{aligned}
$$

where

$$
\mathbf{A} d\boldsymbol{x} = \left.\frac{d}{d \eta}r_\eta\right|_{\eta=0} = \left( c_d \mathbf{K}\left(\tilde{\boldsymbol{d}}+c_{d}\,\boldsymbol{x}_{n+1}\right) + c_v \mathbf{C} + c_a \mathbf{M}\right) d \boldsymbol{x}
$$


### Displacements and Velocity as Unknown


$$
\boldsymbol{y}_{n+1} =  \mathsf{B}_x \, \boldsymbol{y}_n + \mathsf{c}_{yx}\, \boldsymbol{x}_{n+1}
$$

$$
c_{xy} = \begin{pmatrix}
1 & \frac{\gamma}{\beta \Delta t} & \frac{1}{\beta \Delta t^2} \\
\frac{\Delta t \beta}{\gamma} & 1 & \frac{1}{\gamma \Delta t}  \\
\beta \Delta t^2 & \gamma \Delta t & 1
\end{pmatrix}
$$

$$
\begin{pmatrix}
\tilde{ d} \\
\tilde{ v} \\
\tilde{ a}
\end{pmatrix}_a =
\begin{pmatrix}
1 & \Delta t & \Delta t^2 \left(\tfrac{1}{2} - \beta\right) \\ 
& 1 & \Delta t \left(1 - \gamma\right) \\
&& 0
\end{pmatrix}
\begin{pmatrix}
 d \\ 
 v \\ 
 a
\end{pmatrix}_n
$$


$$
\begin{pmatrix}
\tilde{ d} \\
\tilde{ v} \\
\tilde{ a}
\end{pmatrix}_v =
\begin{pmatrix}
1 & - \Delta t\frac{\beta}{\gamma}\left(1 - \frac{\gamma}{\beta}\right) & \Delta t^2 \frac{\beta}{\gamma}\left(\frac{\gamma}{2\beta} - 1 \right) \\ 
 & 0 & \\
&\frac{-1}{\gamma \Delta t} & 1 - \frac{1}{\gamma}
\end{pmatrix}
\begin{pmatrix}
 d \\ 
 v \\ 
 a
\end{pmatrix}_n
$$


$$
\begin{pmatrix}
\tilde{ d} \\
\tilde{ v} \\
\tilde{ a}
\end{pmatrix}_d =
\begin{pmatrix}
0 \\ 
-\frac{\gamma}{\beta\Delta t} & 1 - \frac{\gamma}{\beta} & \Delta t \left(1 - \frac{\gamma}{2\beta}\right) \\
\frac{1}{\beta \Delta t^2} & \frac{-1}{\beta\Delta t} & \quad \left(1 + \tfrac{1}{2\beta}\right)
\end{pmatrix}
\begin{pmatrix}
 d \\ 
 v \\ 
 a
\end{pmatrix}_n
$$



--------------------

<!--

## Parameters

By properly choosing the parameters, we can recover HHT, Newmark, or WBZ methods. 

Parameters can be specified by taking $\rho_\infty\in [0,1]$ as a design
variable to control the numerical dissipation above the normal frequency
$\frac{h}{T}$, where $T$ is the period associated with the highest frequency of
interest. The following relationships are used to obtain a good algorithm that
are accurate and preserve low-frequency modes

$$\begin{aligned}
\alpha_{m} &= \frac{2\rho_\infty-1}{\rho_\infty+1}\\
\alpha_{f} &= \frac{\rho_\infty}{\rho_\infty+1}   \\
\gamma &= \tfrac{1}{2} - \alpha_{m} + \alpha_{f}   \\
\beta  &= \tfrac{1}{4} (1-\alpha_{m}+\alpha_{f})^2
\end{aligned}$$

-->

<!--

## Algorithm

$$
\begin{aligned}
 q_{n+1}&= q_n+h \dot{ q}_n+h^2(0.5-\beta)  a \\
\dot{ q}_{n+1}&=\dot{ q}_n+h(1-\gamma)  a \\
 a&=1 /\left(1-\alpha_m\right)\left(\alpha_f \ddot{ q}_n-\alpha_m  a\right) \\
 q_{n+1}&= q_{n+1}+h^2 \beta  a \\
\dot{ q}_{n+1}&=\dot{ q}_{n+1}+h \gamma  a \\
\ddot{ q}_{n+1}&= 0
\end{aligned}
$$

$$
\begin{aligned}
& \text { for } i=1 \text { to } i_{\max } \text { do } \\
& {\left[\begin{array}{c}
\Delta  q \\
\Delta \boldsymbol{\lambda}
\end{array}\right]:=- S_t^{-1}\left[\begin{array}{l}
 r^{ q} \\
 r^{\boldsymbol{\lambda}}
\end{array}\right]} \\
&  q_{n+1}:= q_{n+1}+\Delta  q \\
& \dot{ q}_{n+1}:=\dot{ q}_{n+1}+\gamma^{\prime} \Delta  q \\
& \ddot{ q}_{n+1}:=\ddot{ q}_{n+1}+\beta^{\prime} \Delta  q \\
& \boldsymbol{\lambda}_{n+1}:=\boldsymbol{\lambda}_{n+1}+\Delta \boldsymbol{\lambda} \\
& \text { end for } \\
&  a:= a+\left(1-\alpha_f\right) /\left(1-\alpha_m\right) \ddot{ q}_{n+1} \\
&
\end{aligned}
$$

-->

## References

- J. Chung, G.M.Hubert. "A Time Integration Algorithm for Structural Dynamics with Improved Numerical Dissipation: The Generalized-$\alpha$ Method" ASME Journal of Applied Mechanics, 60, 371:375, 1993.

