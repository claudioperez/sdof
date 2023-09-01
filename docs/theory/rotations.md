1. Define the updated configuration via the exponential map as:

   $$
   \boldsymbol{\Lambda}_{n+1}=\mathbf{\Lambda}_n \operatorname{Exp} \mathbf{\Theta}
   $$

   where $\Theta \in \mathbb{R}^3$ is the convected relative (incremental) rotation vector.

2. Define the convected relative (incremental) rotation $\Theta \in \mathbb{R}^3$ in terms of $\left(\mathbf{W}_n, \mathbf{A}_n, \mathbf{A}_{n+1}\right)$ by the formula

   $$
   \boldsymbol{\Theta}=h \mathbf{W}_n+h^2\left[\left(\frac{1}{2}-\beta\right) \mathbf{A}_n+\beta \mathbf{A}_{n+1}\right]
   $$

   where $\beta \in\left[0, \frac{1}{2}\right]$ is a parameter with identical significance as in the classical Newmark algorithm.

3. Define the updated convected angular velocity by the formula

   $$
   \mathbf{W}_{n+1}=\mathbf{W}_n+h\left[(1-\gamma) \mathbf{A}_n+\gamma \mathbf{A}_{n+1}\right]
   $$

   where $\gamma \in[0,1]$ is a parameter with identical significance as in the classical Newmark algorithm.

4. Enforce rate of momentum balance at $t_n+\gamma h$ :

   $$
   \boldsymbol{J} \mathbf{A}_{n+\gamma}+\mathbf{W}_{n+\gamma} \times \boldsymbol{J} \mathbf{W}_{n+\gamma}=\boldsymbol{\Lambda}_{n+\gamma}^{\mathrm{t}} \overline{\mathbf{m}}_{n+\gamma}
   $$

