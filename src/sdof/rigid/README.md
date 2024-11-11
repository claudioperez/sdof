$$
\newcommand{\dt}{\Delta t \,}
\newcommand{\VelocL}{\bm{w}}
\newcommand{\velmat}{\bm{w}}
\newcommand{\velspa}{\bm{\omega}}
\newcommand{\VelocR}{\bm{\omega}}
\newcommand{\DisplR}{\bm{\Lambda}}
\newcommand{\matmom}{\bm{\Pi}}
\newcommand{\momspa}{\bm{\pi}}
\newcommand{\accspa}{\bm{\alpha}}
\newcommand{\accmat}{\bm{a}}
\newcommand{\inrmat}{\mathbb{J}}
% \newcommand{\inrmat}{\mathbf{I}}
$$

$$
\begin{aligned}
\dot{\pi} & =\frac{\mathrm{d}}{\mathrm{~d} t}\left[\mathbf{I}_t \velspa\right] \\
& =\dot{\mathbf{I}}_t \velspa+\mathbf{I}_{\mathbf{t}} \dot{\velspa} \\
& =\velspa \times \mathbb{I}_t \velspa+\mathbf{I}_t \dot{\velspa}=\overline{\mathbf{m}}
\end{aligned}
$$

Convected description

$$
\begin{aligned}
\inrmat \dot{\velmat} + \velmat \times \inrmat \velmat & = \boldsymbol{\Lambda}^{\mathrm{t}} \overline{\mathbf{m}} \\
\dot{\boldsymbol{\Lambda}} & =\boldsymbol{\Lambda} \velmat^{\times} \\
\end{aligned}
\qquad
\begin{aligned}
& \mathbf{I}_{\mathbf{t}} \dot{\velspa} + \velspa \times \mathbf{I}_{\mathbf{t}} \velspa = \overline{\mathbf{m}} = \dot{\momspa} \\
& \dot{\mathbf{\Lambda}}= \velspa^{\times} \mathbf{\Lambda} \\
\end{aligned}
$$


$$
\momspa = \mathbf{\Lambda} \matmom
=\mathbf{\Lambda} \mathbf{J} \velmat
$$

$$
\begin{aligned} 
& \dot{\matmom}
=-\left[\mathbf{I}^{-1} \matmom \right] \times \matmom +\mathbf{T} \\ 
& \dot{\bm{\Psi}}=\left(d \exp _{-\Psi}\right)^{-1} \mathbf{I}^{-1} \matmom
\end{aligned}
$$

## Explicit

![alt text](./images/image-4.png)

![alt text](./images/image-6.png)

![alt text](./images/image-7.png)

### SVQ-VLN (ELN)

$$
\left\{\begin{aligned}
\boldsymbol{\Omega}_{k+\frac{1}{2}} & =\boldsymbol{\Omega}_k+\frac{h}{2} \mathbf{J}^{-1}\left(\mathbf{J} \boldsymbol{\Omega}_k \times \boldsymbol{\Omega}_k+\boldsymbol{\tau}\left(\mathbf{Q}_k\right)\right) \\
\mathbf{Q}_{k+1} & =\mathbf{Q}_k \operatorname{cay}\left(\dt \boldsymbol{\Omega}_{k+\frac{1}{2}}\right), \\
\boldsymbol{\Omega}_{k+1} & =\boldsymbol{\Omega}_{k+\frac{1}{2}}+\tfrac{1}{2} \dt \mathbf{J}^{-1}\left(\mathbf{J} \boldsymbol{\Omega}_{k+1} \times \boldsymbol{\Omega}_{k+1}+\tau\left(\mathbf{Q}_{k+1}\right)\right)
\end{aligned}\right.
$$

### SVQ-TLN

$$
\left\{\begin{array}{l}
\left(I+\frac{h}{2} \widehat{\boldsymbol{\Omega}}_{k+1}\right) \mathbf{J} \boldsymbol{\Omega}_{k+1}-\frac{h}{2} \boldsymbol{\tau}\left(\mathbf{Q}_{k+1}\right)=\left(I-\frac{h}{2} \widehat{\boldsymbol{\Omega}}_k\right) \mathbf{J} \boldsymbol{\Omega}_k+\frac{h}{2} \boldsymbol{\tau}\left(\mathbf{Q}_k\right), \\
\mathbf{Q}_{k+1}=\mathbf{Q}_k \operatorname{cay}\left(h \frac{\boldsymbol{\Omega}_k+\boldsymbol{\Omega}_{k+1}}{2}\right) .
\end{array}\right.
$$

### TRAP

$$
\begin{aligned}
\matmom_{t+\Delta t} 
&=\matmom_t+\frac{1}{2}\Delta t \left(
    -\left[\mathbf{I}^{-1} \matmom_t           \right]^{\times} \matmom_t+\mathbf{T}_t
    -\left[\mathbf{I}^{-1} \matmom_{t+\Delta t}\right]^{\times} \matmom_{t+\Delta t}+\mathbf{T}_{t+\Delta t}
\right) \\
\boldsymbol{\Psi}_{t+\Delta t / 2} &=\frac{\Delta t}{2} \mathbf{I}^{-1} \matmom_t, \\
\boldsymbol{\Psi}_{t+\Delta t}     &=\frac{\Delta t}{2} \mathbf{I}^{-1} \matmom_{t+\Delta t}
\end{aligned}
$$

$$
\begin{aligned}
\matmom_{t+\Delta t} 
&=\matmom_t+\frac{1}{2}\Delta t \left(
    -\operatorname{skew}\left[\mathbf{I}^{-1} \matmom_t           \right] \matmom_t+\mathbf{T}_t
    -\operatorname{skew}\left[\mathbf{I}^{-1} \matmom_{t+\Delta t}\right] \matmom_{t+\Delta t}+\mathbf{T}_{t+\Delta t}
\right) \\
\boldsymbol{\Psi}_{t+\Delta t / 2} &=\frac{\Delta t}{2} \mathbf{I}^{-1} \matmom_t, \\
\boldsymbol{\Psi}_{t+\Delta t}     &=\frac{\Delta t}{2} \mathbf{I}^{-1} \matmom_{t+\Delta t}
\end{aligned}
$$

![alt text](images/image-3.png)

![alt text](images/SW.png)

![alt text](images/image.png)

![alt text](images/image-1.png)

![](images/image-1.png)

### `RKMK`

![alt text](images/image-7.png)

### `LIEMID`

![alt text](images/image-6.png)

### `LIEMID[I]` (`IMIDM`)

![alt text](images/image-4.png)
![alt text](images/image-5.png)

RKMK with tableau:
$$
\begin{array}{c|c}
1 / 2 & 1 / 2 \\
\hline & 1
\end{array}
$$

---------

![alt text](./images/LIEMIDE.png)

![alt text](./images/image-3.png)

![alt text](./images/image-2.png)

### `RKMK-TRAP`

![alt text](images/image-2.png)

$$
\left\{\begin{array}{l}
\boldsymbol{\omega}_{k+\frac{1}{2}}=\boldsymbol{\omega}_k+\frac{h}{2} \boldsymbol{J}^{-1}\left(\tau\left(\mathbf{R}_k\right)-\boldsymbol{\omega}_{k+\frac{1}{2}} \times \mathbf{J} \boldsymbol{\omega}_{k+\frac{1}{2}}\right) \\
\mathbf{R}_{k+1}=\mathbf{R}_k \operatorname{cay}\left(\dt  \boldsymbol{\omega}_{k+\frac{1}{2}}\right) \\
\boldsymbol{\omega}_{k+1}=\boldsymbol{\omega}_{k+\frac{1}{2}}+\frac{h}{2} \mathbf{J}^{-1}\left(\tau\left(\mathbf{R}_{k+1}\right)-\boldsymbol{\omega}_{k+\frac{1}{2}} \times \mathbf{J} \boldsymbol{\omega}_{k+\frac{1}{2}}\right)
\end{array}\right.
$$

### `VLV`

$$
\left\{\begin{array}{l}
\boldsymbol{\omega}_{k+\frac{1}{2}}=\boldsymbol{\omega}_k+\frac{h}{2} \mathbf{J}^{-1}\left[\mathbf{J} \boldsymbol{\omega}_{k+\frac{1}{2}} \times \boldsymbol{\omega}_{k+\frac{1}{2}}-\frac{h}{2}\left(\boldsymbol{\omega}_{k+\frac{1}{2}}^T \mathbf{J} \boldsymbol{\omega}_{k+\frac{1}{2}}\right) \boldsymbol{\omega}_{k+\frac{1}{2}}+\tau\left(\mathbf{R}_k\right)\right] \\
\mathbf{R}_{k+1}=\mathbf{R}_k \operatorname{cay}\left(h \boldsymbol{\omega}_{k+\frac{1}{2}}\right) \\
\boldsymbol{\omega}_{k+1}=\boldsymbol{\omega}_{k+\frac{1}{2}}+\frac{h}{2} \mathbf{J}^{-1}\left[\mathbf{J} \boldsymbol{\omega}_{k+\frac{1}{2}} \times \boldsymbol{\omega}_{k+\frac{1}{2}}+\frac{h}{2}\left(\boldsymbol{\omega}_{k+\frac{1}{2}}^T \mathbf{J} \boldsymbol{\omega}_{k+\frac{1}{2}}\right) \boldsymbol{\omega}_{k+\frac{1}{2}}+\tau\left(\mathbf{R}_{k+1}\right)\right]
\end{array}\right.
$$
Solve:
$$
f(\boldsymbol{\omega})=-\boldsymbol{\omega} + \boldsymbol{\omega}_{k+\frac{1}{2}}+\tfrac{1}{2} \dt \left[(\mathbf{J} \boldsymbol{\omega}) \times \boldsymbol{\omega} - \tfrac{1}{2} \dt \left(\boldsymbol{\omega} \cdot \mathbf{J} \boldsymbol{\omega}\right) \boldsymbol{\omega}+\tau(\mathbf{R})\right]
$$
with Jacobian:
$$
J_f(\boldsymbol{\omega})=-\mathbf{I}+\frac{h}{2} \mathbf{J}^{-1}\left[(\mathbf{J} \boldsymbol{\omega})^{\wedge}-\hat{\boldsymbol{\omega}} \mathbf{J}-h \xi \boldsymbol{\omega}^T \mathbf{J}-\frac{h}{2}\left(\boldsymbol{\omega}^T \mathbf{J} \boldsymbol{\omega}\right) \mathbf{I}\right]
$$

### `LIEMID[EA]`

![alt text](images/LIEMIDEA.png)

$$
\left\{\begin{aligned}
\boldsymbol{\Theta}_{k+\frac{1}{2}} & =\tfrac{1}{2} \dt \mathbf{J}^{-1} \exp \left(-\frac{1}{2} \boldsymbol{\Theta}_{k+\frac{1}{2}}\right)\left(\mathbf{J} \boldsymbol{\Omega}_k+\frac{h}{2} \boldsymbol{\tau}\left(\mathbf{Q}_k\right)\right) \\
\mathbf{Q}_{k+\frac{1}{2}} & =\mathbf{Q}_k \exp \left(\boldsymbol{\Theta}_{k+\frac{1}{2}}\right), \\
\boldsymbol{\Omega}_{k+\frac{1}{2}} & =\mathbf{J}^{-1} \exp \left(-\boldsymbol{\Theta}_{k+\frac{1}{2}}\right)\left(\mathbf{J} \boldsymbol{\Omega}_k+\frac{h}{2} \boldsymbol{\tau}\left(\mathbf{Q}_k\right)\right), \\
\boldsymbol{\Theta}_{k+1} & =\frac{h}{2} \mathbf{J}^{-1} \exp \left(-\frac{1}{2} \boldsymbol{\Theta}_{k+1}\right) \mathbf{J} \boldsymbol{\Omega}_{k+\frac{1}{2}}, \\
\mathbf{Q}_{k+1} & =\mathbf{Q}_{k+\frac{1}{2}} \exp \left(\boldsymbol{\Theta}_{k+1}\right), \\
\boldsymbol{\Omega}_{k+1} & =\mathbf{J}^{-1}\left(\exp \left(-\boldsymbol{\Theta}_{k+1}\right) \mathbf{J} \boldsymbol{\Omega}_{k+\frac{1}{2}}+\frac{h}{2} \boldsymbol{\tau}\left(\mathbf{Q}_{k+1}\right)\right),
\end{aligned}\right.
$$

---------------


![alt text](./images/HSTAGM.png)


bb_rkmk_trap
bb_rkmk_trap_wdexp
dyneq_imid
dyneq_imidm
dyneq_trap
dyneq_trapm
dyneq_trapm_zeta
IncrSO3_AKW
LIEMIDEA
liemid_Newmark
mleok
rotint_nmb
simo_wong_algo_c1
simo_wong_algo_c2
SVQ
VLV
