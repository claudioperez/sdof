import sdof
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(3,1,sharex=True)

n = 25

for alpha in [0.8, 1.0]:
    for beta in [0.25, 0.3, 0.5]:
        for gamma in [1.0]:
            u,v,a = sdof.integrate(
                    1., 1.e3, 0.,
                    -np.sin(np.linspace(0, 8*np.pi, n)),
                    8*np.pi/n,
                    alpha_f = alpha,
                    alpha_m = alpha,
                    beta    = beta,
                    gamma   = gamma
            )

            ax[0].plot(a)
            ax[1].plot(v)
            ax[2].plot(u)

ax[0].set_ylabel("a")
ax[1].set_ylabel("v")
ax[2].set_ylabel("u")

plt.show()
