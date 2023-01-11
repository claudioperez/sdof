import sdof
import numpy as np


print(sdof.peaks(
        1., 1., 0.,
        np.loadtxt("data/elCentro.txt"), 0.01
).max_accel)



sdof.integrate(
        1., 1., 0.,
        np.loadtxt("data/elCentro.txt"), 0.01
)


import matplotlib.pyplot as plt
u,v,a = sdof.integrate(
        1., 1., 0.,
        np.sin(np.linspace(0, 5*np.pi, 200)), 5*np.pi/200
)
fig, ax = plt.subplots()
ax.plot(u)
ax.plot(v)
ax.plot(a)
plt.show()
