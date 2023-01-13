import sdof
import numpy as np


print(sdof.peaks(
        1., 1., 0.,
        np.loadtxt("data/elCentro.txt"), 0.01
).max_accel)

f =     np.array([
         0.0000,
         5.0000,
         8.6603,
        10.0000,
         8.6603,
         5.0000,
         0.0000,
         0.0000,
         0.0000,
         0.0000,
         0.0000])

from pandas import DataFrame
U = sdof.integrate(0.2533, 0.1592, 10., f, 0.10)
print(DataFrame(U.T))

import matplotlib.pyplot as plt

u,v,a = sdof.integrate(0.2533, 0.0, 10., np.zeros(100), 0.01, u0=1.0)

fig, ax = plt.subplots()
# ax.plot(a)
# ax.plot(v)
ax.plot(u)
plt.show()

u,v,a = sdof.integrate(
        1., 1.e3, 0.,
        -np.sin(np.linspace(0, 5*np.pi, 200)), 5*np.pi/200
)

fig, ax = plt.subplots()
ax.plot(a)
ax.plot(v)
ax.plot(u)
plt.show()

f = np.loadtxt("data/elCentro.txt")
u,v,a = sdof.integrate(
        1., 1., 0.,
        f, 0.01
)
fig, ax = plt.subplots()
ax.plot(a)
ax.plot(v)
ax.plot(u)
ax.plot(f)
plt.show()

