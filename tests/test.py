import sdof
import numpy as np


print(sdof.peaks(
        1., 1., 0.,
        np.loadtxt("data/elCentro.txt"), 0.01
).max_accel)



print(sdof.integrate(
        1., 1., 0.,
        np.loadtxt("data/elCentro.txt"), 0.01
))

