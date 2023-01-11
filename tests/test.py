import sdof
import numpy as np


print(sdof.peaks(
        1., 1., 0.,
        np.loadtxt("data/elCentro.txt"), 0.01
))

