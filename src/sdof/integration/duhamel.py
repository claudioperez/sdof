

import numpy as np
#
# DUHAMEL
#

def duhamel(motion, dt, period, zeta):
    """
    single_elastic_response
    Perform Duhamels integral to get the displacement.
    http://www.civil.utah.edu/~bartlett/CVEEN7330/Duhamel%27s_integral.pdf
    http://www1.aucegypt.edu/faculty/mharafa/MENG%20475/Forced%20Vibration.pdf

    :param motion: acceleration
    :param zeta: damping ratio
    """
    w_n = (2.0 * np.pi) / period
    w_d = w_n * np.sqrt(1 - zeta ** 2)
    zeta_wn = zeta * w_n
    nt = len(motion)

    time = dt * np.arange(nt + 1)
    u = np.zeros(nt)

    p = motion * dt / w_d

    for i in range(nt):
        dtn = time[:-i - 1]
        u[i:] = p[i] * np.exp(-zeta_wn*dtn) * np.sin(w_d * dtn)

    return u

