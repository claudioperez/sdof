
from sdof.rigid import *

class E61:
    def __init__(self, J, w0, R0, C1, C2, tstar):
        self.C1 = C1 
        self.C2 = C2 
        self.tstar = tstar

    def torque(self, t, R):
        if t < self.tstar:
            return [self.C1, 0, 0]

if __name__ == "__main__":
    import numpy as np
    import time
    import sys
    import matplotlib.pyplot as plt

    vm  = np.array([1, 0.0, 1])
    vm  = vm/np.linalg.norm(vm)
    Gm  = ExpSO3(2.5*vm); # attraction point
    alpha = 0.3  # magnetic potential

    def potential(G, Gm):
        # Potential Energy
        dm = np.sqrt(2*(3 - np.trace(Gm@G)))
        de = np.sqrt(2*(3 - np.trace(G)))
        V = (de - 1.0)**2 - alpha/dm
        return V

    def potentialtorque(t, G):
        """
        Left Trivialized Potential Torque

        Parameters:
        t (float): time
        G (np.array): 3x3 matrix
        Loading (object): object with attributes alpha and Gm

        Returns:
        np.array: 3x1 torque vector
        """

        dm = np.sqrt(2 * (3 - np.trace(Gm.T@G)))
        de = np.sqrt(2 * (3 - np.trace(G)))

        tau = 2.0 * (1.0 - 1.0 / de) * np.array([
            G[1, 2] - G[2, 1],
            G[2, 0] - G[0, 2],
            G[0, 1] - G[1, 0]
        ]) + alpha * 1.0 / dm**3 * np.array([
            G[0, 2] * Gm[0, 1] - G[0, 1] * Gm[0, 2] - G[1, 1] * Gm[1, 2] + G[1, 2] * Gm[1, 1] - G[2, 1] * Gm[2, 2] + G[2, 2] * Gm[2, 1],
            G[0, 0] * Gm[0, 2] - G[0, 2] * Gm[0, 0] + G[1, 0] * Gm[1, 2] - G[1, 2] * Gm[1, 0] + G[2, 0] * Gm[2, 2] - G[2, 2] * Gm[2, 0],
            G[0, 1] * Gm[0, 0] - G[0, 0] * Gm[0, 1] - G[1, 0] * Gm[1, 1] + G[1, 1] * Gm[1, 0] - G[2, 0] * Gm[2, 1] + G[2, 1] * Gm[2, 0]
        ])

        return tau

    # Parameters and initial conditions
    J = np.diag([2.0, 2.0, 4.0])

    dt = 0.25  # timestep
    hc = 1     # coarse timestep
    T = 10000./1  # total simulation time

    # Loading
    Wi0 = np.array([0.0, 0.0, 0.625])  # initial body angular velocity
    Gi0 = ExpSO3([0.0, 0.7227, 0.0])   # initial configuration (identity matrix)

    Loading = Source(Torque = potentialtorque)

    # Methods to be analyzed
    Methods = [
#       ("dyneq_trap",  "incrso3_trap"),

#       ("incrso3_trapm",), #  "incrso3_trapm_old"),

#       ("incrso3_svq", "incrso3_svq_alt"), # "rotint_nmb"),

#       ("incrso3_imid", "kb_e"),
        ("incrso3_mleok",),

#       ("bb_rkmk_trap", "incrso3_akw"),
        (
            "incrso3_vlv",
            #
            # 

            # "IncrSO3_imidm", 
        ),
        (

          "IncrSO3_LIEMIDEA",

            # Broken
    #       "bb_rkmk_trap_wdexp",
    #       "liemid_newmark",
    #       "rotint_nmb",
            # "incrso3_swc2",
        )
    ]

    for methods in Methods:
        markers = iter("-ox.:")
        for method in methods:
            marker = next(markers)
            if True:
                step = getattr(sys.modules[__name__], method.lower())

                start_time = time.time()
                print(method)
                t, _, _, _, E = integrate(J, dt, hc, T, Wi0, Gi0, Gm, step,
                                                Loading, potential)
                end_time = time.time()
                print(f"Elapsed time: {end_time - start_time} seconds")

                plt.plot(t, E - E[0],  marker, label=method)
            try:
                pass
            except Exception as e:
                print(e)

    plt.legend()
    plt.show()

