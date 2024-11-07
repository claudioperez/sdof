import numpy as np
import matplotlib.pyplot as plt
from gamma import average_gamma, newmark_gamma, collocate
from wilson import Wilson_theta
from integrate import GammaScheme, integrate

if __name__ == "__main__":
    import sdof
    import sys 

    f = np.array([
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

    m, k, c = 0.2533, 0.1592, 0.0 #, 10.

    m, k, c = 1, 2, 0 #10
    T = 2 * np.pi * np.sqrt(m/k)

    dt = 10*T #0.10 # size of time step
    nt = 40 # number of time steps
    f = np.zeros(nt)
    v0 = 0 #1/m
    u0 = 1
    Ey = "a"
    if len(sys.argv) > 2:
        Ey = sys.argv[2]

    beta  = 0.3025
    gamma = 0.6
#   m /= 2

    model = (k, c, m)

    Gn, gn = newmark_gamma(Ey, dt, beta, gamma)


    #
    scheme = GammaScheme(Ey, Gn, gn, Gn, gn, predictor=Ey)

    u, v, a = integrate(scheme, model, f, v0=v0, u0=u0)

    plt.plot(u, "o", label="newmark")

    #

#   beta   = 0.2
#   gamma  = 1/2
#   theta  = 1.0
#   u, v, a = Wilson_theta((k, c, m),dt,f,u0=u0,v0=v0,beta=beta,gamma=gamma,theta=theta)
#   plt.plot(u, "-", label=f"$\\{theta = }$")
#   plt.legend()
#   plt.show()
#   sys.exit()

    #
    alpha = 1 + (-0.1)
    alpha_a = 1
    gamma = 1/2 + alpha_a - alpha
    beta = (1 + alpha_a - alpha)**2/4
    Ga, ga = average_gamma([alpha, 1, alpha_a], (Gn, gn))
    scheme = GammaScheme(Ey, Gn, gn, Ga, ga, alpha, predictor=Ey)

    u, v, a = integrate(scheme, model, f, v0=v0, u0=u0)

    plt.plot(u, "x", label=f"$\\{alpha = }$")

    #
    beta   = 1/6 #0.24
    gamma  = 1/2
    beta  = 0.3025
    gamma = 0.6

    theta  = 1.021712

    Gn, gn = newmark_gamma(Ey, dt, beta, gamma)
    Gc, gc = collocate(Ey, dt, theta, beta, gamma)
    scheme = GammaScheme(Ey, Gn, gn, Gc, gc, theta, predictor=Ey)

    u, v, a = integrate(scheme, model, f, v0=v0, u0=u0)

    plt.plot(u, "-", label=f"$\\{theta = }$")


#   # Control
#   u, v, a = sdof.integrate(f, dt, k, c, m,
#                            beta=beta, gamma=gamma, u0=u0, v0=v0
#                            # const avg accel
#                  beta=1/6,  gamma=0.5, # linear accel
#   )
#   plt.plot(u, ".")

    plt.legend()

    plt.show()
