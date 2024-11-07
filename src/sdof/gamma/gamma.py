"""

newmark_gamma()
average_gamma()
collocate_gamma()


gsss = average(Ey, alpha,
               newmark(Ey,...),
               newmark(Ey,...))

alph = average(alpha,
               newmark(Ey,...))

wils = collocate(Ey, theta, (),
                 lambda E,dt: newmark)
"""
import numpy as np


def average_gamma(alpha, scheme_a, scheme_b=None):
    G, g = scheme_a
    A = np.diag(alpha)
    Ga = A@(G - np.eye(3)) + np.eye(3)
    ga = A@g
    return Ga, ga


def collocate(Ey, dt, theta, *params):
    E = Ey
#   Ey = "a"
    Gc, gc = newmark_gamma("a", dt*theta, *params)
    Gy, gy = newmark_gamma("a", dt, *params)

    iy =  "uva".index(Ey)
    ia = 2

    G = [[None for j in range(3)] for i in range(3)]
    g = [None for i in range(3)]
    for i in range(3):
        for j in range(3):
            G[i][j] = Gc[i][j] + gc[i]*((1-theta)*float(j == ia and i!=ia) - theta/gy[iy]*Gy[iy][j])
        g[i] = gc[i]*theta/gy[iy]

    return G, g #convert((G, g), "a", Ey)


def newmark_gamma(Ey, dt, beta, gamma, betas=None, betass=None, gammas=None):
    if gammas is None:
        gammas = (1 -  gamma)

    if betass is None:
        betass = 1

    if betas is None:
        betas = (1 - 2*beta)/2

    Ga = [[1, betass*dt,  betas*dt**2],
          [0,   1,           gammas*dt],
          [0,   0,                   0]]

    ga = [beta*dt**2, gamma*dt, 1]

    if Ey != "a":
        return convert((Ga, ga), "a", Ey)
    else:
        return Ga, ga


def convert(scheme, Ey, Ez, other=None):
    if Ey == Ez:
        return scheme

    Gy, gy = scheme
    iy = "uva".index(Ey)
    iz = "uva".index(Ez)
    iw = "uva".index("uva".replace(Ey,"").replace(Ez,""))

    mu = 1/gy[iz]
    gz = [gi*mu for gi in gy]
    Gz = [

            [((Gy[iw][j] if i!=iy else 0)) - mu*gy[i]*Gj for j,Gj in enumerate(Gy[iz]) ] if i!=iz else
            [0, 0, 0]

            for i in range(3)
    ]
    return Gz, gz

