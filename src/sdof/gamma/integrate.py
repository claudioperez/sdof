import numpy as np
class GammaScheme:
    def __init__(self, Ey, Gn, gn, Ga=None, ga=None, lam=1, predictor=None):
        self.Ey = Ey
        self.iy = "uva".index(Ey)

        self.Gn = Gn
        self.gn = gn

        self.Ga = Ga
        self.ga = ga
        self.lam = lam


def integrate(self, model, f, u0=0, v0=0, maxiter=10):
    nt = len(f)
    x = np.zeros((nt, 3))

    x[0] = initialize(model, u0, v0)

    ip = self.iy

    for i in range(nt-1):
        fa = increment(self, i, f)

        xa = predict(self, x[i], self.iy, ip)

        # Solve
        for _ in range(maxiter):
            p, (k, c, m) = evaluate(model, xa)
            ga = p - fa
            Dg = self.ga[2]*m + self.ga[1]*c + self.ga[0]*k
            dy = -ga / Dg

            # Equation (8)
            for r in range(3):
                xa[r] += self.ga[r]*dy


        x[i+1][:] = advance(self, x[i], xa)

    return x.T


def initialize(model, u0, v0):

    p, (k, c, m) = evaluate(model, [u0, v0, 0])
    
    return u0, v0, -p/m
    
def increment(self, i, f):
    fa = (1-self.lam)*f[i] + self.lam*f[i+1]
    return fa


def evaluate(model, xa):
    k, c, m = model
    p = m*xa[2] + c*xa[1] + k*xa[0]
    return p, model


def predict(self, xo, ip=None, Gp=None, gp=None):
    if ip is None:
        ip = self.iy

    if Gp is None:
        Gp = self.Ga
        gp = self.ga

    return [
            sum((self.Ga[r][s] + self.ga[r]/gp[ip]*(int(ip == s) - Gp[ip][s]))*xo[s] for s in range(3))
            for r in range(3)
    ]


def correct(self, xo, xa, state):
    pass


def advance(self, xo, xa):
    iy = self.iy
    x = [0.0, 0.0, 0.0]
    for s in range(3):
        for r in range(3):
            x[s] += (self.Gn[s][r] - self.gn[r]/self.ga[iy]*self.Ga[iy][r])*xo[r]

        x[s] += xa[iy]*self.gn[s]/self.ga[iy]
    return x

