from collections import namedtuple
from numpy import zeros
import numpy as np
from scipy.optimize import root

def alpha_a(dt,n,model,Force, u0, v0, a0):
    Mass, Damp, Stiff, pi, pe = model

    # Generalized-alpha method
    # M u_tt + C u_t +K u = F

    alpha_m = 1/2
    alpha_f = 1/2
    gamma   = 1/2 - alpha_m + alpha_f
    beta    = 1/4*(1 - alpha_m + alpha_f)**2

    b1 = dt
    b2 = dt**2/2*(1-2*beta)
    b3 = dt*(1-gamma)

    c1 = beta*dt**2
    c2 = gamma*dt
    c3 = 1.0

    def residual(an,uc,vc,ac,dt,tc,Force):
        # M a_(n+1-alpha_m) + C v_(n+1-alpha_f) + K d_(n+1-alpha_f) = F_(n+1-alpha_f)
        # x -- an, the accelerations at the next step

        un = uc + b1*vc +  b2*ac  + c1*an
        vn =         vc +  b3*ac  + c2*an

        tn = tc + dt

        uf = (1-alpha_f)*un + alpha_f*uc
        vf = (1-alpha_f)*vn + alpha_f*vc
        am = (1-alpha_m)*an + alpha_m*ac
        tf = (1-alpha_f)*tn + alpha_f*tc

        M  = Mass(uf,vf,am)
        C  = Damp(uf,vf,am)
        K  = Stiff(uf,vf,am)
        Ff = Force(tf)

        return (M@am.T + C@vf.T + K@uf.T - (Ff[:]).T).flatten()

    ndf = len(u0)
    dM  = np.zeros((n,ndf))
    vM  = np.zeros((n,ndf))
    aM  = np.zeros((n,ndf))
    dM[0,:] = u0
    vM[0,:] = v0
    aM[0,:] = a0

    for i in range(1,n):
        uc = dM[i-1,:]
        vc = vM[i-1,:]
        ac = aM[i-1,:]
        tc = (i-1)*dt
        fun = lambda x: residual(x,uc,vc,ac,dt,tc,Stiff,Force)

        sol = root(fun,ac)
        assert sol.success

        an = sol.x
        un = uc + dt*vc + b2*ac + c1*an
        vn =         vc + b3*ac + c2*an

        dM[i,:] = un
        vM[i,:] = vn
        aM[i,:] = an

    return dM,vM,aM


def newmark(M, C, K, F, dt, nt=None, n=0, u0=0.0, v0=0.0):
# 0 - the constant-average accelaration method (stable)
# 1 - the linear accelaration method (conditionally stable)
# 2 - the central difference method (conditionally stable) --- This need more explanation
# 3 - the Galerkin method (stable)
# 4 - the backward difference method (stable)
# b = 2 beta
    if n == 0:
        a = 0.5
        b = 0.5     # 1/4
    elif n == 1:
        a = 0.5
        b = 1.0/3.0 # 1/6
    elif n == 2:
        a = 0.5
        b = 0.0
    elif n == 3:
        a = 3.0/2.0
        b = 8.0/5.0
    elif n == 4:
        a = 3.0/2.0
        b = 2.0

    beta  = b/2
    gamma = a

    c1 =  1.0
    c2 =  1.0*gamma/(beta*dt)
    c3 =  1.0/(beta*dt**2)

    b1 =  1.0 - gamma/beta
    b2 = dt*(1.0 - gamma/(2*beta))
    b3 = -1.0/(beta*dt)
    b4 =  1.0 + 1.0/(beta*2)

    if nt is None:
        nt = len(F)

    u,v,a = np.zeros((3,nt))
    u[0] = u0
    v[0] = v0
    a[0] = (F[0] - C*v0 - K*u0)/M
    Ke = c3*M + c2*C + K

    for i in range(1,nt):
        if i == len(F):
            dF = -F[-1]
        elif i > len(F):
            dF = 0.0
        else:
            dF = F[i] - F[i-1]

        Fe = dF + ((2.0/(b * dt))*M + (2.0*gamma/b)*C)*v[i-1] + (1.0/(2*beta)*M  + b2*C)*a[i-1]
        du = Fe/Ke

        u[i] = u[i-1] + du
        v[i] = b1*v[i-1] + b2*a[i-1] + c2*du
        a[i] = b4*a[i-1] + b3*v[i-1] + c3*du

    return u, v, a




if __name__ == "__main__":
    # example
    # u1_tt + 4 u1 = 0
    # u2_tt + 9 u2 = 0
    # initial condition : u1(0) = u2(0) = 0, u1_t(0) = 2, u2_t(0) = 3
    # the exact solution is u1 = sin(2t), u2 = sin(3t)
    d0 = [0,0]
    v0 = [2,3]
    a0 = [0,0]
    dt = 1e-3
    n  = 9001
    mass = lambda u,v,a: [[1,0],
                          [0,1]]
    damp = lambda u,v,a: np.zeros((2,2))
    stif = lambda u,v,a: [[4,0],
                          [0,9]]

    Model = namedtuple("Model", "mass damp stif pi pe")

    model = Model(mass, damp, stif, None, None)
    force = lambda t:    np.zeros((2,1))
    dM,vM,aM = alpha_a(dt,n,model,force,d0,v0,a0);

    # check
    tspan = np.arange(n)*dt
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(tspan,dM[:,0],'b',linewidth=2.5, label="alpha")
    ax.plot(tspan,np.sin(2*tspan),'r', linewidth=1.5, label="exact")
    ax.set_xlabel('t')
    ax.set_ylabel('u_1')
    fig.legend()

    #
    fig, ax = plt.subplots()
    ax.plot(tspan,dM[:,1],'b', linewidth=2.5, label="alpha")
    ax.plot(tspan,np.sin(3*tspan),'r', linewidth=1.5, label="exact")
    fig.legend()
    ax.set_xlabel('t')
    ax.set_ylabel('u_2')
    plt.show()


