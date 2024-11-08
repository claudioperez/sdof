import numpy as np

def Wilson_theta(model,dt,F,u0=0,v0=0,beta=0.25,gamma=0.5,theta=1.4):
    """
    beta,gamma,theta: parameters.
    u0,v0,a0: initial state.
    T: time list with uniform interval.
    F: list of time-dependent force vectors.
    """
    dt_ = theta*dt
    b1 = 1/(beta*dt_**2)
    b2 = -1/(beta*dt_)
    b3 = (1/2-beta)/beta
    b4 = gamma*dt_*b1
    b5 = 1+gamma*dt_*b2
    b6 = dt_*(1+gamma*b3-gamma)
    K, C, M = model
    K_= K + b4*C + b1*M
    u = [u0]
    v = [v0]
    a = [(-C*v0 - K*u0)/M]
    R_=[F[0]]
    for i in range(len(F)-1):
        R_.append(F[i]+theta*(F[i+1]-F[i]))

    for ft in R_:
        ft_ = ft + M*(b1*u[-1]-b2*v[-1]-b3*a[-1]) \
                 + C*(b4*u[-1]-b5*v[-1]-b6*a[-1])

        ut_ = ft_ / K_
        vt_ = b4*(ut_-u[-1])+b5*v[-1]+b6*a[-1]
        at_ = b1*(ut_-u[-1])+b2*v[-1]+b3*a[-1]

        at=a[-1]+1/theta*(at_-a[-1])
        vt=v[-1]+((1-gamma)*a[-1]+gamma*at)*dt
        ut=u[-1]+v[-1]*dt + (1/2-beta)*a[-1]*dt**2+beta*at*dt**2

        u.append(ut)
        v.append(vt)
        a.append(at)

    return u, v, a
