#===----------------------------------------------------------------------===#
#
#         STAIRLab -- STructural Artificial Intelligence Laboratory
#
#===----------------------------------------------------------------------===#
#
# Claudio Perez
#
"""
[1] Simo JC, Wong KK. Unconditionally stable algorithms for the orthogonal group that exactly preserve energy and momentum.
    International Journal for Numerical Methods in Engineering 1991;
[2] Iserles A, Munthe-Kaas HZ, Nørsett SP, Zanna A. Lie-group methods. Acta Numerica 2000; 9:215-365.
[3] Krysl P. Dynamically equivalent implicit algorithms for the integration of rigid body rotations
[4] Hairer E, Nørsett SP, Wanner G. Solving Ordinary Differential Equations I. Nonstiff Problems (revised 2nd edn).
    Springer Series in Computational Mathematics, vol 8. Springer: Berlin, 1993.
[5] Munthe-Kaas H. Runge-Kutta methods on Lie groups. British Information Technology 1998; 38(1):92-111.

[6] Krysl P, Endres L Explicit Newmark/Verlet algorithm for time integration of the rotational dynamics of rigid bodies
    INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN ENGINEERING 62 (15): 2154-2177 APR 21 2005
"""
import numpy as np
from tqdm import tqdm
from collections import namedtuple
import warnings

from .util import ExpSO3
from .simo import *
from .rkmk import *
from .cgco import *
from .gcim import *
from .other import *


State = namedtuple("State", ["Alpha", "T"])
Source = namedtuple("Source", ["Torque"])


def kb_e              (t0, p0, w0, R0, dt, J, Loading, state):
    """
    """
    DR = ExpSO3(0.5 * dt * w0)
    Ra = R0 @ DR
    Ma = Loading.Torque(t0 + 0.5 * dt, Ra)
    wa = np.linalg.solve(J, DR.T @ J @ w0 + 0.5 * dt * Ma)

    w1 = w0 + dt * np.linalg.solve(J, Ma - np.cross(wa, J @ wa))

    DR2 = ExpSO3(0.5 * dt * w1)
    R1 = Ra @ DR2
    w1 = np.linalg.solve(J, DR2.T @ (DR.T @ J @ w0 + dt * Ma))
    P1 = J @ w1

    return P1, w1, R1, None


def integrate(J, dt, hc, T, wi, R0, Gm, step, Loading, potential):
    """
    Perform a simulation.

    Parameters:
    J (np.array): inertia tensor
    dt (float): time step size for the simulation
    hc (float): coarse time-stepsize for saving data
    T (float): total length of the simulation
    wi (np.array): initial body angular velocity (rotating coordinates)
    R0 (np.array): initial configuration
    alpha (float): strength of magnetic field
    Gm (np.array): Gm matrix
    step (function): step function to be used in the simulation
    potential (function): potential function
    Loading (np.array): Loading

    Returns:
    tuple: (t, g, p, U, E)
    """
    # Setup
    t = np.arange(0, T, dt)
    g = np.zeros((len(t), 3, 3))
    p = np.zeros((len(t), 3))
    U = np.zeros(len(t))
    E = np.zeros(len(t))

    # Initialize algorithm
    w0 = wi;         # body angular velocity at time 0
    p0 = J@w0;       # body angular momentum at time 0

    # Initial conditions
    g[0] = R0
    p[0] = potential(R0,Gm)
    U[0] = 0.5 * np.dot(w0, J@w0)
    E[0] = 0.5*w0.dot(p[0]) + U[0]


    state = State(Alpha=[], T=[])
    j = 0
    # Simulation loop
    for i in tqdm(range(1, len(t))):
        pn, wn, Rn, state = step((i-1)*dt, p0, w0, R0, dt, J, Loading, state)

        if True : # i % int(hc//dt) == 0:
            j += 1
            t[j] = t[j-1] + int(hc//dt)*dt
            U[j] = potential(Rn, Gm)
            E[j] = U[j] + 0.5 * wn.dot(pn)
            p[j][:] = pn
            g[j][:,:] = Rn

        p0 = pn
        w0 = wn
        R0 = Rn
    return t, g, p, U, E
