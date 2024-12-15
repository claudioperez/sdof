import numpy as np 
import warnings
from .util import dExpSO3, dLogSO3, ExpSO3, HatSO3, EPS


def bb_rkmk_trap_wdexp(t0, p0, w0, R0, dt, J, Loading, state, maxiter=20):
    """
    Momentum and energy conserving integrator from Bottasso,
    Borri (1998), Integrating finite rotations. Trapezoidal Runge-Kutta
    (Munthe-Kaas) algorithm. With the dexp included.
    """
    tol = 100 * EPS

    # calculate effective loads at t_n-1
    m0 = R0 @ Loading.Torque(t0, R0)

    w1 = w0
    Nun = dt * w1
    iter = 0
    pNu = np.zeros_like(w0)

    while True:
        Nun = dt / 2 * (w0 + dLogSO3(-Nun) @ w1)
        R1 = R0 @ ExpSO3(Nun)
        m1 = R1 @ Loading.Torque(t0 + dt, R1)
        w1 = np.linalg.solve(J, R1.T @ (R0 @ J @ w0 + dt / 2 * (m1 + m0)))
        iter += 1
        if iter > maxiter:
            if np.max(np.abs(Nun - pNu)) > 1e-9:
                warnings.warn(f"Failed to converge: ||residual||={np.linalg.norm(Nun - pNu)}")
            break
        if np.linalg.norm(Nun - pNu) < tol:
            break
        pNu = Nun

    P1 = None  # TODO! Placeholder, as P1 is not defined in the provided code
    return J@w1, w1, R1, None

def bb_rkmk_trap      (t0, p0, w0, R0, dt, J, Loading, state, maxiter = 150):
    """
    Momentum and energy conserving integrator from Bottasso,
    Borri (1998), Integrating finite rotations. 
    Trapezoidal Runge-Kutta (Munthe-Kaas) algorithm.
    """
    tol = 100 * EPS

    # calculate effective loads at t_n-1
    tn1 = R0 @ Loading.Torque(t0, R0)

    w1 = w0
    iter = 0
    pTheta = np.zeros_like(w0)

    while True:
        Theta = dt*0.5 * (w1 + w0)
        DRn = ExpSO3(Theta)
        Rn = R0 @ DRn
        tn = Rn @ Loading.Torque(t0 + dt, Rn)
        w1 = np.linalg.solve(J, Rn.T@(R0@J@w0 + dt / 2 * (tn + tn1)))
        iter += 1
        if iter > maxiter:
            if np.max(np.abs(Theta - pTheta)) > 1e-9:
                print(f"Warning: bb_rkmk_trap - Failed to converge: ||residual||={np.linalg.norm(Theta - pTheta)}")
            break
        if np.linalg.norm(Theta - pTheta) < tol:
            break
        pTheta = Theta

    P1 = None  # Placeholder, as P1 is not defined in the provided code
    return J@w1, w1, Rn, None

def dyneq_trap        (t0, p0, w0, R0, dt, J, Loading, state, maxiter = 10):
    """
    Trapezoidal rule integrator.
    Described as TRAP in [3].
    """
    Ptol = 100 * EPS

    t  = t0
    M0 = Loading.Torque(t, R0)
    H0 = J @ w0
    Hn = H0 - dt * np.cross(w0, H0) + dt * M0  # predictor
    iter = 0

    while True:
        pPin = Hn
        w1 = np.linalg.solve(J, Hn)
        Rn = R0 @ ExpSO3(dt*0.5 * w0) @ ExpSO3(dt*0.5 * w1)
        M1 = Loading.Torque(t + dt, Rn)
        Hn = H0 \
                - dt*0.5 * np.cross(w0, H0) \
                - dt*0.5 * np.cross(w1, Hn) \
                + dt / 2 * (M0 + M1)
        iter += 1
        if iter > maxiter:
            if np.max(np.abs(Hn - pPin)) > 1e-9:
                warnings.warn(f"Failed to converge: ||residual||={np.linalg.norm(Hn)}")
            break
        if np.linalg.norm(Hn - pPin) < Ptol:
            break

    return Hn, w1, Rn, None


def incrso3_imid      (t0, p0, w0, R0, dt, J, Loading, state):
    """
    Implicit mid-point integrator described as IMID in [3].
    This integrator is equivalent to the trapezoidal rule TRAP.
    """

    wa = _solve_imid(J, t0, dt, R0, J @ w0, Loading)
    Rn = R0 @ ExpSO3(dt * wa)
    Ra = R0 @ ExpSO3(dt*0.5 * wa)
    Ma = Loading.Torque(t0 + 0.5*dt, Ra)
    wn = w0 + np.linalg.solve(J, -dt * np.cross(wa, J @ wa) + dt * Ma)
    return J@wn, wn, Rn, None

def _solve_imid(J, t, dt, R0, H0, Loading):
    tol = 100 * EPS
    maxiter = 50

    wa = np.zeros(3)

    iter = 0
    while True:
        pOmegamid = wa
        Ra = R0 @ ExpSO3(dt / 2 * wa)
        Ma = Loading.Torque(t + dt / 2, Ra)
        wa = np.linalg.solve(J, H0 - dt / 2 * np.cross(wa, J @ wa) + dt / 2 * Ma)
        iter += 1
        if iter > maxiter:
            if np.max(np.abs(wa - pOmegamid)) > 1e-12:
                warnings.warn(f"Failed to converge: ||residual||={np.linalg.norm(wa)}")
            break
        if np.linalg.norm(wa - pOmegamid) < tol:
            break

    return wa


def incrso3_trap      (t0, p0, w0, R0, dt, J, Loading, state, maxiter = 10):
    """
    Trapezoidal rule integrator.
    Described as TRAP in [3].
    """
    Ptol = 100 * EPS

    M0 = Loading.Torque(t0, R0)
    H0 = J @ w0
    Hn = H0 \
        - dt * np.cross(w0, H0) \
        + dt * M0  # predictor

    iter = 0
    while True:
        pPin = Hn
        w1 = np.linalg.solve(J, Hn)
        Rn = R0 @ ExpSO3(dt*0.5 * w0) @ ExpSO3(dt*0.5 * w1)
        M1 = Loading.Torque(t0 + dt, Rn)
        Hn = H0 \
            - dt*0.5 * np.cross(w0, H0) \
            - dt*0.5 * np.cross(w1, Hn) \
            + dt*0.5 * (M0 + M1)

        iter += 1
        if iter > maxiter:
            if np.max(np.abs(Hn - pPin)) > 1e-9:
                warnings.warn(f"Failed to converge: ||residual||={np.linalg.norm(Hn)}")
            break

        if np.linalg.norm(Hn - pPin) < Ptol:
            break

    return Hn, w1, Rn, None
