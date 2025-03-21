import numpy as np
import warnings
from .util import dExpSO3, dLogSO3, ExpSO3, HatSO3, CaySO3, EPS 



def incrso3_vlv(t0, p0, w0, R0, dt, J, loading, state):
    """
    Lie-Verlet algorithm proposed in:

    Bou-Rabee, N. and Marsden, J. E. [2008].
    Hamilton-Pontryagin Integrators on Lie Groups.
    FOUNDATIONS OF COMPUTATIONAL MATHEMATICS.
    9: 197-219.

    Last Modified: May 27, 2010
    Authors:
        Nawaf Bou-Rabee (nawaf@cims.nyu.edu)
        Giulia Ortolan (ortolang@dei.unipd.it)
        Alessandro Saccon (asaccon@isr.ist.utl.pt)
    """
    maxiter = 10

    # torque applied at time n
    tau = loading.Torque(t0, R0)

    # Step 1) Newton Solve
    p1 = p0.copy()  # initialization
    w1 = w0.copy()  # initialization

    II = J

    for j in range(maxiter):
        F = w1.dot(p1)
        PxW = 0.5 * dt * np.cross(p1, w1)

        f = p0 - p1 + PxW - 0.25 * dt**2 * F * w1 + 0.5 * dt * tau

        if np.linalg.norm(f) < 1e-10:
            break
        elif j == maxiter - 1:
            print(f"{__file__} - maximum # of iterations exceeded")

        Jac = (-J + 0.5 * dt * np.array([
            [0, -p1[2] + II[1,1]*w1[2], p1[1] - II[2,2]*w1[1]],
            [p1[2] - II[0,0]*w1[2], 0, -p1[0] + II[2,2]*w1[0]],
            [-p1[1] + II[0,0]*w1[1], p1[0] - II[1,1]*w1[0], 0]
        ]) - 0.25 * dt**2 * F * np.eye(3))

        # Solve
        w1 = w1 - np.linalg.solve(Jac, f)
        p1 = J @ w1

    # Step 2) Explicit Configuration Update
    R1 = R0 @ CaySO3(dt * w1)

    # Step 3) Explicit Momentum Update
    tau = loading.Torque(t0 + dt, R1)  # torque applied at time n+1

    p1 = p1 + PxW + 0.25 * dt**2 * np.dot(w1, p1) * w1 + 0.5 * dt * tau  # body angular momentum
    w1 = np.linalg.solve(J, p1)  # body angular velocity

    return p1, w1, R1, state


##
def incrso3_vlv_      (t0, p0, w0, R0, dt, J, Loading, state, maxiter=10):
    """
    Incremental SO(3) Variational Lie Group Integrator (VLV)
    """

    tau = Loading.Torque(t0, R0)  # torque applied at time 0
    w1 = _solve_vlv(J, w0, tau, dt, maxiter, tol=1e-10)
    R1 = R0@ExpSO3(dt * w1)  # update rotation matrix

    return J@w1, w1, R1, None

def _solve_vlv(J, w1, tau, dt, maxiter, tol):
    p0 = J@w1
    p1 = p0
    II = np.diag(J)  # inertia vector
    for j in range(maxiter):
        pxw = 0.5 * dt * np.cross(p1, w1)

        f = p0 - p1 + pxw + 0.5 * dt * tau

        if np.linalg.norm(f) < tol:
            break
        elif j == maxiter - 1:
            warnings.warn("Failed to converge")

        Jac = -J + 0.5 * dt * np.array([
            [0, -p1[2] + II[1] * w1[2], p1[1] - II[2] * w1[1]],
            [p1[2] - II[0] * w1[2], 0, -p1[0] + II[2] * w1[0]],
            [-p1[1] + II[0] * w1[1], p1[0] - II[1] * w1[0], 0]
        ])

        w1 = w1 - np.linalg.solve(Jac, f)  # body angular velocity
        p1 = J@w1  # body angular momentum

    return w1


def incrso3_svq       (t0, p0, W0, R0, dt, J, Loading, state, maxiter=10):
    """

    """

    tol = 1e-10
    w1 = W0  # initialization
    p1 = p0
    tau = Loading.Torque(t0, R0)  # torque applied at time 0
    II = np.diag(J)  # inertia vector

    for j in range(maxiter):
        pxw = 0.5 * dt * np.cross(p1, w1)

        f = p0 - p1 + pxw + 0.5 * dt * tau

        if np.linalg.norm(f) < tol:
            break
        elif j == maxiter - 1:
            warnings.warn("Failed to converge")

        df = -J + 0.5 * dt * np.array([
            [0, -p1[2] + II[1] * w1[2], p1[1] - II[2] * w1[1]],
            [p1[2] - II[0] * w1[2], 0, -p1[0] + II[2] * w1[0]],
            [-p1[1] + II[0] * w1[1], p1[0] - II[1] * w1[0], 0]
        ])

        w1 = w1 - np.linalg.solve(df, f)  # body angular velocity
        p1 = J@w1  # body angular momentum

    R1 = R0@ExpSO3(dt * w1)  # update rotation

    return p1, w1, R1, None


def simo_wong_algo_c1(t0, p0, w0, R0, dt, J, state, Loading):
    """
    Momentum and energy conserving integrator from Simo, Wong.
    The momentum and energy conserving integrator from Simo, Wong.
    """
    gamma = 1
    beta = 1/2  # the Newmark parameters: these are values determined as optimal by SW

    nsteps = 1
    maxiter = 30
    tol = 1000 * EPS
    alph = min(beta / gamma, 1)

    t = t0

    # calculate effective loads at t_n-1
    tn1 = Loading.Torque(t, R0, Loading)
    if not hasattr(state, 'Alpha') or state.Alpha is None:
        # Alphan1 = I\(Rn1'*tn1 - Spin(Omegan1)*I*Omegan1)
        Alphan1 = np.linalg.solve(J, tn1 - np.cross(w0, J @ w0))
        state.Alpha = Alphan1

    Alphan1 = state.Alpha
    for n in range(nsteps):
        w1 = w0.copy()
        iter_count = 0
        pTheta = np.zeros_like(w0)
        while True:
            Theta = (beta / gamma) * dt * (
                w1 + w0 - (2 - gamma / beta) * (w0 + dt / 2 * Alphan1)
            )
            DRn = ExpSO3(Theta)
            Rn = R0 @ DRn
            Ra = R0 @ ExpSO3(alph * Theta)
            ta = Loading.Torque(t + alph * dt, Ra, Loading)
            
            # Omegan = I\(DRn'*I*Omegan1 + dt*Rn'*ta);
            w1 = np.linalg.solve(J, DRn.T @ (J @ w0) + dt * ta)
            
            iter_count += 1
            if iter_count > maxiter:
                if np.max(np.abs(Theta - pTheta)) > 1e-12:
                    print(f"c1 - Failed to converge at t={t}, ||residual||={np.linalg.norm(Theta - pTheta)}")
                break
            elif np.linalg.norm(Theta - pTheta) < tol:
                break
            pTheta = Theta

        # swap variables for next step
        state.Alpha = (w1 - w0) / (dt * gamma) + (1 - 1 / gamma) * Alphan1

    # Return statement placeholder (assuming Rn is defined outside the loop)
    return J@w0, w1, Rn, state

def incrso3_svq_alt   (t0, p0, W0, R0, dt, J, Loading, state, maxiter=10):
    """

    """

    tol = 1e-10
    w1 = W0  # initialization
    p1 = p0
    tau = Loading.Torque(t0, R0)  # torque applied at time 0
    II = np.diag(J)  # inertia vector

    for j in range(maxiter):
        pxw = 0.5 * dt * np.cross(p1, w1)

        f = p0 - p1 + pxw + 0.5 * dt * tau

        if np.linalg.norm(f) < tol:
            break
        elif j == maxiter - 1:
            warnings.warn("Failed to converge")

        df = -J + 0.5 * dt * np.array([
            [0, -p1[2] + II[1] * w1[2], p1[1] - II[2] * w1[1]],
            [p1[2] - II[0] * w1[2], 0, -p1[0] + II[2] * w1[0]],
            [-p1[1] + II[0] * w1[1], p1[0] - II[1] * w1[0], 0]
        ])

        w1 = w1 - np.linalg.solve(df, f)  # body angular velocity
        p1 = J@w1  # body angular momentum

    R1 = R0@ExpSO3(dt * w1)  # update rotation

    return p1, w1, R1, None



def incrso3_swc2      (t0, p0, w0, R0, dt, J, Loading, state):
    """
    Implements the ALGO_C2 algorithm?
    """
    t = t0

    M0 = Loading.Torque(t, R0)
    if state.Alpha is None:
        A0 = np.linalg.solve(J, M0 - np.cross(w0, J @ w0))
        state.Alpha = A0

    Theta = dt * (w0 + dt / 2 * state.Alpha)
    DR = ExpSO3(Theta)
    Rn = R0 @ DR
    M1 = Loading.Torque(t + dt, Rn)
    wn = np.linalg.solve(J, Rn.T @ (R0 @ J @ w0 + dt / 2 * (M0 + M1)))  # trap approx of impulse

    state.Alpha = (wn - w0) / dt

    P1 = J @ wn
    return P1, wn, Rn, state

##
def liemid_newmark    (t0, p0, w0, R0, dt, J, Loading, state, maxiter=10):
    tol = 100 * EPS

    def approx_newton(J, I_inv, vec, tol):
        Wa = np.zeros_like(vec)
        for _ in range(maxiter):
            F = J @ Wa - vec
            if np.linalg.norm(F) < tol:
                break
            dWa = -I_inv @ F
            Wa += dWa
        return Wa

    Omegan1 = w0
    Rn1 = R0

    t = t0
    if not state.T:
        state.T = Loading.Torque(t, Rn1)

    Tn1 = state.T
    I_inv = np.linalg.inv(J)

    # half-step with algorithm 1
    Wa = approx_newton(J, I_inv, J @ Omegan1 + dt*0.5 * Tn1, tol)
    Rn = Rn1 @ ExpSO3(dt*0.5 * Wa)
    Omegan = np.linalg.solve(J, Rn.T@Rn1@(J @ Omegan1 + dt*0.5 * Tn1))

    # half-step with algorithm 2
    Rn1 = Rn
    Omegan1 = Omegan
    Wa = approx_newton(J, I_inv, J @ Omegan1, tol)
    Rn = Rn1 @ ExpSO3(dt*0.5 * Wa)

    return J@Omegan, Omegan, Rn, state

def rotint_nmb        (t0, p0, w0, R0, dt, J, Loading, state):
    """
    See NMB method in [6]

        Body-frame version of the gamma=1/2, beta=0,
        with Newton iteration to solve the eqn of motion.
    """
    tol =  100 * EPS

    a0 = np.linalg.solve(J, Loading.Torque(t0, R0) - np.cross(w0, J@w0))

    # update rotation
    Rn = R0@ExpSO3(dt * w0 + (dt ** 2)*0.5 * a0)

    # calculate effective loads at t_n
    M0 = Loading.Torque(t0 + dt, Rn)
    an = _solve_nmb(J, M0, dt, a0, w0, tol)

    # update velocity
    wn = w0 + dt*0.5 * (a0 + an)
    P1 = J@wn

    return P1, wn, Rn, None

def _solve_nmb(J, M0, dt, a0, w0, tol, maxiter=100):
    """
    Solve for Alpha using Newton's method.

    Parameters:
    J (np.array): inertia matrix
    M0 (np.array): torque at time t+dt
    dt (float): time step size
    Alphan1 (np.array): Alpha at time t
    w0 (np.array): angular velocity at time t
    tol (float): tolerance for convergence

    Returns:
    np.array: Alpha at time t+dt
    """
    an = a0
    for _ in range(maxiter):
        wn = w0 + dt / 2 * an
        f = J@an + np.cross(wn, J@wn)  - M0
        if np.linalg.norm(f) < tol:
            break
        df = J + dt / 2 * HatSO3(J@wn)
        an = an -  np.linalg.solve(df, f)
    return an
