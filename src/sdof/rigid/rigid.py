#===----------------------------------------------------------------------===#
#
#         STAIRLab -- STructural Artificial Intelligence Laboratory
#
#===----------------------------------------------------------------------===#
#
# Claudio Perez
#
"""
[3] Dynamically equivalent implicit integrators
"""
import numpy as np
from tqdm import tqdm
from collections import namedtuple
import warnings 

EPS = np.finfo(float).eps
State = namedtuple('State', ['Alpha', 'T'])
Source = namedtuple('Source', ['Torque'])

def CaySO3(w):
    theta = np.linalg.norm(w)
    if theta < 1e-10:
        return np.eye(3)
    else:
        W = HatSO3(w)
        return np.eye(3) + (np.sin(theta) / theta) * W + ((1 - np.cos(theta)) / (theta ** 2)) * (W @ W)

def ExpSO3(omega):
    """
    Exponential map from so(3) to SO(3).
    
    Parameters:
    omega (np.array): 3x1 vector in the Lie algebra so(3)
    
    Returns:
    np.array: 3x3 rotation matrix in SO(3)
    """
    theta = np.linalg.norm(omega)
    if theta < 1e-10:
        return np.eye(3)
    omega_hat = HatSO3(omega)
    return np.eye(3) + (np.sin(theta) / theta) * omega_hat + ((1 - np.cos(theta)) / (theta ** 2)) * np.dot(omega_hat, omega_hat)

def dExpSO3(omega):
    """
    Compute the differential of the exponential map on SO(3).
    
    Parameters:
    omega (np.array): 3x1 vector in the Lie algebra so(3)
    
    Returns:
    np.array: 3x3 matrix representing the differential of the exponential map
    """
    theta = np.linalg.norm(omega)
    if theta < 1e-10:
        return np.eye(3)
    
    omega_hat = HatSO3(omega)
    omega_hat_sq = omega_hat @ omega_hat
    
    a1 = (1 - np.cos(theta)) / (theta ** 2)
    a2 = (theta - np.sin(theta)) / (theta ** 3)
    
    return np.eye(3) + a1 * omega_hat + a2 * omega_hat_sq

def LogSO3(R):
    """
    Logarithm map from SO(3) to so(3) using Spurrier's algorithm.
    
    Parameters:
    R (np.array): 3x3 rotation matrix in SO(3)
    
    Returns:
    np.array: 3x1 vector in the Lie algebra so(3)
    """
    theta = np.arccos((np.trace(R) - 1) / 2)
    if theta < 1e-10:
        return np.zeros(3)
    omega_hat = (theta / (2 * np.sin(theta))) * (R - R.T)
    return VeeSO3(omega_hat)

def dLogSO3(R):
    """
    Differential of the logarithm map (inverse of the exponential tangent) on SO(3).
    
    Parameters:
    R (np.array): 3x3 rotation matrix in SO(3)
    
    Returns:
    np.array: 3x3 matrix representing the differential of the logarithm map
    """
    omega = LogSO3(R)
    theta = np.linalg.norm(omega)
    if theta < 1e-10:
        return np.eye(3)
    omega_hat = HatSO3(omega)
    A = (1 - np.cos(theta)) / (theta ** 2)
    B = (theta - np.sin(theta)) / (theta ** 3)
    return np.eye(3) - 0.5 * omega_hat + (1 / theta ** 2) * (1 - A / (2 * B)) * np.dot(omega_hat, omega_hat)

def HatSO3(omega):
    """
    Hat operator for so(3).
    
    Parameters:
    omega (np.array): 3x1 vector
    
    Returns:
    np.array: 3x3 skew-symmetric matrix
    """
    return np.array([
        [0, -omega[2], omega[1]],
        [omega[2], 0, -omega[0]],
        [-omega[1], omega[0], 0]
    ])

def VeeSO3(omega_hat):
    """
    Vee operator for so(3).
    
    Parameters:
    omega_hat (np.array): 3x3 skew-symmetric matrix
    
    Returns:
    np.array: 3x1 vector
    """
    return np.array([omega_hat[2, 1], omega_hat[0, 2], omega_hat[1, 0]])

def bb_rkmk_trap_wdexp(t0, p0, w0, R0, dt, J, Loading, state):
    """
    Momentum and energy conserving integrator from Bottasso,
    Borri (1998), Integrating finite rotations. Trapezoidal Runge-Kutta
    (Munthe-Kaas) algorithm. With the dexp included.
    """
    maxiter = 20
    tol = 100 * np.finfo(float).eps

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

    P1 = None  # Placeholder, as P1 is not defined in the provided code
    return P1, w1, R1, None

def bb_rkmk_trap      (t0, p0, w0, R0, dt, J, Loading, state):
    """
    Momentum and energy conserving integrator from Bottasso,
    Borri (1998), Integrating finite rotations. Trapezoidal Runge-Kutta
    (Munthe-Kaas) algorithm.
    """
    tol = 100 * np.finfo(float).eps
    maxiter = 150

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
    return P1, w1, Rn, None

def dyneq_trap        (t0, p0, w0, R0, dt, J, Loading, state):
    """
    Trapezoidal rule integrator.
    Described as TRAP in the paper on dynamically equivalent implicit integrators.
    """
    maxiter = 10
    Ptol = 100 * np.finfo(float).eps

    t = t0
    M0 = Loading.Torque(t, R0)
    Pin1 = J @ w0
    Pin = Pin1 - dt * np.cross(w0, Pin1) + dt * M0  # predictor
    iter = 0

    while True:
        pPin = Pin
        w1 = np.linalg.solve(J, Pin)
        Rn = R0 @ ExpSO3(dt*0.5 * w0) @ ExpSO3(dt*0.5 * w1)
        M1 = Loading.Torque(t + dt, Rn)
        Pin = Pin1 - dt*0.5 * HatSO3(w0) @ Pin1 - dt*0.5 * HatSO3(w1) @ Pin + dt / 2 * (M0 + M1)
        iter += 1
        if iter > maxiter:
            if np.max(np.abs(Pin - pPin)) > 1e-9:
                warnings.warn(f"Warning: IncrSO3_Trap - Failed to converge: ||residual||={np.linalg.norm(Pin)}")
            break
        if np.linalg.norm(Pin - pPin) < Ptol:
            break

    return Pin, w1, Rn, None

def incrso3_akw       (t0, p0, w0, R0, dt, J, Loading, state):
    """
    Implicit Austin, Krishnaprasad, Wang mid-point integrator constructor.
    """
    # akw_mid 

    H0 = J @ w0
    H1 = _solve_akw(J, t0, dt, R0, H0, Loading)
    wa = np.linalg.solve(J, (H1 + H0) / 2)
    R1 = R0 @ CaySO3(dt * wa)
    w1 = np.linalg.solve(J, H1)
    return H1, w1, R1, None

def _solve_akw(J, t0, dt, R0, H0, Loading):
    maxiter = 25
    tol = 100 * np.finfo(float).eps

    H1 = H0
    M0 = Loading.Torque(t0, R0)
    iter = 0

    while True:
        Hx = H1
        wa = np.linalg.solve(J, (H1 + H0) / 2)
        Rn = R0 @ CaySO3(dt * wa)
        M1 = Loading.Torque(t0 + dt, Rn)
        H1 = H0 + dt * np.cross(J@wa,  wa) + dt / 2 * (M1 + M0)
        iter += 1
        if iter > maxiter:
            if np.max(np.abs(H1 - Hx)) > 1e-12:
                warnings.warn(f"Warning: {__file__} - Failed to converge: ||residual||={np.linalg.norm(H1)}")
            break
        if np.linalg.norm(H1 - Hx) < tol:
            break

    return H1

def incrso3_imid      (t0, p0, w0, R0, dt, J, Loading, state):
    """
    Implicit mid-point integrator described as IMID in the paper on the
    dynamic equivalence of implicit integrators. This integrator is
    equivalent to the trapezoidal rule TRAP.
    """
    tol = 100 * np.finfo(float).eps

    wa = _solve_imid(J, t0, dt, R0, J @ w0, Loading, tol)
    Rn = R0 @ ExpSO3(dt * wa)
    Ra = R0 @ ExpSO3(dt*0.5 * wa)
    Ma = Loading.Torque(t0 + 0.5*dt, Ra)
    w1 = w0 + np.linalg.solve(J, -dt * HatSO3(wa)@J@ wa + dt * Ma)
    P1 = J @ w1

    return P1, w1, Rn, None

def _solve_imid(J, t, dt, R0, H0, Loading, ceps):
    maxi = 50

    wa = np.zeros(3)
    iter = 0

    while True:
        pOmegamid = wa
        Ra = R0 @ ExpSO3(dt / 2 * wa)
        Ma = Loading.Torque(t + dt / 2, Ra)
        wa = np.linalg.solve(J, H0 - dt / 2 * HatSO3(wa) @ J @ wa + dt / 2 * Ma)
        iter += 1
        if iter > maxi:
            if np.max(np.abs(wa - pOmegamid)) > 1e-12:
                print(f"Warning: {__file__} - Failed to converge: ||residual||={np.linalg.norm(wa)}")
            break
        if np.linalg.norm(wa - pOmegamid) < ceps:
            break

    return wa

def incrso3_imidm     (t0, p0, w0, R0, dt, J, Loading, state):
    """
    Implicit mid-point Lie integrator constructor. Described as LIEMID[I] in
    the paper on the explicit midpoint Lie algorithm, and as IMIDM in the
    paper on dynamically equivalent implicit algorithms.
    """
    wa = _solve_imidm(J, t0, dt, R0, J @ w0, Loading)
    Rn = R0 @ ExpSO3(dt * wa)
    ma = R0 @ Loading.Torque(t0 + dt*0.5, R0@ExpSO3(dt*0.5 * wa))
    w1 = np.linalg.solve(J, Rn.T @ (R0 @ J @ w0 + dt * ma))

    return J@w1, w1, Rn, None

def _solve_imidm(J, t, dt, R0, H0, Loading):
    maxiter = 50
    tol = 1000 * np.finfo(float).eps


    wa = np.linalg.solve(J, H0)
    ma = R0 @ Loading.Torque(t + dt*0.5, R0@ExpSO3(dt*0.5 * wa))
    B = H0 + dt / 2 * R0.T @ ma
    res = wa - np.linalg.solve(J, ExpSO3(-dt / 2 * wa) @ B)
    invIdt2 = np.linalg.inv(J) * dt**2 / 8

    iter = 0
    while np.max(np.abs(res)) > tol:
        K = np.eye(3) - np.linalg.solve(J, dt / 2 * HatSO3(B)) \
            + invIdt2 @ (HatSO3(wa) @ HatSO3(B) + HatSO3(HatSO3(wa) @ B))
        wa = wa - np.linalg.solve(K, res)
        ma = R0 @ Loading.Torque(t + dt / 2, R0 @ ExpSO3(dt / 2 * wa))
        B = H0 + dt / 2 * R0.T @ ma
        res = wa - np.linalg.solve(J, ExpSO3(-dt / 2 * wa) @ B)
        iter += 1
        if iter > maxiter:
            print(f"Warning: {__file__} - Failed to converge after {maxiter} iterations")
            break

    return wa

def KB_E              (t0, p0, w0, R0, dt, J, Loading, state):
    """
    Implicit integrator described in the provided MATLAB code.
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

def incrso3_mleok     (t0, p0, w0, R0, dt, J, Loading, state):
    """
    Melvin Leok variational integrator.
    A Lie Group Variational Integrator for the Attitude Dynamics of a Rigid Body with Applications to
    the 3D Pendulum (with T.Y. Lee, and N.H. McClamroch), submitted, IEEE Conference on Control Applications, 2005.
    """
    tol = 100 * np.finfo(float).eps
    maxiter = 20

    f = dt * w0
    r = f + 100 * tol
    iter = 1

    while np.linalg.norm(r) > tol:
        nf = np.linalg.norm(f)
        r = dt * J @ w0 \
          - np.sin(nf) / nf * J @ f \
          - (1 - np.cos(nf)) / nf**2 * np.cross(f, J @ f)
        f = f + np.linalg.solve(DG(f, J), r)
        if iter > maxiter:
            if np.max(np.abs(r)) > 1e-9:
                print(f"Warning: {__file__} - Failed to converge: ||residual||={np.linalg.norm(r)}")
            break
        iter += 1

    DR = ExpSO3(f)
    Rn = R0 @ DR
    Tn = Loading.Torque(t0 + dt, Rn)
    Omegan = np.linalg.solve(J, DR.T @ J @ w0 + dt * Tn)

    P1 = J @ Omegan
    return P1, Omegan, Rn, None

def incrso3_liemidea  (t0, p0, w0, R0, dt, J, Loading, state, maxiter = 10):
    """
    Incremental SO(3) Lie Midpoint Method (LIEMIDEA)
    
    Parameters:
    t0 (float): initial time
    p0 (np.array): initial momentum
    W0 (np.array): initial angular velocity
    R0 (np.array): initial rotation matrix
    dt (float): time step size
    J (np.array): inertia matrix
    state (dict): state dictionary
    Loading (object): object with method Torque
    
    Returns:
    tuple: (p1, w1, R1, state)
    """

    TH1 = dt * w0  # initialization
    M0 = Loading.Torque(t0, R0)  # torque applied at time 0
    
    for j in range(maxiter):
        expx = ExpSO3(-TH1*0.5)

        v = 0.5 * dt * np.dot(expx, (p0 + 0.5*dt * M0))
        f = -J@TH1 + v

        if np.linalg.norm(f) < 1e-10:
            break
        elif j == maxiter - 1:
            print('maximum # of iterations exceeded')
        
        Jac = -J + 0.5 * HatSO3(v)@dExpSO3(-TH1 * 0.5)
        
        TH1 = TH1 - np.linalg.solve(Jac, f)
    
    # Step 2) Explicit Configuration Update
    R1 = np.dot(R0, ExpSO3(TH1))
    
    # Update momentum and angular velocity
    p1 = p0 + dt*M0
    w1 = np.dot(np.linalg.inv(J), p1)
    
    return p1, w1, R1, None

def incrso3_svq       (t0, p0, W0, R0, dt, J, Loading, state, maxiter=10):
    """
    Incremental SO(3) Symplectic Variational Integrator (SVQ)
    
    Parameters:
    t0 (float): initial time
    p0 (np.array): initial momentum
    W0 (np.array): initial angular velocity
    R0 (np.array): initial rotation matrix
    dt (float): time step size
    J (np.array): inertia matrix
    state (dict): state dictionary
    Loading (object): object with method Torque
    maxiter (int): maximum number of iterations
    
    Returns:
    tuple: (p1, w1, R1, state)
    """
    tol = 1e-10
    w1 = W0  # initialization
    p1 = p0
    tau = Loading.Torque(t0, R0)  # torque applied at time 0
    II = np.diag(J)  # inertia vector

    for j in range(maxiter):
        pxw = 0.5 * dt * np.array([
            p1[1] * w1[2] - p1[2] * w1[1],
            p1[2] * w1[0] - p1[0] * w1[2],
            p1[0] * w1[1] - p1[1] * w1[0]
        ])  # 0.5 * dt * (p1 x w1)

        f = p0 - p1 + pxw + 0.5 * dt * tau

        if np.linalg.norm(f) < tol:
            break
        elif j == maxiter - 1:
            print('maximum # of iterations exceeded')

        Jac = -J + 0.5 * dt * np.array([
            [0, -p1[2] + II[1] * w1[2], p1[1] - II[2] * w1[1]],
            [p1[2] - II[0] * w1[2], 0, -p1[0] + II[2] * w1[0]],
            [-p1[1] + II[0] * w1[1], p1[0] - II[1] * w1[0], 0]
        ])

        w1 = w1 - np.linalg.solve(Jac, f)  # body angular velocity
        p1 = II * w1  # body angular momentum

    R1 = R0@ExpSO3(dt * w1)  # update rotation

    return p1, w1, R1, None

def incrso3_swc2      (t0, p0, w0, R0, dt, J, Loading, state):
    """
    Implements the IncrSO3_SWC2 algorithm.
    """
    t = t0

    M0 = Loading.Torque(t, R0)
    if state.Alpha is None:
        A0 = np.linalg.solve(J, M0 - HatSO3(w0) @ J @ w0)
        state.Alpha = A0

    Theta = dt * (w0 + dt / 2 * state.Alpha)
    DR = ExpSO3(Theta)
    Rn = R0 @ DR
    M1 = Loading.Torque(t + dt, Rn)
    wn = np.linalg.solve(J, Rn.T @ (R0 @ J @ w0 + dt / 2 * (M0 + M1)))  # trap approx of impulse

    state.Alpha = (wn - w0) / dt

    P1 = J @ wn
    return P1, wn, Rn, state

def incrso3_trap      (t0, p0, w0, R0, dt, J, Loading, state, maxiter = 10):
    """
    Trapezoidal rule integrator.
    Described as TRAP in [3].
    """
    Ptol = 100 * np.finfo(float).eps

    M0 = Loading.Torque(t0, R0)
    Pin1 = J @ w0
    Pin = Pin1 \
        - dt * np.cross(w0, Pin1) \
        + dt * M0  # predictor
    
    iter = 0
    while True:
        pPin = Pin
        w1 = np.linalg.solve(J, Pin)
        Rn = R0 @ ExpSO3(dt*0.5 * w0) @ ExpSO3(dt*0.5 * w1)
        M1 = Loading.Torque(t0 + dt, Rn)
        Pin = Pin1 \
            - dt*0.5 * np.cross(w0, Pin1) \
            - dt*0.5 * np.cross(w1, Pin )\
            + dt*0.5 * (M0 + M1)

        iter += 1
        if iter > maxiter:
            if np.max(np.abs(Pin - pPin)) > 1e-9:
                print(f"Warning: {__file__} - Failed to converge: ||residual||={np.linalg.norm(Pin)}")
            break

        if np.linalg.norm(Pin - pPin) < Ptol:
            break

    return Pin, w1, Rn, None

def incrso3_trapm_zeta(t0, p0, w0, R0, dt, J, Loading, state):
    """
    Trapezoidal rule with momentum conservation integrator constructor.
    Described as TRAPM in the paper on dynamically equivalent implicit integrators.
    Dynamically equivalent to the implicit midpoint Lie integrator.
    """
    maxiter = 150
    zeta = 0.5
    tol = 100 * np.finfo(float).eps

    M0 = Loading.Torque(t0, R0)
    H0 = J @ w0
    DRn = ExpSO3(dt * w0)  # predictor
    Pin = DRn.T @ (H0 + dt*0.5 * M0) + dt*0.5 * M0  # predictor
    iter = 0
    while True:
        Hx = Pin
        Omegan = np.linalg.solve(J, Pin)
        DRn = ExpSO3(zeta * dt * w0) @ ExpSO3((1 - zeta) * dt * Omegan)
        Rn = R0 @ DRn
        M1 = Loading.Torque(t0 + dt, Rn)
        Pin = DRn.T @ H0 \
            + zeta * dt * DRn.T @ M0 \
            + (1 - zeta) * dt * M1
        iter += 1
        if iter > maxiter:
            if np.max(np.abs(Pin - Hx)) > 1e-9:
                print(f"Warning: {__file__} - Failed to converge: ||residual||={np.linalg.norm(Pin)}")
            break
        if np.linalg.norm(Pin - Hx) < tol:
            break

    return Pin, Omegan, Rn, None

def incrso3_trapm     (t0, p0, w0, R0, dt, J, Loading, state):
    """
    Trapezoidal rule with momentum conservation integrator constructor.
    Described as TRAPM in the paper on dynamically equivalent implicit integrators.
    Dynamically equivalent to the implicit midpoint Lie integrator.
    """
    maxiter = 10
    Ptol = 100 * np.finfo(float).eps

    M0 = Loading.Torque(t0, R0)
    H0 = J @ w0
    DRn = ExpSO3(dt / 2 * w0) @ ExpSO3(dt / 2 * w0)  # predictor
    H1 = DRn.T @ (H0 + dt / 2 * M0) + dt / 2 * M0  # predictor
    iter = 0

    while True:
        Hx = H1
        w1 = np.linalg.solve(J, H1)
        DRn = ExpSO3(dt / 2 * w0) @ ExpSO3(dt / 2 * w1)
        Rn = R0 @ DRn
        M1 = Loading.Torque(t0 + dt, Rn)
        H1 = DRn.T @ H0 + dt / 2 * DRn.T @ M0 + dt / 2 * M1
        iter += 1
        if iter > maxiter:
            if np.max(np.abs(H1 - Hx)) > 1e-9:
                print(f"Warning: {__file__} - Failed to converge: ||residual||={np.linalg.norm(H1)}")
            break
        if np.linalg.norm(H1 - Hx) < Ptol:
            break

    return H1, w1, Rn, None

##
def incrso3_vlv       (t0, p0, w0, R0, dt, J, Loading, state, maxiter=10):
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
        pxw = 0.5 * dt * np.array([
            p1[1] * w1[2] - p1[2] * w1[1],
            p1[2] * w1[0] - p1[0] * w1[2],
            p1[0] * w1[1] - p1[1] * w1[0]
        ])  # 0.5 * dt * (p1 x w1)

        f = p0 - p1 + pxw + 0.5 * dt * tau

        if np.linalg.norm(f) < tol:
            break
        elif j == maxiter - 1:
            warnings.warn('maximum # of iterations exceeded')

        Jac = -J + 0.5 * dt * np.array([
            [0, -p1[2] + II[1] * w1[2], p1[1] - II[2] * w1[1]],
            [p1[2] - II[0] * w1[2], 0, -p1[0] + II[2] * w1[0]],
            [-p1[1] + II[0] * w1[1], p1[0] - II[1] * w1[0], 0]
        ])

        w1 = w1 - np.linalg.solve(Jac, f)  # body angular velocity
        p1 = J@w1  # body angular momentum

    return w1

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
    Krysl P, Endres L
    Explicit Newmark/Verlet algorithm for time integration of the rotational dynamics of rigid bodies
    INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN ENGINEERING 62 (15): 2154-2177 APR 21 2005
        Body-frame version of the gamma=1/2, beta=0,
        with Newton iteration to solve the eqn of motion.
    
    
    Parameters:
    t0 (float): initial time
    p0 (np.array): initial momentum
    w0 (np.array): initial angular velocity
    R0 (np.array): initial rotation matrix
    dt (float): time step size
    J (np.array): inertia matrix
    state (dict): state dictionary
    Loading (object): object with method Torque
    
    Returns:
    tuple: (P1, Omegan, Rn, state)
    """
    tol =  100 * np.finfo(float).eps

    Alphan1 = np.linalg.solve(J, Loading.Torque(t0, R0) - np.cross(w0, J@w0))

    # update rotation
    Rn = R0@ExpSO3(dt * w0 + (dt ** 2)*0.5 * Alphan1)
  
    # calculate effective loads at t_n
    M0 = Loading.Torque(t0 + dt, Rn)
    Alphan = solve_Alpha(J, M0, dt, Alphan1, w0, tol)
    
    # update velocity
    Omegan = w0 + dt*0.5 * (Alphan1 + Alphan)
    P1 = np.dot(J, Omegan)
    
    return P1, Omegan, Rn, None

def solve_Alpha(J, M0, dt, Alphan1, w0, tol, maxiter=100):
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
    Alpha = Alphan1
    for _ in range(maxiter):  # maximum 100 iterations
        f = np.dot(J, Alpha) - M0 + np.dot(HatSO3(w0 + dt / 2 * Alpha), np.dot(J, (w0 + dt / 2 * Alpha)))
        if np.linalg.norm(f) < tol:
            break
        df = J + dt / 2 * HatSO3(np.dot(J, (w0 + dt / 2 * Alpha)))
        Alpha = Alpha - np.linalg.solve(df, f)
    return Alpha

def integrate(dt, hc, T, J, wi, R0, Gm, step, Loading, potential):
    """
    Perform the simulation.
    
    Parameters:
    dt (float): time step size for the simulation
    hc (float): coarse time-stepsize for saving data
    T (float): total length of the simulation
    II (np.array): inertia vector
    wi (np.array): initial body angular velocity (rotating coordinates)
    R0 (np.array): initial configuration
    alpha (float): strength of magnetic field
    Gm (np.array): Gm matrix
    step (function): step function to be used in the simulation
    potential (function): potential function
    Loading (np.array): Loading matrix
    
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
    w0 = wi;          # body angular velocity at time 0
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

        if i % int(hc//dt) == 0:
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
        V = (de - 1.0)**2 - alpha/dm; 
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
    T = 10000.  # total simulation time

    # Loading
    Wi0 = np.array([0.0, 0.0, 0.625])  # initial body angular velocity
    Gi0 = ExpSO3([0.0, 0.7227, 0.0])   # initial configuration (identity matrix)

    Loading = Source(Torque = potentialtorque)

    # Methods to be analyzed
    Methods = [
        # "dyneq_trap",
        "IncrSO3_trapm",
        # "IncrSO3_trapm_zeta",
        # "IncrSO3_imidm",
        # "IncrSO3_AKW",
        # "rotint_nmb",
        # "liemid_Newmark",
        # "IncrSO3_LIEMIDEA",
        "IncrSO3_SVQ",
        "IncrSO3_VLV"
    ]

    for method in Methods:
        try:
            step = getattr(sys.modules[__name__], method.lower())

            start_time = time.time()
            print(method)
            ttNMB, _, _, _, HNMB = integrate(dt, hc, T, J, Wi0, Gi0, Gm, step, 
                                             Loading, potential)
            end_time = time.time()
            print(f"Elapsed time: {end_time - start_time} seconds")
        
            plt.plot(ttNMB, HNMB - HNMB[0], label=method[8:])
            pass
        except Exception as e: 
            print(e)
    
    plt.legend()
    plt.show()
