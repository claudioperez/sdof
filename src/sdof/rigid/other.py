import numpy as np
import warnings
from .util import dExpSO3, dLogSO3, ExpSO3, CaySO3, HatSO3, EPS

def incrso3_trapm     (t0, p0, w0, R0, dt, J, Loading, state):
    """
    Trapezoidal rule with momentum conservation.
    Described as TRAPM in [3].
    Dynamically equivalent to the implicit midpoint Lie integrator.
    """
    maxiter = 150
    zeta = 0.5
    tol = 100 * EPS

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
                warnings.warn(f"Failed to converge: ||residual||={np.linalg.norm(Pin)}")
            break
        if np.linalg.norm(Pin - Hx) < tol:
            break

    return Pin, Omegan, Rn, None


def incrso3_trapm_old (t0, p0, w0, R0, dt, J, Loading, state):
    """
    Trapezoidal rule with momentum conservation integrator constructor.
    Described as TRAPM in the paper on dynamically equivalent implicit integrators.
    Dynamically equivalent to the implicit midpoint Lie integrator.
    """
    maxiter = 10
    Ptol = 100 * np.finfo(float).eps

    M0 = Loading.Torque(t0, R0)
    H0 = J@w0
    DRn = ExpSO3(dt / 2 * w0) @ ExpSO3(dt / 2 * w0)  # predictor
    H1 = DRn.T @ (H0 + dt / 2 * M0) + dt / 2 * M0  # predictor

    iter = 0
    while True:
        Hx = H1
        w1 = np.linalg.solve(J, H1)
        DRn = ExpSO3(dt*0.5 * w0) @ ExpSO3(dt*0.5 * w1)
        Rn = R0 @ DRn
        M1 = Loading.Torque(t0 + dt, Rn)
        H1 = DRn.T @ H0 + dt*0.5 * DRn.T @ M0 + dt*0.5 * M1
        iter += 1
        if iter > maxiter:
            if np.max(np.abs(H1 - Hx)) > 1e-9:
                warnings.warn(f"Failed to converge: ||residual||={np.linalg.norm(H1)}")
            break
        if np.linalg.norm(H1 - Hx) < Ptol:
            break

    return H1, w1, Rn, None


def incrso3_imidm     (t0, p0, w0, R0, dt, J, Loading, state):
    """
    Implicit mid-point Lie integrator constructor. Described as LIEMID[I] in
    the paper on the explicit midpoint Lie algorithm, and as IMIDM in [3]
    """
    wa = _solve_imidm(J, t0, dt, R0, J @ w0, Loading)
    Rn = R0 @ ExpSO3(dt * wa)
    ma = R0 @ Loading.Torque(t0 + dt*0.5, R0@ExpSO3(dt*0.5 * wa))
    w1 = np.linalg.solve(J, Rn.T @ (R0 @ J @ w0 + dt * ma))

    return J@w1, w1, Rn, None

def _solve_imidm(J, t, dt, R0, H0, Loading):
    maxiter = 50
    tol = 1000 * EPS

    wa = np.linalg.solve(J, H0)
    ma = R0 @ Loading.Torque(t + dt*0.5, R0@ExpSO3(dt*0.5 * wa))
    B = H0 + dt / 2 * R0.T @ ma
    res = wa - np.linalg.solve(J, ExpSO3(-dt / 2 * wa) @ B)
    invIdt2 = np.linalg.inv(J) * dt**2 / 8

    iter = 0
    while np.max(np.abs(res)) > tol:
        K = np.eye(3) \
          - np.linalg.solve(J, dt / 2 * HatSO3(B)) \
            + invIdt2 @ (HatSO3(wa) @ HatSO3(B) + HatSO3(HatSO3(wa) @ B))
        wa = wa - np.linalg.solve(K, res)
        ma = R0 @ Loading.Torque(t + dt / 2, R0 @ ExpSO3(dt / 2 * wa))
        B  = H0 + dt / 2 * R0.T @ ma
        res = wa - np.linalg.solve(J, ExpSO3(-dt / 2 * wa) @ B)
        iter += 1
        if iter > maxiter:
            warnings.warn(f"Failed to converge after {maxiter} iterations")
            break

    return wa


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
        f += np.linalg.solve(DG(f, J), r)
        if iter > maxiter:
            if np.max(np.abs(r)) > 1e-9:
                warnings.warn(f"Failed to converge: ||residual||={np.linalg.norm(r)}")
            break
        iter += 1

    DR = ExpSO3(f)
    Rn = R0 @ DR
    Tn = Loading.Torque(t0 + dt, Rn)
    Omegan = np.linalg.solve(J, DR.T @ J @ w0 + dt * Tn)

    P1 = J @ Omegan
    return P1, Omegan, Rn, None

def DG(f, J):
    """
    Compute the DG matrix based on vector f and matrix J.
    """
    nf = np.linalg.norm(f)
    cs = np.cos(nf)
    sn = np.sin(nf)

    a1 = sn / nf
    a2 = (1 - cs) / nf**2
    b1 = (cs * nf - sn) / nf**3
    b2 = (sn * nf - 2 * (1 - cs)) / nf**4

    d = (a1 * J
         + a2 * (-HatSO3(J @ f) + HatSO3(f) @ J)
         + b1 * J @ np.outer(f, f)
         + b2 * HatSO3(f) @ J @ np.outer(f, f))
    return d

def incrso3_liemidea  (t0, p0, w0, R0, dt, J, Loading, state, maxiter = 10):

    TH1 = dt * w0  # initialization
    M0 = Loading.Torque(t0, R0)

    for j in range(maxiter):
        expx = ExpSO3(-TH1*0.5)
        v = 0.5 * dt * np.dot(expx, (p0 + 0.5*dt * M0))
        f = -J@TH1 + v

        if np.linalg.norm(f) < 1e-10:
            break
        elif j == maxiter - 1:
            warnings.warn(f"maximum number of iterations exceeded")

        Jac = -J + 0.5 * HatSO3(v)@dExpSO3(-TH1 * 0.5)

        TH1 = TH1 + np.linalg.solve(Jac, f)

    # Step 2) Explicit Configuration Update
    R1 = R0@ExpSO3(TH1)

    # Update momentum and angular velocity
    p1 = p0 + dt*M0
    w1 = np.linalg.solve(J, p1)

    return p1, w1, R1, None

def incrso3_liemidea  (t0, p0, w0, R0, dt, J, Loading, state, maxiter = 10):

    """
    Implements LIEMID[EA] algorithm proposed in  

    Krysl, P. [2005].
      Explicit momentum-conserving integrator for dynamics of rigid bodies 
      approximating the midpoint Lie algorithm
      INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN ENGINEERING 
      63 (15): 2171-2193.

    Last Modified: May 27, 2010
    Authors:       Nawaf Bou-Rabee (nawaf@cims.nyu.edu)
                    Giulia Ortolan (ortolang@dei.unipd.it)
                    Alessandro Saccon (asaccon@isr.ist.utl.pt)
    """
    #==========================================================================
    # dt =                the time-stepsize for the simulation;
    # w0 = wi;          # body angular velocity at time 0
    # p0 = II * w0      # body angular momentum at time 0
    # G0 = Gi;          # configuration at time 0
    #==========================================================================
    maxiter = 10          # maximum # of Newton iterations

    #... Step 1) Newton Solve ............................................#
    
    TH1 = dt * w0  # initialization
    M0 = Loading.Torque(t0, R0)  # torque applied at time 0

    for j in range(maxiter):
        expx = ExpSO3(-TH1 / 2)
        
        v = 0.5 * dt * expx @ (p0 + 0.5 * dt * M0)
        f = -J @ TH1 + v

        if np.linalg.norm(f) < 1e-10:
            break
        elif j == maxiter - 1:
            print('maximum # of iterations exceeded')
  
        Jac = -J + 0.5 * HatSO3(v) @ dExpSO3(-TH1 / 2)
        
        # Solve for delta_TH1
        delta_TH1 = np.linalg.solve(Jac, f)
        TH1 = TH1 - delta_TH1
        
    #... Step 2) Explicit Configuration Update ...........................#
    
    expx = ExpSO3(TH1)
    R1 = R0 @ expx
      
    #... Step 3) Explicit Velocity Update ................................#
    
    P1 = expx.T @ (p0 + 0.5 * dt * M0)
    W1 = np.linalg.solve(J, P1)
    
    #... Step 4) Newton Solve ............................................#
    
    TH1 = W1  # initialization

    for j in range(maxiter):
        v = 0.5 * dt * ExpSO3(-TH1 / 2) @ P1
        f = -J @ TH1 + v

        if np.linalg.norm(f) < 1e-10:
            break
        elif j == maxiter - 1:
            print('maximum # of iterations exceeded')
        
        Jac = -J + 0.5 * HatSO3(v) @ dExpSO3(-TH1 / 2)
        
        # Solve for delta_TH1
        delta_TH1 = np.linalg.solve(Jac, f)
        TH1 = TH1 - delta_TH1
    
    #... Step 5) Explicit Configuration Update ...........................#
    expx = ExpSO3(TH1)
    R1 = R1 @ expx

    #... Note: only one force evaluation per step
    M0 = Loading.Torque(t0 + dt, R1)

    #... Step 6) Explicit Velocity Update ................................#
    p1 = expx.T @ P1 + 0.5 * dt * M0    # body angular momentum
    w1 = np.linalg.solve(J, p1)         # body angular velocity

    return p1, w1, R1, state


def incrso3_akw       (t0, p0, w0, R0, dt, J, Loading, state):
    """
    Implicit Austin, Krishnaprasad, Wang mid-point method
    """
    # akw_mid 
    H0 = J @ w0
    Hn = _solve_akw(J, t0, dt, R0, H0, Loading)
    wn = np.linalg.solve(J,  Hn)
    wa = np.linalg.solve(J, (Hn + H0) / 2)
    R1 = R0@CaySO3(dt * wa)
    return (Hn,
            wn, 
            R1, 
            None)

def _solve_akw(J, t0, dt, R0, H0, Loading):
    maxiter = 25
    tol = 100 * EPS

    Hn = H0
    M0 = Loading.Torque(t0, R0)
    iter = 0

    while True:
        Hx = Hn
        wa = np.linalg.solve(J, (Hn + H0) / 2)
        Rn = R0 @ CaySO3(dt * wa)
        M1 = Loading.Torque(t0 + dt, Rn)
        Hn = H0 + dt * HatSO3(J @ wa) @ wa + dt / 2 * (M1 + M0)
        iter += 1
        if iter > maxiter:
            if np.max(np.abs(Hn - Hx)) > 1e-12:
                warnings.warn(f"Failed to converge: residual = {np.linalg.norm(Hn)}")
            break
        if np.linalg.norm(Hn - Hx) < tol:
            break

    return Hn


def discrete_moser_veselov_step(dt, R, Pi, I):
    """
    Perform a single time step of the Moser–Veselov algorithm for the free rigid body.

    Parameters:
    dt (float): Time step size.
    R (np.ndarray): Rotation matrix (3x3) representing the current orientation.
    Pi (np.ndarray): Angular momentum vector (3,) in the body frame.
    I (np.ndarray): Inertia tensor (3x3), assumed diagonal for simplicity.

    Returns:
    tuple: Updated (R_new, Pi_new).
    """
    # Calculate the angular velocity in the body frame
    omega = np.linalg.solve(I, Pi)  # omega = I^(-1) * Pi

    # Compute the skew-symmetric matrix of omega
    Omega = np.array([
        [0, -omega[2], omega[1]],
        [omega[2], 0, -omega[0]],
        [-omega[1], omega[0], 0]
    ])

    # Calculate the update rotation matrix by using the matrix exponential
    R_new = R @ ExpSO3(dt * Omega)

    # Calculate the midpoint angular momentum (implicit method for stability)
    identity = np.eye(3)
    term = identity + (dt / 2) * Omega
    term_inv = np.linalg.inv(term)

    # Update angular momentum Pi using the implicit midpoint method
    Pi_new = term_inv @ (identity - (dt / 2) * Omega) @ Pi

    return R_new, Pi_new

def shake_rattle_step(pos, R, vel, omega, force, torque, mass, inertia, dt, constraint_length, tol=1e-6, max_iter=100):
    """
    Integrate the position and rotation of a rigid body with SHAKE and RATTLE constraints.

    Parameters:
    - pos: Initial center of mass position of the rigid body.
    - R: Initial rotation matrix (3x3) representing orientation.
    - vel: Initial linear velocity.
    - omega: Initial angular velocity (in body frame).
    - force: External force acting on the center of mass.
    - torque: External torque acting on the center of mass.
    - mass: Mass of the rigid body.
    - inertia: Inertia tensor (3x3 matrix in body frame).
    - dt: Time step.
    - constraint_length: Fixed distance constraint for demonstration.
    - tol: Tolerance for SHAKE/RATTLE constraint satisfaction.
    - max_iter: Maximum number of iterations for constraint correction.

    Returns:
    - pos_new, R_new, vel_new, omega_new: Updated position, rotation, linear and angular velocities.

    Example:
    # Initial conditions
    pos = np.array([0.0, 0.0, 0.0])  # Initial position of the center of mass
    R = np.eye(3)  # Initial rotation matrix (identity for no rotation)
    vel = np.array([0.0, 0.0, 0.0])  # Initial linear velocity
    omega = np.array([0.0, 1.0, 0.0])  # Initial angular velocity in body frame
    force = np.array([0.0, -9.8, 0.0])  # External force (gravity)
    torque = np.array([0.0, 0.0, 0.1])  # External torque
    mass = 1.0  # Mass of the rigid body
    inertia = np.diag([0.1, 0.2, 0.3])  # Arbitrary inertia tensor in body frame
    dt = 0.01  # Time step
    constraint_length = 1.0  # Arbitrary constraint distance for demo

    # Run a single time step
    pos_new, R_new, vel_new, omega_new = integrate_rigid_body_step(
        pos, R, vel, omega, force, torque, mass, inertia, dt, constraint_length
    )

    print("New position:", pos_new)
    print("New rotation matrix:\n", R_new)
    print("New linear velocity:", vel_new)
    print("New angular velocity:", omega_new)

    """
    # Predict new position and velocity
    pos_pred = pos + vel * dt + 0.5 * force / mass * dt**2
    vel_pred = vel + force / mass * dt

    # Predict new rotation and angular velocity
    R_pred = R + R @ HatSO3(omega) * dt  # Update rotation matrix
    inertia_world = R @ inertia @ R.T  # Inertia tensor in world frame
    omega_pred = omega + np.linalg.inv(inertia) @ (torque - np.cross(omega, inertia @ omega)) * dt

    # Apply SHAKE for position constraints
    for _ in range(max_iter):
        delta_pos = pos_pred - pos
        dist = np.linalg.norm(delta_pos)
        if abs(dist - constraint_length) < tol:
            break
        correction = (dist - constraint_length) * (delta_pos / dist) / 2
        pos_pred -= correction

    # Apply SHAKE for rotation constraints
    for _ in range(max_iter):
        rotation_deviation = R_pred.T @ R - np.eye(3)  # Measure deviation from orthogonality
        deviation_norm = np.linalg.norm(rotation_deviation)
        if deviation_norm < tol:
            break
        R_pred -= 0.5 * rotation_deviation  # Apply correction towards orthogonality

    # RATTLE for velocity constraints
    for _ in range(max_iter):
        delta_v = vel_pred - vel
        if np.linalg.norm(delta_v) < tol:
            break
        correction_v = (np.dot(delta_pos, delta_v) / np.dot(delta_pos, delta_pos)) * delta_pos
        vel_pred -= correction_v

    # RATTLE for angular velocity constraints
    for _ in range(max_iter):
        rotation_axis = R_pred[:, 0]  # Use any fixed axis in body frame for constraint
        correction_omega = np.cross(rotation_axis, omega_pred)
        if np.linalg.norm(correction_omega) < tol:
            break
        omega_pred -= correction_omega / np.dot(rotation_axis, rotation_axis)

    return pos_pred, R_pred, vel_pred, omega_pred

def integrate_step(pos, vel, force, mass, dt, constraint_length, tol=1e-6, max_iter=100):
    # Predictor step: Update positions and velocities
    pos_pred = pos + vel * dt + 0.5 * force / mass * dt**2
    vel_pred = vel + force / mass * dt
    
    # Apply SHAKE to enforce position constraints
    for _ in range(max_iter):
        delta = pos_pred[1] - pos_pred[0]
        dist = np.linalg.norm(delta)
        if abs(dist - constraint_length) < tol:
            break
        correction = (dist - constraint_length) * (delta / dist) / 2
        pos_pred[0] += correction
        pos_pred[1] -= correction

    # RATTLE for velocity constraints
    for _ in range(max_iter):
        delta_v = vel_pred[1] - vel_pred[0]
        dot_product = np.dot(delta, delta_v)
        if abs(dot_product) < tol:
            break
        correction = (dot_product / dist**2) * (delta / dist) / 2
        vel_pred[0] += correction
        vel_pred[1] -= correction

    return pos_pred, vel_pred

# Initial conditions
pos = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])  # Initial positions of two points
vel = np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])  # Initial velocities
force = np.array([[0.0, -9.8, 0.0], [0.0, -9.8, 0.0]])  # Gravitational force
mass = 1.0  # Mass of points
dt = 0.01  # Time step
constraint_length = 1.0  # Fixed distance constraint

# Run a single time step
pos_new, vel_new = integrate_step(pos, vel, force, mass, dt, constraint_length)

print("New positions:", pos_new)
print("New velocities:", vel_new)

