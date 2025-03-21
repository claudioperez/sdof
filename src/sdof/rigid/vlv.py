"""
For svq:

% Implements Lie-Newmark algorithm with gamma = 1/2, beta = 0, proposed in  
%
% Simo, J. C. and Vu-Quoc, L. [1988].
% On the dynamics in space of rods undergoing large motions - a
% geometrically exact approach.
% COMPUTER METHODS IN APPLIED MECHANICS AND ENGINEERING. 
% vol. 66, pp. 125-161.
%
%... Last Modified: May 27, 2010
%... Authors:       Nawaf Bou-Rabee (nawaf@cims.nyu.edu)
%                   Giulia Ortolan (ortolang@dei.unipd.it)
%                   Alessandro Saccon (asaccon@isr.ist.utl.pt)
"""

def IncrSO3_LIEMIDEA(t0, p0, W0, R0, dt, J, state, Loading):
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

    TH1 = dt * W0  # initialization
    M0 = Loading.Torque(t0, R0, Loading)  # torque applied at time 0

    for j in range(maxiter):
        expx = ExpSO3(-TH1 / 2)

        v = 0.5 * dt * expx @ (p0 + 0.5 * dt * M0)
        f = -J @ TH1 + v

        if np.linalg.norm(f) < 1e-10:
            break
        elif j == maxiter - 1:
            print('maximum # of iterations exceeded')
  
        Jac = -J + 0.5 * Spin(v) @ dExpSO3(-TH1 / 2)
        
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
    M0 = Loading.Torque(t0 + dt, R1, Loading)

    #... Step 6) Explicit Velocity Update ................................#
    p1 = expx.T @ P1 + 0.5 * dt * M0    # body angular momentum
    w1 = np.linalg.solve(J, p1)         # body angular velocity

    return p1, w1, R1, state

def ExpSO3(phi):
    """
    Exponential map for SO(3)
    """
    phi_norm = np.linalg.norm(phi)
    if phi_norm < 1e-10:
        return np.eye(3) + Spin(phi)
    else:
        axis = phi / phi_norm
        skew_phi = Spin(phi)
        return (
            np.eye(3) +
            (np.sin(phi_norm) / phi_norm) * skew_phi +
            ((1 - np.cos(phi_norm)) / (phi_norm ** 2)) * (skew_phi @ skew_phi)
        )

def Spin(v):
    """
    Skew-symmetric matrix (Spin matrix) for vector v
    """
    return np.array([
        [    0, -v[2],  v[1]],
        [ v[2],     0, -v[0]],
        [-v[1],  v[0],     0]
    ])

def dExpSO3(phi):
    """
    Derivative of the exponential map for SO(3)
    """
    phi_norm = np.linalg.norm(phi)
    if phi_norm < 1e-10:
        return 0.5 * Spin(phi)
    else:
        axis = phi / phi_norm
        skew_phi = Spin(phi)
        a = (1 - np.cos(phi_norm)) / (phi_norm ** 2)
        b = (phi_norm - np.sin(phi_norm)) / (phi_norm ** 3)
        return (
            0.5 * skew_phi +
            a * (skew_phi @ skew_phi) +
            b * (skew_phi @ skew_phi @ skew_phi)
        )

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
         + a2 * (-HatSO3(J @ f) + Spin(f) @ J)
         + b1 * J @ np.outer(f, f)
         + b2 * HatSO3(f) @ J @ np.outer(f, f))
    return d


def simo_wong_algo_c1(t0, p0, w0, R0, dt, J, state, Loading):
    """
    Constructor for momentum and energy conserving integrator from Simo, Wong.
    The momentum and energy conserving integrator from Simo, Wong.
    """
    gamma = 1
    beta = 1/2  # the Newmark parameters: these are values determined as optimal by SW

    nsteps = 1
    maxiter = 30
    tol = 1000 * np.finfo(float).eps
    alph = min(beta / gamma, 1)

    t = t0

    # calculate effective loads at t_n-1
    tn1 = Loading.Torque(t, R0, Loading)
    if not hasattr(state, 'Alpha') or state.Alpha is None:
        # Alphan1 = I\(Rn1'*tn1 - Spin(Omegan1)*I*Omegan1)
        Alphan1 = np.linalg.solve(J, tn1 - Spin(w0) @ J @ w0)
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


def incr_so3_vlv(t0, p0, w0, R0, dt, J, state, loading):
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
    # maximum # of Newton iterations
    maxiter = 10

    # torque applied at time n
    tau = loading.Torque(t0, R0, loading)

    # Step 1) Newton Solve
    P1 = p0.copy()  # initialization
    W1 = w0.copy()  # initialization

    II = np.diag(J)

    for j in range(maxiter):
        F = np.dot(W1.T, P1)
        PxW = 0.5 * dt * np.array([
            P1[1]*W1[2] - P1[2]*W1[1],
            P1[2]*W1[0] - P1[0]*W1[2],
            P1[0]*W1[1] - W1[0]*P1[1]
        ])  # 0.5*dt*(P1 x W1)

        f = p0 - P1 + PxW - 0.25 * dt**2 * F * W1 + 0.5 * dt * tau

        if np.linalg.norm(f) < 1e-10:
            break
        elif j == maxiter - 1:
            print(f"{__file__} - maximum # of iterations exceeded")

        Jac = (-J + 0.5 * dt * np.array([
            [0, -P1[2] + II[1,1]*W1[2], P1[1] - II[2,2]*W1[1]],
            [P1[2] - II[0,0]*W1[2], 0, -P1[0] + II[2,2]*W1[0]],
            [-P1[1] + II[0,0]*W1[1], P1[0] - II[1,1]*W1[0], 0]
        ]) - 0.25 * dt**2 * F * np.eye(3))

        # Solve Jac * delta = f for delta
        delta = np.linalg.solve(Jac, f)
        W1 = W1 - delta
        P1 = J @ W1

    # Step 2) Explicit Configuration Update
    R1 = R0 @ CaySO3(dt * W1)

    # Step 3) Explicit Momentum Update
    tau = loading.Torque(t0 + dt, R1, loading)  # torque applied at time n+1

    p1 = P1 + PxW + 0.25 * dt**2 * (np.dot(W1, P1)) * W1 + 0.5 * dt * tau  # body angular momentum
    w1 = np.linalg.solve(J, p1)  # body angular velocity

    return p1, w1, R1, state
