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
    
    M0 = Loading.Torque(t0, R0, Loading)  # torque applied at time 0
    TH1 = dt * W0  # initialization

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


