import numpy as np

def lsdof(dt, omega, p, zeta=None, u0=None, v0=None):
    """
    Transient response of a linear SDOF system using exact integration
    with piecewise linear excitation.
    
    Parameters:
    Deltat (float): Time step of integration.
    omega (np.ndarray): Row vector of eigenfrequencies.
    p (np.ndarray): Acceleration history (force/mass).
    zeta (np.ndarray, optional): Row vector of damping ratios. Defaults to zero.
    u0 (np.ndarray, optional): Row vector of initial displacements. Defaults to zero.
    udot0 (np.ndarray, optional): Row vector of initial velocities. Defaults to zero.
    
    Returns:
    u (np.ndarray): Displacement history arranged column-wise.
    v (np.ndarray): Velocity history arranged column-wise.
    a (np.ndarray): Acceleration history arranged column-wise.

    Compare with FEDEASLab 'LSDOF_LinearWilson'
    """
    nstep = len(p)
    nfreq = len(omega)
    
    # Default damping ratio
    if zeta is None:
        zeta = np.zeros(nfreq)
    
    omega_d = omega * np.sqrt(1 - zeta ** 2)
    omegaDt = omega_d * dt
    
    s = np.sin(omegaDt)
    c = np.cos(omegaDt)
    
    k = omega ** 2           # stiffness per unit mass
    cd = 2 * zeta * omega    # viscous damping per unit mass
    
    ex = np.exp(-zeta * omega * dt)
    C1 = zeta / np.sqrt(1 - zeta ** 2)
    C2 = 2 * zeta / (omega * dt)
    C3 = omega / np.sqrt(1 - zeta**2)
    C4 = (1 - 2 * zeta ** 2) / omegaDt

    A = ex * (C1 * s + c)
    B = ex * s / omega_d
    C = (C2 + ex * ((C4 - C1) * s - (1 + C2) * c)) / k
    D = (1 - C2 + ex * (-C4 * s + C2 * c)) / k
    
    Ap = -ex * C3 * s
    Bp = ex * (c - C1 * s)
    Cp = (-1 + ex * ((C3 * dt + C1) * s + c)) / (k * dt)
    Dp = (1 - ex * (C1 * s + c)) / (k * dt)
    
    u = np.zeros((nstep, nfreq))
    v = np.zeros((nstep, nfreq))
    a = np.zeros((nstep, nfreq))
    
    if u0 is not None and v0 is not None:
        if len(u0) != nfreq or len(v0) != nfreq:
            raise ValueError("Number of initial condition values must match number of modes")
        u[0, :] = u0
        v[0, :] = v0
    
    a[0, :] = p[0] - cd * v[0, :] - k * u[0, :]
    
    for i in range(nstep - 1):
        u[i + 1, :] = A * u[i, :] + B * v[i, :] + C * p[i] + D * p[i + 1]
        v[i + 1, :] = Ap * u[i, :] + Bp * v[i, :] + Cp * p[i] + Dp * p[i + 1]
        a[i + 1, :] = p[i + 1] - cd * v[i + 1, :] - k * u[i + 1, :]
    
    return u, v, a


def _euler(m=1, c=.1, k=1, u0=1, v0=0, n=8, dt=0.05):
    """
    Free response using Euler's method to perform numerical integration.

    Parameters
    ----------
    m, c, k: float
        Mass, damping and stiffness.
    u0, v0: float
        Initial conditions
    n: int
        The number of steps
    dt: float
        The step size.

    Returns
    -------
    t, x, v: array
        Time, displacement, and velocity

    Examples
    --------
    >>> import vibration_toolbox as vtb
    >>> vtb._euler(m=1, c=.1, k=1, u0=1, v0=0, n=8, dt=0.05)  # doctest: +SKIP
    (array([0.  , 0.05, 0.1 , 0.15, 0.2 , 0.25, 0.3 , 0.35, 0.4 ]),
    array([[ 1.  ,  0.  ],
           [ 1.  , -0.05],
           [ 1.  , -0.1 ],
           [ 0.99, -0.15],
           [ 0.99, -0.2 ],
           [ 0.98, -0.25],
           [ 0.96, -0.29],
           [ 0.95, -0.34],
           [ 0.93, -0.39]]))
    # creates the state space matrix
    """
    A = np.array([[0, 1],
                  [-k / m, -c / m]])
    # creates the x array and set the first line according to the initial
    # conditions
    x = np.zeros((n + 1, 2))
    x[0] = u0, v0

    for i in range(0, n):
        x[i + 1] = x[i] + dt * A@x[i]

    t = np.linspace(0, n * dt, n + 1)

    return x


def LIDA(dt,xgtt,omega,ksi,u0,ut0,rinf):
#
# Linear Implicit Dynamic Analysis (LIDA)
# https://viewer.mathworks.com/?viewer=plain_code&url=https%3A%2F%2Fjp.mathworks.com%2Fmatlabcentral%2Fmlc-downloads%2Fdownloads%2Fsubmissions%2F44835%2Fversions%2F11%2Fcontents%2FRESPONSE_SPECTRA%2FLIDA.m&embed=web
#
#     Linear implicit direct time integration of second order differential
#     equation of motion of dynamic response of linear elastic SDOF systems
#
# Description
#     The General Single Step Single Solve (GSSSS) family of algorithms
#     published by X.Zhou & K.K.Tamma (2004) is employed for direct time
#     integration of the general linear or nonlinear structural Single
#     Degree of Freedom (SDOF) dynamic problem. The optimal numerical
#     dissipation and dispersion zero order displacement zero order
#     velocity algorithm designed according to the above journal article,
#     is used in this routine. This algorithm encompasses the scope of
#     Linear Multi-Step (LMS) methods and is limited by the Dahlquist
#     barrier theorem (Dahlquist,1963). The force - displacement - velocity
#     relation of the SDOF structure is linear.
#
# Input parameters
#     #dt# (scalar): time step
#     #xgtt# ([#nstep# x 1]): column vector of the acceleration history of
#         the excitation imposed at the base. #nstep# is the number of time
#         steps of the dynamic response.
#     #omega# (scalar): eigenfrequency of the structure in rad/sec
#     #ksi# (scalar): ratio of critical damping of the SDOF system
#     #u0# (scalar): initial displacement of the SDOF system
#     #ut0# (scalar): initial velocity of the SDOF system
#     #rinf# (scalar): minimum absolute value of the eigenvalues of the
#         amplification matrix. For the amplification matrix see eq.(61) in
#         Zhou & Tamma (2004).
#
# Output parameters
#     #u# ([#nstep# x 1]): time-history of displacement
#     #ut# ([#nstep# x 1]): time-history of velocity
#     #utt# ([#nstep# x 1]): time-history of acceleration
#
# Example (Figure 6.6.1 in Chopra, Tn=1sec)
#     dt=0.02
#     fid=fopen('elcentro.dat','r')
#     text=textscan(fid,'%f %f')
#     fclose(fid)
#     xgtt=9.81*text{1,2}
#     Tn=1
#     omega=2*pi/Tn
#     ksi=0.02
#     u0=0
#     ut0=0
#     rinf=1
#     [u,ut,utt] = Linear_implicit_dynamic_analysis(dt,xgtt,omega,ksi,u0,ut0,rinf)
#     D=max(abs(u))/0.0254
#
#__________________________________________________________________________
# Copyright (c) 13-Sep-2015
#     George Papazafeiropoulos
#     First Lieutenant, Infrastructure Engineer, Hellenic Air Force
#     Civil Engineer, M.Sc., Ph.D. candidate, NTUA
#     Email: gpapazafeiropoulos@yahoo.gr
#     Website: http://users.ntua.gr/gpapazaf/
# _________________________________________________________________________
#% Integration constants
# zero-order displacement & velocity overshooting behavior and
# optimal numerical dissipation and dispersion
    W1, W1L1, W2L2, W3L3, W1L4, W2L5, W1L6, l1, l2, l3, l4, l5 = parameters(rinf)


    #% Transfer function denominator
    Omega=omega*dt
    D = W1L6 + 2*W2L5*ksi*Omega + W3L3*Omega**2
    A31 = -Omega**2/D
    A32 = -1/D*(2*ksi*Omega + W1L1*Omega**2)
    A33 = 1-1/D*(1+2*W1L4*ksi*Omega + W2L2*Omega**2)
    A11 = 1+l3*A31
    A12 = l1+l3*A32
    A13 = l2-l3*(1-A33)
    A21 = l5*A31
    A22 = 1+l5*A32
    A23 = l4-l5*(1-A33)
    # Amplification matrix
    A=[[A11, A12, A13],
       [A21, A22, A23],
       [A31, A32, A33]]
    # Amplification matrix invariants
    A1=A[0,0]+A[1,1]+A[2,2]
    A2=A[0,0]*A[1,1]-A(1,2)*A(2,1)+A[0,0]*A[2,2]-A(1,3)*A(3,1)+A[1,1]*A[2,2]-A(2,3)*A(3,2)
    A3=A[0,0]*A[1,1]*A[2,2]-A[0,0]*A(2,3)*A(3,2)-A(1,2)*A(2,1)*A[2,2]+A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)-A(1,3)*A[1,1]*A(3,1)
    
    # Transfer function denominator
    a = [1, -A1, A2, -A3]
    #% Transfer function nominator
    B1=1/D*dt**2*l3*W1
    B2=1/D*dt**2*(l3*(1-W1)-(A22+A33)*l3*W1+A12*l5*W1+A13*W1)
    B3=1/D*dt**2*(-(A22+A33)*l3*(1-W1)+A12*l5*(1-W1)+A13*(1-W1)+(A22*A33-A23*A32)*l3*W1-(A12*A33-A13*A32)*l5*W1+(A12*A23-A13*A22)*W1)
    B4=1/D*dt**2*((A22*A33-A23*A32)*l3*(1-W1)-(A12*A33-A13*A32)*l5*(1-W1)+(A12*A23-A13*A22)*(1-W1))
    b=[B1,B2,B3,B4]
    #% form initial conditions for filter function
    # equivalent external force
    f=-xgtt
    # stiffness
    k=omega**2
    # damping constants
    c=2*omega*ksi
    # initial acceleration
    utt0=-f(1)-(k*u0+c*ut0)
    U_1 = np.linalg.solve(A, [u0, dt*ut0, dt**2*utt0])
    u_1 = U_1[0]
    U_2 = np.linalg.solve(A,U_1)
    u_2 = U_2(1)
    ypast=[u0,u_1,u_2]
    vinit = np.zeros((1,3))
    vinit(3:-1:1) = filter(-a(4:-1:2), 1, ypast)

    #% main dynamic analysis
    u = filter(b,a,f,vinit)

    #% calculate velocity from the following system of equations:
    # 1st: the first scalar equation of the matrix equation (60) in [1]
    # 2nd: equation of motion (eq.6.12.3 in Chopra:Dynamics of Structures)
    C_u=omega**2*A(1,3)*dt**2-A[0,0]
    C_f=-A(1,3)*dt**2
    C_ut=A(1,2)*dt - A(1,3)*dt**2*2*ksi*omega
    L=1/D*l3*dt**2*((1 - W1)*[0, f(1:end-1)] + W1*f)
    ut = (u + C_u@[u0, u[:-2]] + C_f*[0, f(1:end-1)] - L)/C_ut
    #ut=[ut0; diff(u)/dt]
    #% calculate acceleration from equation of motion
    utt = -omega**2*u - 2*ksi*omega*ut

    return [u,ut,utt]

#
# STATE-SPACE
#
def _compute_a_and_b(zeta, w, dt):
    """
    From the paper by Nigam and Jennings (1968), computes the two matrices.

    :param zeta: damping ratio
    :param w: angular frequencies
    :param dt: time step
    :return: matrices A and B
    """

    zeta2 = zeta * zeta  # D2
    w2 = w ** 2  # W2
    one_ov_w2 = 1. / w2  # A7
    sqrt_b2   = np.sqrt(1. - zeta2)
    w_sqrt_b2 = w * sqrt_b2  # A1

    exp_b       = np.exp(-zeta * w * dt)  # A0
    two_b_ov_w2 = (2 * zeta ** 2 - 1) / (w ** 2 * dt)
    two_b_ov_w3 = 2 * zeta / (w ** 3 * dt)

    sin_wsqrt = np.sin(w_sqrt_b2 * dt)  # A2
    cos_wsqrt = np.cos(w_sqrt_b2 * dt)  # A3

    # Transition matrix
    a_11 = exp_b * (zeta / sqrt_b2 * sin_wsqrt + cos_wsqrt)  # Eq 2.7d(1)
    a_12 = exp_b / (w * sqrt_b2) * sin_wsqrt  # Eq 2.7d(2)
    a_21 = -w / sqrt_b2 * exp_b * sin_wsqrt    # Eq 2.7d(3)
    a_22 = exp_b * (cos_wsqrt - zeta / sqrt_b2 * sin_wsqrt)  # Eq 2.7d(4)

    A = np.array([[a_11, a_12],
                  [a_21, a_22]])

    # B matrix
    bsqrd_ov_w2_p_zeta_ov_w = two_b_ov_w2 + zeta / w
    sin_ov_wsqrt = sin_wsqrt / w_sqrt_b2
    xwcos = zeta * w * cos_wsqrt
    wsqrtsin = w_sqrt_b2 * sin_wsqrt

    # Eq 2.7e
    b_11 =  exp_b * (bsqrd_ov_w2_p_zeta_ov_w * sin_ov_wsqrt + (two_b_ov_w3 + one_ov_w2) * cos_wsqrt) - two_b_ov_w3
    b_12 = -exp_b * (two_b_ov_w2 * sin_ov_wsqrt + two_b_ov_w3 * cos_wsqrt) - one_ov_w2 + two_b_ov_w3
    b_21 =  exp_b * (bsqrd_ov_w2_p_zeta_ov_w * (cos_wsqrt - zeta / sqrt_b2 * sin_wsqrt)
                  - (two_b_ov_w3 + one_ov_w2) * (wsqrtsin + xwcos)) + one_ov_w2 / dt
    b_22 = -exp_b * (two_b_ov_w2 * (cos_wsqrt - zeta / sqrt_b2 * sin_wsqrt) - two_b_ov_w3 * (wsqrtsin + xwcos)) - one_ov_w2 / dt

    B = np.array([[b_11, b_12],
                  [b_21, b_22]])

    return A, B


def exponential(acc, dt, model):
    """

    Implementation of the response spectrum calculation.

    See Nigam and Jennings (1968).

    Ref: Nigam, N. C., Jennings, P. C. (1968) Digital calculation of response spectra from strong-motion earthquake
    records. National Science Foundation.

    :param acc: acceleration in m/s2
    :param periods: response periods of interest
    :param dt: time step of the acceleration time series
    :param zeta: damping ratio
    :return: response displacement, response velocity, response acceleration

    from eqsig.sdof
    """

    acc     = -np.array(acc, dtype=float)
#   periods =  np.array(periods, dtype=float)
#   if periods[0] == 0:
#       s = 1
#   else:
#       s = 0

    w    = model.w # 6.2831853 / periods[s:]

    zeta = model.zeta

    # TODO: delta_t should be less than period / 20

    a, b = _compute_a_and_b(zeta, w, dt)
#   resp_u = np.zeros([len(periods), len(acc)], dtype=float)
#   resp_v = np.zeros([len(periods), len(acc)], dtype=float)
    x = np.zeros((len(acc), 2))
    for i in range(len(acc) - 1):  # possibly speed up using scipy.signal.lfilter
        # x_i+1 = A cross (u, v) + B cross (acc_i, acc_i+1)  # Eq 2.7a
        x[i+1] = a@x[i] + b@[acc[i], acc[i+1]]
#       resp_u[s:, i + 1] = a[0][0]*resp_u[s:, i] + a[0][1]*resp_v[s:, i] + b[0][0]*acc[i] + b[0][1]*acc[i + 1]
#       resp_v[s:, i + 1] = a[1][0]*resp_u[s:, i] + a[1][1]*resp_v[s:, i] + b[1][0]*acc[i] + b[1][1]*acc[i + 1]
    return x.T
    w2 = w ** 2

    if s:
        sdof_acc = np.zeros_like(resp_u, dtype=float)
        sdof_acc[s:] = -2 * zeta * w[:, np.newaxis] * resp_v[s:] - w2[:, np.newaxis] * resp_u[s:]
        sdof_acc[0] = acc
    else:
        sdof_acc = -2 * zeta * w[:, np.newaxis] * resp_v[s:] - w2[:, np.newaxis] * resp_u[s:]

    return resp_u, resp_v, sdof_acc


def tma1(m,c,k,dt,b2,p,d0,v0,nsteps,lvs,dtf):
    # **************************************************************************
    #                            ***** TMA1 *****
    #                   Prepared by D.Bernal - February, 1999.
    # **************************************************************************

    # Program computes the response of linear MDOF systems using the transition
    # matrix approach. The solution is exact for the assumed variation of the
    # load within the step. Two load variations are possible
    # 1) Linear (default)
    # 2) Constant (Requires setting flag ='cst')
    # ***************************************************************************

    # INPUT
    # m,c,k = mass dampig and stiffness matrices for the system.
    # dt = time step
    # b2 = load influence matrix. Column j contains the spatial distribution of input j.
    # p = inputs. Row j contains the time history of input j.
    # d0,v0 = initial conditions
    # nsteps = number of time steps in the solution.  
    # lvs -- if lvs ='' the solution is carried out assuming a linear variation of
    # the load within the step. Set lv ='cst' for constant load within the step.

    # OUTPUT
    # u = matrix of displacements
    # ud = matrix of velocities
    # udd = matrix of accelerations (if flag =7)
    # (each column contains one DOF)

    # Make sure that the length of specified loading is adequate
    dum, ns = p.shape

    if nsteps*dt>ns*dtf:
       raise ValueError('nsteps*dt is larger than the time over which the loading is defined')


    # Interpolate or extrapolate the applied load to match the integration time step.
    if dtf/dt != 1:
      x = np.arange(0,(ns-1)*dtf, dtf)
      xi = np.arange(0, nsteps*dt, dt)
      F = interp1(x,p.T,xi,'spline').T
    else:
      F = p

    p = F

    # Basic matrices.
    dof,dum = m.shape
    A = [[ eye(dof)*0,  eye(dof)],
         [-inv(m)*k  , -inv(m)*c]]

    B = [[eye(dof)*0,  eye(dof)*0],
         [eye(dof)*0,   inv(m)   ]]

    y = _transition(A,B, p, dt, [d0, v0])

    u  = y[:dof,:].T
    u  = [d0.T, u]
    ud = y[dof:2*dof,:].T
    ud = [v0.T,ud]

    # Calculate the acceleration vector if requested.
    if True:
      X = np.linalg.inv(m)
      for i in range(nsteps+1):
        uddi = X*(b2*p[:,i] - c*ud[i,:].T - k*u[i,:].T)
        udd = [udd, uddi]

      udd = udd.T

    return u,ud,udd

def _transition(A,B, b2, p, dt, y0, lvs=None):
  # Compute the matrices used in the marching algorithm.
  AA = A*dt
  Ad = scipy.linalg.expm(AA)
  X  = scipy.linalg.inv(A)
  P1 = X*(Ad-np.eye(2*dof))*B
  P2 = X*(B - P1/dt)

  # Set P2 to zero if the load is constant within the step
  if lvs =='cst':
    P2 = P2*0

  D0 = np.zeros((dof,1))

# # Perform the integration.
# for i in range(nsteps):
#   y1[:,i] = (P1+P2)*[D0;b2*p[:,i]] - P2*[D0;b2*p[:,i+1]] + Ad@y0
#   y0 = y1[:,i]

  return y1

