"""
See also:
https://www.mathworks.com/matlabcentral/fileexchange/32182-wilson-tita-time-integration-method
"""
import numpy as np

def _force_function(f, dt, theta=1)->callable:

    if callable(f):
        return f 

    def f(t):
        pass


def wilson_gem(M, C, K, F, dt, t_end, theta=1.4, u0=0, v0=0):
    """
    Implements the Wilson-Theta method for structural dynamics.

    Args:
        F: Force vector (function of time)
        u0: Initial displacement vector
        v0: Initial velocity vector
        dt: Time step
        t_end: End time
        theta: Wilson-Theta parameter (default=1.4)

    Returns:
        t: Time array
        u: Displacement array
        v: Velocity array
        a: Acceleration array
    """

    # Time array
    t  = np.arange(0, t_end + dt, dt)
    nt = len(t)

    # Initialize displacement, velocity, and acceleration arrays
    u = np.zeros((nt, len(u0)))
    v = np.zeros((nt, len(v0)))
    a = np.zeros((nt, len(v0)))

    # Set initial conditions
    u[0] = u0
    v[0] = v0
    a[0] = np.linalg.solve(M, F(0) - C @ v0 - K @ u0)

    # Effective stiffness matrix
    K_eff = K + (theta / (dt ** 2)) * M + (theta / dt) * C

    for i in range(len(t) - 1):
        # Calculate effective force vector
        F_eff = (
            F(t[i] + theta * dt)
            + theta * (F(t[i + 1]) - F(t[i]))
            + M @ (
                (theta ** 2 / (dt ** 2)) * u[i]
                + (theta / dt) * v[i]
                + (theta - 1) * a[i]
            )
            + C @ ((theta / dt) * u[i] + (theta - 1) * v[i] + (theta / 2 - 1) * dt * a[i])
        )

        # Solve for displacement increment
        du = np.linalg.solve(K_eff, F_eff)

        # Calculate displacement, velocity, and acceleration at theta * dt
        u_theta = u[i] + theta * du
        v_theta = v[i] + (theta * dt) * a[i] + (theta ** 2 / (2 * dt)) * du
        a_theta = (
            (theta / (dt ** 2)) * du - (theta / dt) * v[i] - ((theta - 2) / 2) * a[i]
        )

        # Calculate displacement, velocity, and acceleration at t[i+1]
        u[i + 1] = u[i] + du
        v[i + 1] = v[i] + dt * ((1 - theta) * a[i] + theta * a_theta)
        a[i + 1] = (1 / (theta * dt ** 2)) * du - (1 / (theta * dt)) * v[i] - ((1 - theta) / (2 * theta)) * a[i]

    return u, v, a


def wilson_gpt(model, f, dt, theta, t_end,  u0=0, v0=0, a0=None):
    """
    Wilson-theta method for SDOF system.

    Parameters:
    m : float : Mass of the system.
    c : float : Damping coefficient.
    k : float : Stiffness coefficient.
    f : function : External force as a function of time (f(t)).
    u0 : float : Initial displacement.
    v0 : float : Initial velocity.
    a0 : float : Initial acceleration.
    dt : float : Time step.
    theta : float : Wilson-theta parameter (typically 1.4 for stability).
    t_end : float : End time.

    Returns:
    t : array : Time points.
    u : array : Displacements at each time step.
    v : array : Velocities at each time step.
    a : array : Accelerations at each time step.
    """
    k, c, m = model

    # Time steps and initialize arrays
    n_steps = int(t_end / dt) + 1
    t = np.linspace(0, t_end, n_steps)
    u = np.zeros(n_steps)
    v = np.zeros(n_steps)
    a = np.zeros(n_steps)

    # Initial conditions
    u[0], v[0], a[0] = u0, v0, a0

    # Effective time step and constants
    dt_eff = theta * dt
    k_eff = m / (dt_eff ** 2) + c / (2 * dt_eff)
    a_eff = m / (dt_eff ** 2) - c / (2 * dt_eff)
    b_eff = c / dt_eff - k

    # Time-stepping loop
    for i in range(1, n_steps):
        # Effective force
        p_eff = f(t[i]) + a_eff * u[i-1] + b_eff * v[i-1] + m * a[i-1]

        # Solve for displacement
        u_theta = p_eff / k_eff
        u[i] = u[i-1] + theta * (u_theta - u[i-1])

        # Calculate velocity and acceleration
        v[i] = (u[i] - u[i-1]) / dt_eff
        a[i] = (v[i] - v[i-1]) / dt

    return t, u, v, a

# Define external force function (e.g., sinusoidal load)
def external_force(t):
    return np.sin(2 * np.pi * t)

# Run Wilson-theta method
u, v, a = wilson_gpt((k, c, m), external_force, dt=dt, theta=theta, t_end=5)


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
