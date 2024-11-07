import numpy as np 

def _rk4(m, c, k, u0=1, v0=0, n=8, dt=0.05):
    """Runge-Kutta solution of underdamped system.

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
    >>> vtb._rk4(m=1, c=.1, k=1, u0=1, v0=0, n=8, dt=0.05)
    (array([0.  , 0.05, 0.1 , 0.15, 0.2 , 0.25, 0.3 , 0.35, 0.4 ]),
    array([[ 1.,  0.],
           [ 1.  , -0.05],
           [ 1.  , -0.1 ],
           [ 0.99, -0.15],
           [ 0.98, -0.2 ],
           [ 0.97, -0.24],
           [ 0.96, -0.29],
           [ 0.94, -0.34],
           [ 0.92, -0.38]]))

    """
    x = np.zeros((n + 1, 2))
    x[0, :] = u0, v0
    A = np.array([[0, 1],
                  [-k / m, -c / m]])

    def f(x_): return A@x_

    for i in range(n):
        k1 = dt * f(x[i])
        k2 = dt * f(x[i] + k1 / 2)
        k3 = dt * f(x[i] + k2 / 2)
        k4 = dt * f(x[i] + k3)
        x[i + 1] = x[i] + (k1 + 2.0 * (k2 + k3) + k4) / 6.0

    return x

