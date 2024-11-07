"""
Deconvolution
"""
import numpy as np


def time_from_circ(freq):
    pass 

def freq_from_time(time):
    pass

def fourier(m, c, k, f, dt, n=None):
    """
    Integration in the frequency domain
    """
    from numpy.fft import fft, fftfreq, ifft

    if n is None:
        n = len(f)
        f_pad = f
    else:
        assert n > len(f)
        f_pad = np.zeros(n)
        f_pad[:len(f)] = f

    F  = fft(f_pad)

    # Convert frequency (Hertz) to circular frequency
    w  = 2.0*np.pi*fftfreq(n, dt)

    # Conpute array of transfer function values H(w)
    Hw = 1/(k - m*w**2 + c*w*1j)

    # Element-wise multiplication of two 1D arrays
    U  = Hw*F

    return (
           np.real(ifft(U))[:len(f)]/dt,
           np.real(ifft((0+w*1j)*U))[:len(f)],
           np.real(ifft((0+w*1j)**2*U))[:len(f)]
    )



def fourier_series(dat, t, n):
    """Fourier series approximation to a function.

    returns Fourier coefficients of a function.
    The coefficients are numerical approximations of the true
    coefficients.

    Parameters
    ----------
    dat: array
        Array of data representing the function.
    t: array
        Corresponding time array.
    n: int
        The desired number of terms to use in the
        Fourier series.

    Returns
    -------
    a, b: tuple
        Tuple containing arrays with the Fourier coefficients
        as float arrays.
        The function also produces a plot of the approximation.

    Examples
    --------
    >>> import vibration_toolbox as vtb
    >>> f = np.hstack((np.arange(-1, 1, .04), np.arange(1, -1, -.04)))
    >>> f += 1
    >>> t = np.arange(0, len(f))/len(f)
    >>> a, b = vtb.fourier_series(f, t, 5)

    """
    len_ = len(dat) / 2
    fs = (np.fft.fft(dat)) / len_
    a0 = fs[0]
    a = np.real(np.hstack((a0, fs[1:len(fs / 2)])))
    b = -np.imag(fs[1:len(fs / 2)])
    len_ *= 2
    dt = 2 * np.pi / len_
    tp = np.arange(0, 2 * np.pi, dt)
    dataapprox = a[0] / 2 + np.zeros_like(dat)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(t, dat)

    for i in range(1, n):
        newdat = a[i] * np.cos(tp * i) + b[i] * np.sin(tp * i)
        dataapprox += newdat
        if i == n - 1:
            ax1.plot(t, newdat)

    ax1.plot(t, dataapprox)

    return a, b


def fourier_approximation(a0, aodd, aeven, bodd, beven, N, T):
    r"""Plot the Fourier series defined by coefficient falues.

    Coefficients are defined by Inman [1]_.

    :math:`a_0=\frac{2}{T}\int_0^T F(t)dt`

    :math:`a_n=\frac{2}{T}\int_0^T F(t) \cos(n\omega_T t)dt`

    :math:`b_n=\frac{2}{T}\int_0^T F(t) \sin(n\omega_T t)dt`

    Parameters
    ----------
    a0: float or function
        :math:`a_0`-  Fourier coefficient.
    aodd: float or function
        :math:`a_n`-  Fourier coefficient for n odd.
    aeven: float or function
        :math:`a_n`-  Fourier coefficient for n even.
    bodd: float or function
        :math:`b_n`-  Fourier coefficient for n odd
    beven: float or function
        :math:`b_n`-  Fourier coefficient for n even

    Returns
    -------
    t, F: tuple
        Tuple with time and F(t). It also returns
        a plot with the Fourier approximation.

    References
    ----------
    .. [1] Daniel J. Inman, "Engineering Vibration", 4th ed., Prentice Hall,
           2013.

    Examples
    --------
    >>> # Square wave
    >>> import vibration_toolbox as vtb
    >>> bodd_square = lambda n: -3*(-1+(-1)**n)/n/np.pi
    >>> beven_square = lambda n: -3*(-1+(-1)**n)/n/np.pi
    >>> t, F = vtb.fourier_approximation(-1, 0, 0, bodd_square, beven_square, 20, 2)
    >>> # Triangular wave
    >>> aeven_triangle = lambda n: -8/np.pi**2/n**2
    >>> t, F = vtb.fourier_approximation(0,aeven_triangle,0,0,0,20,10)

    """
    def make_const(value):
        return lambda t: value

    def fourier_contribution(a, b, n):
        return a(n) * np.cos(n * 2 * np.pi * t / T) + \
               b(n) * np.sin(n * 2 * np.pi * t / T)

    a0, aodd, aeven, bodd, beven = (arg if callable(arg)
                                    else make_const(arg)
                                    for arg in (a0, aodd, aeven, bodd, beven))

    dt = min(T / 400, T / 10 * N)
    t = np.arange(0, T * 3, dt)
    F = 0 * t + a0(0) / 2

    # n odd
    for n in range(1, N, 2):
        F = F + fourier_contribution(aodd, bodd, n)
    # n even
    for n in range(2, N, 2):
        F = F + fourier_contribution(aeven, beven, n)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Time, t')
    ax1.set_ylabel('F(t)')
    ax1.plot(t, F)
    return t, F


if __name__ == "__main__":
    from sdof import sdof

    k = 10
    c = 0.05*k
    m = 2

    # Make frequency array, w
    nw  = 500
    dw = 0.2
    wmax = (nw-1)*dw
    w  = np.arange(0, wmax, dw)

    # Make time array, t
    dt = np.pi/wmax
    tmax = 2*np.pi/dw
    t = np.arange(0, tmax, dt)

    print(f"{len(t) = }\n{len(w) = }")

    # Create the SDOF
    model = sdof(k=k, m=m, c=c)

    Hw = np.zeros(2*nw-1, dtype=complex)
    Hw[:nw] = model.transfer(w)

    import numpy.fft as fft
    input = [1] + [0]*(len(t)-1)
    # Note: F is all ones
    F  = fft.fft(input)
    U  = Hw*F
    u = 2.0*np.real(fft.ifft(U))/dt

    print(f"{w = }")

    import matplotlib.pyplot as plt

    _, ax = plt.subplots()
    ax.plot(t, u[:len(t)], 'o')
    ax.plot(t, model.impulse(t)[0], 'x')
    ax.plot(t, fourier(m,c,k,input,dt, 2*len(input)-1)[0], '.')
    # ax.plot(t, model.homogeneous(t, v0=1/m)[0], '.')

#   from integration.transition import exponential
#   ax.plot(t, exponential(np.sin(t), dt, model)[0], '.')
    plt.show()

