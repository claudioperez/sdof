#
# BSD 2-Clause License
#
# Copyright (c) 2022-2023, Claudio M. Perez
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
import os
import pathlib
import ctypes
from ctypes import c_double, c_int, c_bool, CFUNCTYPE, POINTER
import numpy as np
from numpy.ctypeslib import ndpointer

if os.name == "nt":
    so_ext = ".pyd"

else:
    import   distutils.ccompiler
    so_ext = distutils.ccompiler.new_compiler().shared_lib_extension

class sdof_peaks_t(ctypes.Structure):
    _fields_ = [
        ("max_displ",      c_double),
        ("max_veloc",      c_double),
        ("max_accel",      c_double)
    ]

class _sdof_config(ctypes.Structure):
    _fields_ = [
        ("alpha_m",      c_double),
        ("alpha_f",      c_double),
        ("beta",         c_double),
        ("gamma",        c_double)
    ]

# TODO
libfile = str(next(pathlib.Path(__file__).parents[0].glob("_integrate.*"+so_ext)))
lib = ctypes.cdll.LoadLibrary(libfile)

_sdof_integrate_peaks = lib.sdof_integrate_peaks
_sdof_integrate_peaks.restype = c_int
_sdof_integrate_peaks.argtypes = (
    POINTER(_sdof_config),
    c_double,  c_double,  c_double,
    c_double, c_int, POINTER(c_double), c_double,
    POINTER(sdof_peaks_t)
)

_sdof_integrate_peaks_2 = lib.sdof_integrate_peaks_2
_sdof_integrate_peaks_2.restype  = c_int
_sdof_integrate_peaks_2.argtypes = (
    POINTER(_sdof_config),
    c_double,  c_double,  c_double,
    c_double, c_int, POINTER(c_double), c_double,
    POINTER(sdof_peaks_t)
)

_sdof_integrate = lib.sdof_integrate
_sdof_integrate.restype  = c_int
_sdof_integrate.argtypes = (
    POINTER(_sdof_config),
    c_double,  c_double,  c_double,
    c_double, c_int, POINTER(c_double), c_double,
    ndpointer(c_double, flags="C_CONTIGUOUS")
)

_sdof_integrate_plastic = lib.sdof_integrate_plastic
_sdof_integrate_plastic.restype  = c_int
_sdof_integrate_plastic.argtypes = (
    POINTER(_sdof_config),
    c_double,  c_double,  c_double,
    c_double, c_int, POINTER(c_double), c_double,
    ndpointer(c_double, flags="C_CONTIGUOUS")
)

_sdof_integrate_unrolled = lib.sdof_integrate_unrolled
_sdof_integrate_unrolled.restype  = c_int
_sdof_integrate_unrolled.argtypes = (
    POINTER(_sdof_config),
    c_double,  c_double,  c_double,
    c_double, c_int, POINTER(c_double), c_double,
    ndpointer(c_double, flags="C_CONTIGUOUS")
)

# _sdof_integrate_0 = lib.sdof_integrate_0
# _sdof_integrate_0.restype = c_int
# _sdof_integrate_0.argtypes = (
#     POINTER(_sdof_config),
#     c_double,  c_double,  c_double,
#     c_double, c_int, POINTER(c_double), c_double,
#     ndpointer(c_double, flags="C_CONTIGUOUS")
# )

CONFIG = _sdof_config()


try:
    # Try importing the threaded _spectrum library

    # TODO
    libfile = str(next(pathlib.Path(__file__).parents[0].glob("_spectrum*"+so_ext)))
    lib = ctypes.cdll.LoadLibrary(libfile)
    _sdof_spectrum = lib.sdof_spectrum
    _sdof_spectrum.restype  = c_int
    _sdof_spectrum.argtypes = (
        POINTER(_sdof_config),
        POINTER(c_double), c_int, c_double,
        c_double,  c_double,  c_int, # p_min, p_max, np
        c_double,                    # damping
        c_int,
        POINTER(sdof_peaks_t)
    )

except Exception as e:
    import warnings
    warnings.warn(f"Failed to load native library: {e}")


def integrate_0(f, dt,
              k      : float,
              c      : float = 0.0,
              m      : float = 1.0,
              u0=0.0, v0=0.0,
              out  =  None,
              alpha_m: float = 1.0,
              alpha_f: float = 1.0,
              beta   : float = 0.25,
              gamma  : float = 0.50
    ):

    if out is None:
        output = np.empty((3,len(f)))
    else:
        output = out

    output[:2,0] = u0, v0

    config = _sdof_config(
                alpha_m = alpha_m,
                alpha_f = alpha_f,
                beta    = beta,
                gamma   = gamma
    )

    _sdof_integrate_0(config, m, c, k, 1.0, len(f), np.asarray(f).ctypes.data_as(POINTER(c_double)), dt, output)
    return output

def integrate(
              f,
              dt: float,
              k:  float,
              c:  float,
              m:  float,
              u0: float = 0.0,
              v0: float = 0.0,
              out  =  None,
              fy: float = None,
              beta   : float = 0.25,
              gamma  : float = 0.50,
              alpha_m: float = 1.0,
              alpha_f: float = 1.0
    ):
    r"""
    This function integrates scalar differential equations of the form

    .. math::

       m \ddot{u} + c \dot{u} + k u = f(t)

    for constant coefficients :math:`m`, :math:`c` and :math:`k`.

    Parameters:
          f       (list)  : Loading
          dt      (float) : Time step
          k       (float) : Stiffness
          c       (float) : Damping
          m       (float) : Mass
          u0      (float, optional): Initial displacement
          v0      (float, optional): Initial velocity
          out     (array, optional): Array to store outputs
          alpha_m (float, optional): Integrator :math:`\alpha_m` parameter
          alpha_f (float, optional): Integrator :math:`\alpha_f` parameter
          beta    (float, optional): Newmark :math:`\beta` parameter
          gamma   (float, optional): Newmark :math:`\gamma` parameter


    Integration is carried out using a Generalized - :math:`\alpha`
    scheme by calling the C function :c:func:`sdof_integrate`.
    """
    if out is None:
        output = np.empty((len(f),3))
    else:
        output = out

    # Store initial velocity and displacement
    output[0,:2] = u0, v0

    config = _sdof_config(
                alpha_m = alpha_m,
                alpha_f = alpha_f,
                beta    = beta,
                gamma   = gamma
    )
    if fy is not None:
        _sdof_integrate_plastic(config, m, c, k, 1.0, len(f),
                                np.asarray(f).ctypes.data_as(POINTER(c_double)), dt, output)
    else:
        _sdof_integrate_unrolled(config, m, c, k, 1.0, len(f),
                                 np.asarray(f).ctypes.data_as(POINTER(c_double)), dt, output)
    return output.T



def _spectrum_cthreads(response, accel, dt, damping, periods: tuple, config: _sdof_config, n_threads=8):
    Sd, Sv, Sa = response
    output = (sdof_peaks_t * periods[2])()
    for i,damp in enumerate(damping):
        _sdof_spectrum(config, np.asarray(accel).ctypes.data_as(POINTER(c_double)), len(accel), dt,
                        periods[0], periods[1], periods[2],
                        damp, n_threads, output)

        for j in range(periods[2]):
            Sd[i,j] = output[j].max_displ
            Sv[i,j] = output[j].max_veloc
            Sa[i,j] = output[j].max_accel

def _spectrum_pythreads(n_threads=1):
    import threading
    threads = []
    for i in range(n_threads):
        threads.append(
                threading.Thread(
                    target=some_ctypes_func, args=(arg, )
                )
        )

    for thread in threads:
        thread.start()

    for thread in threads:
        thread.join()


def spectrum(accel, dt, damping, periods=None, interp=None, threads:int=None, **kwds):
    """
    This function is a wrapper around the C function :c:func:`sdof_spectrum`.

    Parameters:
        accel (list), dt (float): Excitation array and time step
    """

    # Setup default parameters

    if isinstance(damping, float):
        damping = [damping]

    if periods is None:
        periods = np.linspace(0.02, 3.0, 200)
    elif isinstance(periods, tuple):
        periods = np.linspace(*periods)
    elif isinstance(periods, list):
        periods = np.array(periods)

    if interp is None:
        from scipy.interpolate import interp1d as interp

    config = _sdof_config(
                alpha_m = kwds.get("alpha_m", 1.00),
                alpha_f = kwds.get("alpha_f", 1.00),
                beta    = kwds.get("beta",    0.25),
                gamma   = kwds.get("gamma",   0.50)
    )

    # Initialize some memory

    Sd, Sv, Sa = np.zeros((3,1+len(damping), len(periods)))
    Sd[0,:] = periods[:]
    Sv[0,:] = periods[:]
    Sa[0,:] = periods[:]

    if threads is not None:
        Sd[0,:] = Sv[0,:] = Sa[0,:] = np.linspace(periods[0], periods[-1], len(periods))
        _spectrum_cthreads((Sd[1:,:], Sv[1:,:], Sa[1:,:]), accel, dt, damping,
                           periods=(periods[0], periods[-1], len(periods)),
                           config=config, n_threads=threads)
        return Sd,Sv,Sa

    pi = np.pi
    mass = 1.0
    numdata = len(accel)
    t       = np.arange(0,(numdata)*dt,dt)
    t_max   = max(t)

    u0, v0 = 0.0, 0.0

    for i,dmp in enumerate(damping):
        for j,period in enumerate(periods):
            if dt/period > 0.02:
                dtp  = period*0.02
                dtpx = np.arange(0,t_max,dtp)
                dtpx = dtpx
                accfrni = interp(t, accel)(dtpx)
                accfrn  = accfrni[1:len(accfrni)-1]
                numdatan = len(accfrn)
            else:
                dtp = dt
                accfrn = accel
                numdatan = numdata

            p = -mass*accfrn
            k = 4*pi**2*mass/period**2
            c = 2*dmp*np.sqrt(k/mass)

            d, v, a = integrate(mass, c, k, p, dt)

            # resp = peaks(m, c, k, p, dt)
            # resp.max_displ
            # resp.max_accel

            Sd[1+i,j] = abs(max(d, key=abs))
            Sv[1+i,j] = abs(max(v, key=abs))
            Sa[1+i,j] = abs(max(a+accfrn, key=abs))

    return Sd,Sv,Sa


def peaks(m,c,k, f, dt):
    response = sdof_peaks_t()
    config = _sdof_config(
                alpha_m = kwds.get("alpha_m", 1.00),
                alpha_f = kwds.get("alpha_f", 1.00),
                beta    = kwds.get("beta",    0.25),
                gamma   = kwds.get("gamma",   0.50)
    )
    _sdof_integrate_peaks(config, m, c, k, 1.0, len(f),
                          f.ctypes.data_as(POINTER(c_double)), dt, response)
    return response

