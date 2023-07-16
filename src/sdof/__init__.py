import os
import pathlib
import ctypes
from ctypes import c_double, c_int, c_bool, CFUNCTYPE, POINTER
import numpy as np
from numpy.ctypeslib import ndpointer

if os.name == 'nt':
    so_ext = ".pyd"
else:
    import distutils.ccompiler
    so_ext = distutils.ccompiler.new_compiler().shared_lib_extension

class _sdof_peaks(ctypes.Structure):
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


try:
    libfile = str(next(pathlib.Path(__file__).parents[0].glob("_fsdof.*"+so_ext)))
    lib = ctypes.cdll.LoadLibrary(libfile)
    # conf = lib.CONF
    _fsdof_peaks = lib.fsdof_peaks
    _fsdof_peaks.restype = c_int
    _fsdof_peaks.argtypes = (
        POINTER(_sdof_config),
        c_double,  c_double,  c_double,
        c_double, c_int, POINTER(c_double), c_double,
        POINTER(_sdof_peaks)
    )

    _fsdof_integrate = lib.fsdof_integrate
    _fsdof_integrate.restype = c_int
    _fsdof_integrate.argtypes = (
        POINTER(_sdof_config),
        c_double,  c_double,  c_double,
        c_double, c_int, POINTER(c_double), c_double,
        ndpointer(c_double, flags="C_CONTIGUOUS")
    )

    _fsdof_integrate2 = lib.fsdof_integrate2
    _fsdof_integrate2.restype = c_int
    _fsdof_integrate2.argtypes = (
        POINTER(_sdof_config),
        c_double,  c_double,  c_double,
        c_double, c_int, POINTER(c_double), c_double,
        ndpointer(c_double, flags="C_CONTIGUOUS")
    )

    CONFIG = _sdof_config()

except:
    raise

#   elastic_sdof()
#   plastic


def integrate(m,c,k,f,dt, u0=0.0, v0=0.0,
              out  =  None,
              alpha_m: float = 1.0,
              alpha_f: float = 1.0,
              beta   : float = 0.25,
              gamma  : float = 0.50
    ):

    import numpy as np
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

    _fsdof_integrate(config, m, c, k, 1.0, len(f), np.asarray(f).ctypes.data_as(POINTER(c_double)), dt, output)
    return output

def integrate2(m,c,k,f,dt, u0=0.0, v0=0.0,
              out  =  None,
              alpha_m: float = 1.0,
              alpha_f: float = 1.0,
              beta   : float = 0.25,
              gamma  : float = 0.50
    ):
    import numpy as np
    if out is None:
        output = np.empty((3,len(f)))
    else:
        output = out
    output[0,:2] = u0, v0

    config = _sdof_config(
                alpha_m=alpha_m,
                alpha_f=alpha_f,
                beta   =beta,
                gamma  =gamma
    )

    _fsdof_integrate2(CONFIG, m, c, k, 1.0, len(f), np.asarray(f).ctypes.data_as(POINTER(c_double)), dt, output)
    return output.T


def _thread_spectrum(n_threads=1):
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

def spectrum(accel, dt, damping, periods=None, interp=None):
    if interp is None:
        from scipy.interpolate import interp1d as interp

    if isinstance(damping, float):
        damping = [damping]

    if periods is None:
        periods = np.linspace(0.02, 3.0, 200)
    elif isinstance(periods, tuple):
        periods = np.linspace(*periods)
    elif isinstance(periods, list):
        periods = np.array(periods)


    pi = np.pi
    mass = 1.0
    numdata = len(accel)
    t       = np.arange(0,(numdata)*dt,dt)
    t_max   = max(t)

    u0, v0 = 0.0, 0.0

    Sd, Sv, Sa = np.zeros((3,1+len(damping), len(periods)))
    Sd[0,:] = periods[:]
    Sv[0,:] = periods[:]
    Sa[0,:] = periods[:]

    for i,dmp in enumerate(damping):
        for j,period in enumerate(periods):
            if dt/period > 0.02:
                dtp = period*0.02
                dtpx = np.arange(0,t_max,dtp)
                dtpx = dtpx
                accfrni = interp(t, accel)(dtpx)
                accfrn = accfrni[1:len(accfrni)-1]
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
    response = _sdof_peaks()
    _fsdof_peaks(CONFIG, m, c, k, 1.0, len(f), f.ctypes.data_as(POINTER(c_double)), dt, response)
    return response

