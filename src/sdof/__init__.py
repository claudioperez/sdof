import pathlib
import ctypes
from ctypes import c_double, c_int, c_bool, CFUNCTYPE, POINTER

class SDOF_Peaks(ctypes.Structure):
    _fields_ = [
        ("max_displ",      c_double),
        ("max_accel",      c_double),
        ("time_max_accel", c_double)
    ]

class SDOF_Config(ctypes.Structure):
    _fields_ = [
        ("alpha_m",      c_double),
        ("alpha_f",      c_double),
        ("beta",         c_double),
        ("gamma",        c_double)
    ]


try:
    lib = ctypes.cdll.LoadLibrary(pathlib.Path(__file__).parents[0]/"_fsdof.so")
    conf = lib.CONF
    _fsdof_peaks = lib.fsdof_peaks
    _fsdof_peaks.restype = c_int
    _fsdof_peaks.argtypes = (
        POINTER(SDOF_Config),
        c_double,  c_double,  c_double,
        c_double, c_int, POINTER(c_double), c_double,
        POINTER(SDOF_Peaks)
    )

    CONFIG = SDOF_Config()

except:
    raise


def peaks(m,c,k, dt,f):
    response = SDOF_Peaks()
    _fsdof_peaks(CONFIG, m, c, k, 1.0, len(f), f.ctypes.data_as(POINTER(c_double)), dt, response)
    return response



# import numpy as np
# r = fsdof(1.,1.,1., 0.01, np.sin(np.linspace(0, 5, 100)))
# print(r.max_displ)

def install_me(install_opt=None):
    import os
    import sys
    import subprocess
    if install_opt == "dependencies":
        subprocess.check_call([
            sys.executable, "-m", "pip", "install", *REQUIREMENTS.strip().split("\n")
        ])
        sys.exit()
    try:
        from setuptools import setup
    except ImportError:
        from distutils.core import setup

#   sys.argv = sys.argv[:1] + sys.argv[2:]

    setup(name = "sdof",
          version = "0.0.1",
          description = "Fast SDOF Solver",
          long_description = "",
          author = "",
          author_email = "",
          url = "",
          py_modules = ["sdof"],
          scripts = ["sdof.py"],
          license = "",
          install_requires = [] #*REQUIREMENTS.strip().split("\n")],
    )

if __name__ == "__main__":
    install_me()




