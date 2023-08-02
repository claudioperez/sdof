import timeit
import numpy as np
import sdof
from sdof import CONFIG, c_double, POINTER


# def integrate():

m, k = 1.0, 10.0
N = 50000
dt = 0.01
f = np.sin(np.linspace(0, 5*np.pi, N)).ctypes.data_as(POINTER(c_double))
C = np.linspace(0.001, 0.10, 1)

output = np.zeros((3,N))

print(timeit.timeit("""
for c in C:
    sdof._sdof_integrate(CONFIG, m, c, k, 1.0, N, f, dt, output)
""",globals=globals(), number=2000))

print(timeit.timeit("""
for c in C:
    sdof._sdof_integrate_unrolled(CONFIG, m, c, k, 1.0, N, f, dt, output)
""",globals=globals(), number=2000))

print(timeit.timeit("""
for c in C:
    sdof._sdof_integrate(CONFIG, m, c, k, 1.0, N, f, dt, output)
""",globals=globals(), number=2000))

if False:
    output = np.zeros((N,3))
    print(timeit.timeit("""
    for c in C:
        sdof._sdof_integrate2(CONFIG, m, c, k, 1.0, N, f, dt, output)
    """,globals=globals(), number=1000))




if False:
    response = sdof.SDOF_Peaks()

    print(timeit.timeit("""
    for c in C:
        sdof._sdof_peaks(CONFIG, m, c, k, 1.0, N, f, dt, response)
    """,globals=globals(), number=1000))


