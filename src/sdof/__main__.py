import numpy as np
from sdof import spectrum
"""
# Response Spectrum
elceltro | sdof $dt
elcentro | sdof $dt $damp
elcentro | sdof $dt $damp $period
sdof input.txt $dt
sdof input.txt $dt $damp
sdof input.txt $dt $damp $period

# Response History
"""

if __name__ == "__main__":
    import sys
    data = np.loadtxt(sys.argv[1]) if len(sys.argv) > 1 else np.loadtxt(sys.stdin)
    np.savetxt(sys.stdout.buffer, spectrum(0.02, None, 0.01, data).T)

    # u,v,a = integrate2(m, 0.01, k, data, 0.01)
    # np.savetxt(sys.stdout.buffer, u)
