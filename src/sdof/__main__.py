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

def parse_args(args):
    options = {
        "include_columns": ("periods",  "accel"), #"pseudo_accel",
        "return_mode": "spectra"
    }
    return options

COLUMNS = {
    "periods":      lambda Sd, *_:  Sd[0,:][:,None],
    "pseudo_accel": lambda Sd, *_: (Sd[1:,:]*(2*np.pi/Sd[0,:])**2).T,
    "accel":        lambda Sd, Sv, Sa: Sa[1:,:].T
}


if __name__ == "__main__":
    import sys
    options = parse_args(sys.argv)

    accel  = np.loadtxt(sys.argv[1]) if len(sys.argv) > 1 else np.loadtxt(sys.stdin)
    dt     = float(sys.argv[-1])

    # Response Spectra mode

    Sd, Sv, Sa = spectrum(accel, dt, [0.02, 0.05], None)
    stride = 2
    output = np.zeros((len(Sd[0,:]),stride*(len(options["include_columns"])-1)+1))

    output[:,0] = Sd[0,:]

    for i, key in enumerate(options["include_columns"][1:]):
        resp   = COLUMNS[key](Sd, Sv, Sa)
        # stride = resp.shape[1] if len(resp.shape) == 2 else 1

        output[:,i*stride+1:(i+1)*stride+1] = resp #(Sd[1:,:]*(2*np.pi/Sd[0,:])**2).T


    np.savetxt(sys.stdout.buffer, output)


