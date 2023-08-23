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

import numpy as np
from sdof import spectrum

def parse_args(args):
    options = {
        "include_columns": ("periods",  "accel"), #"pseudo_accel",
        "return_mode": "spectra",
        "time_step": None,
        "num_threads": 4,
        "file": None
    }
    argi = iter(args)
    for arg in argi:
        if arg[:2] == "-j":
            options["num_threads"] = int(next(argi))

        elif options["file"] is None:
            options["file"] = arg

        elif options["time_step"] is None:
            options["time_step"] = float(arg)

    if options["file"] is None:
        options["file"] = sys.stdin

    return options


COLUMNS = {
    "periods":      lambda Sd, *_:  Sd[0,:][:,None],
    "pseudo_accel": lambda Sd, *_: (Sd[1:,:]*(2*np.pi/Sd[0,:])**2).T,
    "accel":        lambda Sd, Sv, Sa: Sa[1:,:].T
}


if __name__ == "__main__":
    import sys
    options = parse_args(sys.argv[1:])

    accel  = np.loadtxt(options["file"])
    dt     = options["time_step"]

    # Response Spectra mode

    Sd, Sv, Sa = spectrum(accel, dt, [0.02], None, threads=options["num_threads"])
    stride = 1
    output = np.zeros((len(Sd[0,:]),stride*(len(options["include_columns"])-1)+1))

    output[:,0] = Sd[0,:]

    for i, key in enumerate(options["include_columns"][1:]):
        resp   = COLUMNS[key](Sd, Sv, Sa)
        # stride = resp.shape[1] if len(resp.shape) == 2 else 1

        output[:,i*stride+1:(i+1)*stride+1] = resp #(Sd[1:,:]*(2*np.pi/Sd[0,:])**2).T


    np.savetxt(sys.stdout.buffer, output)

