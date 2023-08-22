import matplotlib.pyplot as plt

try:
    import scienceplots
    plt.style.use(['science','ieee'])
except:
    pass


import numpy as np
from sdof import spectrum


def parse_args(args):
    options = {
        "include_columns": ("periods",  "pseudo_accel"),
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

    accel  = np.loadtxt("data/elCentro.txt")
    dt     = 0.02

    # Response Spectra mode

    Sd, Sv, Sa = spectrum(accel, dt, [0.02, 0.05, 0.1], None)
    stride = 3
    output = np.zeros((len(Sd[0,:]),stride*(len(options["include_columns"])-1)+1))

    output[:,0] = Sd[0,:]

    for i, key in enumerate(options["include_columns"][1:]):
        resp   = COLUMNS[key](Sd, Sv, Sa)
        # stride = resp.shape[1] if len(resp.shape) == 2 else 1

        output[:,i*stride+1:(i+1)*stride+1] = resp #(Sd[1:,:]*(2*np.pi/Sd[0,:])**2).T


    fig, ax = plt.subplots()
    for y in output.T[1:]:
        ax.plot(Sd[0,:], y, linewidth=0.2)
    ax.axis("off")
#   fig.savefig("spectrum.svg")
    plt.show()




