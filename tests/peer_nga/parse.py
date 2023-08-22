import numpy as np

def parse_spectra(lines):
    "-- Unscaled Horizontal & Vertical Spectra, as recorded --"


def parse(filename):
    print(np.loadtxt(filename, skiprows=0, delimiter=","))


if __name__ == "__main__":
    import sys
    print(parse(sys.argv[1]))
