import numpy as np
import sys
from sys import argv
from lightcurve import *

def main():
    # --------------------------- #
    # ----- Setting (begin) ----- #
    # --------------------------- #
    if not len(sys.argv)==4:
        print('Error: Number of inputs is wrong(0).(must be {1})'\
            .format(len(sys.argv), 4))
        sys.exit()
    # Filename of input light curve
    name_infits=argv[1]
    # Filename of output light curve
    name_outfits=argv[2]
    # Sampling interval [s]
    dt_rebin=float(argv[3])
    # --------------------------- #
    # ----- Setting (end)   ----- #
    # --------------------------- #

    # ------------------------------------ #
    # ----- Read light curve (begin) ----- #
    # ------------------------------------ #
    lc=LightCurve()
    lc.read_data(name_file=name_infits)
    # ---------------------------------- #
    # ----- Read light curve (end) ----- #
    # ---------------------------------- #

    # ------------------------------------- #
    # ----- Print information (begin) ----- #
    # ------------------------------------- #
    lc.print_obs_info()
    # ------------------------------------- #
    # ----- Print information (end)   ----- #
    # ------------------------------------- #

    # ------------------------------------- #
    # ----- Rebin light curve (begin) ----- #
    # ------------------------------------- #
    lc.rebin_lightcurve(dt_rebin=dt_rebin)
    # ------------------------------------- #
    # ----- Rebin light curve (end)   ----- #
    # ------------------------------------- #

    # ------------------------------------- #
    # ----- Write light curve (begin) ----- #
    # ------------------------------------- #
    lc.write_lightcurve(name_file=name_outfits)
    # ------------------------------------- #
    # ----- Write light curve (end)   ----- #
    # ------------------------------------- #

    print('Light curve was rebinned successfully.')
    print('Results are stored in {0}.'.format(name_outfits))

if __name__=='__main__':
    main()
