import numpy as np
import sys
from sys import argv
from lightcurve import *

def main():
    # --------------------------- #
    # ----- Setting (begin) ----- #
    # --------------------------- #
    if not len(sys.argv)==3:
        print('Error: Number of inputs is wrong(0).(must be {1})'\
            .format(len(sys.argv), 3))
        sys.exit()
    # Filename of input light curve
    name_inlc=argv[1]
    # Filename of output light curve
    name_outlc=argv[2]
    # --------------------------- #
    # ----- Setting (end)   ----- #
    # --------------------------- #

    # ------------------------------------ #
    # ----- Read light curve (begin) ----- #
    # ------------------------------------ #
    lc=LightCurve()
    lc.read_data(name_file=name_inlc)
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
    lc.convert_counts2rate()
    # ------------------------------------- #
    # ----- Rebin light curve (end)   ----- #
    # ------------------------------------- #

    # ------------------------------------- #
    # ----- Write light curve (begin) ----- #
    # ------------------------------------- #
    lc.write_lightcurve(name_file=name_outlc)
    # ------------------------------------- #
    # ----- Write light curve (end)   ----- #
    # ------------------------------------- #

    print('COUNTS was converted into RATE successfully.')
    print('Results are stored in {0}.'.format(name_outlc))

if __name__=='__main__':
    main()
