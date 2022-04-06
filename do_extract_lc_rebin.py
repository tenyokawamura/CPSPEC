import numpy as np
import sys
from sys import argv
from lightcurve import *
# ---------------------------------------------------------- #
# Extract data which can be compared between LE, ME and HE   #
# due to common sampling times.                              #
# ---------------------------------------------------------- #
def main():
    # --------------------------- #
    # ----- Setting (begin) ----- #
    # --------------------------- #
    if not len(sys.argv)==4:
        print('Error: Number of inputs is wrong(0).(must be {1})'\
            .format(len(sys.argv), 4))
        sys.exit()
    # Filename of input light curve
    name_inlc =argv[1]
    # Filename of input common indices
    name_intxt=argv[2]
    # Filename of output light curve
    name_outlc=argv[3]
    # --------------------------- #
    # ----- Setting (end)   ----- #
    # --------------------------- #

    # --------------------------------------- #
    # ----- Read common indices (begin) ----- #
    # --------------------------------------- #
    with open(name_intxt, 'r') as fin:
        lines=fin.readlines()
    if len(lines)==1:
        print('No common data (skipped).')
        sys.exit()
    is_min, is_max\
        =np.loadtxt(fname=name_intxt, \
            dtype='int', \
            skiprows=1, \
            unpack=True)
    is_min=np.ravel(np.array([is_min]))
    is_max=np.ravel(np.array([is_max]))
    # --------------------------------------- #
    # ----- Read common indices (end)   ----- #
    # --------------------------------------- #

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

    # --------------------------------------- #
    # ----- Extract light curve (begin) ----- #
    # --------------------------------------- #
    if len(lc.data[0])==4:
        first=True
        for i_min, i_max in zip(is_min, is_max):
            times_extract   =lc.data['TIME'   ][i_min:i_max+1]
            rates_extract   =lc.data['RATE'   ][i_min:i_max+1]
            errors_extract  =lc.data['ERROR'  ][i_min:i_max+1]
            fracexps_extract=lc.data['FRACEXP'][i_min:i_max+1]
            if first==True:
                first=False
                lc.times_rebin   =times_extract
                lc.rates_rebin   =rates_extract
                lc.errors_rebin  =errors_extract
                lc.fracexps_rebin=fracexps_extract
            else:
                lc.times_rebin   =np.append(lc.times_rebin,    times_extract)
                lc.rates_rebin   =np.append(lc.rates_rebin,    rates_extract)
                lc.errors_rebin  =np.append(lc.errors_rebin,   errors_extract)
                lc.fracexps_rebin=np.append(lc.fracexps_rebin, fracexps_extract)
    elif len(lc.data[0])==5:
        first=True
        for i_min, i_max in zip(is_min, is_max):
            times_extract   =lc.data['TIME'   ][i_min:i_max+1]
            rates_extract   =lc.data['RATE'   ][i_min:i_max+1]
            errors_extract  =lc.data['ERROR'  ][i_min:i_max+1]
            fracexps_extract=lc.data['FRACEXP'][i_min:i_max+1]
            deadcs_extract  =lc.data['DEADC'  ][i_min:i_max+1]
            if first==True:
                first=False
                lc.times_rebin   =times_extract
                lc.rates_rebin   =rates_extract
                lc.errors_rebin  =errors_extract
                lc.fracexps_rebin=fracexps_extract
                lc.deadcs_rebin  =deadcs_extract
            else:
                lc.times_rebin   =np.append(lc.times_rebin,    times_extract)
                lc.rates_rebin   =np.append(lc.rates_rebin,    rates_extract)
                lc.errors_rebin  =np.append(lc.errors_rebin,   errors_extract)
                lc.fracexps_rebin=np.append(lc.fracexps_rebin, fracexps_extract)
                lc.deadcs_rebin  =np.append(lc.deadcs_rebin,   deadcs_extract)
    # --------------------------------------- #
    # ----- Extract light curve (end)   ----- #
    # --------------------------------------- #

    # ------------------------------------- #
    # ----- Write light curve (begin) ----- #
    # ------------------------------------- #
    lc.write_lightcurve(name_file=name_outlc)
    # ------------------------------------- #
    # ----- Write light curve (end)   ----- #
    # ------------------------------------- #

    print('Light curve was extracted successfully.')
    print('Results are stored in {0}.'.format(name_outlc))

if __name__=='__main__':
    main()
