import numpy as np
import sys
from sys import argv
from lightcurve import *
from fourier import *

def main():
    # --------------------------- #
    # ----- Setting (begin) ----- #
    # --------------------------- #
    if not len(sys.argv)==9:
        print('Error: Number of inputs is wrong(0).(must be {1})'\
            .format(len(sys.argv), 9))
        sys.exit()
    ## Minimum channel
    #ch_min=int(argv[1])
    ## Maximum channel
    #ch_max=int(argv[2])
    # Minimum channel
    ch_min=int(float(argv[1]))
    # Maximum channel
    ch_max=int(float(argv[2]))
    # Number of bins per interval
    n_bin=int(argv[3])
    # Number of intervals
    n_int=int(argv[4])
    # Maximize intervals? (True: maximize, False: stop at n_int)
    maximize=argv[5]
    # Fraction of acceptable absent bins (data gaps)
    frac_gap=float(argv[6])
    # Filename of light curve
    name_infits=argv[7]
    # Filename of fourier transform
    name_outfits=argv[8]
    # --------------------------- #
    # ----- Setting (end)   ----- #
    # --------------------------- #

    # ------------------------------------ #
    # ----- Read light curve (begin) ----- #
    # ------------------------------------ #
    lc=LightCurve()
    lc.set_par(n_bin=n_bin, n_int=n_int, maximize=maximize, frac_gap=frac_gap)
    lc.read_data(name_file=name_infits)
    # ---------------------------------- #
    # ----- Read light curve (end) ----- #
    # ---------------------------------- #

    # ------------------------------------------------- #
    # ----- Prepare for Fourier transform (begin) ----- #
    # ------------------------------------------------- #
    lc.prepare_ft()
    # ------------------------------------------------- #
    # ----- Prepare for Fourier transform (end)   ----- #
    # ------------------------------------------------- #

    # ------------------------------------- #
    # ----- Print information (begin) ----- #
    # ------------------------------------- #
    lc.print_obs_info()
    lc.print_ana_info()
    # ------------------------------------- #
    # ----- Print information (end)   ----- #
    # ------------------------------------- #

    # ---------------------------------- #
    # ----- Call FFT class (begin) ----- #
    # ---------------------------------- #
    i_int=0 # Index of interval
    while True:
        t_start, t_end, rates, drates=lc.extract_interval(i_int=i_int)
        #print(rates)
        if t_start==None:
            break
        ft=Fourier()
        ft.t_start=t_start
        ft.t_end=t_end
        ft.f_min=lc.f_min
        ft.f_max=lc.f_max
        ft.dt=lc.dt
        ft.fft_calc(xs=rates, dxs=drates)
        ft.i_int=i_int
        ft.write_ft(\
            ch_min=ch_min,\
            ch_max=ch_max,\
            name_file=name_outfits,\
            telescope=lc.telescope,\
            instrument=lc.instrument,\
            source=lc.source,\
            exposure=lc.exposure)
        i_int+=1
        if maximize==False and n_int==i_int:
            break
    n_int_fin=i_int

    print('\n{0} intervals were processed successfully.'.format(n_int_fin))
    print('Results are stored in {0}.'.format(name_outfits))

if __name__=='__main__':
    main()
