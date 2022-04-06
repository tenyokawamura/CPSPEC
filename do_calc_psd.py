import numpy as np
import sys
from sys import argv
from fourier import *
from powerspectrum import *

def main():
    # --------------------------- #
    # ----- Setting (begin) ----- #
    # --------------------------- #
    if not len(sys.argv)==4:
        print('Error: Number of inputs is wrong(0).(must be {1})'\
            .format(len(sys.argv), 4))
        sys.exit()
    # Number of bins per interval
    rebin=float(argv[1])
    # Filename of Fourier transform
    name_infits=argv[2]
    # Filename of PSD
    name_outfits=argv[3]
    # --------------------------- #
    # ----- Setting (end)   ----- #
    # --------------------------- #

    # -------------------------------------------- #
    # ----- Read overall information (begin) ----- #
    # -------------------------------------------- #
    ft=Fourier()
    ft.read_info(name_file=name_infits)
    ft.print_info()
    # -------------------------------------------- #
    # ----- Read overall information (end)   ----- #
    # -------------------------------------------- #

    # ----------------------------------------- #
    # ----- Calculate periodogram (begin) ----- #
    # ----------------------------------------- #
    ps=PowerSpectrum()
    ps.dt=ft.dt
    ps.n_bin=ft.n_bin
    for i_int in range(ft.n_int):
        print('(Interval {0})'.format(i_int+1))
        fs, bs, rate_mea, rate_var_mea=ft.read_ft(i_int=i_int)
        ps.calc_psd(bs=bs, rate_mea=rate_mea, rate_var_mea=rate_var_mea)
    # ----------------------------------------- #
    # ----- Calculate periodogram (end)   ----- #
    # ----------------------------------------- #

    # --------------------------------- #
    # ----- Calculate PSD (begin) ----- #
    # --------------------------------- #
    ps.fs=ft.fs
    if rebin<=1:
        print('Error: parameter \'rebin\' must be greater than 1.')
        sys.exit()
    ps.rebin=rebin
    ps.average_psd()
    # --------------------------------- #
    # ----- Calculate PSD (end)   ----- #
    # --------------------------------- #

    # ----------------------------- #
    # ----- Write PSD (begin) ----- #
    # ----------------------------- #
    ps.write_psd(\
        ch_min=ft.ch_min,\
        ch_max=ft.ch_max,\
        name_fits=name_outfits,\
        telescope=ft.telescope,\
        instrument=ft.instrument,\
        source=ft.source,\
        exposure=ft.exposure)
    # ----------------------------- #
    # ----- Write PSD (end)   ----- #
    # ----------------------------- #
    print('PSD was calculated successfully.')
    print('Results are stored in {0}.'.format(name_outfits))

if __name__=='__main__':
    main()
