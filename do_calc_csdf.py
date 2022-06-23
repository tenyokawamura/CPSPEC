import numpy as np
import sys
from sys import argv
from fourier import *
from powerspectrum import *
from crossspectrum import *
##########################################################################
##### Positive lag if the reference band lags behind the energy band #####
##########################################################################
def main():
    # --------------------------- #
    # ----- Setting (begin) ----- #
    # --------------------------- #
    if not len(sys.argv)==8:
        print('Error: Number of inputs is wrong({0}).(must be {1})'\
            .format(len(sys.argv), 8))
        sys.exit()
    # Number of bins per interval
    rebin              =float(argv[1])
    # Whether enery band is overlaped
    overlap            =argv[2]
    # Filename of Fourier transform
    name_infits_fft    =argv[3]
    # Filename of Fourier transform for reference
    name_infits_fft_ref=argv[4]
    # Filename of PSD
    name_infits_psd    =argv[5]
    # Filename of PSD for reference
    name_infits_psd_ref=argv[6]
    # Filename of CSD
    name_outfits       =argv[7]
    if overlap=='True':
        overlap=True
    else:
        overlap=False
    # --------------------------- #
    # ----- Setting (end)   ----- #
    # --------------------------- #

    # -------------------------------------------- #
    # ----- Read overall information (begin) ----- #
    # -------------------------------------------- #
    # Energy band
    ft=Fourier()
    ft.read_info(name_file=name_infits_fft)
    ft.print_info()

    # Reference band
    ft_ref=Fourier()
    ft_ref.read_info(name_file=name_infits_fft_ref)
    ft_ref.print_info()
    # -------------------------------------------- #
    # ----- Read overall information (end)   ----- #
    # -------------------------------------------- #

    # --------------------------------- #
    # ----- Calculate CSD (begin) ----- #
    # --------------------------------- #
    cs=CrossSpectrum()
    cs.dt=ft.dt
    cs.n_bin=ft.n_bin
    cs.overlap=overlap
    for i_int in range(ft.n_int):
        print('(Interval {0})'.format(i_int+1))
        fs,     bs,     rate_mea,     rate_var_mea    =ft.read_ft(i_int=i_int)
        fs_ref, bs_ref, rate_ref_mea, rate_ref_var_mea=ft_ref.read_ft(i_int=i_int)
        cs.calc_csdf(bs=bs,         rate_mea=rate_mea,         rate_var_mea=rate_var_mea,\
                     bs_ref=bs_ref, rate_ref_mea=rate_ref_mea, rate_ref_var_mea=rate_ref_var_mea)
    # --------------------------------- #
    # ----- Calculate CSD (end)   ----- #
    # --------------------------------- #

    # Necessary for calculating errors of CSD.
    # ---------------------------- #
    # ----- Read PSD (begin) ----- #
    # ---------------------------- #
    # Energy band
    ps=PowerSpectrum()
    fs, psds_raw_mea, psds_raw_sig, psds_noi, ns_fbin, ns_int\
        =ps.read_psd(name_file=name_infits_psd)

    # Reference band
    ps_ref=PowerSpectrum()
    fs_ref, psds_ref_raw_mea, psds_ref_raw_sig, psds_ref_noi, ns_ref_fbin, ns_ref_int\
        =ps_ref.read_psd(name_file=name_infits_psd_ref)
    # ---------------------------- #
    # ----- Read PSD (end)   ----- #
    # ---------------------------- #

    # ----------------------------- #
    # ----- Merge CSD (begin) ----- #
    # ----------------------------- #
    cs.fs=ft.fs
    if rebin<=1:
        print('Error: parameter \'rebin\' must be greater than 1.')
        sys.exit()
    cs.rebin=rebin
    cs.average_csdf()
    cs.calc_csdf_error(\
        psds_raw_mea    =psds_raw_mea,\
        psds_noi        =psds_noi,\
        ns_fbin         =ns_fbin,\
        ns_int          =ns_int,\
        psds_ref_raw_mea=psds_ref_raw_mea,\
        psds_ref_noi    =psds_ref_noi,\
        ns_ref_fbin     =ns_ref_fbin,\
        ns_ref_int      =ns_ref_int)
    # ----------------------------- #
    # ----- Merge CSD (end)   ----- #
    # ----------------------------- #

    # ----------------------------- #
    # ----- Write PSD (begin) ----- #
    # ----------------------------- #
    cs.write_csdf(\
        ch_min=ft.ch_min,\
        ch_max=ft.ch_max,\
        ch_ref_min=ft_ref.ch_min,\
        ch_ref_max=ft_ref.ch_max,\
        name_fits=name_outfits,\
        telescope=ft.telescope,\
        instrument=ft.instrument,\
        source=ft.source,\
        exposure=ft.exposure)
    # ----------------------------- #
    # ----- Write PSD (end)   ----- #
    # ----------------------------- #
    print('CSD was calculated successfully.')
    print('Results are stored in {0}.'.format(name_outfits))

if __name__=='__main__':
    main()
