import os
import subprocess
import shlex
import time
import sys
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from include import inputs_nicer
def main():
    # -------------------------- #
    # ---------- Main ---------- #
    # -------------------------- #
    '''
    ##########################################
    ########## Extract light curves ##########
    ##########################################
    #multi=True
    multi=False
    time_start=time.time()
    if multi==False:
        # ----- Single process ----- #
        for ch_min, ch_max in zip(inputs_nicer.chs_min, inputs_nicer.chs_max):
            pars=[inputs_nicer.name_inevt, inputs_nicer.session, ch_min, ch_max, inputs_nicer.binsize]
            command_do_extract_lc(pars=pars)
    elif multi==True:
        # ----- Multi process ----- #
        first=True
        for ch_min, ch_max in zip(inputs_nicer.chs_min, inputs_nicer.chs_max):
            pars=[inputs_nicer.name_inevt, inputs_nicer.session, ch_min, ch_max, inputs_nicer.binsize]
            if first==True:
                first=False
                parss=[pars]
            else:
                parss.append(pars)
        with ProcessPoolExecutor() as executor:
            executor.map(command_do_extract_lc, parss)
    time_end=time.time()
    time_process=time_end-time_start
    print('------------------------------')
    print('Run time: {0:.1f} s'.format(time_process))

    #########################
    ########## FFT ##########
    #########################
    #multi=True
    multi=False
    time_start=time.time()
    if multi==False:
        # ----- Single process ----- #
        for ch_min, ch_max in zip(inputs_nicer.chs_min, inputs_nicer.chs_max):
            name_inlc=inputs_nicer.name_inevt.replace('.evt', '_{0:04}_{1:04}.lc'.format(int(ch_min), int(ch_max)))
            name_outfits=name_inlc.replace('.lc', '_fft.fits')
            pars=[name_inlc, name_outfits, ch_min, ch_max, inputs_nicer.n_bin, inputs_nicer.n_int, inputs_nicer.maximize, inputs_nicer.frac_gap]
            command_do_fft(pars=pars)
    elif multi==True:
        # ----- Multi process ----- #
        first=True
        for ch_min, ch_max in zip(inputs_nicer.chs_min, inputs_nicer.chs_max):
            name_inlc=inputs_nicer.name_inevt.replace('.evt', '_{0:04}_{1:04}.lc'.format(int(ch_min), int(ch_max)))
            name_outfits=name_inlc.replace('.lc', '_fft.fits')
            pars=[name_inlc, name_outfits, ch_min, ch_max, inputs_nicer.n_bin, inputs_nicer.n_int, inputs_nicer.maximize, inputs_nicer.frac_gap]
            if first==True:
                first=False
                parss=[pars]
            else:
                parss.append(pars)
        with ProcessPoolExecutor() as executor:
            executor.map(command_do_fft, parss)
    time_end=time.time()
    time_process=time_end-time_start
    print('------------------------------')
    print('Run time: {0:.1f} s'.format(time_process))

    #########################
    ########## PSD ##########
    #########################
    #multi=True
    multi=False
    time_start=time.time()
    if multi==False:
        # ----- Single process ----- #
        for ch_min, ch_max in zip(inputs_nicer.chs_min, inputs_nicer.chs_max):
            name_infits =inputs_nicer.name_inevt.replace('.evt', '_{0:04}_{1:04}_fft.fits'.format(int(ch_min), int(ch_max)))
            name_outfits=inputs_nicer.name_inevt.replace('.evt', '_{0:04}_{1:04}_psd.fits'.format(int(ch_min), int(ch_max)))
            pars=[name_infits, name_outfits, ch_min, ch_max, inputs_nicer.rebin]
            command_do_calc_psd(pars=pars)
    elif multi==True:
        # ----- Multi process ----- #
        first=True
        for ch_min, ch_max in zip(inputs_nicer.chs_min, inputs_nicer.chs_max):
            name_infits =inputs_nicer.name_inevt.replace('.evt', '_{0:04}_{1:04}_fft.fits'.format(int(ch_min), int(ch_max)))
            name_outfits=inputs_nicer.name_inevt.replace('.evt', '_{0:04}_{1:04}_psd.fits'.format(int(ch_min), int(ch_max)))
            pars=[name_infits, name_outfits, ch_min, ch_max, inputs_nicer.rebin]
            if first==True:
                first=False
                parss=[pars]
            else:
                parss.append(pars)
        with ProcessPoolExecutor() as executor:
            executor.map(command_do_fft, parss)
    time_end=time.time()
    time_process=time_end-time_start
    print('------------------------------')
    print('Run time: {0:.1f} s'.format(time_process))
    '''

    #########################
    ########## CSD ##########
    #########################
    #multi=True
    multi=False
    time_start=time.time()
    if multi==False:
        # ----- Single process ----- #
        name_infits_fft_ref=inputs_nicer.name_inevt.replace('.evt', '_{0:04}_{1:04}_fft.fits'.format(int(inputs_nicer.ch_ref_min), int(inputs_nicer.ch_ref_max)))
        name_infits_psd_ref=inputs_nicer.name_inevt.replace('.evt', '_{0:04}_{1:04}_psd.fits'.format(int(inputs_nicer.ch_ref_min), int(inputs_nicer.ch_ref_max)))
        for ch_min, ch_max in zip(inputs_nicer.chs_min, inputs_nicer.chs_max):
            name_infits_fft=inputs_nicer.name_inevt.replace('.evt', '_{0:04}_{1:04}_fft.fits'.format(int(ch_min), int(ch_max)))
            name_infits_psd=inputs_nicer.name_inevt.replace('.evt', '_{0:04}_{1:04}_psd.fits'.format(int(ch_min), int(ch_max)))
            name_outfits=inputs_nicer.name_inevt.replace(\
                '.evt',\
                '_{0:04}_{1:04}_{2:04}_{3:04}_csdf.fits'.format(int(ch_min), int(ch_max), int(inputs_nicer.ch_ref_min), int(inputs_nicer.ch_ref_max)))
            pars=[name_infits_fft, name_infits_psd, name_infits_fft_ref, name_infits_psd_ref, name_outfits,\
                  ch_min, ch_max, inputs_nicer.ch_ref_min, inputs_nicer.ch_ref_max, inputs_nicer.rebin]
            command_do_calc_csdf(pars=pars)
    elif multi==True:
        # ----- Multi process ----- #
        name_infits_fft_ref=inputs_nicer.name_inevt.replace('.evt', '_{0:04}_{1:04}_fft.fits'.format(int(inputs_nicer.ch_ref_min), int(inputs_nicer.ch_ref_max)))
        name_infits_psd_ref=inputs_nicer.name_inevt.replace('.evt', '_{0:04}_{1:04}_psd.fits'.format(int(inputs_nicer.ch_ref_min), int(inputs_nicer.ch_ref_max)))
        first=True
        for ch_min, ch_max in zip(inputs_nicer.chs_min, inputs_nicer.chs_max):
            name_infits_fft=inputs_nicer.name_inevt.replace('.evt', '_{0:04}_{1:04}_fft.fits'.format(int(ch_min), int(ch_max)))
            name_infits_psd=inputs_nicer.name_inevt.replace('.evt', '_{0:04}_{1:04}_psd.fits'.format(int(ch_min), int(ch_max)))
            name_outfits=inputs_nicer.name_inevt.replace(\
                '.evt',\
                '_{0:04}_{1:04}_{2:04}_{3:04}_csdf.fits'.format(int(ch_min), int(ch_max), int(inputs_nicer.ch_ref_min), int(inputs_nicer.ch_ref_max)))
            pars=[name_infits_fft, name_infits_psd, name_infits_fft_ref, name_infits_psd_ref, name_outfits,\
                  ch_min, ch_max, inputs_nicer.ch_ref_min, inputs_nicer.ch_ref_max, inputs_nicer.rebin]
            command_do_calc_csdf(pars=pars)
            if first==True:
                first=False
                parss=[pars]
            else:
                parss.append(pars)
        with ProcessPoolExecutor() as executor:
            executor.map(command_do_fft, parss)
    time_end=time.time()
    time_process=time_end-time_start
    print('------------------------------')
    print('Run time: {0:.1f} s'.format(time_process))

def command_do_extract_lc(pars):
    name_inevt  =pars[0]
    session     =pars[1]
    ch_min      =int(pars[2])
    ch_max      =int(pars[3])
    binsize     =pars[4]
    name_outlc  =name_inevt.replace('.evt', '_{0:04}_{1:04}.lc'.format(int(ch_min), int(ch_max))) # No extension

    cmd='sh xselect_lc.sh {0} {1} {2} {3} {4} {5}'\
        .format(name_inevt, name_outlc, session, ch_min, ch_max, binsize)
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

def command_do_fft(pars):
    name_inlc   =pars[0]
    name_outfits=pars[1]
    ch_min      =pars[2]
    ch_max      =pars[3]
    n_bin       =pars[4]
    n_int       =pars[5]
    maximize    =pars[6]
    frac_gap    =pars[7]

    cmd='python do_fft.py {0} {1} {2} {3} {4} {5} {6} {7}'\
        .format(ch_min, ch_max, n_bin, n_int, maximize, frac_gap, name_inlc, name_outfits)
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

def command_do_calc_psd(pars):
    name_infits =pars[0]
    name_outfits=pars[1]
    ch_min      =pars[2]
    ch_max      =pars[3]
    rebin       =pars[4]

    cmd='python do_calc_psd.py {0} {1} {2}'\
        .format(rebin, name_infits, name_outfits)
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

def command_do_calc_csdf(pars):
    name_infits_fft    =pars[0]
    name_infits_psd    =pars[1]
    name_infits_fft_ref=pars[2]
    name_infits_psd_ref=pars[3]
    name_outfits       =pars[4]
    ch_min             =pars[5]
    ch_max             =pars[6]
    ch_ref_min         =pars[7]
    ch_ref_max         =pars[8]
    rebin              =pars[9]

    # Assume that (ch_min<=ch_ref_min) & (ch_ref_max<=ch_max) = False
    if ((ch_ref_min<=ch_min) & (ch_max<=ch_ref_max))==True:
        overlap=True
    else:
        overlap=False
    cmd='python do_calc_csdf.py {0} {1} {2} {3} {4} {5} {6}'\
        .format(rebin, overlap, name_infits_fft, name_infits_fft_ref, name_infits_psd, name_infits_psd_ref, name_outfits)
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

if __name__=='__main__':
    main()
