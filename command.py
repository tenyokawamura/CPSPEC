import os
import subprocess
import shlex
import time
import sys
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
def main():
    # ----------------------------- #
    # ---------- Setting ---------- #
    # ----------------------------- #
    #name_inevt='ni1200120104_0mpu7_cl_bary_fpm.evt'
    name_inevt='ni1200120104_0mpu7_cl_fpm.evt'

    # --- #
    #n_dir=12
    #names_dir_raw=['./P0114661003_lc/P0114661003{0:02}/'.format(i+1) for i in range(n_dir)]
    #names_dir_pro=['./P0114661003{0:02}_product/'.format(i+1) for i in range(n_dir)]
    ## Common data extraction
    #e_min_le=2.6
    #e_max_le=4.8
    #e_min_me=7
    #e_max_me=11
    #e_min_he=25
    #e_max_he=35
    ## Energy bands
    #es_min       =[2.6,   4.8, 7,    7,    11,   23,   25,   35,   48,   67,   100,  150]
    #es_max       =[4.8,   7,   11,   11,   23,   35,   35,   48,   67,   100,  150,  200]
    #instruments  =['LE', 'LE', 'LE', 'ME', 'ME', 'ME', 'HE', 'HE', 'HE', 'HE', 'HE', 'HE']
    ## Reference band
    #e_ref_min      =2.6
    #e_ref_max      =4.8
    #instrument_ref ='LE'
    #### Light curve rebining ###
    #dt_rebin=1./128. #[s]
    #### FFT ###
    #n_bin=2**15
    #n_int=100
    #maximize=True
    #frac_gap=0
    #### PSD, CSD ###
    #rebin=1.2
    # --- #

    # -------------------------- #
    # ---------- Main ---------- #
    # -------------------------- #
    ###################################################
    ########## Generate energy spectra (PHA) ##########
    ###################################################
    # ----- PHA ----- #
    cmd='sh xselect_pha.sh {0}'.format(name_inevt)
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    ##########################################
    ########## Light curve rebining ##########
    ##########################################
    # ----- Single process ----- #
    #for name_dir_raw, name_dir_pro in zip(names_dir_raw, names_dir_pro):
    #    for e_min, e_max, instrument in zip(es_min, es_max, instruments):
    #        pars=[name_dir_raw, name_dir_pro, e_min, e_max, instrument, dt_rebin]
    #        command_do_rebin_lc(pars=pars)

    # ----- Multi process ----- #
    #time_start=time.time()
    #first=True
    #for name_dir_raw, name_dir_pro in zip(names_dir_raw, names_dir_pro):
    #    for e_min, e_max, instrument in zip(es_min, es_max, instruments):
    #        pars=[name_dir_raw, name_dir_pro, e_min, e_max, instrument, dt_rebin]
    #        if first==True:
    #            first=False
    #            parss=[pars]
    #        else:
    #            parss.append(pars)
    #with ProcessPoolExecutor() as executor:
    #    executor.map(command_do_rebin_lc, parss)
    #time_end=time.time()
    #time_process=time_end-time_start
    #print('------------------------------')
    #print('(Process efficiency for rebinning light curves)')
    #print('Run time: {0:.1f} s'.format(time_process))

    ###########################################################
    ########## Extract common (rebinned) light curve ##########
    ###########################################################
    # ----- Single process ---- #
    #for name_dir_pro in names_dir_pro:
    #    pars=[name_dir_pro, e_min_le, e_max_le, e_min_me, e_max_me, e_min_he, e_max_he]
    #    command_do_extract_index_rebin(pars=pars)
    #for name_dir_pro in names_dir_pro:
    #    for e_min, e_max, instrument in zip(es_min, es_max, instruments):
    #        pars=[name_dir_pro, e_min, e_max, instrument]
    #    command_do_extract_lc_rebin(pars=pars)

    # ----- Multi process ---- #
    #time_start=time.time()
    #first=True
    #for name_dir_pro in names_dir_pro:
    #    pars=[name_dir_pro, e_min_le, e_max_le, e_min_me, e_max_me, e_min_he, e_max_he]
    #    if first==True:
    #        first=False
    #        parss=[pars]
    #    else:
    #        parss.append(pars)
    #with ProcessPoolExecutor() as executor:
    #    executor.map(command_do_extract_index_rebin, parss)
    #time_end=time.time()
    #time_process=time_end-time_start
    #print('------------------------------')
    #print('(Process efficiency for extracting common indices)')
    #print('Run time: {0:.1f} s'.format(time_process))

    #time_start=time.time()
    #first=True
    #for name_dir_pro in names_dir_pro:
    #    for e_min, e_max, instrument in zip(es_min, es_max, instruments):
    #        pars=[name_dir_pro, e_min, e_max, instrument]
    #        if first==True:
    #            first=False
    #            parss=[pars]
    #        else:
    #            parss.append(pars)
    #with ProcessPoolExecutor() as executor:
    #    executor.map(command_do_extract_lc_rebin, parss)
    #time_end=time.time()
    #time_process=time_end-time_start
    #print('------------------------------')
    #print('(Process efficiency for extracting common light curves)')
    #print('Run time: {0:.1f} s'.format(time_process))

    #########################
    ########## FFT ##########
    #########################
    # ----- Single process ----- #
    #for name_dir_raw, name_dir_pro in zip(names_dir_raw, names_dir_pro):
    #    for e_min, e_max, instrument in zip(es_min, es_max, instruments):
    #        pars=[name_dir_raw, name_dir_pro, e_min, e_max, instrument, n_bin, n_int, maximize, frac_gap]
    #        command_do_fft(pars=pars)

    # ----- Multi process ----- #
    #time_start=time.time()
    #first=True
    #for name_dir_raw, name_dir_pro in zip(names_dir_raw, names_dir_pro):
    #    for e_min, e_max, instrument in zip(es_min, es_max, instruments):
    #        pars=[name_dir_raw, name_dir_pro, e_min, e_max, instrument, n_bin, n_int, maximize, frac_gap]
    #        if first==True:
    #            first=False
    #            parss=[pars]
    #        else:
    #            parss.append(pars)
    #with ProcessPoolExecutor() as executor:
    #    executor.map(command_do_fft, parss)
    #time_end=time.time()
    #time_process=time_end-time_start
    #print('------------------------------')
    #print('(Process efficiency for FFT)')
    #print('Run time: {0:.1f} s'.format(time_process))

    #########################
    ########## PSD ##########
    #########################
    # ----- Single process ----- #
    #for name_dir_raw, name_dir_pro in zip(names_dir_raw, names_dir_pro):
    #    for e_min, e_max, instrument in zip(es_min, es_max, instruments):
    #        pars=[name_dir_pro, e_min, e_max, instrument, rebin]
    #        command_do_calc_psd(pars=pars)

    # ----- Multi process ----- #
    #time_start=time.time()
    #first=True
    #for name_dir_raw, name_dir_pro in zip(names_dir_raw, names_dir_pro):
    #    for e_min, e_max, instrument in zip(es_min, es_max, instruments):
    #        pars=[name_dir_pro, e_min, e_max, instrument, rebin]
    #        if first==True:
    #            first=False
    #            parss=[pars]
    #        else:
    #            parss.append(pars)
    #with ProcessPoolExecutor() as executor:
    #    executor.map(command_do_calc_psd, parss)
    #time_end=time.time()
    #time_process=time_end-time_start
    #print('------------------------------')
    #print('(Process efficiency for calculating PSD)')
    #print('Run time: {0:.1f} s'.format(time_process))

    #########################
    ########## CSD ##########
    #########################
    # ----- Single process ----- #
    #for name_dir_raw, name_dir_pro in zip(names_dir_raw, names_dir_pro):
    #    for e_min, e_max, instrument in zip(es_min, es_max, instruments):
    #        pars=[name_dir_pro,\
    #            e_min,     e_max,     instrument,\
    #            e_ref_min, e_ref_max, instrument_ref,\
    #            rebin]
    #        command_do_calc_csdf(pars=pars)

    # ----- Multi process ----- #
    #time_start=time.time()
    #first=True
    #for name_dir_raw, name_dir_pro in zip(names_dir_raw, names_dir_pro):
    #    for e_min, e_max, instrument in zip(es_min, es_max, instruments):
    #        pars=[name_dir_pro,\
    #            e_min,     e_max,     instrument,\
    #            e_ref_min, e_ref_max, instrument_ref,\
    #            rebin]
    #        if first==True:
    #            first=False
    #            parss=[pars]
    #        else:
    #            parss.append(pars)
    #with ProcessPoolExecutor() as executor:
    #    executor.map(command_do_calc_csdf, parss)
    #time_end=time.time()
    #time_process=time_end-time_start
    #print('------------------------------')
    #print('(Process efficiency for calculating CSD)')
    #print('Run time: {0:.1f} s'.format(time_process))

def command_do_convert_counts2rate(pars):
    name_dir_raw=pars[0]
    name_dir_pro=pars[1]
    e_min       =pars[2]
    e_max       =pars[3]
    instrument  =pars[4]

    if instrument=='LE':
        name_inlc =name_dir_raw+'{0}_{1}_{2}_g0_0-94.lc'.format(instrument, e_min, e_max)
    elif instrument=='ME':
        name_inlc =name_dir_raw+'{0}_{1}_{2}_g0_0-53.lc'.format(instrument, e_min, e_max)
    elif instrument=='HE':
        name_inlc =name_dir_raw+'{0}_{1}_{2}_g0_0-17.lc'.format(instrument, e_min, e_max)
    else:
        print('Warning: unexpected instrument {0}.'.format(instrument))
        return
    name_outlc=name_dir_pro+'{0}_{1:04}_{2:04}.lc'.format(instrument, int(e_min), int(e_max))

    exist=os.path.exists(name_inlc)
    if exist==False:
        print('Error: {0} does not exist.'.format(name_inlc))
        return
    exist=os.path.exists(name_dir_pro)
    if exist==False:
        print('Error: {0} does not exist.'.format(name_dir_pro))
        return
    cmd='python do_convert_counts2rate.py {0} {1}'\
        .format(name_inlc, name_outlc)
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

def command_do_rebin_lc(pars):
    name_dir_raw=pars[0]
    name_dir_pro=pars[1]
    e_min       =pars[2]
    e_max       =pars[3]
    instrument  =pars[4]
    dt_rebin    =pars[5]

    name_inlc =name_dir_pro+'{0}_net_{1:04}_{2:04}.lc'.format(instrument, int(e_min), int(e_max))
    name_outlc=name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin.lc'.format(instrument, int(e_min), int(e_max))
    # File (directory) existence
    exist=os.path.exists(name_inlc)
    if exist==False:
        print('Error: {0} does not exist.'.format(name_inlc))
        return
    exist=os.path.exists(name_dir_pro)
    if exist==False:
        print('Error: {0} does not exist.'.format(name_dir_pro))
        return
    cmd='python do_rebin_lc.py {0} {1} {2}'\
        .format(name_inlc, name_outlc, dt_rebin)
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

def command_do_extract_index_rebin(pars):
    name_dir_pro=pars[0]
    e_min_le    =pars[1]
    e_max_le    =pars[2]
    e_min_me    =pars[3]
    e_max_me    =pars[4]
    e_min_he    =pars[5]
    e_max_he    =pars[6]

    name_inlc_le=name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin.lc'.format('LE', int(e_min_le), int(e_max_le))
    name_inlc_me=name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin.lc'.format('ME', int(e_min_me), int(e_max_me))
    name_inlc_he=name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin.lc'.format('HE', int(e_min_he), int(e_max_he))

    exist=os.path.exists(name_inlc_le)
    if exist==False:
        print('Error: {0} does not exist.'.format(name_inlc_le))
        return
    exist=os.path.exists(name_inlc_me)
    if exist==False:
        print('Error: {0} does not exist.'.format(name_inlc_me))
        return
    exist=os.path.exists(name_inlc_he)
    if exist==False:
        print('Error: {0} does not exist.'.format(name_inlc_he))
        return

    name_outtxt_le=name_dir_pro+'{0}_net_rebin_common_index.txt'.format('LE')
    name_outtxt_me=name_dir_pro+'{0}_net_rebin_common_index.txt'.format('ME')
    name_outtxt_he=name_dir_pro+'{0}_net_rebin_common_index.txt'.format('HE')

    cmd='python do_extract_index_rebin.py {0} {1} {2} {3} {4} {5}'\
        .format(name_inlc_le, name_inlc_me, name_inlc_he,\
            name_outtxt_le, name_outtxt_me, name_outtxt_he)
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

def command_do_extract_lc_rebin(pars):
    name_dir_pro=pars[0]
    e_min       =pars[1]
    e_max       =pars[2]
    instrument  =pars[3]

    name_inlc=name_dir_pro+'{0}_{1:04}_{2:04}_rebin.lc'.format(instrument, int(e_min), int(e_max))
    exist=os.path.exists(name_inlc)
    if exist==False:
        print('Warning: {0} does not exist.'.format(name_inlc))
        return

    if instrument=='LE':
        name_intxt=name_dir_pro+'{0}_net_rebin_common_index.txt'.format('LE')
    elif instrument=='ME':
        name_intxt=name_dir_pro+'{0}_net_rebin_common_index.txt'.format('ME')
    elif instrument=='HE':
        name_intxt=name_dir_pro+'{0}_net_rebin_common_index.txt'.format('HE')
    else:
        print('Warning: unexpected instrument {0}.'.format(instrument))
        return

    name_outlc=name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin_common.lc'.format(instrument, int(e_min), int(e_max))

    cmd='python do_extract_lc_rebin.py {0} {1} {2}'\
        .format(name_inlc, name_intxt, name_outlc)
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

def command_do_fft(pars):
    name_dir_raw=pars[0]
    name_dir_pro=pars[1]
    e_min       =pars[2]
    e_max       =pars[3]
    instrument  =pars[4]
    n_bin       =pars[5]
    n_int       =pars[6]
    maximize    =pars[7]
    frac_gap    =pars[8]

    # For calculation of power spectrum
    #name_inlc   =name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin.lc'.format(instrument, int(e_min), int(e_max))
    #name_outfits=name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin_fft.fits'.format(instrument, int(e_min), int(e_max))

    # For calculation of cross spectrum
    name_inlc   =name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin_common.lc'.format(instrument, int(e_min), int(e_max))
    name_outfits=name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin_common_fft.fits'.format(instrument, int(e_min), int(e_max))

    # File (directory) existence
    exist=os.path.exists(name_inlc)
    if exist==False:
        print('Warning: {0} does not exist.'.format(name_inlc))
        return
    cmd='python do_fft.py {0} {1} {2} {3} {4} {5} {6} {7}'\
        .format(e_min, e_max, n_bin, n_int, maximize, frac_gap, name_inlc, name_outfits)
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

def command_do_calc_psd(pars):
    name_dir_pro=pars[0]
    e_min       =pars[1]
    e_max       =pars[2]
    instrument  =pars[3]
    rebin       =pars[4]

    # For calculation of power spectrum
    #name_infits =name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin_fft.fits'.format(instrument, int(e_min), int(e_max))
    #name_outfits=name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin_psd.fits'.format(instrument, int(e_min), int(e_max))

    # For calculation of cross spectrum
    name_infits =name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin_common_fft.fits'.format(instrument, int(e_min), int(e_max))
    name_outfits=name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin_common_psd.fits'.format(instrument, int(e_min), int(e_max))

    # File (directory) existence
    exist=os.path.exists(name_infits)
    if exist==False:
        print('Warning: {0} does not exist.'.format(name_infits))
        return
    cmd='python do_calc_psd.py {0} {1} {2}'\
        .format(rebin, name_infits, name_outfits)
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

def command_do_calc_csdf(pars):
    name_dir_pro  =pars[0]
    e_min         =pars[1]
    e_max         =pars[2]
    instrument    =pars[3]
    e_ref_min     =pars[4]
    e_ref_max     =pars[5]
    instrument_ref=pars[6]
    rebin         =pars[7]

    name_infits_fft_ref=name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin_common_fft.fits'.format(instrument_ref, int(e_ref_min), int(e_ref_max))
    name_infits_psd_ref=name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin_common_psd.fits'.format(instrument_ref, int(e_ref_min), int(e_ref_max))
    exist=os.path.exists(name_infits_fft_ref)
    if exist==False:
        print('Warning: {0} does not exist.'.format(name_infits_fft_ref))
        return
    exist=os.path.exists(name_infits_psd_ref)
    if exist==False:
        print('Warning: {0} does not exist.'.format(name_infits_psd_ref))
        return

    name_infits_fft=name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin_common_fft.fits'.format(instrument, int(e_min), int(e_max))
    name_infits_psd=name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin_common_psd.fits'.format(instrument, int(e_min), int(e_max))
    exist=os.path.exists(name_infits_fft)
    if exist==False:
        print('Warning: {0} does not exist.'.format(name_infits_fft))
        return
    exist=os.path.exists(name_infits_psd)
    if exist==False:
        print('Warning: {0} does not exist.'.format(name_infits_psd))
        return
    name_outfits=name_dir_pro+'{0}_net_{1:04}_{2:04}_{3}_net_{4:04}_{5:04}_rebin_common_csdf.fits'\
        .format(instrument, int(e_min), int(e_max), instrument_ref, int(e_ref_min), int(e_ref_max))

    # Assume that (e_min<=e_ref_min) & (e_ref_max<=e_max) & (instrument_ref==instrument) = False
    if ((e_ref_min<=e_min) & (e_max<=e_ref_max) & (instrument==instrument_ref)):
        overlap=True
    else:
        overlap=False
    cmd='python do_calc_csdf.py {0} {1} {2} {3} {4} {5} {6}'\
        .format(rebin, overlap, name_infits_fft, name_infits_fft_ref, name_infits_psd, name_infits_psd_ref, name_outfits)
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

if __name__=='__main__':
    main()
