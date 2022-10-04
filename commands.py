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
    # --- #
    obsid='P0114661036'
    instruments=['LE', 'ME', 'HE']
    #chs_min={'LE': [0,    106, 296, 556, 816],  'ME': [0,    137, 342], 'HE': [0,   14, 24, 37, 60]}
    #chs_max={'LE': [1535, 295, 555, 815, 1288], 'ME': [1023, 341, 546], 'HE': [255, 23, 36, 59, 93]}
    chs_min={'LE': [106, 296, 556, 816],  'ME': [137, 342], 'HE': [14, 24, 37, 60]}
    chs_max={'LE': [295, 555, 815, 1288], 'ME': [341, 546], 'HE': [23, 36, 59, 93]}
    instrument_ref='LE'
    ch_ref_min=296
    ch_ref_max=555
    # FFT
    n_bin=2**15
    n_int=100
    maximize=True
    frac_gap=0
    # PSD/CSD
    rebin=1.1
    # Operations
    operations={'FFT': False, 'PSD': False, 'CSDF': True}
    multi=True
    # --- #

    # -------------------------- #
    # ---------- Main ---------- #
    # -------------------------- #
    #########################
    ########## FFT ##########
    #########################
    if operations['FFT']==True:
        if multi==False:
            # ----- Single process ----- #
            for instrument in instruments:
                for ch_min, ch_max in zip(chs_min[instrument], chs_max[instrument]):
                    name_inlc   ='HXMT_{0}_{1}_{2:04}_{3:04}_net_common.lc'.format(obsid, instrument, int(ch_min), int(ch_max))
                    name_outfits='HXMT_{0}_{1}_{2:04}_{3:04}_net_common_fft.fits'.format(obsid, instrument, int(ch_min), int(ch_max))
                    pars=[name_inlc, name_outfits, ch_min, ch_max, instrument, n_bin, n_int, maximize, frac_gap]
                    command_do_fft(pars=pars)

        elif multi==True:
            # ----- Multi process ----- #
            time_start=time.time()
            first=True
            for instrument in instruments:
                for ch_min, ch_max in zip(chs_min[instrument], chs_max[instrument]):
                    name_inlc   ='HXMT_{0}_{1}_{2:04}_{3:04}_net_common.lc'      .format(obsid, instrument, int(ch_min), int(ch_max))
                    name_outfits='HXMT_{0}_{1}_{2:04}_{3:04}_net_common_fft.fits'.format(obsid, instrument, int(ch_min), int(ch_max))
                    pars=[name_inlc, name_outfits, ch_min, ch_max, instrument, n_bin, n_int, maximize, frac_gap]
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
            print('(Process efficiency for FFT)')
            print('Run time: {0:.1f} s'.format(time_process))

    #########################
    ########## PSD ##########
    #########################
    if operations['PSD']==True:
        if multi==False:
            # ----- Single process ----- #
            for instrument in instruments:
                for ch_min, ch_max in zip(chs_min[instrument], chs_max[instrument]):
                    name_infits ='HXMT_{0}_{1}_{2:04}_{3:04}_net_common_fft.fits'.format(obsid, instrument, int(ch_min), int(ch_max))
                    name_outfits='HXMT_{0}_{1}_{2:04}_{3:04}_net_common_psd.fits'.format(obsid, instrument, int(ch_min), int(ch_max))
                    pars=[name_infits, name_outfits, ch_min, ch_max, instrument, rebin]
                    command_do_calc_psd(pars=pars)

        elif multi==True:
            # ----- Multi process ----- #
            time_start=time.time()
            first=True
            for instrument in instruments:
                for ch_min, ch_max in zip(chs_min[instrument], chs_max[instrument]):
                    name_infits ='HXMT_{0}_{1}_{2:04}_{3:04}_net_common_fft.fits'.format(obsid, instrument, int(ch_min), int(ch_max))
                    name_outfits='HXMT_{0}_{1}_{2:04}_{3:04}_net_common_psd.fits'.format(obsid, instrument, int(ch_min), int(ch_max))
                    pars=[name_infits, name_outfits, ch_min, ch_max, instrument, rebin]
                    if first==True:
                        first=False
                        parss=[pars]
                    else:
                        parss.append(pars)
            with ProcessPoolExecutor() as executor:
                executor.map(command_do_calc_psd, parss)
            time_end=time.time()
            time_process=time_end-time_start
            print('------------------------------')
            print('(Process efficiency for calculating PSD)')
            print('Run time: {0:.1f} s'.format(time_process))

    #########################
    ########## CSD ##########
    #########################
    if operations['CSDF']==True:
        if multi==False:
            # ----- Single process ----- #
            for instrument in instruments:
                for ch_min, ch_max in zip(chs_min[instrument], chs_max[instrument]):
                    name_infits_ref    ='HXMT_{0}_{1}_{2:04}_{3:04}_net_common_fft.fits' .format(obsid, instrument_ref, int(ch_ref_min), int(ch_ref_max))
                    name_infits        ='HXMT_{0}_{1}_{2:04}_{3:04}_net_common_fft.fits' .format(obsid, instrument,     int(ch_min),     int(ch_max))
                    name_infits_psd_ref='HXMT_{0}_{1}_{2:04}_{3:04}_net_common_psd.fits' .format(obsid, instrument_ref, int(ch_ref_min), int(ch_ref_max))
                    name_infits_psd    ='HXMT_{0}_{1}_{2:04}_{3:04}_net_common_psd.fits' .format(obsid, instrument,     int(ch_min),     int(ch_max))
                    name_outfits       ='HXMT_{0}_{1}_{2:04}_{3:04}_{4}_{5:04}_{6:04}_net_common_csdf.fits'.format(obsid, instrument_ref, int(ch_ref_min), int(ch_ref_max), instrument, int(ch_min), int(ch_max))
                    if ((ch_ref_min<=ch_min) & (ch_max<=ch_ref_max) & (instrument==instrument_ref)):
                        overlap=True
                    else:
                        overlap=False
                    pars=[name_infits_ref, name_infits, name_infits_psd_ref, name_infits_psd, name_outfits, overlap, rebin]
                    command_do_calc_csdf(pars=pars)

        elif multi==True:
            # ----- Multi process ----- #
            time_start=time.time()
            first=True
            for instrument in instruments:
                for ch_min, ch_max in zip(chs_min[instrument], chs_max[instrument]):
                    name_infits_ref    ='HXMT_{0}_{1}_{2:04}_{3:04}_net_common_fft.fits' .format(obsid, instrument_ref, int(ch_ref_min), int(ch_ref_max))
                    name_infits        ='HXMT_{0}_{1}_{2:04}_{3:04}_net_common_fft.fits' .format(obsid, instrument,     int(ch_min),     int(ch_max))
                    name_infits_psd_ref='HXMT_{0}_{1}_{2:04}_{3:04}_net_common_psd.fits' .format(obsid, instrument_ref, int(ch_ref_min), int(ch_ref_max))
                    name_infits_psd    ='HXMT_{0}_{1}_{2:04}_{3:04}_net_common_psd.fits' .format(obsid, instrument,     int(ch_min),     int(ch_max))
                    name_outfits       ='HXMT_{0}_{1}_{2:04}_{3:04}_{4}_{5:04}_{6:04}_net_common_csdf.fits'.format(obsid, instrument_ref, int(ch_ref_min), int(ch_ref_max), instrument, int(ch_min), int(ch_max))
                    if ((ch_ref_min<=ch_min) & (ch_max<=ch_ref_max) & (instrument==instrument_ref)):
                        overlap=True
                    else:
                        overlap=False
                    pars=[name_infits_ref, name_infits, name_infits_psd_ref, name_infits_psd, name_outfits, overlap, rebin]
                    if first==True:
                        first=False
                        parss=[pars]
                    else:
                        parss.append(pars)
            with ProcessPoolExecutor() as executor:
                executor.map(command_do_calc_csdf, parss)
            time_end=time.time()
            time_process=time_end-time_start
            print('------------------------------')
            print('(Process efficiency for calculating CSD)')
            print('Run time: {0:.1f} s'.format(time_process))

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
    name_inlc   =pars[0]
    name_outfits=pars[1]
    ch_min      =pars[2]
    ch_max      =pars[3]
    instrument  =pars[4]
    n_bin       =pars[5]
    n_int       =pars[6]
    maximize    =pars[7]
    frac_gap    =pars[8]

    # File (directory) existence
    exist=os.path.exists(name_inlc)
    if exist==False:
        print('Warning: {0} does not exist.'.format(name_inlc))
        return
    cmd='python do_fft.py {0} {1} {2} {3} {4} {5} {6} {7}'\
        .format(ch_min, ch_max, n_bin, n_int, maximize, frac_gap, name_inlc, name_outfits)
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

def command_do_calc_psd(pars):
    name_infits =pars[0]
    name_outfits=pars[1]
    ch_min      =pars[2]
    ch_max      =pars[3]
    instrument  =pars[4]
    rebin       =pars[5]

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
    name_infits_ref    =pars[0]
    name_infits        =pars[1]
    name_infits_psd_ref=pars[2]
    name_infits_psd    =pars[3]
    name_outfits       =pars[4]
    overlap            =pars[5]
    rebin              =pars[6]

    # File (directory) existence
    if os.path.exists(name_infits_ref)==False:
        print('Warning: {0} does not exist.'.format(name_infits_ref))
        return
    if os.path.exists(name_infits)==False:
        print('Warning: {0} does not exist.'.format(name_infits))
        return
    if os.path.exists(name_infits_psd_ref)==False:
        print('Warning: {0} does not exist.'.format(name_infits_psd_ref))
        return
    if os.path.exists(name_infits_psd)==False:
        print('Warning: {0} does not exist.'.format(name_infits_psd))
        return
    cmd='python do_calc_csdf.py {0} {1} {2} {3} {4} {5} {6}'\
        .format(rebin, overlap, name_infits, name_infits_ref, name_infits_psd, name_infits_psd_ref, name_outfits)
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

#def command_do_calc_csdf(pars):
#    name_dir_pro  =pars[0]
#    e_min         =pars[1]
#    e_max         =pars[2]
#    instrument    =pars[3]
#    e_ref_min     =pars[4]
#    e_ref_max     =pars[5]
#    instrument_ref=pars[6]
#    rebin         =pars[7]
#
#    name_infits_fft_ref=name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin_common_fft.fits'.format(instrument_ref, int(e_ref_min), int(e_ref_max))
#    name_infits_psd_ref=name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin_common_psd.fits'.format(instrument_ref, int(e_ref_min), int(e_ref_max))
#    exist=os.path.exists(name_infits_fft_ref)
#    if exist==False:
#        print('Warning: {0} does not exist.'.format(name_infits_fft_ref))
#        return
#    exist=os.path.exists(name_infits_psd_ref)
#    if exist==False:
#        print('Warning: {0} does not exist.'.format(name_infits_psd_ref))
#        return
#
#    name_infits_fft=name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin_common_fft.fits'.format(instrument, int(e_min), int(e_max))
#    name_infits_psd=name_dir_pro+'{0}_net_{1:04}_{2:04}_rebin_common_psd.fits'.format(instrument, int(e_min), int(e_max))
#    exist=os.path.exists(name_infits_fft)
#    if exist==False:
#        print('Warning: {0} does not exist.'.format(name_infits_fft))
#        return
#    exist=os.path.exists(name_infits_psd)
#    if exist==False:
#        print('Warning: {0} does not exist.'.format(name_infits_psd))
#        return
#    name_outfits=name_dir_pro+'{0}_net_{1:04}_{2:04}_{3}_net_{4:04}_{5:04}_rebin_common_csdf.fits'\
#        .format(instrument, int(e_min), int(e_max), instrument_ref, int(e_ref_min), int(e_ref_max))
#
#    # Assume that (e_min<=e_ref_min) & (e_ref_max<=e_max) & (instrument_ref==instrument) = False
#    if ((e_ref_min<=e_min) & (e_max<=e_ref_max) & (instrument==instrument_ref)):
#        overlap=True
#    else:
#        overlap=False
#    cmd='python do_calc_csdf.py {0} {1} {2} {3} {4} {5} {6}'\
#        .format(rebin, overlap, name_infits_fft, name_infits_fft_ref, name_infits_psd, name_infits_psd_ref, name_outfits)
#    tokens=shlex.split(cmd)
#    subprocess.run(tokens)

if __name__=='__main__':
    main()
