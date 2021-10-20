from parameters import *
from constants import *
import subprocess
import shlex
def main():
    '''
    #####################################
    ########## Energy spectrum ##########
    #####################################
    # ----- PHA ----- #
    cmd='sh xselect_pha.sh {0}'.format(NAME_INEVT)
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    # ----- Count rate ----- #
    name_inpha=NAME_INEVT.replace(EXTENSION_EVT, EXTENSION_PHA)
    name_outqdp=NAME_INEVT.replace(EXTENSION_EVT, '_cr')
    cmd='sh xspec_cr.sh {0} {1} {2} {3}'.format(name_inpha,\
                                                name_outqdp,\
                                                CH_MIN,\
                                                CH_MAX)
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    ################################################
    ########## Light curve/Power spectrum ##########
    ################################################
    ##### Extract light curves #####
    ### Reference ###
    ref_min_str='{:04}'.format(REF_MIN)
    ref_max_str='{:04}'.format(REF_MAX)
    cmd='sh xselect_lc.sh {0} {1} {2} {3} {4} {5} {6}'.format(NAME_INEVT,\
                                                              NAME_SESSION,\
                                                              REF_MIN,\
                                                              REF_MAX,\
                                                              ref_min_str,\
                                                              ref_max_str,\
                                                              BINSIZE)
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    ### Channel-of-Interest ###
    for coi_min, coi_max in zip(COIS_MIN, COIS_MAX):
        coi_min_str='{:04}'.format(coi_min)
        coi_max_str='{:04}'.format(coi_max)
        cmd='sh xselect_lc.sh {0} {1} {2} {3} {4} {5} {6}'.format(NAME_INEVT,\
                                                                  NAME_SESSION,\
                                                                  coi_min,\
                                                                  coi_max,\
                                                                  coi_min_str,\
                                                                  coi_max_str,\
                                                                  BINSIZE)
        tokens=shlex.split(cmd)
        subprocess.run(tokens)

    ########## FFT ##########
    ### Reference ###
    name_inlc=NAME_INEVT.replace(EXTENSION_EVT, '_{0:04}_{1:04}{2}'.format(REF_MIN, REF_MAX, EXTENSION_LC))
    cmd='python lc_fft.py {0}'.format(name_inlc)
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    ### Channel-of-Interest ###
    for coi_min, coi_max in zip(COIS_MIN, COIS_MAX):
        name_inlc=NAME_INEVT.replace(EXTENSION_EVT, '_{0:04}_{1:04}{2}'.format(coi_min, coi_max, EXTENSION_LC))
        cmd='python lc_fft.py {0}'.format(name_inlc)
        tokens=shlex.split(cmd)
        subprocess.run(tokens)

    ########## Power spectra ##########
    ### Reference ###
    name_intxt=NAME_INEVT.replace(EXTENSION_EVT, '_{0:04}_{1:04}_fft{2}'.format(REF_MIN, REF_MAX, EXTENSION_TXT))
    cmd='python psd_calc_ave_ns.py {0}'.format(name_intxt)
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    ### Channel-of-Interest ###
    for coi_min, coi_max in zip(COIS_MIN, COIS_MAX):
        name_intxt=NAME_INEVT.replace(EXTENSION_EVT, '_{0:04}_{1:04}_fft{2}'.format(coi_min, coi_max, EXTENSION_TXT))
        cmd='python psd_calc_ave_ns.py {0}'.format(name_intxt)
        tokens=shlex.split(cmd)
        subprocess.run(tokens)
    '''

    ########## Frequency-dependent cross spectra ##########
    for coi_min, coi_max in zip(COIS_MIN, COIS_MAX):
        cmd='python csd_f_calc.py {0} {1} {2} {3} {4}'.format(NAME_INEVT,\
                                                              coi_min,\
                                                              coi_max,\
                                                              REF_MIN,\
                                                              REF_MAX)
        tokens=shlex.split(cmd)
        subprocess.run(tokens)

    ########## Energy-dependent cross spectra ##########
    cmd='python csd_e_calc.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    ########## RMS spectra ##########
    cmd='python rms_calc.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    ########## Covariance spectra ##########
    cmd='python cov_calc.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

if __name__=='__main__':
    main()
