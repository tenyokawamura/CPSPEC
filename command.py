import subprocess
import shlex
def main():
    # ----------------------------- #
    # ---------- Setting ---------- #
    # ----------------------------- #
    ### General ###
    name_inevt='ni1200120110_0mpu7_cl_bary_fpm.evt'
    chs_min=[51]
    chs_max=[150]
    ### Light curve ###
    dt=5.e-2 #[s]
    ### FFT ###
    n_bin=16384
    n_int=40
    maximize=True
    frac_gap=0
    ### PSD ###
    rebin=1.2

    # -------------------------- #
    # ---------- Main ---------- #
    # -------------------------- #
    #####################################
    ########## Energy spectrum ##########
    #####################################
    cmd='sh xselect_pha.sh {0}'.format(name_inevt, 'session')
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    #################################
    ########## Light curve ##########
    #################################
    for ch_min, ch_max in zip(chs_min, chs_max):
        ch_min_str='{:04}'.format(ch_min)
        ch_max_str='{:04}'.format(ch_max)
        cmd='sh xselect_lc.sh {0} {1} {2} {3} {4} {5} {6}'\
            .format(name_inevt, 'session', ch_min, ch_max, ch_min_str, ch_max_str, dt)
        tokens=shlex.split(cmd)
        subprocess.run(tokens)

    #########################
    ########## FFT ##########
    #########################
    for ch_min, ch_max in zip(chs_min, chs_max):
        name_infits =name_inevt.replace('.evt', '_{0:04}_{1:04}{2}'    .format(ch_min, ch_max, '.lc'))
        name_outfits=name_inevt.replace('.evt', '_{0:04}_{1:04}_fft{2}'.format(ch_min, ch_max, '.fits'))
        cmd='python do_fft.py {0} {1} {2} {3} {4} {5} {6} {7}'\
            .format(ch_min, ch_max, n_bin, n_int, maximize, frac_gap, name_infits, name_outfits)
        tokens=shlex.split(cmd)
        subprocess.run(tokens)

    #########################
    ########## PSD ##########
    #########################
    for ch_min, ch_max in zip(chs_min, chs_max):
        name_infits =name_inevt.replace('.evt', '_{0:04}_{1:04}_fft{2}'.format(ch_min, ch_max, '.fits'))
        name_outfits=name_inevt.replace('.evt', '_{0:04}_{1:04}_psd{2}'.format(ch_min, ch_max, '.fits'))
        cmd='python do_calc_psd.py {0} {1} {2}'\
            .format(rebin, name_infits, name_outfits)
        tokens=shlex.split(cmd)
        subprocess.run(tokens)

if __name__=='__main__':
    main()
