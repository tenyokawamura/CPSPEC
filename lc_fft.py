from parameters import *
from constants import *
from functions import *
import numpy as np
from numpy import random
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages                            
import sys
from sys import argv
import astropy.io.fits as fits
def main():
    if not len(sys.argv)==2:
        print('Error: The number of inputs is wrong(0).(must be {1})'.format(len(sys.argv), 2))
        sys.exit()
    name_inlc=argv[1]

    ##### Setting (Start) #####
    dt=NEWBIN
    n_d=NEWBINS_PER_INTERVAL
    n_i_exp=INTERVALS_PER_FRAME
    r_g=R_G
    name_outtxt=name_inlc.replace(EXTENSION_LC, '_fft{0}'.format(EXTENSION_TXT))
    display=True
    ##### Setting (End) #####

    ##### Check (Start) #####
    check_existence(name_inlc)
    check_extension(name_inlc, EXTENSION_LC)
    check_extension(name_outtxt, EXTENSION_TXT)
    ##### Check (End) #####

    ##### Read file (Start) #####
    fin=fits.open(name_inlc)
    header=fin[1].header
    data=fin[1].data
    telescope=header['TELESCOP']
    instrument=header['INSTRUME']
    source=header['OBJECT']
    exposure=header['EXPOSURE']
    resol_t=header['TIMEDEL']
    n_row=header['NAXIS2']
    print_obs_info(telescope=telescope,\
                   instrument=instrument,\
                   source=source,\
                   exposure=exposure,\
                   resol_t=resol_t,\
                   n_row=n_row)
    ##### Read file (End) #####

    ##### Preparation (Start) #####
    n_i_max=n_row//n_d
    lar_t=(n_d-1)*dt
    f_min=1/lar_t
    f_max=1/(2.*dt)
    print_ana_info(dt=dt,\
                   n_d=n_d,\
                   lar_t=lar_t,\
                   n_i_max=n_i_max,\
                   n_i_exp=n_i_exp,\
                   f_min=f_min,\
                   f_max=f_max)
    ##### Preparation (End) #####

    ##### Write information (Start) #####
    write_info(name_inlc=name_inlc,\
               name_outtxt=name_outtxt,\
               telescope=telescope,\
               instrument=instrument,\
               source=source,\
               exposure=exposure,\
               resol_t=resol_t,\
               n_row=n_row,\
               dt=dt,\
               n_d=n_d,\
               lar_t=lar_t,\
               n_i_max=n_i_max,\
               n_i_exp=n_i_exp,\
               f_min=f_min,\
               f_max=f_max)
    ##### Write information (End) #####

    ##### Check (Start) #####
    if not abs(dt-resol_t)<10**(-6):
        print('Error: Under construction')
        sys.exit()
    ##### Check (End) #####

    ##### FFT (Start) #####
    n_g_lim=int(n_d*r_g)
    n_i_suc=0
    n_g_c=0 #Nubmer of gaps in one segment
    is_g=np.empty(0, dtype=np.int8) #index of gaps in one segment
    ns_g=np.empty(0, dtype=np.int8) #index of gaps in one segment
    i_d_s=0
    n_d_c=0 #Current number of data in one segment
    print('(No. {0})'.format(n_i_suc+1), end=' ', flush=True)
    for i_d in range(len(data)):
        t_d=data[i_d][0]
        if i_d==i_d_s:
            t_s=t_d
        t_exp=t_s+dt*n_d_c
        if abs(t_d-t_exp)<dt:
            n_d_c=n_d_c+1
            if n_d_c==n_d:
                i_d_e=i_d+1
                t_e=t_d
                data_grp=data[i_d_s:i_d_e]

                n_d_real=len(data_grp)
                xs=np.zeros(n_d_real)
                for i_d_real in range(n_d_real):
                    xs[i_d_real]=data_grp[i_d_real][1]
                x_mea, x_var=mea_var_calc_1d(xs=xs)
                x_sig=np.sqrt(x_var)

                print('\nNumber of actual data: {0}'.format(n_d_real))
                print('Number of generated data: {0}'.format(n_g_c))
                
                for i_g, n_g in zip(is_g[::-1], ns_g[::-1]):
                    xs_g=np.zeros(n_g)
                    for x_g in xs_g:
                        x_g=random.normal(x_mea, x_sig)
                    i_x=i_g-i_d_s
                    xs=np.insert(xs, i_x, xs_g)

                if not len(xs)==n_d:
                    print('Error')
                    sys.exit()

                ##### FT (Start) #####
                #f, b=ft_calc(data_grp,\
                #             n_d,\
                #             dt,\
                #             lar_t)
                #f, b, cr_mean=fft_calc(data=data_grp,\
                #                       n_d=n_d,\
                #                       dt=dt,\
                #                       lar_t=lar_t)
                f, b, cr_mean=fft_calc_v2(xs=xs,\
                                          n_d=n_d,\
                                          dt=dt)
                del data_grp
                ##### FT (End) #####

                ##### Write PSD (Start) #####
                write_ft(name_outtxt,\
                          n_i_suc,\
                          t_s,\
                          t_e,\
                          cr_mean,\
                          f,\
                          b)
                del f
                del b
                ##### Write PSD (Start) #####

                n_d_c=0
                n_g_c=0
                is_g=np.empty(0, dtype=np.int8) 
                ns_g=np.empty(0, dtype=np.int8) 
                i_d_s=i_d+1
                n_i_suc=n_i_suc+1
                if n_i_suc==n_i_exp:
                    break
                print('\n(No. {0})'.format(n_i_suc+1), end=' ', flush=True)
        else:
            t_g=t_d-t_exp
            n_g=int(round((t_d-t_exp)/dt-1.))
            n_g_c=n_g_c+n_g
            if n_g_c>n_g_lim:
                print('.', end='', flush=True)
                n_d_c=0
                n_g_c=0
                is_g=np.empty(0, dtype=np.int8) 
                ns_g=np.empty(0, dtype=np.int8) 
                i_d_s=i_d+1
            else:
                n_d_c=n_d_c+n_g+1
                is_g=np.append(is_g, i_d)
                ns_g=np.append(ns_g, n_g)

    print('')
    del data

    if n_i_suc>0:
        print('{0} intervals were processed successfully.'.format(n_i_suc))
    elif n_i_suc==0:
        print('Error: No interval was processed successfully.')
    ##### FFT (End) #####

if __name__=='__main__':
    main()
