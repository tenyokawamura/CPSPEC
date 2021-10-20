from parameters import *
from constants import *
from functions import *
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages                            
import sys
from sys import argv
def main():
    ##### Setting (Start) #####
    rebin=REBIN
    #name_intxt='ni1200120106_0mpu7_cl_0051_0100_fft.txt'
    #name_intxt='ni1200120106_0mpu7_cl_0101_1000_fft.txt'
    name_intxt='ni1200120189_0mpu7_cl_0051_1200_fft.txt'
    name_outtxt=name_intxt.replace('_fft', '_psd_ns_ave')
    name_outpdf=name_outtxt.replace(EXTENSION_TXT, EXTENSION_PDF)
    display=True
    ##### Setting (End) #####

    ##### Check (Start) #####
    check_extension(name_intxt, EXTENSION_TXT)
    check_extension(name_outtxt, EXTENSION_TXT)
    check_extension(name_outpdf, EXTENSION_PDF)
    ##### Check (End) #####

    ##### Read file (End) #####
    telescope,\
    instrument,\
    source,\
    exposure,\
    resol_t,\
    n_row\
    =read_obs_info(name_intxt=name_intxt)
    print_obs_info(telescope=telescope,\
                   instrument=instrument,\
                   source=source,\
                   exposure=exposure,\
                   resol_t=resol_t,\
                   n_row=n_row)

    dt,\
    n_d,\
    lar_t,\
    n_i_max,\
    n_i_exp,\
    f_min,\
    f_max\
    =read_ana_info(name_intxt=name_intxt)
    print_ana_info(dt=dt,\
                   n_d=n_d,\
                   lar_t=lar_t,\
                   n_i_max=n_i_max,\
                   n_i_exp=n_i_exp,\
                   f_min=f_min,\
                   f_max=f_max)
    ##### Read file (End) #####

    ##### Calculate PSD (Start) #####
    i_i=0
    while True: 
        print('(No. {0})'.format(i_i+1))
        f,\
        b,\
        t_s,\
        t_e,\
        cr_mean,\
        flag_end\
        =read_ft(name_intxt=name_intxt,\
                 no_i=i_i+1)
        psd=psd_calc(bs=b,\
                     x_mean=cr_mean,\
                     n_d=n_d,\
                     dt=dt)
        if i_i==0:
            f_buf=f
            psds=psd
        else:
            psds=np.vstack((psds, psd))

        i_i=i_i+1
        if flag_end==True:
            break

    n_p=i_i
    f=f_buf
    if n_p>1:
        print('{0} periodograms are calculated successfully.'.format(n_p))
        f_fin,\
        psd_fin_mean,\
        psd_fin_sigma,\
        n_mer_p,\
        n_mer_f\
        =psd_unify(f,\
                   psds,\
                   n_p,\
                   rebin)

    elif n_p==1:
        print('Error: Under construction.')
        sys.exit()

    else:
        print('Error: No periodogram was calculated successfully.')
        sys.exit()
    ##### Calculate PSD (End) #####

    ##### Write PSD (Start) #####
    write_psd(name_outtxt=name_outtxt,\
              fs=f_fin,\
              psds_mean=psd_fin_mean,\
              psds_sigma=psd_fin_sigma,\
              ns_mer_p=n_mer_p,\
              ns_mer_f=n_mer_f)
    ##### Write PSD (Start) #####

    ##### Plot PSD (Start) #####
    f_plx,\
    f_ply,\
    psd_mean_pl,\
    f_errx,\
    f_erry,\
    psd_mean_err,\
    psd_sigma_err\
    =params_adjust_plot(f_fin,\
                        psd_fin_mean,\
                        psd_fin_sigma)

    out_pdf=PdfPages(name_outpdf)
    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(f_plx,\
            f_ply*psd_mean_pl,\
            drawstyle='steps-post',\
            color='black',\
            linestyle='solid',\
            linewidth=1.2,\
            alpha=0.8)
    ax.errorbar(f_errx,\
                f_erry*psd_mean_err,\
                yerr=f_erry*psd_sigma_err,\
                capsize=0,\
                color='black',\
                linestyle='None',\
                linewidth=1.2,\
                alpha=0.8)
    #ax.set_xlim(0, int(num_tot-1))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    ax.set_ylim(10**(-5), 10**(-1))
    #ax.set_yticks(np.linspace(0, 700, 8))
    #ax.set_yticks(np.linspace(0, 700, 15), minor=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Frequency (Hz)', fontsize=18)
    ax.set_ylabel("νP(ν) ((rms/mean)$^2$)", fontsize=18)
    ax.tick_params(which='major',\
                   direction="in",\
                   length=10,\
                   width=1.0,\
                   labelsize=16,\
                   top=True,\
                   bottom=True,\
                   right=True,\
                   left=True)
    ax.tick_params(which='minor',\
                   direction="in",\
                   length=5,\
                   width=1.0,\
                   labelsize=16,\
                   top=True,\
                   bottom=True,\
                   right=True,\
                   left=True)
    #ax.legend(bbox_to_anchor=(0, 1),\
    #          loc="upper left",\
    #          borderaxespad=1.,\
    #          fancybox=0,\
    #          edgecolor='black',\
    #          fontsize=16)
    ax.set_title('Number of periodograms: {0}'.format(n_p), fontsize=18)
    #plt.show()

    out_pdf.savefig()
    out_pdf.close()
    ########## Plot PSD (Start) ##########

if __name__=='__main__':
    main()
