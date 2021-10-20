from parameters import *
from constants import *
from functions import *
import numpy as np
import cmath
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages                            
import sys
from sys import argv
def main():
    ##### Setting (Start) #####
    rebin=REBIN
    ch_min_1=51
    ch_max_1=151
    ch_min_2=51
    ch_max_2=1200
    name_intxt_psd_def='ni1200120106_0mpu7_cl_XXXX_YYYY_psd_ave_ns.txt'
    name_intxt_csd_def='ni1200120106_0mpu7_cl_XXXX_YYYY_ZZZZ_WWWW_csd.txt'
    display=True
    ##### Setting (End) #####

    ##### Check (Start) #####
    name_intxt_psd_1=name_intxt_psd_def.replace('XXXX', '{:04}'.format(ch_min_1))
    name_intxt_psd_1=name_intxt_psd_1.replace('YYYY', '{:04}'.format(ch_max_1))
    name_intxt_psd_2=name_intxt_psd_def.replace('XXXX', '{:04}'.format(ch_min_2))
    name_intxt_psd_2=name_intxt_psd_2.replace('YYYY', '{:04}'.format(ch_max_2))
    name_intxt_csd=name_intxt_csd_def.replace('XXXX', '{:04}'.format(ch_min_1))
    name_intxt_csd=name_intxt_csd.replace('YYYY', '{:04}'.format(ch_max_1))
    name_intxt_csd=name_intxt_csd.replace('ZZZZ', '{:04}'.format(ch_min_2))
    name_intxt_csd=name_intxt_csd.replace('WWWW', '{:04}'.format(ch_max_2))
    name_outtxt=name_intxt_csd.replace('_csd', '_coh')
    name_outpdf=name_outtxt.replace(EXTENSION_TXT, EXTENSION_PDF)
    check_existence(name_intxt_psd_1)
    check_existence(name_intxt_psd_2)
    check_existence(name_intxt_csd)
    check_extension(name_intxt_psd_1, EXTENSION_TXT)
    check_extension(name_intxt_psd_2, EXTENSION_TXT)
    check_extension(name_intxt_csd, EXTENSION_TXT)
    check_extension(name_outtxt, EXTENSION_TXT)
    check_extension(name_outpdf, EXTENSION_PDF)
    ##### Check (End) #####

    ##### Read PSD (Start) #####
    f_psd_1,\
    psd_raw_mean_1,\
    psd_raw_sigma_1,\
    psd_noise_1,\
    n_mer_p_1,\
    n_mer_f_psd_1,\
    cr_mean_1,\
    exposure_1\
    =read_psd_ave_ns(name_intxt_psd_1)

    f_psd_2,\
    psd_raw_mean_2,\
    psd_raw_sigma_2,\
    psd_noise_2,\
    n_mer_p_2,\
    n_mer_f_psd_2,\
    cr_mean_2,\
    exposure_2\
    =read_psd_ave_ns(name_intxt_psd_2)
    ##### Read PSD (End) #####

    ##### Read CSD (Start) #####
    f_csd,\
    csd_mean,\
    csd_sigma,\
    n_mer_c,\
    n_mer_f_csd\
    =read_csd(name_intxt=name_intxt_csd)
    ##### Read CSD (End) #####

    ##### Check (Start) #####
    if len(f_psd_1)!=len(f_psd_2) or len(f_psd_1)!=len(f_csd):
        print('Error')
        sys.exit()
    if np.any(n_mer_p_1!=n_mer_p_2) or np.any(n_mer_p_1!=n_mer_c):
        print('Error')
        sys.exit()
    if np.any(n_mer_f_psd_1!=n_mer_f_psd_2) or np.any(n_mer_f_psd_1!=n_mer_f_csd):
        print('Error')
        sys.exit()

    f=f_psd_1
    n_mer_i=n_mer_p_1
    n_mer_f=n_mer_f_psd_1
    ##### Check (End) #####

    ##### Calculate coherence (Start) #####
    coh=coh_calc(psd_raw_1=psd_raw_mean_1,\
                 psd_noise_1=psd_noise_1,\
                 psd_raw_2=psd_raw_mean_2,\
                 psd_noise_2=psd_noise_2,\
                 csd=csd_mean,\
                 n_mer_i=n_mer_i,\
                 n_mer_f=n_mer_f)
    ##### Calculate coherence (End) #####

    ##### Write coherence (Start) #####
    write_coh(name_outtxt=name_outtxt,\
              fs=f,\
              cohs=coh,\
              ns_mer_i=n_mer_i,\
              ns_mer_f=n_mer_f)
    ##### Write coherence (End) #####

    ##### Plot coherence (Start) #####
    coh_sigma=np.zeros(len(coh)) #how to evaluate error of coherence?
    f_plx,\
    f_ply,\
    coh_pl,\
    f_errx,\
    f_erry,\
    coh_err,\
    coh_err\
    =params_adjust_plot(f,\
                        coh,\
                        coh_sigma)

    out_pdf=PdfPages(name_outpdf)
    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(f_plx,\
            coh_pl,\
            drawstyle='steps-post',\
            color='black',\
            linestyle='solid',\
            linewidth=1.2,\
            alpha=0.8)
    '''
    ax.errorbar(f_errx,\
                f_erry*csd_mean_err,\
                yerr=f_erry*csd_sigma_err,\
                capsize=0,\
                color='black',\
                linestyle='None',\
                linewidth=1.2,\
                alpha=0.8)
    '''
    #ax.set_xlim(0, int(num_tot-1))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    #ax.set_ylim(10**(-5), 10**(-1))
    #ax.set_yticks(np.linspace(0, 700, 8))
    #ax.set_yticks(np.linspace(0, 700, 15), minor=True)
    ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.set_xlabel('Frequency (Hz)', fontsize=18)
    ax.set_ylabel('Coherence (-)', fontsize=18)
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
    #ax.set_title('Number of cross spectra: {0}'.format(n_c), fontsize=18)
    #plt.show()

    out_pdf.savefig()
    out_pdf.close()
    ########## Plot CSD (Start) ##########

if __name__=='__main__':
    main()
