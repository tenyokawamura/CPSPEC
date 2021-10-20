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
    coi_min=101
    coi_max=1000
    ref_min=51
    ref_max=100
    name_intxt_def='ni1200120106_0mpu7_cl_XXXX_YYYY_ZZZZ_WWWW_csd_f.txt'
    ##### Setting (End) #####

    ##### Check (Start) #####
    name_intxt=name_intxt_def.replace('XXXX', '{:04}'.format(coi_min))
    name_intxt=name_intxt.replace('YYYY', '{:04}'.format(coi_max))
    name_intxt=name_intxt.replace('ZZZZ', '{:04}'.format(ref_min))
    name_intxt=name_intxt.replace('WWWW', '{:04}'.format(ref_max))

    name_outpdf=name_intxt.replace(EXTENSION_TXT, '_draw{0}'.format(EXTENSION_PDF))
    check_existence(name_intxt)
    check_extension(name_intxt, EXTENSION_TXT)
    check_extension(name_outpdf, EXTENSION_PDF)
    out_pdf=PdfPages(name_outpdf)
    ##### Check (End) #####

    ##### Read CSD (Start) #####
    fs_mea,\
    csds_re_mea,\
    csds_im_mea,\
    csds_ab_mea,\
    csds_ph_mea,\
    csds_ti_mea,\
    csds_re_sig,\
    csds_im_sig,\
    csds_ab_sig,\
    csds_ph_sig,\
    csds_ti_sig,\
    gamma2s,\
    n_mer_p,\
    n_mer_f\
    =read_csd_f(name_intxt=name_intxt)
    n_c=int(n_mer_p[0])
    ##### Read CSD (End) #####

    ##### Plot CSD (Start) #####
    ### Re ###
    f_plx,\
    f_ply,\
    csd_mea_pl,\
    f_errx,\
    f_erry,\
    csd_mea_err,\
    csd_sig_err\
    =params_adjust_plot(fs_mea,\
                        csds_re_mea,\
                        csds_re_sig)

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(f_plx,\
            f_ply*csd_mea_pl,\
            drawstyle='steps-post',\
            color='red',\
            linestyle='solid',\
            linewidth=1.2,\
            alpha=0.8)
    ax.errorbar(f_errx,\
                f_erry*csd_mea_err,\
                yerr=f_erry*csd_sig_err,\
                capsize=0,\
                color='red',\
                linestyle='None',\
                linewidth=1.2,\
                alpha=0.8)
    #ax.set_xlim(0, int(num_tot-1))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    ax.set_ylim(1.*10**(-7), 1.*10**(-1))
    #ax.set_yticks(np.linspace(0, 700, 8))
    #ax.set_yticks(np.linspace(0, 700, 15), minor=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Frequency (Hz)', fontsize=18)
    ax.set_ylabel('Re[νC(ν)] ((rms/mean)$^2$)', fontsize=18)
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
    ax.set_title('Number of cross spectra: {0}'.format(n_c), fontsize=18)
    #plt.show()
    out_pdf.savefig()
    plt.close()

    ### Im ###
    f_plx,\
    f_ply,\
    csd_mea_pl,\
    f_errx,\
    f_erry,\
    csd_mea_err,\
    csd_sig_err\
    =params_adjust_plot(fs_mea,\
                        csds_im_mea,\
                        csds_im_sig)

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(f_plx,\
            f_ply*csd_mea_pl,\
            drawstyle='steps-post',\
            color='red',\
            linestyle='solid',\
            linewidth=1.2,\
            alpha=0.8)
    ax.errorbar(f_errx,\
                f_erry*csd_mea_err,\
                yerr=f_erry*csd_sig_err,\
                capsize=0,\
                color='red',\
                linestyle='None',\
                linewidth=1.2,\
                alpha=0.8)
    #ax.set_xlim(0, int(num_tot-1))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    ax.set_ylim(1.*10**(-7), 1.*10**(-1))
    #ax.set_yticks(np.linspace(0, 700, 8))
    #ax.set_yticks(np.linspace(0, 700, 15), minor=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Frequency (Hz)', fontsize=18)
    ax.set_ylabel('Im[νC(ν)] ((rms/mean)$^2$)', fontsize=18)
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
    ax.set_title('Number of cross spectra: {0}'.format(n_c), fontsize=18)
    #plt.show()
    out_pdf.savefig()
    plt.close()

    ### Abs ###
    f_plx,\
    f_ply,\
    csd_mea_pl,\
    f_errx,\
    f_erry,\
    csd_mea_err,\
    csd_sig_err\
    =params_adjust_plot(fs_mea,\
                        csds_ab_mea,\
                        csds_ab_sig) #needed to be fixed

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(f_plx,\
            f_ply*csd_mea_pl,\
            drawstyle='steps-post',\
            color='red',\
            linestyle='solid',\
            linewidth=1.2,\
            alpha=0.8)
    ax.errorbar(f_errx,\
                f_erry*csd_mea_err,\
                yerr=f_erry*csd_sig_err,\
                capsize=0,\
                color='red',\
                linestyle='None',\
                linewidth=1.2,\
                alpha=0.8)
    #ax.set_xlim(0, int(num_tot-1))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    ax.set_ylim(1.*10**(-7), 1.*10**(-1))
    #ax.set_yticks(np.linspace(0, 700, 8))
    #ax.set_yticks(np.linspace(0, 700, 15), minor=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Frequency (Hz)', fontsize=18)
    ax.set_ylabel('ν|C(ν)| ((rms/mean)$^2$)', fontsize=18)
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
    ax.set_title('Number of cross spectra: {0}'.format(n_c), fontsize=18)
    #plt.show()
    out_pdf.savefig()
    plt.close()

    ### Phase ###
    f_plx,\
    f_ply,\
    csd_mea_pl,\
    f_errx,\
    f_erry,\
    csd_mea_err,\
    csd_sig_err\
    =params_adjust_plot(fs_mea,\
                        csds_ph_mea,\
                        csds_ph_sig) 

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(f_plx,\
            csd_mea_pl,\
            drawstyle='steps-post',\
            color='red',\
            linestyle='solid',\
            linewidth=1.2,\
            alpha=0.8)
    ax.errorbar(f_errx,\
                csd_mea_err,\
                yerr=csd_sig_err,\
                capsize=0,\
                color='red',\
                linestyle='None',\
                linewidth=1.2,\
                alpha=0.8)
    #ax.set_xlim(0, int(num_tot-1))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    #ax.set_ylim(10**(-5), 10**(-1))
    #ax.set_yticks(np.linspace(0, 700, 8))
    #ax.set_yticks(np.linspace(0, 700, 15), minor=True)
    ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.set_xlabel('Frequency (Hz)', fontsize=18)
    ax.set_ylabel('Phase lag (radian)', fontsize=18)
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
    ax.set_title('Number of cross spectra: {0}'.format(n_c), fontsize=18)
    #plt.show()
    out_pdf.savefig()
    plt.close()

    ### Time ###
    f_plx,\
    f_ply,\
    csd_mea_pl,\
    f_errx,\
    f_erry,\
    csd_mea_err,\
    csd_sig_err\
    =params_adjust_plot(fs_mea,\
                        csds_ti_mea,\
                        csds_ti_sig) 
    zeros=np.zeros(len(f_plx))

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(f_plx,\
            zeros,\
            color='gray',\
            linestyle='dashed',\
            linewidth=1.5,\
            alpha=0.8)
    ax.plot(f_plx,\
            csd_mea_pl*ONE2MILI,\
            drawstyle='steps-post',\
            color='red',\
            linestyle='solid',\
            linewidth=1.2,\
            alpha=0.8)
    ax.errorbar(f_errx,\
                csd_mea_err*ONE2MILI,\
                yerr=csd_sig_err*ONE2MILI,\
                capsize=0,\
                color='red',\
                linestyle='None',\
                linewidth=1.2,\
                alpha=0.8)
    ax.set_xlim(6.0*10**(-1), 1.0*10**(2))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    ax.set_ylim(-5.0, 15.0)
    #ax.set_yticks(np.linspace(-1, 3, 5))
    #ax.set_yticks(np.linspace(-1.6, 3.0, 24), minor=True)
    ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.set_xlabel('Frequency (Hz)', fontsize=18)
    ax.set_ylabel('Time lag (ms)', fontsize=18)
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
    ax.set_title('Number of cross spectra: {0}'.format(n_c), fontsize=18)
    #plt.show()
    out_pdf.savefig()
    plt.close()

    ### Intrinsic coherence ###
    gamma2s_sig=np.zeros(len(gamma2s))
    f_plx,\
    f_ply,\
    csd_mea_pl,\
    f_errx,\
    f_erry,\
    csd_mea_err,\
    csd_sig_err\
    =params_adjust_plot(fs_mea,\
                        gamma2s,\
                        gamma2s_sig) #needed to be modified
    ones=np.ones(len(f_plx))

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(f_plx,\
            ones,\
            color='gray',\
            linestyle='dashed',\
            linewidth=1.5,\
            alpha=0.8)
    ax.plot(f_plx,\
            csd_mea_pl,\
            drawstyle='steps-post',\
            color='red',\
            linestyle='solid',\
            linewidth=1.2,\
            alpha=0.8)
    ax.errorbar(f_errx,\
                csd_mea_err,\
                yerr=csd_sig_err,\
                capsize=0,\
                color='red',\
                linestyle='None',\
                linewidth=1.2,\
                alpha=0.8)
    #ax.set_xlim(1.5*10**(0), 8.5*10**(1))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    ax.set_ylim(0, 2.0)
    #ax.set_yticks(np.linspace(-1, 3, 5))
    #ax.set_yticks(np.linspace(-1.6, 3.0, 24), minor=True)
    ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.set_xlabel('Frequency (Hz)', fontsize=18)
    ax.set_ylabel('Intrinsic coherence (-)', fontsize=18)
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
    ax.set_title('Number of cross spectra: {0}'.format(n_c), fontsize=18)
    #plt.show()
    out_pdf.savefig()
    out_pdf.close()
    ########## Plot CSD (End) ##########

if __name__=='__main__':
    main()
