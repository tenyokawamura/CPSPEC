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
    name_inevt=NAME_INEVT
    name_intxt_ref=name_inevt.replace(EXTENSION_EVT,\
                                      '_{0:04}_{1:04}_psd_ave_ns{2}'\
                                      .format(REF_MIN, REF_MAX, EXTENSION_TXT))
    name_outtxt=name_inevt.replace(EXTENSION_EVT,\
                                   '_rms{0}'\
                                   .format(EXTENSION_TXT))
    name_outpdf=name_outtxt.replace(EXTENSION_TXT, EXTENSION_PDF)
    ##### Setting (End) #####

    ##### Check (Start) #####
    check_extension(name_outtxt, EXTENSION_TXT)
    check_extension(name_outpdf, EXTENSION_PDF)
    out_pdf=PdfPages(name_outpdf)
    ##### Check (End) #####

    n_b=len(COIS_MIN)
    kevs_bin=np.zeros(n_b)
    kevs_mea=np.zeros(n_b)
    kevs_sig=np.zeros(n_b)
    rmses_per_mea=np.zeros(n_b)
    rmses_per_sig=np.zeros(n_b)
    rmses_abs_mea=np.zeros(n_b)
    rmses_abs_sig=np.zeros(n_b)
    crs_mea=np.zeros(n_b)
    crs_sig=np.zeros(n_b)
    for i_b, (coi_min, coi_max) in enumerate(zip(COIS_MIN, COIS_MAX)):
        ##### Filename (Start) #####
        name_intxt_coi=name_inevt.replace(EXTENSION_EVT,\
                                      '_{0:04}_{1:04}_psd_ave_ns{2}'\
                                      .format(coi_min, coi_max, EXTENSION_TXT))
        name_intxt_csd=name_inevt.replace(EXTENSION_EVT,\
                                          '_{0:04}_{1:04}_{2:04}_{3:04}_csd_f{4}'\
                                          .format(coi_min,\
                                                  coi_max,\
                                                  REF_MIN,\
                                                  REF_MAX,\
                                                  EXTENSION_TXT))
        check_existence(name_intxt_coi)
        check_existence(name_intxt_csd)
        check_extension(name_intxt_coi, EXTENSION_TXT)
        check_extension(name_intxt_csd, EXTENSION_TXT)
        ##### Filename (End) #####

        ##### Convert CH into energy (Start) #####
        kev_min=ch2kev(ch=coi_min)
        kev_max=ch2kev(ch=coi_max)
        kev_bin=kev_max-kev_min
        kev_mea=(kev_min+kev_max)/2.
        kev_sig=(kev_max-kev_min)/2.
        kevs_bin[i_b]=kev_bin
        kevs_mea[i_b]=kev_mea
        kevs_sig[i_b]=kev_sig
        ##### Convert CH into energy (End) #####

        ##### Read PSD (Start) #####
        f_psd_ref,\
        psd_raw_mea_ref,\
        psd_raw_sig_ref,\
        psd_noi_ref,\
        n_mer_i_ref,\
        n_mer_f_ref,\
        cr_mea_ref,\
        exposure_ref\
        =read_psd_ave_ns(name_intxt_ref)

        f_psd_coi,\
        psd_raw_mea_coi,\
        psd_raw_sig_coi,\
        psd_noi_coi,\
        n_mer_i_coi,\
        n_mer_f_coi,\
        cr_mea_coi,\
        exposure_coi\
        =read_psd_ave_ns(name_intxt_coi)

        cr_sig_coi=cr_sig_calc(cr_mea=cr_mea_coi,\
                               t=exposure_coi)
        crs_mea[i_b]=cr_mea_coi
        crs_sig[i_b]=cr_sig_coi
        ##### Read PSD (End) #####

        ##### Read CSD (Start) #####
        f_csd,\
        csd_re_mea,\
        csd_im_mea,\
        csd_ab_mea,\
        csd_ph_mea,\
        csd_ti_mea,\
        csd_re_sig,\
        csd_im_sig,\
        csd_ab_sig,\
        csd_ph_sig,\
        csd_ti_sig,\
        gamma2,\
        n_mer_i_csd,\
        n_mer_f_csd\
        =read_csd_f(name_intxt=name_intxt_csd)
        csd_raw_mea=csd_re_mea+1j*csd_im_mea
        ##### Read CSD (End) #####

        ##### Cut (Start) #####
        n_mer_i_cut,\
        n_mer_f_cut,\
        f_cut,\
        psd_raw_ref_cut,\
        psd_noi_ref_cut,\
        psd_raw_coi_cut,\
        psd_noi_coi_cut,\
        csd_raw_cut\
        =cut_psd_csd(f_ref=f_psd_ref,\
                     psd_raw_ref=psd_raw_mea_ref,\
                     psd_noi_ref=psd_noi_ref,\
                     n_mer_i_ref=n_mer_i_ref,\
                     n_mer_f_ref=n_mer_f_ref,\
                     f_coi=f_psd_coi,\
                     psd_raw_coi=psd_raw_mea_coi,\
                     psd_noi_coi=psd_noi_coi,\
                     n_mer_i_coi=n_mer_i_coi,\
                     n_mer_f_coi=n_mer_f_coi,\
                     f_csd=f_csd,\
                     csd_raw=csd_raw_mea,\
                     n_mer_i_csd=n_mer_i_csd,\
                     n_mer_f_csd=n_mer_f_csd,\
                     f_min=F_MIN_RMS,\
                     f_max=F_MAX_RMS)
        psd_coi_cut=psd_raw_coi_cut-psd_noi_coi_cut
        ##### Cut (End) #####

        ##### Calculate rms (Start) #####
        rms_per_mea,\
        rms_per_sig,\
        rms_abs_mea,\
        rms_abs_sig\
        =rms_calc_strict(fs=f_cut,\
                         psds_raw_ref=psd_raw_ref_cut,\
                         psds_noi_ref=psd_noi_ref_cut,\
                         psds_raw_coi=psd_raw_coi_cut,\
                         psds_noi_coi=psd_noi_coi_cut,\
                         csds=csd_raw_cut,\
                         ns_mer_i=n_mer_i_cut,\
                         ns_mer_f=n_mer_f_cut,\
                         cr_mea_coi=cr_mea_coi)
        '''
        rms_per_mea,\
        rms_per_sig,\
        rms_abs_mea,\
        rms_abs_sig\
        =rms_calc(fs=f_cut,\
                  psds_raw=psd_raw_cut,\
                  psds_sig=psd_sig_cut,\
                  psds_noi=psd_noi_cut,\
                  cr_mea=cr_mea,\
                  ns_mer_i=n_mer_i_cut,\
                  ns_mer_f=n_mer_f_cut)
        '''
        rmses_per_mea[i_b]=rms_per_mea
        rmses_per_sig[i_b]=rms_per_sig
        rmses_abs_mea[i_b]=rms_abs_mea
        rmses_abs_sig[i_b]=rms_abs_sig
        ##### Calculate rms (Start) #####

    ##### Write rms (Start) #####
    write_rms(name_outtxt=name_outtxt,\
              cois_min=COIS_MIN,\
              cois_max=COIS_MAX,\
              f_min=F_MIN_RMS,\
              f_max=F_MAX_RMS,\
              rmses_per_mea=rmses_per_mea,\
              rmses_per_sig=rmses_per_sig,\
              rmses_abs_mea=rmses_abs_mea,\
              rmses_abs_sig=rmses_abs_sig,\
              crs_mea=crs_mea,\
              crs_sig=crs_sig,\
              exposure=exposure_coi)
    ##### Write rms (End) #####

    ##### Plot rms spectrum (Start) #####
    rmses_abs_mea_kev=rmses_abs_mea/kevs_bin
    rmses_abs_sig_kev=rmses_abs_sig/kevs_bin

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.errorbar(x=kevs_mea,\
                y=rmses_per_mea*100.,\
                xerr=kevs_sig,\
                yerr=rmses_per_sig*100.,\
                marker='o',\
                markersize=10.0,\
                capsize=0,\
                color='black',\
                linestyle='None',\
                linewidth=1.5,\
                alpha=0.8)
    #ax.set_xlim(0, int(num_tot-1))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    #ax.set_ylim(10**(-5), 10**(-1))
    #ax.set_yticks(np.linspace(0, 700, 8))
    #ax.set_yticks(np.linspace(0, 700, 15), minor=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Energy (keV)', fontsize=18)
    ax.set_ylabel('Fractional rms (%)', fontsize=18)
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
    ax.set_title('Frequency: {0:.1e} Hz - {1:.1e} Hz'.format(F_MIN_RMS, F_MAX_RMS), fontsize=18)
    #plt.show()
    out_pdf.savefig()
    plt.close()

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.errorbar(x=kevs_mea,\
                y=rmses_abs_mea_kev,\
                xerr=kevs_sig,\
                yerr=rmses_abs_sig_kev,\
                marker='o',\
                markersize=10.0,\
                capsize=0,\
                color='black',\
                linestyle='None',\
                linewidth=1.5,\
                alpha=0.8)
    #ax.set_xlim(0, int(num_tot-1))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    #ax.set_ylim(10**(-5), 10**(-1))
    #ax.set_yticks(np.linspace(0, 700, 8))
    #ax.set_yticks(np.linspace(0, 700, 15), minor=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Energy (keV)', fontsize=18)
    ax.set_ylabel('Absolute rms of flux (s$^{-1}$ keV$^{-1}$)', fontsize=18)
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
    ax.set_title('Frequency: {0:.1e} Hz - {1:.1e} Hz'.format(F_MIN_RMS, F_MAX_RMS), fontsize=18)
    #plt.show()
    out_pdf.savefig()
    out_pdf.close()
    ##### Plot rms spectrum (Start) #####

if __name__=='__main__':
    main()
