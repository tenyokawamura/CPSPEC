from constants import *
from parameters import *
from functions import *
import numpy as np
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
                                   '_csd_e{0}'\
                                   .format(EXTENSION_TXT))
    name_outpdf=name_outtxt.replace(EXTENSION_TXT, EXTENSION_PDF)
    ##### Setting (End) #####

    ##### Check (Start) #####
    check_extension(name_outtxt, EXTENSION_TXT)
    check_extension(name_outpdf, EXTENSION_PDF)
    out_pdf=PdfPages(name_outpdf)
    ##### Check (End) #####

    n_b=len(COIS_MIN)
    kevs_mea=np.zeros(n_b)
    kevs_sig=np.zeros(n_b)
    csds_re_mea=np.zeros(n_b)
    csds_im_mea=np.zeros(n_b)
    csds_ab_mea=np.zeros(n_b)
    csds_ph_mea=np.zeros(n_b)
    csds_ti_mea=np.zeros(n_b)
    csds_re_sig=np.zeros(n_b)
    csds_im_sig=np.zeros(n_b)
    csds_ab_sig=np.zeros(n_b)
    csds_ph_sig=np.zeros(n_b)
    csds_ti_sig=np.zeros(n_b)
    gamma2s=np.zeros(n_b)
    ns_mer_i=np.zeros(n_b)
    ns_mer_f=np.zeros(n_b)
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
        check_existence(name_intxt_ref)
        check_existence(name_intxt_coi)
        check_existence(name_intxt_csd)
        check_extension(name_intxt_ref, EXTENSION_TXT)
        check_extension(name_intxt_coi, EXTENSION_TXT)
        check_extension(name_intxt_csd, EXTENSION_TXT)
        ##### Filename (End) #####

        ##### Convert CH into energy (Start) #####
        kev_min=ch2kev(ch=coi_min)
        kev_max=ch2kev(ch=coi_max)
        kev_mea=(kev_min+kev_max)/2.
        kev_sig=(kev_max-kev_min)/2.
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
        cr_mean_ref,\
        exposure_ref\
        =read_psd_ave_ns(name_intxt_ref)

        f_psd_coi,\
        psd_raw_mea_coi,\
        psd_raw_sig_coi,\
        psd_noi_coi,\
        n_mer_i_coi,\
        n_mer_f_coi,\
        cr_mean_coi,\
        exposure_coi\
        =read_psd_ave_ns(name_intxt_coi)
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

        ### Cut (Start) ###
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
                     f_min=F_MIN_LAGE,\
                     f_max=F_MAX_LAGE)
        ### Cut (End) ###

        ### Calculate CSD (Energy dependency) (Start) ###
        f_mea,\
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
        n_mer_i,\
        n_mer_f\
        =csd_e_process(fs=f_cut,\
                       psds_raw_ref=psd_raw_ref_cut,\
                       psds_noi_ref=psd_noi_ref_cut,\
                       psds_raw_coi=psd_raw_coi_cut,\
                       psds_noi_coi=psd_noi_coi_cut,\
                       csds=csd_raw_cut,\
                       ns_mer_i=n_mer_i_cut,\
                       ns_mer_f=n_mer_f_cut)
        csds_re_mea[i_b]=csd_re_mea
        csds_im_mea[i_b]=csd_im_mea
        csds_ab_mea[i_b]=csd_ab_mea
        csds_ph_mea[i_b]=csd_ph_mea
        csds_ti_mea[i_b]=csd_ti_mea
        csds_re_sig[i_b]=csd_re_sig
        csds_im_sig[i_b]=csd_im_sig
        csds_ab_sig[i_b]=csd_ab_sig
        csds_ph_sig[i_b]=csd_ph_sig
        csds_ti_sig[i_b]=csd_ti_sig
        gamma2s[i_b]=gamma2
        ns_mer_i[i_b]=n_mer_i
        ns_mer_f[i_b]=n_mer_f
        ### Calculate CSD (Energy dependency) (End) ###

    ##### Write lag-frequency spectrum (Start) #####
    write_csd_e(name_outtxt=name_outtxt,\
                ref_min=REF_MIN,\
                ref_max=REF_MAX,\
                cois_min=COIS_MIN,\
                cois_max=COIS_MAX,\
                f_min=F_MIN_LAGE,\
                f_max=F_MAX_LAGE,\
                csds_re_mea=csds_re_mea,\
                csds_im_mea=csds_im_mea,\
                csds_ab_mea=csds_ab_mea,\
                csds_ph_mea=csds_ph_mea,\
                csds_ti_mea=csds_ti_mea,\
                csds_re_sig=csds_re_sig,\
                csds_im_sig=csds_im_sig,\
                csds_ab_sig=csds_ab_sig,\
                csds_ph_sig=csds_ph_sig,\
                csds_ti_sig=csds_ti_sig,\
                gamma2s=gamma2s,\
                ns_mer_i=ns_mer_i,\
                ns_mer_f=ns_mer_f)
    ##### Write lag-frequency spectrum (End) #####

    print(f_mea)
    ##### Plot Re[CSD]-energy spectrum (Start) #####
    #f_plx,\
    #f_ply,\
    #tau_mean_pl,\
    #f_errx,\
    #f_erry,\
    #tau_mean_err,\
    #tau_sigma_err\
    #=params_adjust_plot(f,\
    #                    tau_mean,\
    #                    tau_sigma)

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    #fig=plt.figure(figsize=(6, 6))
    ax=fig.add_subplot(1, 1, 1)
    #ax.plot(f_plx,\
    #        tau_mean_pl*ONE2MILI,\
    #        drawstyle='steps-post',\
    #        color='red',\
    #        linestyle='solid',\
    #        linewidth=1.5,\
    #        alpha=0.8)
    ax.errorbar(x=kevs_mea,\
                y=csds_re_mea,\
                xerr=kevs_sig,\
                yerr=csds_re_sig,\
                marker='o',\
                markersize=10.0,\
                capsize=0,\
                color='black',\
                linestyle='None',\
                linewidth=1.5,\
                alpha=0.8)
    #ax.set_xlim(1., 90)
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    #ax.set_ylim(-4, 4.0)
    #ax.set_yticks(np.linspace(0, 700, 8))
    #ax.set_yticks(np.linspace(0, 700, 15), minor=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Energy (keV)', fontsize=18)
    ax.set_ylabel('Re[CSD] ((rms/mean)$^2$)/Hz', fontsize=18)
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
    ax.set_title('Frequency: {0:.1e} Hz - {1:.1e} Hz'.format(F_MIN_LAGE, F_MAX_LAGE), fontsize=18)
    out_pdf.savefig()
    plt.close()
    ##### Plot Re[CSD]-energy spectrum (End) #####

    ##### Plot Im[CSD]-energy spectrum (Start) #####
    #f_plx,\
    #f_ply,\
    #tau_mean_pl,\
    #f_errx,\
    #f_erry,\
    #tau_mean_err,\
    #tau_sigma_err\
    #=params_adjust_plot(f,\
    #                    tau_mean,\
    #                    tau_sigma)

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    #fig=plt.figure(figsize=(6, 6))
    ax=fig.add_subplot(1, 1, 1)
    #ax.plot(f_plx,\
    #        tau_mean_pl*ONE2MILI,\
    #        drawstyle='steps-post',\
    #        color='red',\
    #        linestyle='solid',\
    #        linewidth=1.5,\
    #        alpha=0.8)
    ax.errorbar(x=kevs_mea,\
                y=csds_im_mea,\
                xerr=kevs_sig,\
                yerr=csds_im_sig,\
                marker='o',\
                markersize=10.0,\
                capsize=0,\
                color='black',\
                linestyle='None',\
                linewidth=1.5,\
                alpha=0.8)
    #ax.set_xlim(1., 90)
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    #ax.set_ylim(-4, 4.0)
    #ax.set_yticks(np.linspace(0, 700, 8))
    #ax.set_yticks(np.linspace(0, 700, 15), minor=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Energy (keV)', fontsize=18)
    ax.set_ylabel('Im[CSD] ((rms/mean)$^2$)/Hz', fontsize=18)
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
    ax.set_title('Frequency: {0:.1e} Hz - {1:.1e} Hz'.format(F_MIN_LAGE, F_MAX_LAGE), fontsize=18)
    out_pdf.savefig()
    plt.close()
    ##### Plot Im[CSD]-energy spectrum (End) #####

    ##### Plot Abs[CSD]-energy spectrum (Start) #####
    #f_plx,\
    #f_ply,\
    #tau_mean_pl,\
    #f_errx,\
    #f_erry,\
    #tau_mean_err,\
    #tau_sigma_err\
    #=params_adjust_plot(f,\
    #                    tau_mean,\
    #                    tau_sigma)

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    #fig=plt.figure(figsize=(6, 6))
    ax=fig.add_subplot(1, 1, 1)
    #ax.plot(f_plx,\
    #        tau_mean_pl*ONE2MILI,\
    #        drawstyle='steps-post',\
    #        color='red',\
    #        linestyle='solid',\
    #        linewidth=1.5,\
    #        alpha=0.8)
    ax.errorbar(x=kevs_mea,\
                y=csds_ab_mea,\
                xerr=kevs_sig,\
                yerr=csds_ab_sig,\
                marker='o',\
                markersize=10.0,\
                capsize=0,\
                color='black',\
                linestyle='None',\
                linewidth=1.5,\
                alpha=0.8)
    #ax.set_xlim(1., 90)
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    #ax.set_ylim(-4, 4.0)
    #ax.set_yticks(np.linspace(0, 700, 8))
    #ax.set_yticks(np.linspace(0, 700, 15), minor=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Energy (keV)', fontsize=18)
    ax.set_ylabel('Abs[CSD] ((rms/mean)$^2$)/Hz', fontsize=18)
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
    ax.set_title('Frequency: {0:.1e} Hz - {1:.1e} Hz'.format(F_MIN_LAGE, F_MAX_LAGE), fontsize=18)
    out_pdf.savefig()
    plt.close()
    ##### Plot Abs[CSD]-energy spectrum (End) #####

    ##### Plot phase lag-energy spectrum (Start) #####
    '''
    f_plx,\
    f_ply,\
    pha_mean_pl,\
    f_errx,\
    f_erry,\
    pha_mean_err,\
    pha_sigma_err\
    =params_adjust_plot(f,\
                        pha_mean,\
                        pha_sigma)
    '''

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    #ax.plot(f_plx,\
    #        pha_mean_pl,\
    #        drawstyle='steps-post',\
    #        color='black',\
    #        linestyle='solid',\
    #        linewidth=1.2,\
    #        alpha=0.8)
    ax.errorbar(x=kevs_mea,\
                y=csds_ph_mea,\
                xerr=kevs_sig,\
                yerr=csds_ph_sig,\
                marker='o',\
                markersize=10.0,\
                capsize=0,\
                color='black',\
                linestyle='None',\
                linewidth=1.5,\
                alpha=0.8)
    #ax.set_xlim(1, 80)
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    #ax.set_ylim(-1, 1)
    #ax.set_yticks(np.linspace(0, 700, 8))
    #ax.set_yticks(np.linspace(0, 700, 15), minor=True)
    ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.set_xlabel('Energy (keV)', fontsize=18)
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
    ax.set_title('Frequency: {0:.1e} Hz - {1:.1e} Hz'.format(F_MIN_LAGE, F_MAX_LAGE), fontsize=18)
    out_pdf.savefig()
    plt.close()
    ##### Plot phase lag-energy spectrum (End) #####

    ##### Plot time lag-energy spectrum (Start) #####
    #f_plx,\
    #f_ply,\
    #tau_mean_pl,\
    #f_errx,\
    #f_erry,\
    #tau_mean_err,\
    #tau_sigma_err\
    #=params_adjust_plot(f,\
    #                    tau_mean,\
    #                    tau_sigma)

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    #fig=plt.figure(figsize=(6, 6))
    ax=fig.add_subplot(1, 1, 1)
    #ax.plot(f_plx,\
    #        tau_mean_pl*ONE2MILI,\
    #        drawstyle='steps-post',\
    #        color='red',\
    #        linestyle='solid',\
    #        linewidth=1.5,\
    #        alpha=0.8)
    ax.errorbar(x=kevs_mea,\
                y=csds_ti_mea*ONE2MILI,\
                xerr=kevs_sig,\
                yerr=csds_ti_sig*ONE2MILI,\
                marker='o',\
                markersize=10.0,\
                capsize=0,\
                color='black',\
                linestyle='None',\
                linewidth=1.5,\
                alpha=0.8)
    #ax.set_xlim(1., 90)
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    #ax.set_ylim(-4, 4.0)
    #ax.set_yticks(np.linspace(0, 700, 8))
    #ax.set_yticks(np.linspace(0, 700, 15), minor=True)
    ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.set_xlabel('Energy (keV)', fontsize=18)
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
    ax.set_title('Frequency: {0:.1e} Hz - {1:.1e} Hz'.format(F_MIN_LAGE, F_MAX_LAGE), fontsize=18)
    out_pdf.savefig()
    plt.close()
    ##### Plot time lag-energy spectrum (End) #####

    ##### Plot Intrinsic coherence spectrum (Start) #####
    #f_plx,\
    #f_ply,\
    #tau_mean_pl,\
    #f_errx,\
    #f_erry,\
    #tau_mean_err,\
    #tau_sigma_err\
    #=params_adjust_plot(f,\
    #                    tau_mean,\
    #                    tau_sigma)

    gamma2s_sig=np.zeros(len(gamma2s))

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    #fig=plt.figure(figsize=(6, 6))
    ax=fig.add_subplot(1, 1, 1)
    #ax.plot(f_plx,\
    #        tau_mean_pl*ONE2MILI,\
    #        drawstyle='steps-post',\
    #        color='red',\
    #        linestyle='solid',\
    #        linewidth=1.5,\
    #        alpha=0.8)
    ax.errorbar(x=kevs_mea,\
                y=gamma2s,\
                xerr=kevs_sig,\
                yerr=gamma2s_sig,\
                marker='o',\
                markersize=10.0,\
                capsize=0,\
                color='black',\
                linestyle='None',\
                linewidth=1.5,\
                alpha=0.8)
    #ax.set_xlim(1., 90)
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    #ax.set_ylim(-4, 4.0)
    #ax.set_yticks(np.linspace(0, 700, 8))
    #ax.set_yticks(np.linspace(0, 700, 15), minor=True)
    ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.set_xlabel('Energy (keV)', fontsize=18)
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
    ax.set_title('Frequency: {0:.1e} Hz - {1:.1e} Hz'.format(F_MIN_LAGE, F_MAX_LAGE), fontsize=18)
    out_pdf.savefig()
    plt.close()
    out_pdf.close()
    ##### Plot Intrinsic coherence spectrum (End) #####

if __name__=='__main__':
    main()
