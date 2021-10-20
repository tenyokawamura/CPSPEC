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
    if not len(sys.argv)==6:
        print('Error: The number of inputs is wrong(0).(must be {1})'.format(len(sys.argv), 7))
        sys.exit()
    name_inevt=argv[1]
    # 1: Channel-of-interest
    # 2: Reference
    ch_min_coi=int(argv[2])
    ch_max_coi=int(argv[3])
    ch_min_ref=int(argv[4])
    ch_max_ref=int(argv[5])
    name_intxt_coi=name_inevt.replace(EXTENSION_EVT,\
                                    '_{0:04}_{1:04}_fft{2}'.format(ch_min_coi, ch_max_coi, EXTENSION_TXT))
    name_intxt_psd_coi=name_inevt.replace(EXTENSION_EVT,\
                                        '_{0:04}_{1:04}_psd_ave_ns{2}'.format(ch_min_coi, ch_max_coi, EXTENSION_TXT))
    name_intxt_ref=name_inevt.replace(EXTENSION_EVT,\
                                    '_{0:04}_{1:04}_fft{2}'.format(ch_min_ref, ch_max_ref, EXTENSION_TXT))
    name_intxt_psd_ref=name_inevt.replace(EXTENSION_EVT,\
                                        '_{0:04}_{1:04}_psd_ave_ns{2}'.format(ch_min_ref, ch_max_ref, EXTENSION_TXT))
    name_outtxt=name_inevt.replace(EXTENSION_EVT,\
                                   '_{0:04}_{1:04}_{2:04}_{3:04}_csd_f{4}'\
                                   .format(ch_min_coi,\
                                           ch_max_coi,\
                                           ch_min_ref,\
                                           ch_max_ref,\
                                           EXTENSION_TXT))
    name_outpdf=name_outtxt.replace(EXTENSION_TXT, EXTENSION_PDF)
    rebin=REBIN
    if ch_min_ref<=ch_min_coi and ch_max_coi<=ch_max_ref:
        overlap=True
    else:
        overlap=False
    ##### Setting (End) #####

    ##### Check (Start) #####
    check_existence(name_intxt_coi)
    check_existence(name_intxt_psd_coi)
    check_existence(name_intxt_ref)
    check_existence(name_intxt_psd_ref)
    check_extension(name_intxt_coi, EXTENSION_TXT)
    check_extension(name_intxt_ref, EXTENSION_TXT)
    check_extension(name_outtxt, EXTENSION_TXT)
    check_extension(name_outpdf, EXTENSION_PDF)
    out_pdf=PdfPages(name_outpdf)
    ##### Check (End) #####

    ##### Read file (End) #####
    telescope,\
    instrument,\
    source,\
    exposure,\
    resol_t,\
    n_row\
    =read_obs_info(name_intxt=name_intxt_coi)
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
    =read_ana_info(name_intxt=name_intxt_coi)
    print_ana_info(dt=dt,\
                   n_d=n_d,\
                   lar_t=lar_t,\
                   n_i_max=n_i_max,\
                   n_i_exp=n_i_exp,\
                   f_min=f_min,\
                   f_max=f_max)
    ##### Read file (End) #####

    ##### Calculate CSD (Start) #####
    i_i=0
    while True: 
        print('(No. {0})'.format(i_i+1))
        f_coi,\
        b_coi,\
        t_s_coi,\
        t_e_coi,\
        cr_mean_coi,\
        flag_end_coi\
        =read_ft(name_intxt=name_intxt_coi,\
                 no_i=i_i+1)

        f_ref,\
        b_ref,\
        t_s_ref,\
        t_e_ref,\
        cr_mean_ref,\
        flag_end_ref\
        =read_ft(name_intxt=name_intxt_ref,\
                 no_i=i_i+1)

        csd=csd_calc(bs_coi=b_coi,\
                     bs_ref=b_ref,\
                     x_mean_coi=cr_mean_coi,\
                     x_mean_ref=cr_mean_ref,\
                     n_d=n_d,\
                     dt=dt,\
                     overlap=overlap)
        if i_i==0:
            f_buf=f_coi
            csds=csd
        else:
            csds=np.vstack((csds, csd))

        i_i=i_i+1
        if flag_end_coi==True and flag_end_ref==True:
            break
    ##### Calculate CSD (End) #####

    ##### Post-processing (Start) #####
    n_c=i_i
    f=f_buf
    if n_c>1:
        print('{0} cross spectra are calculated successfully.'.format(n_c))
        ### Read power spectra (Start) ###
        f_psd_coi,\
        psd_raw_mea_coi,\
        psd_raw_sig_coi,\
        psd_noi_coi,\
        n_mer_p_coi,\
        n_mer_f_psd_coi,\
        cr_mean_coi,\
        exposure_coi\
        =read_psd_ave_ns(name_intxt_psd_coi)

        f_psd_ref,\
        psd_raw_mea_ref,\
        psd_raw_sig_ref,\
        psd_noi_ref,\
        n_mer_p_ref,\
        n_mer_f_psd_ref,\
        cr_mean_ref,\
        exposure_ref\
        =read_psd_ave_ns(name_intxt_psd_ref)

        f_psd=f_psd_coi
        ### Read power spectra (End) ###

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
        ns_mer_c,\
        ns_mer_f\
        =csd_f_process(fs=f,\
                       csds=csds,\
                       n_c=n_c,\
                       rebin=rebin,\
                       fs_psd=f_psd,\
                       psds_raw_ref=psd_raw_mea_coi,\
                       psds_noi_ref=psd_noi_coi,\
                       psds_raw_coi=psd_raw_mea_ref,\
                       psds_noi_coi=psd_noi_ref)

    elif n_c==1:
        print('Error: Under construction.')
        sys.exit()

    else:
        print('Error: No cross spectrum was calculated successfully.')
        sys.exit()
    ##### Post-processing (End) #####

    ##### Write CSD (Start) #####
    write_csd_f(name_outtxt=name_outtxt,\
                fs=fs_mea,\
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
                ns_mer_c=ns_mer_c,\
                ns_mer_f=ns_mer_f)
    ##### Write CSD (Start) #####

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
                        csds_re_sig) #needed to be fixed

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(f_plx,\
            csd_mea_pl,\
            drawstyle='steps-post',\
            color='black',\
            linestyle='solid',\
            linewidth=1.2,\
            alpha=0.8)
    ax.errorbar(f_errx,\
                csd_mea_err,\
                yerr=csd_sig_err,\
                capsize=0,\
                color='black',\
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
    ax.set_yscale('log')
    ax.set_xlabel('Frequency (Hz)', fontsize=18)
    ax.set_ylabel('Re[C(ν)] ((rms/mean)$^2$)/Hz', fontsize=18)
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
                        csds_im_sig) #needed to be fixed

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(f_plx,\
            csd_mea_pl,\
            drawstyle='steps-post',\
            color='black',\
            linestyle='solid',\
            linewidth=1.2,\
            alpha=0.8)
    ax.errorbar(f_errx,\
                csd_mea_err,\
                yerr=csd_sig_err,\
                capsize=0,\
                color='black',\
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
    ax.set_yscale('log')
    ax.set_xlabel('Frequency (Hz)', fontsize=18)
    ax.set_ylabel('Im[C(ν)] ((rms/mean)$^2$/Hz)', fontsize=18)
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
            color='black',\
            linestyle='solid',\
            linewidth=1.2,\
            alpha=0.8)
    ax.errorbar(f_errx,\
                f_erry*csd_mea_err,\
                yerr=f_erry*csd_sig_err,\
                capsize=0,\
                color='black',\
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
                        csds_ph_sig) #needed to be fixed

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(f_plx,\
            csd_mea_pl,\
            drawstyle='steps-post',\
            color='black',\
            linestyle='solid',\
            linewidth=1.2,\
            alpha=0.8)
    ax.errorbar(f_errx,\
                csd_mea_err,\
                yerr=csd_sig_err,\
                capsize=0,\
                color='black',\
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
    taus_mea=csds_ph_mea/(2.*np.pi*fs_mea)
    taus_sig=csds_ph_sig/(2.*np.pi*fs_mea)

    f_plx,\
    f_ply,\
    csd_mea_pl,\
    f_errx,\
    f_erry,\
    csd_mea_err,\
    csd_sig_err\
    =params_adjust_plot(fs_mea,\
                        taus_mea,\
                        taus_sig) #needed to be fixed

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(f_plx,\
            csd_mea_pl,\
            drawstyle='steps-post',\
            color='black',\
            linestyle='solid',\
            linewidth=1.2,\
            alpha=0.8)
    ax.errorbar(f_errx,\
                csd_mea_err,\
                yerr=csd_sig_err,\
                capsize=0,\
                color='black',\
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
    ax.set_ylabel('Time lag (s)', fontsize=18)
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

    ### Coherence ###
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
                        gamma2s_sig)

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(f_plx,\
            csd_mea_pl,\
            drawstyle='steps-post',\
            color='black',\
            linestyle='solid',\
            linewidth=1.2,\
            alpha=0.8)
    ax.errorbar(f_errx,\
                csd_mea_err,\
                yerr=csd_sig_err,\
                capsize=0,\
                color='black',\
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
