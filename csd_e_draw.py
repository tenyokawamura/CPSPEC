from constants import *
from parameters import *
from functions import *
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
import sys 
from sys import argv
import numpy as np
def main():
    ##### Display #####
    print_filename(argv[0])

    ##### Setting #####
    name_intxt='ni1200120106_0mpu7_cl_csd_e.txt'
    name_outpdf='ni1200120106_0mpu7_cl_csd_e_draw.pdf'
    check_existence(name_intxt)
    check_extension(name_intxt, EXTENSION_TXT)
    check_extension(name_outpdf, EXTENSION_PDF)
    out_pdf=PdfPages(name_outpdf)

    ##### Draw power spectra #####
    ##### Read data #####
    ref_min,\
    ref_max,\
    f_min,\
    f_max,\
    cois_min,\
    cois_max,\
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
    ns_mer_i,\
    ns_mer_f\
    =read_csd_e(name_intxt=name_intxt)

    ##### Convert CH into energy (Start) #####
    kevs_min=ch2kev(ch=cois_min)
    kevs_max=ch2kev(ch=cois_max)
    kevs_mea=(kevs_min+kevs_max)/2.
    kevs_sig=(kevs_max-kevs_min)/2.
    ##### Convert CH into energy (End) #####

    kevs_pl=np.append(kevs_min, kevs_max[-1])
    csds_re_mea_pl=np.append(csds_re_mea, csds_re_mea[-1])
    csds_im_mea_pl=np.append(csds_im_mea, csds_im_mea[-1])
    csds_ab_mea_pl=np.append(csds_ab_mea, csds_ab_mea[-1])
    csds_ph_mea_pl=np.append(csds_ph_mea, csds_ph_mea[-1])
    csds_ti_mea_pl=np.append(csds_ti_mea, csds_ti_mea[-1])
    gamma2s_pl=np.append(gamma2s, gamma2s[-1])

    ##### Plot Re[CSD]-energy spectra (Start) #####
    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(kevs_pl,\
            csds_re_mea_pl,\
            drawstyle='steps-post',\
            color='red',\
            linestyle='solid',\
            linewidth=2.0,\
            alpha=0.8)
    ax.errorbar(kevs_mea,\
                csds_re_mea,\
                yerr=csds_re_sig,\
                capsize=0,\
                color='red',\
                linestyle='None',\
                elinewidth=2.0,\
                linewidth=1.2,\
                alpha=0.8)
    #ax.set_xlim(1, 80)
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    #ax.set_ylim(1.*10**(-4), 2.*10**(-3))
    ax.set_ylim(1.*10**(-6), 2.*10**(-3))
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
    ##### Plot phase lag-energy spectra (End) #####

    ##### Plot Im[CSD]-energy spectra (Start) #####
    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(kevs_pl,\
            csds_im_mea_pl,\
            drawstyle='steps-post',\
            color='red',\
            linestyle='solid',\
            linewidth=2.0,\
            alpha=0.8)
    ax.errorbar(kevs_mea,\
                csds_im_mea,\
                yerr=csds_im_sig,\
                capsize=0,\
                color='red',\
                linestyle='None',\
                elinewidth=2.0,\
                linewidth=1.2,\
                alpha=0.8)
    #ax.set_xlim(1, 80)
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    #ax.set_ylim(1.*10**(-6), 2.*10**(-4))
    ax.set_ylim(1.*10**(-6), 2.*10**(-3))
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
    ##### Plot phase lag-energy spectra (End) #####

    ##### Plot Abs[CSD]-energy spectra (Start) #####
    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(kevs_pl,\
            csds_ab_mea_pl,\
            drawstyle='steps-post',\
            color='red',\
            linestyle='solid',\
            linewidth=2.0,\
            alpha=0.8)
    ax.errorbar(kevs_mea,\
                csds_ab_mea,\
                yerr=csds_ab_sig,\
                capsize=0,\
                color='red',\
                linestyle='None',\
                elinewidth=2.0,\
                linewidth=1.2,\
                alpha=0.8)
    #ax.set_xlim(1, 80)
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    #ax.set_ylim(1.*10**(-4), 2.*10**(-3))
    ax.set_ylim(1.*10**(-6), 2.*10**(-3))
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
    ##### Plot phase lag-energy spectra (End) #####

    ##### Plot phase lag-energy spectra (Start) #####
    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(kevs_pl,\
            csds_ph_mea_pl,\
            drawstyle='steps-post',\
            color='red',\
            linestyle='solid',\
            linewidth=2.0,\
            alpha=0.8)
    ax.errorbar(kevs_mea,\
                csds_ph_mea,\
                yerr=csds_ph_sig,\
                capsize=0,\
                color='red',\
                linestyle='None',\
                elinewidth=2.0,\
                linewidth=1.2,\
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
    ##### Plot phase lag-energy spectra (End) #####

    ##### Plot time lag-energy spectra (Start) #####
    zeros=np.zeros(len(kevs_pl))

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(kevs_pl,\
            zeros,\
            color='gray',\
            linestyle='dashed',\
            linewidth=1.5,\
            alpha=0.8)
    ax.plot(kevs_pl,\
            csds_ti_mea_pl*ONE2MILI,\
            drawstyle='steps-post',\
            color='red',\
            linestyle='solid',\
            linewidth=2.0,\
            alpha=0.8)
    ax.errorbar(kevs_mea,\
                csds_ti_mea*ONE2MILI,\
                yerr=csds_ti_sig*ONE2MILI,\
                capsize=0,\
                color='red',\
                linestyle='None',\
                elinewidth=2.0,\
                linewidth=1.2,\
                alpha=0.8)
    #ax.set_xlim(0.5, 10)
    #ax.set_xticks(np.linspace(0.5, 10, 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    ax.set_ylim(-1.0, 3.0)
    #ax.set_yticks(np.linspace(-4, 6, 6))
    #ax.set_yticks(np.linspace(-5.5, 7.5, 27), minor=True)
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
    ##### Plot time lag-energy spectra (End) #####

    ##### Plot instrinsic coherence spectra (Start) #####
    ones=np.ones(len(kevs_pl))

    plt.rcParams['font.family'] = 'Arial'
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    ax.plot(kevs_pl,\
            ones,\
            color='gray',\
            linestyle='dashed',\
            linewidth=1.5,\
            alpha=0.8)
    ax.plot(kevs_pl,\
            gamma2s_pl,\
            drawstyle='steps-post',\
            color='red',\
            linestyle='solid',\
            linewidth=2.0,\
            alpha=0.8)
    #ax.errorbar(kevs_mea,\
    #            csds_ph_mea,\
    #            yerr=csds_ph_sig,\
    #            capsize=0,\
    #            color='red',\
    #            linestyle='None',\
    #            elinewidth=2.0,\
    #            linewidth=1.2,\
    #            alpha=0.8)
    #ax.set_xlim(1, 80)
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    ax.set_ylim(0, 2)
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
    ##### Plot instrinsic coherence spectra (End) #####

if __name__=='__main__':
    main()
