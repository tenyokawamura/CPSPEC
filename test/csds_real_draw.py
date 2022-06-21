import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ptick
import matplotlib.gridspec as gridspec
import astropy.io.fits as fits
import numpy as np 
import sys 
from sys import argv
def main():
    # ----- Setting ----- #
    names_infits=[\
        './ni1200120147_0mpu7_cl_bary_fpm_0101_0260_0261_0480_csdf.fits',\
        './ni1200120147_0mpu7_cl_bary_fpm_0261_0480_0261_0480_csdf.fits',\
        './ni1200120147_0mpu7_cl_bary_fpm_0481_0700_0261_0480_csdf.fits',\
        './ni1200120147_0mpu7_cl_bary_fpm_0701_1100_0261_0480_csdf.fits'\
        ]
    ch_ref_min=261
    ch_ref_max=480
    chs_min   =[101, ch_ref_min, 481, 701 ]
    chs_max   =[260, ch_ref_max, 700, 1100]
    name_outpdf='./ni1200120147_0mpu7_cl_bary_fpm_csds_real.pdf'
    colors=['red', 'black', 'green', 'blue']
    xs_pl=[1.e-3, 2.e2]
    ys_pl=[2.e-4, 1.e-1]
    out_pdf=PdfPages(name_outpdf)

    # ----- Main ----- #
    # --------------------- #
    # ----- Read data ----- #
    # --------------------- #
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['figure.subplot.bottom']=0.15
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)

    e_ref_min=1.e-2*ch_ref_min
    e_ref_max=1.e-2*ch_ref_max
    for i_n, name_infits in enumerate(names_infits):
        hdus=fits.open(name_infits)
        fs=hdus[1].data['F']
        dfs=hdus[1].data['DF']
        csds_re_raw=hdus[1].data['CSDRM']
        dcsds_re=hdus[1].data['CSDRS']
        csds_noise=hdus[1].data['CSDN']
        csds_re=csds_re_raw-csds_noise

        e_min=1.e-2*chs_min[i_n]
        e_max=1.e-2*chs_max[i_n]

        ax.errorbar(x=fs,\
                    y=fs*csds_re,\
                    xerr=dfs,\
                    yerr=fs*dcsds_re,\
                    capsize=0.,\
                    color=colors[i_n],\
                    marker='None',\
                    markersize=6.0,\
                    drawstyle='steps-mid',\
                    linestyle='solid',\
                    linewidth=2.4,\
                    alpha=0.8,\
                    label='{0:.1f} - {1:.1f} keV'.format(e_min, e_max))

    ax.set_xlim(xs_pl[0], xs_pl[1])
    #ax.set_xlim(2.*10**(-3), 8.*10**(1))
    #ax.set_xlim(0.1, 10)
    #ax.set_xticks(np.linspace(0, 20, 21))
    #ax.set_xticks(np.linspace(0, 20, 41), minor=True)
    ax.set_ylim(ys_pl[0], ys_pl[1])
    #ax.set_yticks(np.linspace(-1200, 1200, 13))
    #ax.set_yticks(np.linspace(-1200, 1200, 25), minor=True)
    ax.set_xlabel('$f$ (Hz)', fontsize=20)
    ax.set_ylabel('$\mathrm{Re}[f\,C(f)]$ ($(\sigma / \mu)^2$)', fontsize=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    #ax.xaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))
    #ax.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    #ax.xaxis.offsetText.set_fontsize(16)
    #ax.yaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))
    #ax.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
    #ax.yaxis.offsetText.set_fontsize(16)
    ax.spines['top'].set_linewidth(2.0)
    ax.spines['left'].set_linewidth(2.0)
    ax.spines['bottom'].set_linewidth(2.0)
    ax.spines['right'].set_linewidth(2.0)
    ax.legend(bbox_to_anchor=(0, 1),\
              loc='upper left',\
              borderaxespad=1.,\
              fancybox=0,\
              frameon=False,\
              edgecolor='black',\
              fontsize=16).get_frame().set_linewidth(1.2)
    ax.tick_params(which='major',\
                   direction="in",\
                   length=10,\
                   width=2.0,\
                   labelsize=18,\
                   pad=5.0,\
                   top=True,\
                   bottom=True,\
                   right=True,\
                   left=True)
    ax.tick_params(which='minor',\
                   direction="in",\
                   length=5,\
                   width=2.0,\
                   labelsize=18,\
                   pad=5.0,\
                   top=True,\
                   bottom=True,\
                   right=True,\
                   left=True)
    ax.set_title('Reference: {0:.1f} - {1:.1f} keV'.format(e_ref_min, e_ref_max), fontsize=20)
    out_pdf.savefig(transparent=True)
    plt.close()
    out_pdf.close()

if __name__=='__main__':
    main()
