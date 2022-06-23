import numpy as np 
import astropy.io.fits as fits
import sys 
from sys import argv
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ptick
def main():
    # ----- Setting ----- #
    names_infits=[\
        'ni1200120130_0mpu7_cl_fpm_0101_0260_0261_0480_csdf.fits'\
        ]
    names_intxt=[\
        'ni1200120130_0mpu7_cl_fpm_0101_0260_0261_0480_csd_abs.txt'\
        ]
    name_outpdf='./ni1200120130_0mpu7_cl_fpm_csd_abs_comp_xspec.pdf'
    xs_pl=[1.e-3, 2.e2]
    ys_pl=[2.e-4, 1.e-1]
    colors=['black', 'red', 'green', 'blue', 'cyan', 'magenta']
    out_pdf=PdfPages(name_outpdf)

    # ----- Main ----- #
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['figure.subplot.bottom']=0.15
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)

    for i_n, (name_infits, name_intxt) in enumerate(zip(names_infits, names_intxt)):
        # -------------------- #
        # ----- Original ----- #
        # -------------------- #
        hdus=fits.open(name_infits)
        fs=hdus[1].data['F']
        dfs=hdus[1].data['DF']
        csds_ab=hdus[1].data['CSDAM'] # unused
        dcsds_ab=hdus[1].data['CSDAS']

        csds_re_raw=hdus[1].data['CSDRM']
        csds_im=hdus[1].data['CSDIM']
        csds_noise=hdus[1].data['CSDN']
        csds_re=csds_re_raw-csds_noise
        csds_ab=np.sqrt(csds_re**2+csds_im**2)

        ch_min=hdus[1].header['CHMIN']
        ch_max=hdus[1].header['CHMAX']
        ch_ref_min=hdus[1].header['REFMIN']
        ch_ref_max=hdus[1].header['REFMAX']
        e_ref_min=ch_ref_min*1.e-2
        e_ref_max=ch_ref_max*1.e-2
        e_min=ch_min*1.e-2
        e_max=ch_max*1.e-2

        ax.errorbar(x=fs,\
                    y=fs*csds_ab,\
                    xerr=dfs,\
                    yerr=fs*dcsds_ab,\
                    capsize=0.,\
                    color='black',\
                    marker='None',\
                    markersize=6.0,\
                    drawstyle='steps-mid',\
                    linestyle='solid',\
                    linewidth=2.4,\
                    alpha=0.8,\
                    label='Original')

        # ----------------- #
        # ----- XSPEC ----- #
        # ----------------- #
        fs, dfs, csds_ab, dcsds_ab\
        =np.loadtxt(fname=name_intxt,\
                    dtype='float',\
                    skiprows=0,\
                    unpack=True)

        ax.errorbar(x=fs,\
                    y=fs*csds_ab,\
                    xerr=dfs,\
                    yerr=fs*dcsds_ab,\
                    capsize=0.,\
                    color='red',\
                    marker='None',\
                    markersize=6.0,\
                    drawstyle='steps-mid',\
                    linestyle='solid',\
                    linewidth=2.4,\
                    alpha=0.8,\
                    label='XSPEC')

        ax.set_xlim(xs_pl[0], xs_pl[1])
        #ax.set_xlim(2.*10**(-3), 8.*10**(1))
        #ax.set_xlim(0.1, 10)
        #ax.set_xticks(np.linspace(0, 20, 21))
        #ax.set_xticks(np.linspace(0, 20, 41), minor=True)
        ax.set_ylim(ys_pl[0], ys_pl[1])
        #ax.set_yticks(np.linspace(-1200, 1200, 13))
        #ax.set_yticks(np.linspace(-1200, 1200, 25), minor=True)
        ax.set_xlabel('$f$ (Hz)', fontsize=20)
        ax.set_ylabel('$f\,|C(f)|$ ($(\sigma / \mu)^2$)', fontsize=20)
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
                  fontsize=14).get_frame().set_linewidth(1.2)
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
        ax.set_title('{0:.1f}-{1:.1f} keV vs {2:.1f}-{3:.1f} keV'.format(e_min, e_max, e_ref_min, e_ref_max), fontsize=20)
        out_pdf.savefig(transparent=True)
        plt.close()
    out_pdf.close()

if __name__=='__main__':
    main()