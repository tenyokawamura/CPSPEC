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
    name_infits='./ni1200120147_0mpu7_cl_bary_fpm_0101_0260_psd.fits'
    name_outpdf='./ni1200120147_0mpu7_cl_bary_fpm_0101_0260_psd.pdf'
    colors=['blue', 'red', 'magenta', 'orange', 'cyan', 'brown']
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

    hdus=fits.open(name_infits)
    fs=hdus[1].data['F']
    dfs=hdus[1].data['DF']
    psds_raw=hdus[1].data['PSDRM']
    dpsds_raw=hdus[1].data['PSDRS']
    psds_noi=hdus[1].data['PSDN']
    psds=psds_raw-psds_noi
    dpsds=dpsds_raw

    ax.errorbar(x=fs,\
                y=fs*psds,\
                xerr=dfs,\
                yerr=fs*dpsds,\
                capsize=0.,\
                color='black',\
                marker='None',\
                markersize=6.0,\
                drawstyle='steps-mid',\
                linestyle='solid',\
                linewidth=2.4,\
                alpha=0.8)

    ax.set_xlim(xs_pl[0], xs_pl[1])
    #ax.set_xlim(2.*10**(-3), 8.*10**(1))
    #ax.set_xlim(0.1, 10)
    #ax.set_xticks(np.linspace(0, 20, 21))
    #ax.set_xticks(np.linspace(0, 20, 41), minor=True)
    ax.set_ylim(ys_pl[0], ys_pl[1])
    #ax.set_yticks(np.linspace(-1200, 1200, 13))
    #ax.set_yticks(np.linspace(-1200, 1200, 25), minor=True)
    ax.set_xlabel('$f$ (Hz)', fontsize=20)
    ax.set_ylabel('$f\,P(f)$ ($(\sigma / \mu)^2$)', fontsize=20)
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
    #ax.legend(bbox_to_anchor=(1, 0),\
    #          loc='lower right',\
    #          borderaxespad=1.,\
    #          fancybox=0,\
    #          edgecolor='black',\
    #          fontsize=16).get_frame().set_linewidth(1.2)
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
    ax.set_title('Power spectrum', fontsize=20)
    out_pdf.savefig(transparent=True)
    plt.close()
    out_pdf.close()

def read_psd_ave_ns(name_intxt):
    with open(name_intxt, 'r') as fin:
        line=fin.readline()
        line=fin.readline()
        line_str=line.split()
        cr_mean=float(line_str[0])
        exposure=float(line_str[1])

    f,\
    psd_mean,\
    psd_sigma,\
    psd_noise,\
    n_mer_p,\
    n_mer_f\
    =np.loadtxt(fname=name_intxt,\
                dtype='float',\
                skiprows=3,\
                unpack=True)
    return f,\
           psd_mean,\
           psd_sigma,\
           psd_noise,\
           n_mer_p,\
           n_mer_f,\
           cr_mean,\
           exposure

if __name__=='__main__':
    main()
