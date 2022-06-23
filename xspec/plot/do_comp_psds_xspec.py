import numpy as np 
import astropy.io.fits as fits
import sys 
from sys import argv
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ptick
from include import inputs_xspec_plot
def main():
    # ----- Setting ----- #
    name_outpdf=inputs_xspec_plot.name_inevt.replace('.evt', '_psds_xspec.pdf')
    xs_pl=[1.e-3, 2.e2]
    ys_pl=[2.e-4, 1.e-1]
    out_pdf=PdfPages(name_outpdf)

    # ----- Main ----- #
    for i_n, (ch_min, ch_max) in enumerate(zip(inputs_xspec_plot.chs_min, inputs_xspec_plot.chs_max)):
        # -------------------- #
        # ----- Original ----- #
        # -------------------- #
        name_infits=inputs_xspec_plot.name_inevt.replace('.evt', '_{0:04}_{1:04}_psd.fits'.format(int(ch_min), int(ch_max)))
        hdus=fits.open(name_infits)
        fs=hdus[1].data['F']
        dfs=hdus[1].data['DF']
        psds_raw=hdus[1].data['PSDRM']
        dpsds_raw=hdus[1].data['PSDRS']
        psds_noi=hdus[1].data['PSDN']
        psds=psds_raw-psds_noi
        dpsds=dpsds_raw

        ch_min=hdus[1].header['CHMIN']
        ch_max=hdus[1].header['CHMAX']
        e_min=ch_min*1.e-2
        e_max=ch_max*1.e-2

        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['figure.subplot.bottom']=0.15
        fig=plt.figure(figsize=(9, 6))
        ax=fig.add_subplot(1, 1, 1)

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
                    alpha=0.8,\
                    label='Original')

        # ----------------- #
        # ----- XSPEC ----- #
        # ----------------- #
        name_intxt=inputs_xspec_plot.name_inevt.replace('.evt', '_{0:04}_{1:04}_psd_xspec.txt'.format(int(ch_min), int(ch_max)))
        fs, dfs, psds, dpsds\
        =np.loadtxt(fname=name_intxt,\
                    dtype='float',\
                    skiprows=0,\
                    unpack=True)

        ax.errorbar(x=fs,\
                    y=fs*psds,\
                    xerr=dfs,\
                    yerr=fs*dpsds,\
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
        ax.set_ylabel('$f\,P(f)$ ($(\sigma/\mu)^2$)', fontsize=20)
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
        ax.set_title('{0:.1f}-{1:.1f} keV'.format(e_min, e_max), fontsize=20)
        out_pdf.savefig(transparent=True)
        plt.close()
    out_pdf.close()

if __name__=='__main__':
    main()
