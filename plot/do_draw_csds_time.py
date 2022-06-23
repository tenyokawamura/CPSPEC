import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ptick
import matplotlib.gridspec as gridspec
import astropy.io.fits as fits
import numpy as np 
import sys 
from sys import argv
from include import inputs_plot
def main():
    # ----- Setting ----- #
    name_outpdf=inputs_plot.name_inevt.replace('.evt', '_csds_time_pos.pdf')
    colors=['red', 'black', 'green', 'blue']
    xs_pl=[1.e-3, 2.e2]
    ys_pl=[1.e-3, 1.e1]
    out_pdf=PdfPages(name_outpdf)
    # --- #

    # ----- Main ----- #
    # --------------------- #
    # ----- Read data ----- #
    # --------------------- #
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['figure.subplot.bottom']=0.15
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)

    e_ref_min=1.e-2*inputs_plot.ch_ref_min
    e_ref_max=1.e-2*inputs_plot.ch_ref_max
    for i_n, (ch_min, ch_max) in enumerate(zip(inputs_plot.chs_min, inputs_plot.chs_max)):
        name_infits=inputs_plot.name_inevt.replace('.evt', '_{0:04}_{1:04}_{2:04}_{3:04}_csdf.fits'.format(int(ch_min), int(ch_max), int(inputs_plot.ch_ref_min), int(inputs_plot.ch_ref_max)))
        hdus=fits.open(name_infits)
        fs=hdus[1].data['F']
        dfs=hdus[1].data['DF']
        csds_ti=hdus[1].data['CSDTM']
        dcsds_ti=hdus[1].data['CSDTS']

        e_min=1.e-2*ch_min
        e_max=1.e-2*ch_max

        ax.errorbar(x=fs,\
                    y=csds_ti,\
                    xerr=dfs,\
                    yerr=dcsds_ti,\
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
    ax.set_ylabel('$\\tau (f)$ (s)', fontsize=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    #ax.set_yscale('symlog', linthresh=1.e-3)
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
    ax.legend(bbox_to_anchor=(1, 1),\
              loc='upper right',\
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
