import numpy as np 
import astropy.io.fits as fits
import sys 
from sys import argv
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ptick
def main():
    # ----- Setting ----- #
    names_intxt=[\
        'ni1200120110_csd_phase_0051_0100_0101_0260_xspec.txt',\
        'ni1200120110_csd_phase_0261_0480_0101_0260_xspec.txt',\
        'ni1200120110_csd_phase_0481_0700_0101_0260_xspec.txt',\
        'ni1200120110_csd_phase_0701_1100_0101_0260_xspec.txt'\
        ]
    chs_min=[51,  261, 481, 701]
    chs_max=[100, 480, 700, 1100]
    name_outpdf='./csd_phase_nicer_xspec.pdf'
    xs_pl=[6.e-4, 1.e1]
    ys_pl=[-0.2, 0.3]
    colors=['black', 'red', 'green', 'blue', 'cyan', 'magenta']
    out_pdf=PdfPages(name_outpdf)

    # ----- Main ----- #
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['figure.subplot.bottom']=0.15
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)

    for i_n, name_intxt in enumerate(names_intxt):
        e_min=chs_min[i_n]*1.e-2
        e_max=chs_max[i_n]*1.e-2

        fs, dfs, phs, dphs\
        =np.loadtxt(fname=name_intxt,\
                    dtype='float',\
                    skiprows=0,\
                    unpack=True)

        ax.errorbar(x=fs,\
                    y=phs,\
                    xerr=dfs,\
                    yerr=dphs,\
                    capsize=0.,\
                    color=colors[i_n],\
                    marker='None',\
                    markersize=6.0,\
                    drawstyle='steps-mid',\
                    linestyle='solid',\
                    linewidth=2.4,\
                    alpha=0.8,\
                    label='{0:.1f}-{1:.1f} keV'\
                        .format(e_min, e_max))

    ax.set_xlim(xs_pl[0], xs_pl[1])
    #ax.set_xlim(2.*10**(-3), 8.*10**(1))
    #ax.set_xlim(0.1, 10)
    #ax.set_xticks(np.linspace(0, 20, 21))
    #ax.set_xticks(np.linspace(0, 20, 41), minor=True)
    ax.set_ylim(ys_pl[0], ys_pl[1])
    #ax.set_yticks(np.linspace(-1200, 1200, 13))
    #ax.set_yticks(np.linspace(-1200, 1200, 25), minor=True)
    ax.set_xlabel('$f$ (Hz)', fontsize=20)
    ax.set_ylabel('$\phi (f)$ (rad)', fontsize=20)
    ax.set_xscale('log')
    #ax.set_yscale('log')
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
    #ax.set_title('{0:.2f}-{1:.2f} keV'\
    #             .format(e_min, e_max), fontsize=18)
    out_pdf.savefig(transparent=True)
    plt.close()
    out_pdf.close()

if __name__=='__main__':
    main()
