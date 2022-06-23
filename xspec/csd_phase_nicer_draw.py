import numpy as np 
import astropy.io.fits as fits
import sys 
from sys import argv
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ptick
def main():
    # ----- Setting ----- #
    # --- #
    #name_infits='./ni1200120110_csd.fits'
    #name_outpdf='./csd_phase_nicer.pdf'
    #i_f_min=0
    #i_f_max=None
    #xs_pl=[6.e-4, 1.e1]
    #ys_pl=[-0.2, 0.3]
    #colors=['black', 'red', 'green', 'blue', 'cyan', 'magenta']
    #all_coi=True
    #if all_coi==True:
    #    hdus=fits.open(name_infits)
    #    n_coi=len(hdus)-1
    #    is_b=np.arange(n_coi)
    #    hdus.close()
    #else:
    #    is_b=[0]
    #    #is_band=[9, 10, 11]
    #    #is_band=[0]
    # --- #
    name_infits='./ni1200120110_csd.fits'
    name_outpdf='./csd_phase_nicerb.pdf'
    i_f_min=0
    i_f_max=None
    xs_pl=[6.e-4, 1.e1]
    ys_pl=[-0.1, 0.2]
    colors=['black', 'red', 'green', 'blue', 'cyan', 'magenta']
    all_coi=False
    if all_coi==True:
        hdus=fits.open(name_infits)
        n_coi=len(hdus)-1
        is_b=np.arange(n_coi)
        hdus.close()
    else:
        is_b=[1, 2, 3]
        #is_band=[9, 10, 11]
        #is_band=[0]
    # --- #
    out_pdf=PdfPages(name_outpdf)

    # ----- Main ----- #
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['figure.subplot.bottom']=0.15
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)

    hdus=fits.open(name_infits)
    for i_n, i_b in enumerate(is_b):
        # --------------------- #
        # ----- Read data ----- #
        # --------------------- #
        ref_min=hdus[i_b+1].header['REFMIN'] #[keV]
        ref_max=hdus[i_b+1].header['REFMAX'] #[keV]
        ch_min=hdus[i_b+1].header['CHMIN'] #[keV]
        ch_max=hdus[i_b+1].header['CHMAX'] #[keV]
        e_ref_min=ref_min*1.e-2
        e_ref_max=ref_max*1.e-2
        e_min=ch_min*1.e-2
        e_max=ch_max*1.e-2
        fs_min_data=hdus[i_b+1].data['FMIN']
        fs_max_data=hdus[i_b+1].data['FMAX']
        phs=hdus[i_b+1].data['CSDPM']
        dphs=hdus[i_b+1].data['CSDPS']
        fs_mid_data=(fs_min_data+fs_max_data)/2.
        dfs_data=(fs_max_data-fs_min_data)/2.

        # ----- Limit range ----- #
        fs_min_data=fs_min_data[i_f_min:i_f_max]
        fs_max_data=fs_max_data[i_f_min:i_f_max]
        dfs_data   =dfs_data   [i_f_min:i_f_max]
        phs        =phs        [i_f_min:i_f_max]
        dphs       =dphs       [i_f_min:i_f_max]
        phs*=-1 # Reference band leads energy band --> positive lag

        ax.errorbar(x=fs_mid_data,\
                    y=phs,\
                    xerr=dfs_data,\
                    yerr=dphs,\
                    capsize=0.,\
                    color=colors[i_n],\
                    marker='None',\
                    markersize=6.0,\
                    drawstyle='steps-mid',\
                    linestyle='solid',\
                    linewidth=2.4,\
                    alpha=0.8,\
                    label='{0:.1f}-{1:.1f} keV vs {2:.1f}-{3:.1f} keV'\
                        .format(e_min, e_max, e_ref_min, e_ref_max))

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
