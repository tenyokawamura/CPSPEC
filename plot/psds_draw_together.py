import os
import astropy.io.fits as fits
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
def main():
    # --- #
    obsid='P0114661036'
    instruments=['LE', 'ME', 'HE']
    n_seg=3
    chs_min={'LE': [0   ], 'ME': [0   ], 'HE': [0  ]}
    chs_max={'LE': [1535], 'ME': [1023], 'HE': [255]}
    # Plot
    colors=['red', 'green', 'blue']
    xs_pl=[1.e-3, 2.e2]
    ys_pl=[2.e-4, 1.e-1]
    name_outpdf='P0114661036_net_common_psd_cpspec_together.pdf'
    # --- #
    out_pdf=PdfPages(name_outpdf)

    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['figure.subplot.bottom']=0.15
    fig=plt.figure(figsize=(9, 6))
    ax=fig.add_subplot(1, 1, 1)
    
    i_pl=0
    for instrument in instruments:
        for ch_min, ch_max in zip(chs_min[instrument], chs_max[instrument]):
            if instrument=='LE':
                e_min=0.1+13.*ch_min/1536.
                e_max=0.1+13.*ch_max/1536.
            elif instrument=='ME':
                e_min=3.+60.*ch_min/1024.
                e_max=3.+60.*ch_max/1024.
            elif instrument=='HE':
                e_min=15.+370.*ch_min/256.
                e_max=15.+370.*ch_max/256.

            name_infits='HXMT_{0}_{1}_{2:04}_{3:04}_net_common_psd.fits'.format(obsid, instrument, int(ch_min), int(ch_max))
            exist=os.path.exists(name_infits)
            if exist==False:
                print('Warning: {0} does not exist.'.format(name_infits))
                continue
            hdus=fits.open(name_infits)
            fs=hdus[1].data['F']
            dfs=hdus[1].data['DF']
            psds_raw=hdus[1].data['PSDRM']
            dpsds_raw=hdus[1].data['PSDRS']
            psds_noi=hdus[1].data['PSDN']
            psds=psds_raw-psds_noi
            dpsds=dpsds_raw
            hdus.close()

            ax.errorbar(x=fs,\
                        y=fs*psds,\
                        xerr=dfs,\
                        yerr=fs*dpsds,\
                        capsize=0.,\
                        color=colors[i_pl],\
                        marker='None',\
                        markersize=6.0,\
                        drawstyle='steps-mid',\
                        linestyle='solid',\
                        linewidth=2.4,\
                        alpha=0.8,\
                        label='{0}, {1}, {2}-{3}'.format(obsid, instrument, ch_min, ch_max))
            i_pl+=1

    ax.set_xlim(xs_pl[0], xs_pl[1])
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 6))
    #ax.set_xticks(np.linspace(0, int(num_tot-1), 11), minor=True)
    #ax.set_ylim(-15.0, 35.0)
    ax.set_ylim(ys_pl[0], ys_pl[1])
    #ax.set_yticks(np.linspace(-1, 3, 5))
    #ax.set_yticks(np.linspace(-1.6, 3.0, 24), minor=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    #ax.set_yscale('symlog', linthresh=1.e-3)
    ax.set_xlabel('$f$ (Hz)', fontsize=20)
    ax.set_ylabel('$f\,P(f)$ ($(\sigma/\mu)^2$)', fontsize=20)
    ax.spines['top'].set_linewidth(2.0)
    ax.spines['left'].set_linewidth(2.0)
    ax.spines['bottom'].set_linewidth(2.0)
    ax.spines['right'].set_linewidth(2.0)
    ax.tick_params(which='major',\
                   direction="in",\
                   length=12,\
                   width=2.0,\
                   pad=5.0,\
                   labelsize=18,\
                   top=True,\
                   bottom=True,\
                   right=True,\
                   left=True)
    ax.tick_params(which='minor',\
                   direction="in",\
                   length=6,\
                   width=2.0,\
                   pad=5.0,\
                   labelsize=18,\
                   top=True,\
                   bottom=True,\
                   right=True,\
                   left=True)
    ax.legend(bbox_to_anchor=(0, 1),\
              loc='upper left',\
              borderaxespad=1.,\
              fancybox=0,\
              edgecolor='black',\
              fontsize=10).get_frame().set_linewidth(1.2)
    #ax.set_title('{0:.1f}-{1:.1f} keV'.format(e_min, e_max), fontsize=20)
    #plt.show()
    out_pdf.savefig(transparent=True)
    plt.close()
    out_pdf.close()

if __name__=='__main__':
    main()
