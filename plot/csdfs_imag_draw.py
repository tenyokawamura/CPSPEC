import os
import astropy.io.fits as fits 
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
def main():
    # --- #
    obsid='P0114661036'
    instruments=['LE', 'ME', 'HE']
    chs_min={'LE': [106, 296, 556, 816],  'ME': [137, 342], 'HE': [14, 24, 37, 60]}
    chs_max={'LE': [295, 555, 815, 1288], 'ME': [341, 546], 'HE': [23, 36, 59, 93]}
    instrument_ref='LE'
    ch_ref_min=296
    ch_ref_max=555
    # Plot
    colors=['red', 'green', 'blue']
    xs_pl=[1.e-3, 2.e2]
    ys_pl=[2.e-4, 1.e-1]
    name_outpdf='P0114661036_net_common_csdf_imag_cpspec.pdf'
    # --- #
    out_pdf=PdfPages(name_outpdf)

    # Reference band
    if instrument_ref=='LE':
        e_ref_min=0.1+13.*ch_ref_min/1536.
        e_ref_max=0.1+13.*ch_ref_max/1536.
    elif instrument_ref=='ME':
        e_ref_min=3.+60.*ch_ref_min/1024.
        e_ref_max=3.+60.*ch_ref_max/1024.
    elif instrument_ref=='HE':
        e_ref_min=15.+370.*ch_ref_min/256.
        e_ref_max=15.+370.*ch_ref_max/256.
    for instrument in instruments:
        for ch_min, ch_max in zip(chs_min[instrument], chs_max[instrument]):
            plt.rcParams['font.family'] = 'Arial'
            plt.rcParams['figure.subplot.bottom']=0.15
            fig=plt.figure(figsize=(9, 6))
            ax=fig.add_subplot(1, 1, 1)

            # Subject band
            if instrument=='LE':
                e_min=0.1+13.*ch_min/1536.
                e_max=0.1+13.*ch_max/1536.
            elif instrument=='ME':
                e_min=3.+60.*ch_min/1024.
                e_max=3.+60.*ch_max/1024.
            elif instrument=='HE':
                e_min=15.+370.*ch_min/256.
                e_max=15.+370.*ch_max/256.

            name_infits='HXMT_{0}_{1}_{2:04}_{3:04}_{4}_{5:04}_{6:04}_net_common_csdf.fits'.format(obsid, instrument_ref, int(ch_ref_min), int(ch_ref_max), instrument, int(ch_min), int(ch_max))
            exist=os.path.exists(name_infits)
            if exist==False:
                print('Warning: {0} does not exist.'.format(name_infits))
                continue
            hdus=fits.open(name_infits)
            fs=hdus[1].data['F']
            dfs=hdus[1].data['DF']
            csds_im=hdus[1].data['CSDIM']
            dcsds_im=hdus[1].data['CSDIS']
            hdus.close()

            ax.errorbar(x=fs,\
                        y=fs*csds_im,\
                        xerr=dfs,\
                        yerr=fs*dcsds_im,\
                        capsize=0.,\
                        color='black',\
                        marker='None',\
                        markersize=6.0,\
                        drawstyle='steps-mid',\
                        linestyle='solid',\
                        linewidth=2.4,\
                        alpha=0.8,\
                        label='{0}, {1}, {2}-{3} vs {4}, {5}-{6}'.format(obsid, instrument, ch_min, ch_max, instrument_ref, ch_ref_min, ch_ref_max))

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
            ax.set_ylabel('$\mathrm{Im}[f\,C(f)]$ ($(\sigma/\mu)^2$)', fontsize=20)
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
                      fontsize=16).get_frame().set_linewidth(1.2)
            ax.set_title('{0:.1f}-{1:.1f} keV vs {2:.1f}-{3:.1f} keV'.format(e_min, e_max, e_ref_min, e_ref_max), fontsize=20)
            #plt.show()
            out_pdf.savefig(transparent=True)
            plt.close()
    out_pdf.close()

def read_data_psd(name_inqdp):
    f,\
    dt,\
    psd_mean,\
    psd_sigma,\
    col5\
    =np.loadtxt(fname=name_inqdp, \
                dtype='float', \
                skiprows=3, \
                unpack=True)
    return f, dt, psd_mean, psd_sigma, col5

def params_adjust_plot(f,\
                       psd_mean,\
                       psd_sigma):
    ### ax.plot ###
    f_plx=np.empty(0)
    for i_f in range(len(f)):
        if i_f==0:
            f_plx=np.append(f_plx, 0)
        else:
            f_plx=np.append(f_plx, (f[i_f]+f[i_f-1])/2.)
    f_plx[0]=f[0]-(f[1]-f[0])/2. #Assume Delta f=const.
    f_plx=np.append(f_plx, f[-1]+(f[-1]-f[-2])/2.)
    f_ply=np.append(f, f[-1])
    psd_mean_pl=np.append(psd_mean, psd_mean[-1])
    ### ax.errorbar ###
    f_errx=f
    f_erry=f
    psd_mean_err=psd_mean
    psd_sigma_err=psd_sigma
    return f_plx, f_ply, psd_mean_pl, f_errx, f_erry, psd_mean_err, psd_sigma_err

if __name__=='__main__':
    main()
