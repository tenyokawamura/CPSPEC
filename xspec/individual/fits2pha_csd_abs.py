import numpy as np 
import astropy.io.fits as fits
import sys 
from sys import argv
import os
def main():
    # ----- Setting ----- #
    name_infits='./ni1200120110_timing.fits'
    name_outpha_def='./ni1200120110_csd_amp.pha'
    name_rmf='simple_timing.rmf'
    name_arf='const_timing.arf'
    dt=5.e-2 #[s]
    n_data=2**(14) #[-]
    i_f_min=0
    i_f_max=None
    all_coi=False
    if all_coi==True:
        hdus=fits.open(name_infits_obs)
        n_coi=len(hdus)-1
        is_band=np.arange(0, n_coi, 1)
        hdus.close()
    else:
        is_b=[0, 1]
        #is_band=[9, 10, 11]
        #is_band=[0]

    # ----- Set time/frequency ----- #
    lt=n_data*dt #[s]
    ts=np.arange(0, lt, dt)
    if not n_data==len(ts):
        print('Error')
        sys.exit()

    f_min=1./lt #[Hz]
    f_max=1./(2.*dt) #[Hz]
    df=f_min #[Hz]
    fs=np.arange(f_min, f_max+df, df)
    n_ch=len(fs) #Number of channels in XSPEC
    #print(n_ch, 2*n_ch, n_data)

    # ----- Main ----- #
    hdus=fits.open(name_infits)
    for i_b in is_b:
        # --------------------- #
        # ----- Read data ----- #
        # --------------------- #
        ref_min=hdus[i_b+2].header['REFMIN']
        ref_max=hdus[i_b+2].header['REFMAX']
        coi_min=hdus[i_b+2].header['CHMIN']
        coi_max=hdus[i_b+2].header['CHMAX']

        fs_min_data=hdus[i_b+2].data['FMIN']
        fs_max_data=hdus[i_b+2].data['FMAX']
        fs_mid_data=(fs_min_data+fs_max_data)/2.
        dfs_data=(fs_max_data-fs_min_data)/2.
        csds_amp_mea=hdus[i_b+2].data['CSDAM']
        csds_amp_sig=hdus[i_b+2].data['CSDAS']
        dfs=fs_max_data-fs_min_data

        # ----- Limit range ----- #
        fs_min_data=fs_min_data[i_f_min:i_f_max]
        fs_max_data=fs_max_data[i_f_min:i_f_max]
        dfs=dfs[i_f_min:i_f_max]
        csds_amp_mea=csds_amp_mea[i_f_min:i_f_max]
        csds_amp_sig=csds_amp_sig[i_f_min:i_f_max]

        # ----- Conversion from frequency (energy in XSPEC) into channel (PHA) ----- #
        begin=True
        for f_min_data, f_max_data in zip(fs_min_data, fs_max_data):
            chs=np.where((f_min_data<fs) & (fs<=f_max_data))
            ch_rebin_min=chs[0][0]
            ch_rebin_max=chs[0][-1]
            if begin==True:
                chs_rebin_min=ch_rebin_min
                chs_rebin_max=ch_rebin_max
                begin=False
            else:
                chs_rebin_min=np.append(chs_rebin_min, ch_rebin_min)
                chs_rebin_max=np.append(chs_rebin_max, ch_rebin_max)
            
        # -------------------- #
        # ----- Make PHA ----- #
        # -------------------- #
        # ----- PRIMARY ----- #
        hdu_pri=fits.PrimaryHDU(data=None, header=None)

        # ----- SPECTRUM ----- #
        chs=np.arange(n_ch)
        rates=np.zeros(n_ch) # P(f) (power spectrum)
        stats_err=np.zeros(n_ch) # dP(f) (error)
        qualities=np.zeros(n_ch)
        groupings=-np.ones(n_ch)

        for ch_rebin_min, ch_rebin_max, csd_amp_mea, csd_amp_sig, df\
        in zip(chs_rebin_min, chs_rebin_max, csds_amp_mea, csds_amp_sig, dfs):
            # df is necessary for XSPEC to show fP(f) with 'plot euf'.
            rates[ch_rebin_min]=csd_amp_mea*df 
            stats_err[ch_rebin_min]=csd_amp_sig*df
            groupings[ch_rebin_min]=1
        groupings[chs_rebin_max[-1]+1]=1
        for ch in range(chs_rebin_min[0]):
            qualities[ch]=5
        for ch in range(chs_rebin_max[-1]+1, len(chs)):
            qualities[ch]=5
        hdu_spe=\
        fits.BinTableHDU.from_columns([fits.Column(name='CHANNEL',  format='i4', array=chs),\
                                       fits.Column(name='RATE',     format='f4', array=rates),\
                                       fits.Column(name='STAT_ERR', format='f4', array=stats_err),\
                                       fits.Column(name='GROUPING', format='i2', array=groupings),\
                                       fits.Column(name='QUALITY',  format='i2', array=qualities)])

        hdu_spe.header['TUNIT2']='count/s'
        hdu_spe.header['EXTNAME']='SPECTRUM'
        hdu_spe.header['HDUCLASS']='OGIP'
        hdu_spe.header['HDUVERS1']='1.2.0'
        hdu_spe.header['HDUVERS']='1.2.0'
        hdu_spe.header['HDUCLAS3']='RATE'
        hdu_spe.header['TLMIN1']=0
        hdu_spe.header['TLMAX1']=n_ch-1
        hdu_spe.header['TELESCOP']='NICER'
        hdu_spe.header['INSTRUME']='NICER'
        hdu_spe.header['FILTER']='NONE'
        hdu_spe.header['AREASCAL']=1.0
        hdu_spe.header['BACKFILE']='none'
        hdu_spe.header['BACKSCAL']=1.0
        hdu_spe.header['CORRFILE']='none'
        hdu_spe.header['CORRSCAL']=1.0
        hdu_spe.header['RESPFILE']=name_rmf
        hdu_spe.header['ANCRFILE']=name_arf
        hdu_spe.header['DETCHANS']=n_ch
        hdu_spe.header['CHANTYPE']='PHA'
        hdu_spe.header['POISSERR']=False
        hdu_spe.header['STAT_ERR']=0
        hdu_spe.header['SYS_ERR']=0
        hdu_spe.header['GROUPING']=0
        hdu_spe.header['QUALITY']=0
        hdu_spe.header['HDUCLAS1']='SPECTRUM'
        hdu_spe.header['DATAMODE']='PHOTON'
        # Exposure [s]: As long as 'RATE' is used instead of 'COUNT', 
        # any positive value should be fine.
        # But, 1 sec is the safest.
        hdu_spe.header['EXPOSURE']=1 
        hdu_spe.header['OBJECT']='MAXI J1820+070'

        # ----- Write ----- #
        hdus_out=fits.HDUList([hdu_pri, hdu_spe])
        name_outpha=name_outpha_def.replace('.pha',\
                                            '_{0:04}_{1:04}_{2:04}_{3:04}.pha'\
                                            .format(coi_min, coi_max, ref_min, ref_max))
        hdus_out.writeto(name_outpha, overwrite=True)

if __name__=='__main__':
    main()
