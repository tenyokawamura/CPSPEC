import numpy as np 
import astropy.io.fits as fits
import sys 
from sys import argv
import os
def main():
    # ----- Setting ----- #
    name_infits='./ni1200120110_csd.fits'
    name_outpha_def='./ni1200120110_csd_phase.pha'
    name_rmf='response_nicer.rmf'
    name_arf='response_nicer.arf'
    f_min=5.e-3 #[Hz]
    f_max=4.e1  #[Hz]
    n_bin=1000 #[-]
    i_f_min=3
    i_f_max=None
    sys_err=0. # Systematic error [x 100 %]
    all_coi=True
    if all_coi==True:
        hdus=fits.open(name_infits)
        n_coi=len(hdus)-1
        is_b=np.arange(0, n_coi, 1)
        hdus.close()
    else:
        is_b=[0, 1]
        #is_band=[9, 10, 11]
        #is_band=[0]

    # ----- Set time/frequency ----- #
    fs=np.geomspace(f_min, f_max, n_bin+1)

    fs_min=fs
    fs_max=np.roll(fs, -1)
    # Delete the maximum bin
    fs_min=fs_min[:-1]
    fs_max=fs_max[:-1]
    
    fs=(fs_min+fs_max)/2.
    dfs=(fs_max-fs_min)/2.

    n_ch=len(fs_min)

    # ----- Main ----- #
    hdus=fits.open(name_infits)
    for i_b in is_b:
        # --------------------- #
        # ----- Read data ----- #
        # --------------------- #
        ref_min=hdus[i_b+1].header['REFMIN']
        ref_max=hdus[i_b+1].header['REFMAX']
        coi_min=hdus[i_b+1].header['CHMIN']
        coi_max=hdus[i_b+1].header['CHMAX']

        fs_min_data=hdus[i_b+1].data['FMIN']
        fs_max_data=hdus[i_b+1].data['FMAX']
        fs_mid_data=(fs_min_data+fs_max_data)/2.
        dfs_data=(fs_max_data-fs_min_data)/2.
        csds_ph_mea=hdus[i_b+1].data['CSDPM']
        csds_ph_sig=hdus[i_b+1].data['CSDPS']

        # ----- Limit range ----- #
        fs_min_data=fs_min_data[i_f_min:i_f_max]
        fs_max_data=fs_max_data[i_f_min:i_f_max]
        csds_ph_mea=csds_ph_mea[i_f_min:i_f_max]
        csds_ph_sig=csds_ph_sig[i_f_min:i_f_max]
        csds_ph_mea*=-1 # Reference band leads energy band --> Positive lag

        # ----- Conversion from frequency (energy in XSPEC) into channel (PHA) ----- #
        first=True
        for f_min_data, f_max_data in zip(fs_min_data, fs_max_data):
            chs=np.where((f_min_data<fs) & (fs<=f_max_data))[0]
            ch_rebin_min=chs[0]
            ch_rebin_max=chs[-1]
            df_rebin=np.sum(2.*dfs[chs])
            if first==True:
                dfs_rebin=df_rebin
                chs_rebin_min=ch_rebin_min
                chs_rebin_max=ch_rebin_max
                first=False
            else:
                dfs_rebin=np.append(dfs_rebin, df_rebin)
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
        syses_err=np.zeros(n_ch) # dP(f) (error)
        qualities=np.zeros(n_ch)
        groupings=-np.ones(n_ch)

        groupings[0]=1 # asked by XSPEC
        for ch_rebin_min, ch_rebin_max, csd_ph_mea, csd_ph_sig, df_rebin\
        in zip(chs_rebin_min, chs_rebin_max, csds_ph_mea, csds_ph_sig, dfs_rebin):
            rates[ch_rebin_min]=csd_ph_mea*df_rebin
            stats_err[ch_rebin_min]=csd_ph_sig*df_rebin
            syses_err[ch_rebin_min]=sys_err
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
                                       fits.Column(name='SYS_ERR',  format='f4', array=syses_err),\
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
