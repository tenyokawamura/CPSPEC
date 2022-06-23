import astropy.io.fits as fits
import numpy as np 
import os
import sys 
from sys import argv
from include import inputs_xspec

def main():
    # -------------------------- #
    # ---------- Main ---------- #
    # -------------------------- #
    # ----- Set time/frequency ----- #
    fs=np.geomspace(inputs_xspec.f_min, inputs_xspec.f_max, inputs_xspec.n_bin+1)

    fs_min=fs
    fs_max=np.roll(fs, -1)
    # Delete the maximum bin
    fs_min=fs_min[:-1]
    fs_max=fs_max[:-1]

    n_ch=len(fs_min)
    chs=np.arange(n_ch)

    # ---------- Calibration ---------- #
    # ----- PRIMARY ----- #
    hdu_pri=fits.PrimaryHDU(data=None, header=None)

    # ----- EBOUND ----- #
    es_min=fs_min
    es_max=fs_max

    hdu_ebo\
    =fits.BinTableHDU.from_columns([fits.Column(name='CHANNEL', format='i2', array=chs),\
                                    fits.Column(name='E_MIN',   format='f8', array=es_min),\
                                    fits.Column(name='E_MAX',   format='f8', array=es_max)])

    hdu_ebo.header['TUNIT1']=''
    hdu_ebo.header['TUNIT2']='keV'
    hdu_ebo.header['TUNIT3']='keV'
    hdu_ebo.header['EXTNAME']='EBOUNDS'
    hdu_ebo.header['HDUCLASS']='OGIP'
    hdu_ebo.header['HDUCLAS1']='RESPONSE'
    hdu_ebo.header['HDUCLAS2']='EBOUNDS'
    hdu_ebo.header['CHANTYPE']='PHA'
    hdu_ebo.header['HDUVERS']='1.1.0'
    hdu_ebo.header['TELESCOP']='NICER'
    hdu_ebo.header['INSTRUME']='NICER'
    hdu_ebo.header['FILTER']='NONE'
    hdu_ebo.header['DETCHANS']=n_ch

    # ----- MATRIX ----- #
    # --- (Simple) linear binning for frequency (= Energy in XSPEC) --- #
    fs_min_mat=fs_min
    fs_max_mat=fs_max
    # --- Set parameters --- #
    n_f=len(fs_min_mat)
    ns_grp=np.ones(n_f)
    ns_chan=np.ones(n_f)
    first=True
    for f_min_mat, f_max_mat, n_chan in zip(fs_min_mat, fs_max_mat, ns_chan):
        f_mid=(f_min_mat+f_max_mat)/2.
        #ch=int(f_mid/df)-1 # f_min=df corresponds to channel 0.
        ch=((fs_min-f_mid)**2+(fs_max-f_mid)**2).argmin()
        f_chan=ch
        matrix=np.ones(int(n_chan))
        if first==True:
            first=False
            fs_chan=f_chan
            matrices=matrix
        else:
            fs_chan=np.append(fs_chan, f_chan)
            matrices=np.vstack((matrices, matrix))
    es_lo=fs_min_mat
    es_hi=fs_max_mat

    hdu_mat\
    =fits.BinTableHDU.from_columns([fits.Column(name='ENERG_LO', format='f8', array=es_lo),\
                                    fits.Column(name='ENERG_HI', format='f8', array=es_hi),\
                                    fits.Column(name='N_GRP',    format='i2', array=ns_grp),\
                                    fits.Column(name='F_CHAN',   format='i4', array=fs_chan),\
                                    fits.Column(name='N_CHAN',   format='i4', array=ns_chan),\
                                    fits.Column(name='MATRIX',   format='PE(398)', array=matrices)])

    hdu_mat.header['TUNIT1']='keV'
    hdu_mat.header['TUNIT2']='keV'
    hdu_mat.header['TUNIT3']=''
    hdu_mat.header['TUNIT4']=''
    hdu_mat.header['TUNIT5']=''
    hdu_mat.header['TUNIT6']=''
    hdu_mat.header['EXTNAME']='MATRIX'
    hdu_mat.header['HDUCLASS']='OGIP'
    hdu_mat.header['HDUCLAS1']='RESPONSE'
    hdu_mat.header['HDUCLAS2']='RSP_MATRIX'
    hdu_mat.header['HDUCLAS3']='REDIST'
    hdu_mat.header['CHANTYPE']='PHA'
    hdu_mat.header['HDUVERS']='1.3.0'
    hdu_mat.header['TELESCOP']='NICER'
    hdu_mat.header['INSTRUME']='XTI'
    hdu_mat.header['FILTER']='NONE'
    hdu_mat.header['DETCHANS']=n_ch
    hdu_mat.header['TLMIN4']=0

    # ----- Write ----- #
    hdus=fits.HDUList([hdu_pri, hdu_ebo, hdu_mat])
    hdus.writeto(inputs_xspec.name_rmf, overwrite=True)

if __name__=='__main__':
    main()
