import astropy.io.fits as fits
import numpy as np 
import os
import sys 
from sys import argv

def main():
    # ----------------------------- #
    # ---------- Setting ---------- #
    # ----------------------------- #
    # ----- File ----- #
    name_outarf='response_nicer.arf'
    f_min=1.e-3 #[Hz]
    f_max=8.e1  #[Hz]
    n_bin=1000 #[-]

    # -------------------------- #
    # ---------- Main ---------- #
    # -------------------------- #
    # ----- Set time/frequency ----- #
    fs=np.geomspace(f_min, f_max, n_bin+1)

    fs_min=fs
    fs_max=np.roll(fs, -1)
    # Delete the maximum bin
    fs_min=fs_min[:-1]
    fs_max=fs_max[:-1]

    n_ch=len(fs_min)

    # -------------------------- #
    # ---------- Main ---------- #
    # -------------------------- #
    # ----- PRIMARY ----- #
    hdu_pri=fits.PrimaryHDU(data=None, header=None)

    # ----- SPECRESP ----- #
    es_lo=fs_min #[keV]
    es_hi=fs_max #[keV]
    n_e=len(es_lo)
    print(n_e)
    areas=np.ones(n_e)

    hdu_rsp\
    =fits.BinTableHDU.from_columns([fits.Column(name='ENERG_LO', format='f4', array=es_lo),\
                                    fits.Column(name='ENERG_HI', format='f4', array=es_hi),\
                                    fits.Column(name='SPECRESP', format='f4', array=areas)])

    hdu_rsp.header['TUNIT1']='keV'
    hdu_rsp.header['TUNIT2']='keV'
    hdu_rsp.header['TUNIT3']='cm**2'
    hdu_rsp.header['EXTNAME']='SPECRESP'
    hdu_rsp.header['HDUCLASS']='OGIP'
    hdu_rsp.header['HDUCLAS1']='RESPONSE'
    hdu_rsp.header['HDUCLAS2']='SPECRESP'
    hdu_rsp.header['HDUVERS']='1.1.0'
    hdu_rsp.header['TELESCOP']='NICER'
    hdu_rsp.header['INSTRUME']='NICER'
    hdu_rsp.header['FILTER']='NONE'

    # ----- Write ----- #
    hdus=fits.HDUList([hdu_pri, hdu_rsp])
    hdus.writeto(name_outarf, overwrite=True)

if __name__=='__main__':
    main()
