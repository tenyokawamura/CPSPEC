import numpy as np
import astropy.io.fits as fits

class PowerSpectrum:
    def __init__(self):
        self.first_psd=True

    def calc_psd(self, bs, rate_mea, rate_var_mea):
        # Normalization for (rms/mean)^2/Hz
        norm=2.*self.dt/((rate_mea**2)*self.n_bin)
        psds=norm*(np.abs(bs)**2)

        if self.first_psd==True:
            self.rates_mea    =np.array([rate_mea])
            self.rates_var_mea=np.array([rate_var_mea]) # Necessary for Gaussian noise (2022/02/02)
            self.psdss=np.array([psds])
            self.first_psd=False
        else:
            self.rates_mea    =np.append(self.rates_mea,     rate_mea)
            self.rates_var_mea=np.append(self.rates_var_mea, rate_var_mea)
            self.psdss=np.vstack((self.psdss, psds))

    def average_psd(self):
        self.n_int=len(self.psdss)
        self.dfs=(self.fs-np.roll(self.fs, 1))/2.
        self.dfs[0]=(self.fs[1]-self.fs[0])/2.
        self.fs_min=self.fs-self.dfs
        self.fs_max=self.fs+self.dfs
        # ----- Logarithmic rebin (begin) ----- #
        self.fs_min_rebin=np.empty(0)
        self.fs_max_rebin=np.empty(0)
        self.psds_raw_mea=np.empty(0)
        self.psds_raw_sig=np.empty(0)
        self.psds_noi    =np.empty(0)
        self.ns_fbin     =np.empty(0)
        f_mid_min =0.
        f_mid_max =self.fs[0]
        psdss_t=self.psdss.T
        while True:
            is_merge=np.where((f_mid_min<self.fs)&(self.fs<=f_mid_max))[0]
            #print(f_mid_min, f_mid_max, is_merge)
            if len(is_merge)==0:
                if f_mid_min>self.fs[-1]:
                    break
            else:
                n_fbin=len(is_merge)
                fs_min_merge=self.fs_min[is_merge]
                fs_max_merge=self.fs_max[is_merge]
                psds_merge  =np.ravel(psdss_t[is_merge])

                f_min_rebin=fs_min_merge[0]
                f_max_rebin=fs_max_merge[-1]
                psd_raw_mea=np.mean(psds_merge)
                psd_raw_sig=np.sqrt(np.var(psds_merge)/len(psds_merge))

                self.ns_fbin     =np.append(self.ns_fbin     , n_fbin)
                self.fs_min_rebin=np.append(self.fs_min_rebin, f_min_rebin)
                self.fs_max_rebin=np.append(self.fs_max_rebin, f_max_rebin)
                self.psds_raw_mea=np.append(self.psds_raw_mea, psd_raw_mea)
                self.psds_raw_sig=np.append(self.psds_raw_sig, psd_raw_sig)

            f_mid_min=f_mid_max
            f_mid_max*=self.rebin
        # ----- Logarithmic rebin (end)   ----- #
        
        self.fs_rebin =(self.fs_min_rebin+self.fs_max_rebin)/2.
        self.dfs_rebin=(self.fs_max_rebin-self.fs_min_rebin)/2.
        self.ns_int  =self.n_int*np.ones(len(self.fs_min_rebin))
        ######################################
        ### Assume purely Poissonian noise ###
        ######################################
        #self.psds_noi=(2./np.mean(self.rates_mea))*np.ones(len(self.fs_min_rebin))
        #self.psds_noi=np.mean(2./self.rates_mea)*np.ones(len(self.fs_min_rebin)) # This should be more reasonable.
        #######################################################
        ### Assume Gaussian noise (see Vaughan et al. 2003) ###
        #######################################################
        self.psds_noi=np.mean(2.*self.dt*self.rates_var_mea/(self.rates_mea**2))*np.ones(len(self.fs_min_rebin)) 

    def write_psd(\
        self,\
        ch_min,\
        ch_max,\
        name_fits,\
        telescope,\
        instrument,\
        source,\
        exposure):

        self.name_fits=name_fits

        names=['F', 'DF', 'PSDRM', 'PSDRS', 'PSDN', 'NFBIN', 'NINT']
        comments=[\
            'Minumum frequency',\
            'Maximum frequency',\
            'Raw PSD (mean)',\
            'Raw PSD (sigma)',\
            'Noise PSD',\
            'Number of frequency bins merged',\
            'Number of intervals merged']
        units=['Hz', 'Hz', '(rms/mean)^2/Hz', '(rms/mean)^2/Hz', '(rms/mean)^2/Hz', '-', '-']

        hdu_ext=fits.BinTableHDU.from_columns([\
            fits.Column(name=names[0], unit=units[0], format='f4', array=self.fs_rebin),\
            fits.Column(name=names[1], unit=units[1], format='f4', array=self.dfs_rebin),\
            fits.Column(name=names[2], unit=units[2], format='f4', array=self.psds_raw_mea),\
            fits.Column(name=names[3], unit=units[3], format='f4', array=self.psds_raw_sig),\
            fits.Column(name=names[4], unit=units[4], format='f4', array=self.psds_noi),\
            fits.Column(name=names[5], unit=units[5], format='i4', array=self.ns_fbin),\
            fits.Column(name=names[6], unit=units[6], format='i4', array=self.ns_int)])

        for i_n, (name, comment) in enumerate(zip(names, comments)):
            hdu_ext.header['TTYPE{0}'.format(int(i_n+1))]=(name, comment)

        hdu_ext.header['CHMIN']   =(ch_min      , 'Minimum channel')
        hdu_ext.header['CHMAX']   =(ch_max      , 'Minimum channel')
        hdu_ext.header['DT']      =(self.dt     , 'Sampling interval')
        hdu_ext.header['NBIN']    =(self.n_bin  , 'Number of bin per interval')
        hdu_ext.header['REBIN']   =(self.rebin  , 'Logarithmic rebinning')
        #hdu_ext.header['RATEM']   =(self.x_mea  , 'Mean count rate [count s^-1]')
        hdu_ext.header['TELESCOP']=(telescope   , 'Telescope')
        hdu_ext.header['INSTRUME']=(instrument  , 'Instrument')
        hdu_ext.header['OBJECT']  =(source      , 'Object')
        hdu_ext.header['EXPOSURE']=(exposure    , '[s]')

        # --- Primaru=y HDU --- #
        hdu_pri=fits.PrimaryHDU(data=None, header=None)
        hdu=fits.HDUList([hdu_pri, hdu_ext])
        hdu.writeto(name_fits, overwrite=True)

    def read_psd(self, name_file):
        hdus=fits.open(name_file)
        fs          =hdus[1].data['F'] # Index 0 ... Primary HDU
        psds_raw_mea=hdus[1].data['PSDRM']
        psds_raw_sig=hdus[1].data['PSDRS']
        psds_noi    =hdus[1].data['PSDN']
        ns_fbin     =hdus[1].data['NFBIN']
        ns_int      =hdus[1].data['NINT']
        hdus.close()
        return fs, psds_raw_mea, psds_raw_sig, psds_noi, ns_fbin, ns_int

