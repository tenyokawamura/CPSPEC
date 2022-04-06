import numpy as np
import astropy.io.fits as fits
import sys

class CrossSpectrum:
    def __init__(self):
        self.first_csd=True

    def calc_csdf(self, bs, rate_mea, rate_var_mea, bs_ref, rate_ref_mea, rate_ref_var_mea):
        # Normalization for (rms/mean)^2/Hz
        norm=2.*self.dt/(rate_mea*rate_ref_mea*self.n_bin)

        # Noise contribution due to overlapping of energy band
        if self.overlap==True:
            ######################################
            ### Assume purely Poissonian noise ###
            ######################################
            csd_noi_overlap=2./rate_ref_mea
            #######################################################
            ### Assume Gaussian noise (see Vaughan et al. 2003) ###
            #######################################################
            csd_noi_overlap=2.*self.dt*rate_ref_var_mea/(rate_ref_mea**2)
        else:
            csd_noi_overlap=0.

        # Calculate cross spectrum
        #csds=norm*(bs_ref.conjugate())*bs-csd_noi_overlap
        # Noise should not be subtracted at this stage (2022/02/09)
        csds=norm*(bs_ref.conjugate())*bs

        # Store data
        if self.first_csd==True:
            self.rates_mea   =rate_mea
            self.rate_ref_mea=rate_ref_mea
            self.csdss       =csds
            self.csds_noi_overlap=csd_noi_overlap
            self.first_csd   =False
        else:
            self.rates_mea=np.append(self.rates_mea,     rate_mea)
            self.csdss    =np.vstack((self.csdss, csds))
            self.csds_noi_overlap=np.append(self.csds_noi_overlap, csd_noi_overlap)

    # Calculate only mean values of cross spectra (not thier errors)
    def average_csdf(self):
        self.n_int=len(self.csdss)
        self.dfs=(self.fs-np.roll(self.fs, 1))/2.
        self.dfs[0]=(self.fs[1]-self.fs[0])/2.
        self.fs_min=self.fs-self.dfs
        self.fs_max=self.fs+self.dfs
        # ----- Logarithmic rebin (begin) ----- #
        self.fs_min_rebin=np.empty(0)
        self.fs_max_rebin=np.empty(0)
        self.csds_re_mea =np.empty(0)
        self.csds_im_mea =np.empty(0)
        self.csds_ab_mea =np.empty(0)
        self.csds_ph_mea =np.empty(0)
        self.csds_ti_mea =np.empty(0)
        self.ns_fbin     =np.empty(0)
        f_mid_min =0.
        f_mid_max =self.fs[0]
        csdss_t=self.csdss.T
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
                csds_merge  =np.ravel(csdss_t[is_merge])

                f_min_rebin=fs_min_merge[0]
                f_max_rebin=fs_max_merge[-1]
                f_mid_rebin=(f_min_rebin+f_max_rebin)/2.
                csd_mea=np.mean(csds_merge)
                csd_re_mea=np.real(csd_mea)
                csd_im_mea=np.imag(csd_mea)
                csd_ab_mea=np.abs(csd_mea)
                csd_ph_mea=np.angle(csd_mea)
                csd_ti_mea=csd_ph_mea/(2.*np.pi*f_mid_rebin)

                self.ns_fbin     =np.append(self.ns_fbin     , n_fbin)
                self.fs_min_rebin=np.append(self.fs_min_rebin, f_min_rebin)
                self.fs_max_rebin=np.append(self.fs_max_rebin, f_max_rebin)
                self.csds_re_mea =np.append(self.csds_re_mea , csd_re_mea)
                self.csds_im_mea =np.append(self.csds_im_mea , csd_im_mea)
                self.csds_ab_mea =np.append(self.csds_ab_mea , csd_ab_mea)
                self.csds_ph_mea =np.append(self.csds_ph_mea , csd_ph_mea)
                self.csds_ti_mea =np.append(self.csds_ti_mea , csd_ti_mea)

            f_mid_min=f_mid_max
            f_mid_max*=self.rebin
        # ----- Logarithmic rebin (end)   ----- #
        
        self.fs_rebin =(self.fs_min_rebin+self.fs_max_rebin)/2.
        self.dfs_rebin=(self.fs_max_rebin-self.fs_min_rebin)/2.
        self.ns_int   =self.n_int*np.ones(len(self.fs_min_rebin))
        self.csds_noi =np.mean(self.csds_noi_overlap)*np.ones(len(self.fs_min_rebin))

    def calc_csdf_error(\
        self,\
        psds_raw_mea,\
        psds_noi,\
        ns_fbin,\
        ns_int,\
        psds_ref_raw_mea,\
        psds_ref_noi,\
        ns_ref_fbin,\
        ns_ref_int):
        # Check consistency on rebinning
        if ((np.all(ns_fbin==ns_ref_fbin)) & (np.all(ns_ref_fbin==self.ns_fbin)) &\
            (np.all(ns_int==ns_ref_int)) &  (np.all(ns_ref_int==self.ns_int)))==False:
            print('Error: Rebinning inconsistency')
            sys.exit()

        # d(Re[CSD(f)]), d(Im[CSD(f)]), d(|CSD(f)|)
        numes=psds_raw_mea*psds_ref_raw_mea+(self.csds_re_mea**2)-(self.csds_im_mea**2)
        numes*=(numes>=0) # Numerator in \sqrt must not be negative.
        self.csds_re_sig=np.sqrt(numes/(2.*self.ns_fbin*self.ns_int))

        numes=psds_raw_mea*psds_ref_raw_mea-(self.csds_re_mea**2)+(self.csds_im_mea**2)
        numes*=(numes>=0) # Numerator in \sqrt must not be negative.
        self.csds_im_sig=np.sqrt(numes/(2.*self.ns_fbin*self.ns_int))

        self.csds_ab_sig=np.sqrt(psds_raw_mea*psds_ref_raw_mea/(self.ns_fbin*self.ns_int))

        # d(\phi(f))
        first=True
        for psd_raw_mea, psd_noi, psd_ref_raw_mea, psd_ref_noi, csd_ab_mea, n_fbin, n_int \
            in zip(psds_raw_mea, psds_noi, psds_ref_raw_mea, psds_ref_noi, self.csds_ab_mea, ns_fbin, ns_int):

            b2, g2=calc_bias_coherence(\
                psd_raw_mea    =psd_raw_mea,\
                psd_noi        =psd_noi,\
                psd_ref_raw_mea=psd_ref_raw_mea,\
                psd_ref_noi    =psd_ref_noi,\
                csd_ab_mea     =csd_ab_mea,\
                n_fbin         =n_fbin,\
                n_int          =n_int)

            csd_ph_sig=np.sqrt((1.-g2)/(2.*g2*n_fbin*n_int))
            if first==True:
                first=False
                self.csds_ph_sig=csd_ph_sig
                self.b2s        =b2
                self.g2s        =g2
            else:
                self.csds_ph_sig=np.append(self.csds_ph_sig, csd_ph_sig)
                self.b2s        =np.append(self.b2s, b2)
                self.g2s        =np.append(self.g2s, g2)

        # d(\tau(f))
        self.csds_ti_sig=self.csds_ph_sig/(2.*np.pi*self.fs_rebin)

    def write_csdf(\
        self,\
        ch_min,\
        ch_max,\
        name_fits,\
        telescope,\
        instrument,\
        source,\
        exposure):

        self.name_fits=name_fits

        names=[\
            'F'     , 'DF',\
            'CSDRM' , 'CSDRS',\
            'CSDIM' , 'CSDIS',\
            'CSDAM' , 'CSDAS',\
            'CSDPM' , 'CSDPS',\
            'CSDTM' , 'CSDTS',\
            'CSDN'  ,\
            'GAMMA2', 'B2',\
            'NFBIN' , 'NINT']
        comments=[\
            'Minumum frequency',\
            'Maximum frequency',\
            'Re[CSD] (mean)',\
            'Re[CSD] (sigma)',\
            'Im[CSD] (mean)',\
            'Im[CSD] (sigma)',\
            '|CSD|   (mean)',\
            '|CSD|   (sigma)',\
            'Phase   (mean)',\
            'Phase   (sigma)',\
            'Time    (mean)',\
            'Time    (sigma)',\
            'Noise CSD',\
            'Coherence^2',\
            'Bias term',\
            'Number of frequency bins merged',\
            'Number of intervals merged']
        units=[\
            'Hz'             , 'Hz',\
            '(rms/mean)^2/Hz', '(rms/mean)^2/Hz',\
            '(rms/mean)^2/Hz', '(rms/mean)^2/Hz',\
            '(rms/mean)^2/Hz', '(rms/mean)^2/Hz',\
            'radian'         , 'radian',\
            'second'         , 'second',\
            '(rms/mean)^2/Hz',\
            '-'              , '-',\
            '-'              , '-']

        hdu_ext=fits.BinTableHDU.from_columns([\
            fits.Column(name=names[0],  unit=units[0],  format='f4', array=self.fs_rebin),\
            fits.Column(name=names[1],  unit=units[1],  format='f4', array=self.dfs_rebin),\
            fits.Column(name=names[2],  unit=units[2],  format='f4', array=self.csds_re_mea),\
            fits.Column(name=names[3],  unit=units[3],  format='f4', array=self.csds_re_sig),\
            fits.Column(name=names[4],  unit=units[4],  format='f4', array=self.csds_im_mea),\
            fits.Column(name=names[5],  unit=units[5],  format='f4', array=self.csds_im_sig),\
            fits.Column(name=names[6],  unit=units[6],  format='f4', array=self.csds_ab_mea),\
            fits.Column(name=names[7],  unit=units[7],  format='f4', array=self.csds_ab_sig),\
            fits.Column(name=names[8],  unit=units[8],  format='f4', array=self.csds_ph_mea),\
            fits.Column(name=names[9],  unit=units[9],  format='f4', array=self.csds_ph_sig),\
            fits.Column(name=names[10], unit=units[10], format='f4', array=self.csds_ti_mea),\
            fits.Column(name=names[11], unit=units[11], format='f4', array=self.csds_ti_sig),\
            fits.Column(name=names[12], unit=units[12], format='f4', array=self.csds_noi),\
            fits.Column(name=names[13], unit=units[13], format='f4', array=self.g2s),\
            fits.Column(name=names[14], unit=units[14], format='f4', array=self.b2s),\
            fits.Column(name=names[15], unit=units[15], format='i4', array=self.ns_fbin),\
            fits.Column(name=names[16], unit=units[16], format='i4', array=self.ns_int)])

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

def calc_bias_coherence(
    psd_raw_mea,\
    psd_noi,\
    psd_ref_raw_mea,\
    psd_ref_noi,\
    csd_ab_mea,\
    n_fbin,\
    n_int):
    g2=1 # Initial condition of the intrinsic coherence
    num_loop=0
    while True:
        num_loop+=1
        if n_fbin*n_int>=500:
            b2_now=0
        else:
            b2_now=calc_bias(\
                psd_raw_mea    =psd_raw_mea,\
                psd_noi        =psd_noi,\
                psd_ref_raw_mea=psd_ref_raw_mea,\
                psd_ref_noi    =psd_ref_noi,\
                n_fbin         =n_fbin,\
                n_int          =n_int,\
                g2             =g2)
            if csd_ab_mea**2<b2_now:
                b2_now=0

        g2_now=calc_raw_coherence(\
            psd_raw_mea    =psd_raw_mea,\
            psd_ref_raw_mea=psd_ref_raw_mea,\
            csd_ab_mea     =csd_ab_mea,\
            b2             =b2_now)
        g2_now=np.min([g2_now, 1]) # Coherence is not greater than 1.

        # Check convergence
        if np.abs(g2-g2_now)<1.e-3:
            g2=g2_now
            b2=b2_now
            break
        else:
            g2=g2_now
            b2=b2_now

        if num_loop>=100:
            print('Warning: Iteration loop was not finished.')
            b2=0
            g2=calc_raw_coherence(\
                psd_raw_mea    =psd_raw_mea,\
                psd_ref_raw_mea=psd_ref_raw_mea,\
                csd_ab_mea     =csd_ab_mea,\
                b2             =b2_now)
            g2=np.min([g2, 1]) # Coherence is not greater than 1.
            break
    return b2, g2

def calc_bias(\
    psd_raw_mea,\
    psd_noi,\
    psd_ref_raw_mea,\
    psd_ref_noi,\
    n_fbin,\
    n_int,\
    g2):
    b2=(psd_raw_mea*psd_ref_raw_mea-g2*(np.amax([psd_raw_mea-psd_noi, 0]))\
        *(np.amax([psd_ref_raw_mea-psd_ref_noi, 0])))\
        /(n_fbin*n_int)
    return b2

def calc_raw_coherence(\
    psd_raw_mea,\
    psd_ref_raw_mea,\
    csd_ab_mea,\
    b2):
    g2=((csd_ab_mea**2)-b2)/(psd_raw_mea*psd_ref_raw_mea)
    return g2
