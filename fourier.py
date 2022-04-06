import numpy as np
import astropy.io.fits as fits

class Fourier:
    def __init__(self):
        pass

    def fft_calc(self, xs, dxs):
        self.n_bin=len(xs)
        self.x_mea=np.mean(xs)
        xs=xs-self.x_mea

        # 2022/02/02, Necessary for calculating power spectrum of Gaussian noise
        # ------------------------------------------------- #
        # --- Average of variance of count rate (start) --- #
        # ------------------------------------------------- #
        dx2s=dxs**2
        self.dx2_mea=np.mean(dx2s)
        # ------------------------------------------------- #
        # --- Average of variance of count rate (end)   --- #
        # ------------------------------------------------- #

        bs=np.fft.fft(a=xs)
        fs=np.fft.fftfreq(n=self.n_bin, d=self.dt)
        if self.n_bin%2==0:
            i_max=int((self.n_bin/2))
        else:
            i_max=int((self.n_bin-1)/2)
        self.bs=bs[1:i_max]
        self.fs=fs[1:i_max]
        self.f_min=self.fs[0]
        self.f_max=self.fs[-1]

    def write_ft(\
        self,\
        ch_min,\
        ch_max,\
        name_file,\
        telescope,\
        instrument,\
        source,\
        exposure):

        self.name_file=name_file

        names=['FREQ', 'FT']
        comments=['Frequency', 'Fourier transform']
        units=['Hz', 'count s^-1']

        hdu_ext=fits.BinTableHDU.from_columns([\
            fits.Column(name=names[0], unit=units[0], format='f4', array=self.fs),\
            fits.Column(name=names[1], unit=units[1], format='c8', array=self.bs)])

        for i_n, (name, comment) in enumerate(zip(names, comments)):
            hdu_ext.header['TTYPE{0}'.format(int(i_n+1))]=(name, comment)

        hdu_ext.header['CHMIN']   =(ch_min      , 'Minimum channel')
        hdu_ext.header['CHMAX']   =(ch_max      , 'Minimum channel')
        hdu_ext.header['TSTART']  =(self.t_start, '[s]')
        hdu_ext.header['TEND']    =(self.t_end  , '[s]')
        hdu_ext.header['DT']      =(self.dt     , 'Sampling interval')
        hdu_ext.header['NBIN']    =(self.n_bin  , 'Number of bin per interval')
        hdu_ext.header['FMIN']    =(self.f_min  , 'Minimum frequency [Hz]')
        hdu_ext.header['FMAX']    =(self.f_max  , 'Maximum frequency [Hz]')
        hdu_ext.header['INDINT']  =(self.i_int  , 'Index of interval')
        hdu_ext.header['RATEM']   =(self.x_mea  , 'Mean count rate [count s^-1]')
        hdu_ext.header['VRATEM']  =(self.dx2_mea, 'Mean variance of count rate [count^2 s^-2]')
        hdu_ext.header['TELESCOP']=(telescope   , 'Telescope')
        hdu_ext.header['INSTRUME']=(instrument  , 'Instrument')
        hdu_ext.header['OBJECT']  =(source      , 'Object')
        hdu_ext.header['EXPOSURE']=(exposure    , '[s]')

        if self.i_int==0:
            # --- Primaru=y HDU --- #
            hdu_pri=fits.PrimaryHDU(data=None, header=None)
            hdu=fits.HDUList([hdu_pri, hdu_ext])
            hdu.writeto(name_file, overwrite=True)
        else:
            with fits.open(name_file, mode='update', memmap=True) as hdu:
                hdu.append(hdu_ext)
                hdu.flush()

    def read_info(
        self,\
        name_file):

        self.name_file=name_file
        hdus=fits.open(name_file)
        # Number of interval
        self.n_int=len(hdus)-1
        self.fs        =hdus[1].data['FREQ']
        self.ch_min    =hdus[1].header['CHMIN']
        self.ch_max    =hdus[1].header['CHMAX']
        self.telescope =hdus[1].header['TELESCOP']
        self.instrument=hdus[1].header['INSTRUME']
        self.source    =hdus[1].header['OBJECT']
        self.exposure  =hdus[1].header['EXPOSURE']
        self.dt        =hdus[1].header['DT']
        self.n_bin     =hdus[1].header['NBIN']
        self.f_min     =hdus[1].header['FMIN']
        self.f_max     =hdus[1].header['FMAX']

    def print_info(self):
        print('------------------------------')
        print('{0:<32}: {1}'       .format('Filename'                   , self.name_file))
        print('{0:<32}: {1}'       .format('Telescope'                  , self.telescope))
        print('{0:<32}: {1}'       .format('Instrument'                 , self.instrument))
        print('{0:<32}: {1}'       .format('Minimum Channel'            , self.ch_min))
        print('{0:<32}: {1}'       .format('Maximum Channel'            , self.ch_max))
        print('{0:<32}: {1}'       .format('Source'                     , self.source))
        print('{0:<32}: {1:.0f} s' .format('Exposure'                   , self.exposure))
        print('{0:<32}: {1:.2g} s' .format('Sampling interval'          , self.dt))
        print('{0:<32}: {1}'       .format('Number of bins per interval', self.n_bin))
        print('{0:<32}: {1}'       .format('Number of intervals'        , self.n_int))
        print('{0:<32}: {1:.2g} Hz'.format('Minimum frequency'          , self.f_min))
        print('{0:<32}: {1:.2g} Hz'.format('Maximum frequency'          , self.f_max))
        print('')

    def read_ft(self, i_int):
        hdus=fits.open(self.name_file)
        fs=hdus[i_int+1].data['FREQ'] # Index 0 ... Primary HDU
        bs=hdus[i_int+1].data['FT']
        rate_mea=hdus[i_int+1].header['RATEM']
        rate_var_mea=hdus[i_int+1].header['VRATEM']
        hdus.close()
        return fs, bs, rate_mea, rate_var_mea

