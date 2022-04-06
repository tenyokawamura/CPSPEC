import numpy as np
from numpy import random
import astropy.io.fits as fits
import sys

class LightCurve:
    def __init__(self):
        self.i_bin_start=0

    def set_par(self, n_bin, n_int, maximize, frac_gap):
        # Number of bins per interval
        self.n_bin=n_bin
        # Number of intervals
        self.n_int=n_int
        # Maximize intervals or not
        self.maximize=maximize
        # Fraction of acceptable absent bins
        self.frac_gap=frac_gap
        # Number of acceptable absent bins
        self.n_gap=self.n_bin*self.frac_gap

    def read_data(self, name_file):
        self.name_file=name_file
        hdus=fits.open(self.name_file)
        self.header=hdus[1].header
        self.data  =hdus[1].data
        self.telescope =self.header['TELESCOP']
        self.instrument=self.header['INSTRUME']
        self.source    =self.header['OBJECT']
        self.exposure  =self.header['EXPOSURE']
        self.dt        =self.header['TIMEDEL']
        self.n_row     =self.header['NAXIS2']
        hdus.close()

    def print_obs_info(self):
        print('------------------------------')
        print('(Observation information)')
        print('{0:<32}: {1}'      .format('Filename'         , self.name_file))
        print('{0:<32}: {1}'      .format('Telescope'        , self.telescope))
        print('{0:<32}: {1}'      .format('Instrument'       , self.instrument))
        print('{0:<32}: {1}'      .format('Source'           , self.source))
        print('{0:<32}: {1:.0f} s'.format('Exposure'         , self.exposure))
        print('{0:<32}: {1:.2g} s'.format('Sampling interval', self.dt))
        print('{0:<32}: {1}'      .format('Number of rows'   , self.n_row))
        print('')

    def convert_counts2rate(self):
        self.rates_rebin=self.data['COUNTS']/\
            ((np.roll(self.data['TIME'], -1)-self.data['TIME'])*self.data['FRACEXP'])
        self.errors_rebin=self.data['ERROR']/\
            ((np.roll(self.data['TIME'], -1)-self.data['TIME'])*self.data['FRACEXP'])
        self.times_rebin   =self.data['TIME'][:-1]
        self.rates_rebin   =self.rates_rebin[:-1]
        self.errors_rebin  =self.errors_rebin[:-1]
        self.fracexps_rebin=self.data['FRACEXP'][:-1]
        if len(self.data[0])==5:
            self.deadcs_rebin=self.data['DEADC'][:-1]

    def rebin_lightcurve(self, dt_rebin):
        i_start=0
        t_start=self.data['TIME'][i_start]
        if len(self.data[0])==4:
            for i_row in range(self.n_row):
                dt=self.data['TIME'][i_row]-t_start
                if dt>=dt_rebin:
                    i_end=i_row
                    rates   =self.data['RATE'   ][i_start:i_end]
                    errors  =self.data['ERROR'  ][i_start:i_end]
                    fracexps=self.data['FRACEXP'][i_start:i_end]
                    n_rebin=len(rates)
                    time_rebin   =self.data['TIME'][i_start]
                    rate_rebin   =np.mean(rates)
                    error_rebin  =np.sqrt(np.sum(errors**2))/n_rebin # Propagation of error
                    fracexp_rebin=np.mean(fracexps)
                    
                    if i_start==0:
                        self.times_rebin   =time_rebin
                        self.rates_rebin   =rate_rebin
                        self.errors_rebin  =error_rebin
                        self.fracexps_rebin=fracexp_rebin
                    else:
                        self.times_rebin   =np.append(self.times_rebin,    time_rebin)
                        self.rates_rebin   =np.append(self.rates_rebin,    rate_rebin)
                        self.errors_rebin  =np.append(self.errors_rebin,   error_rebin)
                        self.fracexps_rebin=np.append(self.fracexps_rebin, fracexp_rebin)

                    i_start=i_row
                    t_start=self.data['TIME'][i_start]

        elif len(self.data[0])==5:
            for i_row in range(self.n_row):
                dt=self.data['TIME'][i_row]-t_start
                if dt>=dt_rebin:
                    i_end=i_row
                    rates   =self.data['RATE'   ][i_start:i_end]
                    errors  =self.data['ERROR'  ][i_start:i_end]
                    fracexps=self.data['FRACEXP'][i_start:i_end]
                    deadcs  =self.data['DEADC'  ][i_start:i_end]
                    n_rebin=len(rates)
                    time_rebin   =self.data['TIME'][i_start]
                    rate_rebin   =np.mean(rates)
                    error_rebin  =np.sqrt(np.sum(errors**2))/n_rebin # Propagation of error
                    fracexp_rebin=np.mean(fracexps)
                    deadc_rebin  =np.mean(deadcs)
                    
                    if i_start==0:
                        self.times_rebin   =time_rebin
                        self.rates_rebin   =rate_rebin
                        self.errors_rebin  =error_rebin
                        self.fracexps_rebin=fracexp_rebin
                        self.deadcs_rebin  =deadc_rebin
                    else:
                        self.times_rebin   =np.append(self.times_rebin,    time_rebin)
                        self.rates_rebin   =np.append(self.rates_rebin,    rate_rebin)
                        self.errors_rebin  =np.append(self.errors_rebin,   error_rebin)
                        self.fracexps_rebin=np.append(self.fracexps_rebin, fracexp_rebin)
                        self.deadcs_rebin  =np.append(self.deadcs_rebin,   deadc_rebin)

                    i_start=i_row
                    t_start=self.data['TIME'][i_start]

    def write_lightcurve(self, name_file):
        hdus_in=fits.open(self.name_file)
            
        ########## PRIMARY (Column 0) ##########
        hdu0=fits.PrimaryHDU(data=hdus_in[0].data)
        hdu0.header=hdus_in[0].header

        ########## Column 1 ##########
        # D: double precision float (64-bit), E: single precision float (32-bit)
        if len(self.data[0])==4:
            hdu1=fits.BinTableHDU.from_columns([\
                fits.Column(name='TIME'   , format='D', array=self.times_rebin),\
                fits.Column(name='RATE'   , format='E', array=self.rates_rebin),\
                fits.Column(name='ERROR'  , format='E', array=self.errors_rebin),\
                fits.Column(name='FRACEXP', format='E', array=self.fracexps_rebin)])
        elif len(self.data[0])==5:
            hdu1=fits.BinTableHDU.from_columns([\
                fits.Column(name='TIME'   , format='D', array=self.times_rebin),\
                fits.Column(name='RATE'   , format='E', array=self.rates_rebin),\
                fits.Column(name='ERROR'  , format='E', array=self.errors_rebin),\
                fits.Column(name='FRACEXP', format='E', array=self.fracexps_rebin),\
                fits.Column(name='DEADC'  , format='E', array=self.deadcs_rebin)])
        hdu1.header=hdus_in[1].header
        # Time resolution
        hdu1.header['TIMEDEL']=np.amin(self.times_rebin[1:]-self.times_rebin[:-1])

        ########## Column 2 ##########
        hdu2=fits.BinTableHDU(data=hdus_in[2].data)
        hdu2.header=hdus_in[2].header

        ########## Column 3 ##########
        hdu3=fits.BinTableHDU(data=hdus_in[3].data)
        hdu3.header=hdus_in[3].header

        hdus_in.close()

        hdus=fits.HDUList([hdu0, hdu1, hdu2, hdu3])
        hdus.writeto(name_file, overwrite=True)

    def prepare_ft(self):
        self.n_int_max=self.n_row//self.n_bin
        self.c_t      =self.n_bin*self.dt #Capital t = T [s]
        # Minimum frequency [Hz]
        self.f_min    =1/self.c_t
        # Maximum (Nyquist) frequency [Hz]
        self.f_max    =1/(2.*self.dt)

    def print_ana_info(self):
        print('------------------------------')
        print('(Analysis information)')
        print('{0:<32}: {1:.2g} s' .format('Sampling interval'           , self.dt))
        print('{0:<32}: {1}'       .format('Number of bins per interval' , self.n_bin))
        print('{0:<32}: {1}'       .format('Maximum number of intervals' , self.n_int_max))
        print('{0:<32}: {1}'       .format('Proposed number of intervals', self.n_int))
        print('{0:<32}: {1}'       .format('Maximize intervals'          , self.maximize))
        print('{0:<32}: {1:.2g} Hz'.format('Minimum frequency'           , self.f_min))
        print('{0:<32}: {1:.2g} Hz'.format('Maximum frequency'           , self.f_max))
        print('')

    def extract_interval(self, i_int):
        # Nubmer of gaps in one segment (for counting)
        num_gap=0 
        # Index of gaps in one segment
        is_gap=np.empty(0, dtype=np.int8) 
        # Number of gaps in each index
        ns_gap=np.empty(0, dtype=np.int8) 
        # Number of data in one segment (for counting)
        num_bin=0 
        print('(Interval {0})'.format(i_int+1), end=' ', flush=True)
        #for i_d in range(len(data)):
        i_bin=self.i_bin_start
        i_bin_max=len(self.data)-1
        while True:
            if i_bin>i_bin_max:
                return None, None, None, None
            t_bin=self.data[i_bin][0]
            if i_bin==self.i_bin_start:
                t_start=t_bin
            t_exp=t_start+self.dt*num_bin
            if abs(t_bin-t_exp)<self.dt:
                num_bin+=1
                if num_bin==self.n_bin:
                    i_bin_end=i_bin+1
                    t_end=t_bin
                    data_int_real=self.data[self.i_bin_start:i_bin_end]
                    n_bin_real=len(data_int_real)

                    #n_d_real=len(data_int_real)
                    #xs=np.zeros(n_d_real)
                    #for i_d_real in range(n_d_real):
                    #    xs[i_d_real]=data_int_real[i_d_real][1]
                    rates_int =data_int_real['RATE']
                    drates_int=data_int_real['ERROR'] # 2022/02/02, necessary for Gaussian noise

                    # --------------------------------------------------- #
                    # ----- Simulate and insert random data (begin) ----- #
                    # --------------------------------------------------- #
                    ###############################################
                    ### self.n_gap must be zero for now         ###
                    ### because drates_int is not handled here. ###
                    ### (2022/02/02)                            ###
                    ###############################################
                    if self.n_gap!=0:
                        #x_mea, x_var=mea_var_calc_1d(xs=xs)
                        #x_sig=np.sqrt(x_var)
                        rate_mea=np.mean(rates_int)
                        rate_var=np.var(rates_int)
                        rate_sig=np.sqrt(rate_var)
                        for i_gap, n_gap in zip(is_gap[::-1], ns_gap[::-1]):
                            rates_sim=np.zeros(n_gap)
                            for rate_sim in rates_sim:
                                rate_sim=random.normal(rate_mea, rate_sig)
                            i_x=i_gap-i_bin_start
                            rates_int=np.insert(rates_int, i_x, rates_sim)
                    # --------------------------------------------------- #
                    # ----- Simulate and insert random data (end)   ----- #
                    # --------------------------------------------------- #

                    if not len(rates_int)==self.n_bin:
                        print('Error')
                        sys.exit()
                    print('\nNumber of actual data: {0}'.format(n_bin_real))
                    print('Number of simulated data: {0}'.format(num_gap))

                    self.i_bin_start=i_bin+1

                    return t_start, t_end, rates_int, drates_int

            else:
                t_g=t_bin-t_exp
                n_gap=int(round((t_bin-t_exp)/self.dt-1.))
                num_gap=num_gap+n_gap
                if num_gap>self.n_gap:
                    print('.', end='', flush=True)
                    num_bin=0
                    num_gap=0
                    is_gap=np.empty(0, dtype=np.int8) 
                    ns_gap=np.empty(0, dtype=np.int8) 
                    self.i_bin_start=i_bin+1
                else:
                    num_bin=num_bin+n_gap+1
                    is_gap=np.append(is_gap, i_bin)
                    ns_gap=np.append(ns_gap, n_gap)
            i_bin+=1

