import numpy as np
from numpy import random
import astropy.io.fits as fits

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
        hdus=fits.open(name_file)
        self.header=hdus[1].header
        self.data  =hdus[1].data
        self.telescope =self.header['TELESCOP']
        self.instrument=self.header['INSTRUME']
        self.source    =self.header['OBJECT']
        self.exposure  =self.header['EXPOSURE']
        self.dt        =self.header['TIMEDEL']
        self.n_row     =self.header['NAXIS2']

    def print_obs_info(self):
        print('------------------------------')
        print('(Observation information)')
        print('{0:<32}: {1}'      .format('Filename'         , self.name_file))
        print('{0:<32}: {1}'      .format('Telescope'        , self.telescope))
        print('{0:<32}: {1}'      .format('Instrument'       , self.instrument))
        print('{0:<32}: {1}'      .format('Source'           , self.source))
        print('{0:<32}: {1:.0f} s'.format('Exposure'         , self.exposure))
        print('{0:<32}: {1:.2g} s'.format('Sampling interval', self.dt))
        print('{0:<32}: {1}'      .format('Number of rows'   , self.n_row)) # for what?
        print('')

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
                return None, None, None
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
                    rates_int=data_int_real['RATE']

                    # --------------------------------------------------- #
                    # ----- Simulate and insert random data (begin) ----- #
                    # --------------------------------------------------- #
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

                    return t_start, t_end, rates_int

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

