#######################################
########## General parameter ##########
#######################################
#NAME_INEVT='ni1200120189_0mpu7_cl_bary_fpm.evt'
NAME_INEVT='ni1200120189_0mpu7_cl_bary.evt'

#--- XSELECT ---#
##### Session name (arbitrary) #####
NAME_SESSION='MAXI_J1820p070'
##### Binning of time in the extraction of light curve #####
BINSIZE=1.*10**(-3) #[sec]

##### PSD (Start) #####
##### Binnig of time in the calculation of power spectrum #####
NEWBIN=1.*10**(-3) #[sec]

##### Period in the calculation of power spectrum #####
NEWBINS_PER_INTERVAL=262144 #[-]

##### Number of period in the calculation of power spectrum #####
INTERVALS_PER_FRAME=80 #[-]

##### Binning of frequency in the calculation of power spectrum #####
REBIN=-1.2 #[-]

# ----- Acceptable rate of data gaps ----- #
R_G=0 #[-]
##### PSD (End) #####

##### Energy band (Start) #####
REF_MIN=201
REF_MAX=1000
COIS_MIN=[51,  151, 501]
COIS_MAX=[150, 500, 1000]
#####  Energy band (End)  #####

##### Lag-energy spectrum (Start) #####
F_MIN_LAGE=3.*10**(0) #[Hz]
F_MAX_LAGE=3.*10**(1) #[Hz]
##### Lag-energy spectrum (End) #####

##### RMS spectrum (Start) #####
F_MIN_RMS=3.*10**(0) #[Hz]
F_MAX_RMS=3.*10**(1) #[Hz]
##### RMS spectrum (End) #####

##### RMS spectrum (Start) #####
F_MIN_COV=3.*10**(0) #[Hz]
F_MAX_COV=3.*10**(1) #[Hz]
##### RMS spectrum (End) #####
