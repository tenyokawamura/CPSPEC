name_inevt='ni1200120130_0mpu7_cl_fpm.evt'
ch_ref_min=261
ch_ref_max=480
#chs_min   =[101]
#chs_max   =[260]
chs_min   =[101, ch_ref_min, 481, 701 ]
chs_max   =[260, ch_ref_max, 700, 1100]

# Extract light curve
session='MAXI_J1820+070'
binsize=1./128. # [s]

# FFT
n_bin=2**15
n_int=100
maximize=True
frac_gap=0.

# PSD, CSD
rebin=1.2
