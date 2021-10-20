import numpy as np
import cmath
import os
from constants import *

####################################
########## Common (Start) ##########
####################################
########## Display current program ##########
def print_filename(name_file):
    print('#####################################')
    print('Execute: {0}'.format(name_file))
    print('#####################################')

########## Check file extension ##########
def check_extension(name_file, extension):
    if not name_file.endswith(extension):
        print('Error: extension must be {0}. ({1}).'.format(extension, name_file))
        sys.exit()

########## Check file existence ##########
def check_existence(name_file):
    if not os.path.exists(name_file):
        print('Error: {0} does not exist.'.format(name_file))
        sys.exit()

def params_adjust_plot(f,\
                       x_mean,\
                       x_sigma):
    ### ax.plot ###
    f_plx=np.empty(0)
    for i_f in range(len(f)):
        if i_f==0:
            f_plx=np.append(f_plx, 0)
        else:
            f_plx=np.append(f_plx, (f[i_f]+f[i_f-1])/2.)
    f_plx[0]=f[0]-(f[1]-f[0])/2. #Assume Delta f=const.
    f_plx=np.append(f_plx, f[-1]+(f[-1]-f[-2])/2.)
    f_ply=np.append(f, f[-1])
    x_mean_pl=np.append(x_mean, x_mean[-1])
    ### ax.errorbar ###
    f_errx=f
    f_erry=f
    x_mean_err=x_mean
    x_sigma_err=x_sigma
    return f_plx, f_ply, x_mean_pl, f_errx, f_erry, x_mean_err, x_sigma_err

def mean_calc_1d(xs):
    if not xs.ndim==1:
        print('Error')
        sys.exit()
    n_x=len(xs)
    tot=0
    for x in xs:
        tot=tot+x
    x_mean=tot/n_x
    return x_mean

def var_calc_1d(xs):
    if not xs.ndim==1:
        print('Error')
        sys.exit()
    x_mean=mean_calc_1d(xs)
    n_x=len(xs)
    tot=0
    for x in xs:
        tot=tot+(x-x_mean)**2
    x_var=tot/n_x
    return x_var

def mea_var_calc_1d(xs):
    if not xs.ndim==1:
        print('Error')
        sys.exit()
    x_mea=mean_calc_1d(xs)
    n_x=len(xs)
    tot=0
    for x in xs:
        tot=tot+(x-x_mea)**2
    x_var=tot/n_x
    return x_mea, x_var

def mean_calc_2d(xs):
    if not xs.ndim==2:
        print('Error')
        sys.exit()
    n_row=len(xs)
    n_col=len(xs[0])
    x_mean=np.zeros(n_col)
    for i_col in range(n_col):
        x_mean[i_col]=mean_calc_1d(xs.T[i_col])
    return x_mean

def mean_var_calc_2d(xs):
    if not xs.ndim==2:
        print('Error')
        sys.exit()
    n_row=len(xs)
    n_col=len(xs[0])
    x_mean=np.zeros(n_col)
    x_var=np.zeros(n_col)
    for i_col in range(n_col):
        x_mean[i_col]=mean_calc_1d(xs.T[i_col])
        x_var[i_col]=var_calc_1d(xs.T[i_col])
    return x_mean, x_var

def ch2kev(ch):
    kev=0.01*ch
    return kev

def cr_sig_calc(cr_mea,\
                t):
    cr_sig=np.sqrt(cr_mea/t) #sqrt(cr_mea*t)/t, Assume Poisson noise
    return cr_sig

def int_tra_1d(xs,\
               ys):
    if not len(xs)==len(ys):
        print('Error')
        sys.exit()
    n_d=len(xs)
    tot=0
    for i in range(n_d-1):
        tra=(ys[i]+ys[i+1])*(xs[i+1]-xs[i])/2.
        tot=tot+tra
    return tot
##################################
########## Common (End) ##########
##################################

#######################################
########## lc_fft.py (Start) ##########
#######################################
def print_obs_info(telescope,\
                   instrument,\
                   source,\
                   exposure,\
                   resol_t,\
                   n_row):
    print('------------------------------')
    print('(Observation information)')
    print('Telescope: {0}'.format(telescope))
    print('Instrument: {0}'.format(instrument))
    print('Object: {0}'.format(source))
    print('Exposure time: {0:.3f} s'.format(exposure))
    print('Time resolution: {0:.3f} ms'.format(resol_t*ONE2MILI))
    print('Number of rows: {0}'.format(n_row))
    print('')

def print_ana_info(dt,\
                   n_d,\
                   lar_t,\
                   n_i_max,\
                   n_i_exp,\
                   f_min,\
                   f_max):
    print('(Analysis information)')
    print('Time resolution: {0:.3f} ms'.format(dt))
    print('Number of data in one interval: {0}'.format(n_d))
    print('Total time in one interval: {0:.3f} s'.format(lar_t))
    print('Maximum number of intervals: {0}'.format(n_i_max))
    print('Expected number of intervals: {0}'.format(n_i_exp))
    print('Minimum frequency: {0:.3g} mHz'.format(f_min*ONE2MILI))
    print('Maximum frequency: {0:.3g} Hz'.format(f_max))
    print('')

def write_info(name_inlc,\
               name_outtxt,\
               telescope,\
               instrument,\
               source,\
               exposure,\
               resol_t,\
               n_row,\
               dt,\
               n_d,\
               lar_t,\
               n_i_max,\
               n_i_exp,\
               f_min,\
               f_max):
    check_extension(name_inlc, EXTENSION_LC)
    check_extension(name_outtxt, EXTENSION_TXT)

    with open(name_outtxt, 'w') as fout:
        fout.write('--------------------------------------------------------------------------------\n')

        str_info='#Observation information\n'
        fout.write(str_info)
        str_info='{0:<32}: {1}\n'.format('Input file', name_inlc)
        fout.write(str_info)
        str_info='{0:<32}: {1}\n'.format('Output file', name_outtxt)
        fout.write(str_info)
        str_info='{0:<32}: {1}\n'.format('Telescope', telescope)
        fout.write(str_info)
        str_info='{0:<32}: {1}\n'.format('Instrument', instrument)
        fout.write(str_info)
        str_info='{0:<32}: {1}\n'.format('Source', source)
        fout.write(str_info)
        str_info='{0:<32}: {1}\n'.format('Exposure[s]', exposure)
        fout.write(str_info)
        str_info='{0:<32}: {1}\n'.format('Time resolution[s]', resol_t)
        fout.write(str_info)
        str_info='{0:<32}: {1}\n\n'.format('Number of rows', n_row)
        fout.write(str_info)

        str_info='#Analysis information\n'
        fout.write(str_info)
        str_info='{0:<32}: {1:.6g}\n'.format('Time resolution', dt)
        fout.write(str_info)
        str_info='{0:<32}: {1}\n'.format('Number of data in one interval', n_d)
        fout.write(str_info)
        str_info='{0:<32}: {1:.6g}\n'.format('Total time in one interval', lar_t)
        fout.write(str_info)
        str_info='{0:<32}: {1}\n'.format('Maximum number of intervals', n_i_max)
        fout.write(str_info)
        str_info='{0:<32}: {1}\n'.format('Expected number of intervals', n_i_exp)
        fout.write(str_info)
        str_info='{0:<32}: {1:.6g}\n'.format('Minimum frequency[Hz]', f_min)
        fout.write(str_info)
        str_info='{0:<32}: {1:.6g}\n'.format('Maximum frequency[Hz]', f_max)
        fout.write(str_info)
        
        fout.write('--------------------------------------------------------------------------------\n')

def fft_calc(data,\
             n_d,\
             dt,\
             lar_t):
    n_t=len(data)
    if not n_t==n_d:
        print('Error')
        sys.exit()
    xs=np.zeros(n_t)
    for i_t in range(n_t):
        xs[i_t]=data[i_t][1]
    del data
    x_mean=mean_calc_1d(xs)
    xs=xs-x_mean

    bs=np.fft.fft(a=xs,\
                  n=n_t,\
                  axis=0,\
                  norm=None)
    fs=np.fft.fftfreq(n=n_d,\
                      d=dt)
    n_f=n_t//2
    bs=bs[1:n_f]
    fs=fs[1:n_f]
    return fs, bs, x_mean

def fft_calc_v2(xs,\
                n_d,\
                dt):
    n_t=len(xs)
    if not n_t==n_d:
        print('Error')
        sys.exit()
    x_mea=mean_calc_1d(xs)
    xs=xs-x_mea

    bs=np.fft.fft(a=xs,\
                  n=n_t,\
                  axis=0,\
                  norm=None)
    fs=np.fft.fftfreq(n=n_d,\
                      d=dt)
    n_f=n_t//2
    bs=bs[1:n_f]
    fs=fs[1:n_f]
    return fs, bs, x_mea

def write_ft(name_outtxt,\
              i_i,\
              t_s,\
              t_e,\
              cr_mean,\
              fs,\
              bs):
    check_existence(name_outtxt)
    check_extension(name_outtxt, EXTENSION_TXT)
    with open(name_outtxt, 'a') as fout:
        str_header='#Number[-] Start[s] End[s] Mean flux[/s]\n'
        str_data=['{:.0f} ' .format(i_i+1),\
                  '{:.6f} ' .format(t_s),\
                  '{:.6f} ' .format(t_e),\
                  '{:.6f}\n'.format(cr_mean)]
        fout.write(str_header)
        fout.writelines(str_data)

        str_header='#Frequency[Hz] Re[FT[Light curve]][/s] Im[FT[Light curve]][/s] \n'
        fout.write(str_header)
        for f, b in zip(fs, bs):
            str_data=['{:.9g} ' .format(f),\
                      '{:.9g} ' .format(b.real),\
                      '{:.9g}\n'.format(b.imag)]
            fout.writelines(str_data)
#######################################
########## lc_fft.py (Start) ##########
#######################################

################################################
########## psd_calc_ns_ave.py (Start) ##########
################################################
def read_obs_info(name_intxt):
    check_extension(name_intxt, EXTENSION_TXT)
    with open(name_intxt, 'r') as fin:
        while True:
            line=fin.readline()
            if line.startswith('#Observation'):
                break
        while True:
            line=fin.readline()
            if line.startswith('Input'):
                line_str=line.split()
                name_infile=line_str[-1]
            elif line.startswith('Output'):
                line_str=line.split()
                name_outfile=line_str[-1]
            elif line.startswith('Telescope'):
                line_str=line.split()
                telescope=line_str[-1]
            elif line.startswith('Instrument'):
                line_str=line.split()
                instrument=line_str[-1]
            elif line.startswith('Source'):
                line_str=line.split()
                source=line_str[-1]
            elif line.startswith('Exposure'):
                line_str=line.split()
                exposure=float(line_str[-1])
            elif line.startswith('Time'):
                line_str=line.split()
                resol_t=float(line_str[-1])
            elif line.startswith('Number'):
                line_str=line.split()
                n_row=int(line_str[-1])
                break
    return telescope,\
           instrument,\
           source,\
           exposure,\
           resol_t,\
           n_row

def read_ana_info(name_intxt):
    check_extension(name_intxt, EXTENSION_TXT)
    with open(name_intxt, 'r') as fin:
        while True:
            line=fin.readline()
            if line.startswith('#Analysis'):
                break
        while True:
            line=fin.readline()
            if line.startswith('Time'):
                line_str=line.split()
                dt=float(line_str[-1])
            elif line.startswith('Number'):
                line_str=line.split()
                n_d=int(line_str[-1])
            elif line.startswith('Total'):
                line_str=line.split()
                lar_t=float(line_str[-1])
            elif line.startswith('Maximum number'):
                line_str=line.split()
                n_i_max=int(line_str[-1])
            elif line.startswith('Expected'):
                line_str=line.split()
                n_i_exp=int(line_str[-1])
            elif line.startswith('Minimum frequency'):
                line_str=line.split()
                f_min=float(line_str[-1])
            elif line.startswith('Maximum frequency'):
                line_str=line.split()
                f_max=float(line_str[-1])
                break
    return dt,\
           n_d,\
           lar_t,\
           n_i_max,\
           n_i_exp,\
           f_min,\
           f_max

def read_ft(name_intxt,\
            no_i):
    check_extension(name_intxt, EXTENSION_TXT)
    fs=np.empty(0)
    bs=np.empty(0, dtype=np.complex)
    with open(name_intxt, 'r') as fin:
        while True:
            line=fin.readline()
            if line.startswith('#Number[-]'):
                line=fin.readline()
                line_str=line.split()
                if int(line_str[0])==no_i:
                    t_s=float(line_str[1])
                    t_e=float(line_str[2])
                    cr_mean=float(line_str[3])
                    line=fin.readline()
                    break
        while True:
            line=fin.readline()
            if line.startswith('#Number[-]'):
                flag_end=False
                break
            elif not line:
                flag_end=True
                break
            else:
                line_str=line.split()
                f=float(line_str[0])
                b_real=float(line_str[1])
                b_imag=float(line_str[2])
                b=b_real+1j*b_imag
                fs=np.append(fs, f)
                bs=np.append(bs, b)
    return fs,\
           bs,\
           t_s,\
           t_e,\
           cr_mean,\
           flag_end

def psd_calc(bs,\
             x_mean,\
             n_d,\
             dt):
    a=2.*dt/(n_d*(x_mean**2))
    psd_pn=2./x_mean
    psds=a*np.abs(bs)**2 - psd_pn
    return psds

def psd_mean_sigma_calc(f, psd):
    f_mean=mean_calc_1d(f)
    psd_1d=np.ravel(psd)
    n_sample=len(psd_1d)
    psd_mean=mean_calc_1d(psd_1d)
    psd_sigma=np.sqrt(var_calc_1d(psd_1d)/n_sample)
    return f_mean, psd_mean, psd_sigma

def psd_unify(fs,\
              psds,\
              n_p,\
              rebin):
    if not len(psds)==n_p:
        print('Error')
        sys.exit()
    if rebin>-1:
        print('Error: Under construction')
        sys.exit()
    psds_t=psds.T
    alpha=abs(rebin)
    df=fs[1]-fs[0]
    f_s=fs[0]
    i_f_s=0
    f_fin=np.empty(0)
    psd_fin_mean=np.empty(0)
    psd_fin_sigma=np.empty(0)
    ns_mer_p=np.empty(0)
    ns_mer_f=np.empty(0)
    for i_f, f in enumerate(fs):
        #if i_f==0 or f-f_s>df: #more reasonable
        if i_f==0 or i_f==1 or f-f_s>df: #more similar to powspec
            i_f_e=i_f+1
            n_mer_p=n_p
            n_mer_f=i_f_e-i_f_s
            ns_mer_p=np.append(ns_mer_p, n_mer_p)
            ns_mer_f=np.append(ns_mer_f, n_mer_f)
            fs_mer=fs[i_f_s:i_f_e]
            psds_mer=psds_t[i_f_s:i_f_e]
            f_mer, psd_mer_mean, psd_mer_sigma=psd_mean_sigma_calc(fs_mer, psds_mer)
            f_fin=np.append(f_fin, f_mer)
            psd_fin_mean=np.append(psd_fin_mean, psd_mer_mean)
            psd_fin_sigma=np.append(psd_fin_sigma, psd_mer_sigma)

            i_f_s=i_f+1
            f_s=fs[i_f]
            df=df*alpha
    return f_fin,\
           psd_fin_mean,\
           psd_fin_sigma,\
           ns_mer_p,\
           ns_mer_f

def write_psd(name_outtxt,\
              fs,\
              psds_mean,\
              psds_sigma,\
              ns_mer_p,\
              ns_mer_f):
    check_extension(name_outtxt, EXTENSION_TXT)
    str_header='#Frequency[Hz] PSD_mean[(rms/mean)^2 Hz^-1] PSD_sigma[(rms/mean)^2 Hz^-1] \n'
    with open(name_outtxt, 'w') as fout:
        fout.write(str_header)
        for f, psd_mean, psd_sigma, n_mer_p, n_mer_f\
        in zip(fs, psds_mean, psds_sigma, ns_mer_p, ns_mer_f):
            str_data=['{:.9g} '.format(f),\
                      '{:.9g} '.format(psd_mean),\
                      '{:.9g} '.format(psd_sigma),\
                      '{:.0f} '.format(n_mer_p),\
                      '{:.0f}\n'.format(n_mer_f)]
            fout.writelines(str_data)
##############################################
########## psd_calc_ns_ave.py (End) ##########
##############################################

################################################
########## psd_calc_ave_ns.py (Start) ##########
################################################
def psd_raw_calc(bs):
    psds_raw=np.abs(bs)**2
    return psds_raw

def psd_mean_sigma_calc_norm(f,\
                             psd_raw,\
                             cr_mean,\
                             n_d,\
                             dt):
    f_mean=mean_calc_1d(f)
    psd_raw_1d=np.ravel(psd_raw)
    a=2.*dt/(n_d*(cr_mean**2))
    psd_raw_1d_norm=a*psd_raw_1d
    n_sample=len(psd_raw_1d_norm)
    psd_raw_mean=mean_calc_1d(psd_raw_1d_norm)
    psd_raw_sigma=np.sqrt(var_calc_1d(psd_raw_1d_norm)/n_sample)
    return f_mean, psd_raw_mean, psd_raw_sigma

def psd_unify_norm(fs,\
                   psds_raw,\
                   crs_mean,\
                   n_p,\
                   n_d,\
                   dt,\
                   rebin):
    if len(psds_raw)!=n_p or len(crs_mean)!=n_p:
        print('Error')
        sys.exit()
    if rebin>-1:
        print('Error: Under construction')
        sys.exit()

    cr_mean=mean_calc_1d(crs_mean)
    psd_pn=2./cr_mean

    psds_raw_t=psds_raw.T
    alpha=abs(rebin)
    df=fs[1]-fs[0]
    f_s=fs[0]
    i_f_s=0
    f_fin=np.empty(0)
    psd_raw_fin_mean=np.empty(0)
    psd_raw_fin_sigma=np.empty(0)
    ns_mer_p=np.empty(0)
    ns_mer_f=np.empty(0)
    for i_f, f in enumerate(fs):
        #if i_f==0 or f-f_s>df: #more reasonable
        if i_f==0 or i_f==1 or f-f_s>df: #more similar to powspec
            i_f_e=i_f+1
            n_mer_p=n_p
            n_mer_f=i_f_e-i_f_s
            ns_mer_p=np.append(ns_mer_p, n_mer_p)
            ns_mer_f=np.append(ns_mer_f, n_mer_f)
            fs_mer=fs[i_f_s:i_f_e]
            psds_raw_mer=psds_raw_t[i_f_s:i_f_e]
            f_mer,\
            psd_raw_mer_mean,\
            psd_raw_mer_sigma\
            =psd_mean_sigma_calc_norm(f=fs_mer,\
                                      psd_raw=psds_raw_mer,\
                                      cr_mean=cr_mean,\
                                      n_d=n_d,\
                                      dt=dt)
            f_fin=np.append(f_fin, f_mer)
            psd_raw_fin_mean=np.append(psd_raw_fin_mean, psd_raw_mer_mean)
            psd_raw_fin_sigma=np.append(psd_raw_fin_sigma, psd_raw_mer_sigma)

            i_f_s=i_f+1
            f_s=fs[i_f]
            df=df*alpha
    return f_fin,\
           psd_raw_fin_mean,\
           psd_raw_fin_sigma,\
           ns_mer_p,\
           ns_mer_f,\
           psd_pn,\
           cr_mean

def write_psd_ave_ns(name_outtxt,\
                     fs,\
                     psds_raw_mean,\
                     psds_raw_sigma,\
                     psd_noise,\
                     ns_mer_p,\
                     ns_mer_f,\
                     cr_mean,\
                     exposure):
    check_extension(name_outtxt, EXTENSION_TXT)
    with open(name_outtxt, 'w') as fout:
        str_head='#Mean flux[/s] Exposure[s]\n'
        str_data=['{:.9g} '.format(cr_mean),\
                  '{:.9g}\n'.format(exposure)]
        fout.write(str_head)
        fout.writelines(str_data)
        str_head='#Frequency[Hz] PSD_raw_mean[(rms/mean)^2 Hz^-1] PSD_raw_sigma[(rms/mean)^2 Hz^-1] PSD_noise[(rms/mean)^2 Hz^-1] N_p[-] N_f[-]\n'
        fout.write(str_head)
        for f, psd_raw_mean, psd_raw_sigma, n_mer_p, n_mer_f\
        in zip(fs, psds_raw_mean, psds_raw_sigma, ns_mer_p, ns_mer_f):
            str_data=['{:.9g} '.format(f),\
                      '{:.9g} '.format(psd_raw_mean),\
                      '{:.9g} '.format(psd_raw_sigma),\
                      '{:.9g} '.format(psd_noise),\
                      '{:.0f} '.format(n_mer_p),\
                      '{:.0f}\n'.format(n_mer_f)]
            fout.writelines(str_data)
##############################################
########## psd_calc_ave_ns.py (End) ##########
##############################################

#########################################
########## csd_calc.py (Start) ##########
#########################################
def csd_calc(bs_coi,\
             bs_ref,\
             x_mean_coi,\
             x_mean_ref,\
             n_d,\
             dt,\
             overlap):
    if overlap==True:
        nc=2./x_mean_ref
    else:
        nc=0
    a=2.*dt/(n_d*x_mean_coi*x_mean_ref)
    csds=a*(bs_coi.conjugate())*bs_ref-nc
    return csds

def csd_mean_sigma_calc(f, csd):
    f_mean=mean_calc_1d(f)
    csd_1d=np.ravel(csd)
    n_sample=len(csd_1d)
    csd_mean=mean_calc_1d(csd_1d)
    csd_sigma=np.sqrt(var_calc_1d(csd_1d)/n_sample)
    return f_mean, csd_mean, csd_sigma

def csd_unify(fs,\
              csds,\
              n_c,\
              rebin):
    if not len(csds)==n_c:
        print('Error')
        sys.exit()
    if rebin>-1:
        print('Error: Under construction')
        sys.exit()
    csds_t=csds.T
    alpha=abs(rebin)
    df=fs[1]-fs[0]
    f_s=fs[0]
    i_f_s=0
    f_fin=np.empty(0)
    csd_fin_mean=np.empty(0)
    csd_fin_sigma=np.empty(0)
    ns_mer_c=np.empty(0)
    ns_mer_f=np.empty(0)
    for i_f, f in enumerate(fs):
        #if i_f==0 or f-f_s>df: #more reasonable
        if i_f==0 or i_f==1 or f-f_s>df: #more similar to powspec
            i_f_e=i_f+1
            n_mer_c=n_c
            n_mer_f=i_f_e-i_f_s
            ns_mer_c=np.append(ns_mer_c, n_mer_c)
            ns_mer_f=np.append(ns_mer_f, n_mer_f)
            fs_mer=fs[i_f_s:i_f_e]
            csds_mer=csds_t[i_f_s:i_f_e]
            f_mer, csd_mer_mean, csd_mer_sigma=csd_mean_sigma_calc(fs_mer, csds_mer)
            f_fin=np.append(f_fin, f_mer)
            csd_fin_mean=np.append(csd_fin_mean, csd_mer_mean)
            csd_fin_sigma=np.append(csd_fin_sigma, csd_mer_sigma)

            i_f_s=i_f+1
            f_s=fs[i_f]
            df=df*alpha
    return f_fin,\
           csd_fin_mean,\
           csd_fin_sigma,\
           ns_mer_c,\
           ns_mer_f

def csd_mer_mea_calc(fs,\
                     csds):
    f_mea=mean_calc_1d(fs)
    csd_1d=np.ravel(csds)
    csd_mea=mean_calc_1d(csd_1d)
    return f_mea, csd_mea

def csd_e_sig_calc(csd_ab,\
                   psd_raw_ref,\
                   psd_noi_ref,\
                   psd_raw_coi,\
                   psd_noi_coi,\
                   n_i,\
                   n_f):
    b2, gamma2=bias_calc(psd_raw_1=psd_raw_coi,\
                         psd_noi_1=psd_noi_coi,\
                         psd_raw_2=psd_raw_ref,\
                         psd_noi_2=psd_noi_ref,\
                         n_mer_i=n_i,\
                         n_mer_f=n_f,\
                         csd_ab=csd_ab)
    csd_re_sig=np.sqrt(psd_raw_ref*(psd_raw_coi-((csd_ab**2-b2)/(psd_raw_ref-psd_noi_ref)))/(2.*n_i*n_f))
    csd_im_sig=csd_re_sig
    csd_ab_sig=csd_re_sig
    csd_ph_sig=np.sqrt(psd_raw_ref*( (psd_raw_coi/(csd_ab**2-b2)) - (1./(psd_raw_ref-psd_noi_ref)) ) / (2.*n_i*n_f) )

    return csd_re_sig,\
           csd_im_sig,\
           csd_ab_sig,\
           csd_ph_sig,\
           gamma2

def csd_unify_ing(fs,\
                  csds,\
                  n_c,\
                  rebin,\
                  fs_psd,\
                  psds_raw_ref,\
                  psds_noi_ref,\
                  psds_raw_coi,\
                  psds_noi_coi):
    ### Check (Start) ###
    if not len(csds)==n_c:
        print('Error')
        sys.exit()
    if rebin>-1:
        print('Error: Under construction')
        sys.exit()
    ### Check (End) ###

    ### Preparation (Start) ###
    csds_t=csds.T
    alpha=abs(rebin)
    df=fs[1]-fs[0]
    f_s=fs[0]
    i_f_s=0
    fs_mea=np.empty(0)
    csds_re_mea=np.empty(0)
    csds_im_mea=np.empty(0)
    csds_ab_mea=np.empty(0)
    csds_ph_mea=np.empty(0)
    ns_mer_c=np.empty(0)
    ns_mer_f=np.empty(0)
    ### Preparation (End) ###

    ### Merge (Start) ###
    ### Mean ###
    for i_f, f in enumerate(fs):
        #if i_f==0 or f-f_s>df: #more reasonable
        if i_f==0 or i_f==1 or f-f_s>df: #more similar to powspec
            i_f_e=i_f+1
            n_mer_c=n_c
            n_mer_f=i_f_e-i_f_s
            ns_mer_c=np.append(ns_mer_c, n_mer_c)
            ns_mer_f=np.append(ns_mer_f, n_mer_f)
            fs_mer=fs[i_f_s:i_f_e]
            csds_mer=csds_t[i_f_s:i_f_e]
            f_mea, csd_mea=csd_mer_mea_calc(fs=fs_mer,\
                                            csds=csds_mer)
            csd_re_mea=csd_mea.real
            csd_im_mea=csd_mea.imag
            csd_ab_mea=abs(csd_mea)
            csd_ph_mea=cmath.phase(csd_mea)

            fs_mea=np.append(fs_mea, f_mea)
            csds_re_mea=np.append(csds_re_mea, csd_re_mea)
            csds_im_mea=np.append(csds_im_mea, csd_im_mea)
            csds_ab_mea=np.append(csds_ab_mea, csd_ab_mea)
            csds_ph_mea=np.append(csds_ph_mea, csd_ph_mea)

            i_f_s=i_f+1
            f_s=fs[i_f]
            df=df*alpha

    ### Sigma ###
    n_d=len(fs_mea)
    if not n_d==len(fs_psd):
        print('Error')
        sys.exit()
    csds_re_sig=np.zeros(n_d)
    csds_im_sig=np.zeros(n_d)
    csds_ab_sig=np.zeros(n_d)
    csds_ph_sig=np.zeros(n_d)
    gamma2s=np.zeros(n_d)
    for i_d in range(n_d):
        csds_re_sig[i_d],\
        csds_im_sig[i_d],\
        csds_ab_sig[i_d],\
        csds_ph_sig[i_d],\
        gamma2s[i_d]\
        =csd_e_sig_calc(csd_ab=csds_ab_mea[i_d],\
                        psd_raw_ref=psds_raw_ref[i_d],\
                        psd_noi_ref=psds_noi_ref[i_d],\
                        psd_raw_coi=psds_raw_coi[i_d],\
                        psd_noi_coi=psds_noi_coi[i_d],\
                        n_i=ns_mer_c[i_d],\
                        n_f=ns_mer_f[i_d])
    ### Merge (End) ###

    return fs_mea,\
           csds_re_mea,\
           csds_im_mea,\
           csds_ab_mea,\
           csds_ph_mea,\
           csds_re_sig,\
           csds_im_sig,\
           csds_ab_sig,\
           csds_ph_sig,\
           ns_mer_c,\
           ns_mer_f

def csd_f_sig_calc(f,\
                   csd_re,\
                   csd_im,\
                   csd_ab,\
                   psd_raw_ref,\
                   psd_noi_ref,\
                   psd_raw_coi,\
                   psd_noi_coi,\
                   n_i,\
                   n_f):
    b2, gamma2=bias_calc(psd_raw_1=psd_raw_coi,\
                         psd_noi_1=psd_noi_coi,\
                         psd_raw_2=psd_raw_ref,\
                         psd_noi_2=psd_noi_ref,\
                         n_mer_i=n_i,\
                         n_mer_f=n_f,\
                         csd_ab=csd_ab)
    g2=raw_coh_calc(psd_raw_ref=psd_raw_ref,\
                    psd_raw_coi=psd_raw_coi,\
                    csd_ab=csd_ab,\
                    b2=b2)
    csd_re_sig=np.sqrt(((psd_raw_ref*psd_raw_coi) + (csd_re**2) - (csd_im**2))/(2.*n_i*n_f))
    csd_im_sig=np.sqrt(((psd_raw_ref*psd_raw_coi) - (csd_re**2) + (csd_im**2))/(2.*n_i*n_f))
    csd_ab_sig=np.sqrt(psd_raw_ref*psd_raw_coi/(n_i*n_f))
    csd_ph_sig=np.sqrt((1.-g2)/(2.*g2*n_i*n_f))
    csd_ti_sig=csd_ph_sig/(2.*np.pi*f)

    return csd_re_sig,\
           csd_im_sig,\
           csd_ab_sig,\
           csd_ph_sig,\
           csd_ti_sig,\
           gamma2

def csd_f_process(fs,\
                  csds,\
                  n_c,\
                  rebin,\
                  fs_psd,\
                  psds_raw_ref,\
                  psds_noi_ref,\
                  psds_raw_coi,\
                  psds_noi_coi):
    ### Check (Start) ###
    if not len(csds)==n_c:
        print('Error')
        sys.exit()
    if rebin>-1:
        print('Error: Under construction')
        sys.exit()
    ### Check (End) ###

    ### Preparation (Start) ###
    csds_t=csds.T
    alpha=abs(rebin)
    df=fs[1]-fs[0]
    f_s=fs[0]
    i_f_s=0
    fs_mea=np.empty(0)
    csds_re_mea=np.empty(0)
    csds_im_mea=np.empty(0)
    csds_ab_mea=np.empty(0)
    csds_ph_mea=np.empty(0)
    csds_ti_mea=np.empty(0)
    ns_mer_c=np.empty(0)
    ns_mer_f=np.empty(0)
    ### Preparation (End) ###

    ### Merge (Start) ###
    ### Mean ###
    for i_f, f in enumerate(fs):
        #if i_f==0 or f-f_s>df: #more reasonable
        if i_f==0 or i_f==1 or f-f_s>df: #more similar to powspec
            i_f_e=i_f+1
            n_mer_c=n_c
            n_mer_f=i_f_e-i_f_s
            ns_mer_c=np.append(ns_mer_c, n_mer_c)
            ns_mer_f=np.append(ns_mer_f, n_mer_f)
            fs_mer=fs[i_f_s:i_f_e]
            csds_mer=csds_t[i_f_s:i_f_e]
            f_mea, csd_mea=csd_mer_mea_calc(fs=fs_mer,\
                                            csds=csds_mer)
            csd_re_mea=csd_mea.real
            csd_im_mea=csd_mea.imag
            csd_ab_mea=abs(csd_mea)
            csd_ph_mea=cmath.phase(csd_mea)
            csd_ti_mea=csd_ph_mea/(2.*np.pi*f_mea)

            fs_mea=np.append(fs_mea, f_mea)
            csds_re_mea=np.append(csds_re_mea, csd_re_mea)
            csds_im_mea=np.append(csds_im_mea, csd_im_mea)
            csds_ab_mea=np.append(csds_ab_mea, csd_ab_mea)
            csds_ph_mea=np.append(csds_ph_mea, csd_ph_mea)
            csds_ti_mea=np.append(csds_ti_mea, csd_ti_mea)

            i_f_s=i_f+1
            f_s=fs[i_f]
            df=df*alpha

    ### Sigma ###
    n_d=len(fs_mea)
    if not n_d==len(fs_psd):
        print('Error')
        sys.exit()
    csds_re_sig=np.zeros(n_d)
    csds_im_sig=np.zeros(n_d)
    csds_ab_sig=np.zeros(n_d)
    csds_ph_sig=np.zeros(n_d)
    csds_ti_sig=np.zeros(n_d)
    gamma2s=np.zeros(n_d)
    for i_d in range(n_d):
        csds_re_sig[i_d],\
        csds_im_sig[i_d],\
        csds_ab_sig[i_d],\
        csds_ph_sig[i_d],\
        csds_ti_sig[i_d],\
        gamma2s[i_d]\
        =csd_f_sig_calc(f=fs_mea[i_d],\
                        csd_re=csds_re_mea[i_d],\
                        csd_im=csds_im_mea[i_d],\
                        csd_ab=csds_ab_mea[i_d],\
                        psd_raw_ref=psds_raw_ref[i_d],\
                        psd_noi_ref=psds_noi_ref[i_d],\
                        psd_raw_coi=psds_raw_coi[i_d],\
                        psd_noi_coi=psds_noi_coi[i_d],\
                        n_i=ns_mer_c[i_d],\
                        n_f=ns_mer_f[i_d])
    ### Merge (End) ###

    return fs_mea,\
           csds_re_mea,\
           csds_im_mea,\
           csds_ab_mea,\
           csds_ph_mea,\
           csds_ti_mea,\
           csds_re_sig,\
           csds_im_sig,\
           csds_ab_sig,\
           csds_ph_sig,\
           csds_ti_sig,\
           gamma2s,\
           ns_mer_c,\
           ns_mer_f

def write_csd(name_outtxt,\
              fs,\
              csds_mean,\
              csds_sigma,\
              ns_mer_c,\
              ns_mer_f):
    check_extension(name_outtxt, EXTENSION_TXT)
    str_header='#Frequency[Hz] Re[CSD_mean] Im[CSD_mean][(rms/mean)^2 Hz^-1] Re[CSD_sigma] Im[CSD_sigma][(rms/mean)^2 Hz^-1] N_c[-] N_f[-] \n'
    with open(name_outtxt, 'w') as fout:
        fout.write(str_header)
        for f, csd_mean, csd_sigma, n_mer_c, n_mer_f\
        in zip(fs, csds_mean, csds_sigma, ns_mer_c, ns_mer_f):
            csd_mean_re=csd_mean.real
            csd_mean_im=csd_mean.imag
            csd_sigma_re=csd_sigma.real
            csd_sigma_im=csd_sigma.imag
            str_data=['{:.9g} '.format(f),\
                      '{:.9g} '.format(csd_mean_re),\
                      '{:.9g} '.format(csd_mean_im),\
                      '{:.9g} '.format(csd_sigma_re),\
                      '{:.9g} '.format(csd_sigma_im),\
                      '{:.0f} '.format(n_mer_c),\
                      '{:.0f}\n'.format(n_mer_f)]
            fout.writelines(str_data)

def write_csd_ing(name_outtxt,\
                  fs,\
                  csds_re_mea,\
                  csds_im_mea,\
                  csds_ab_mea,\
                  csds_ph_mea,\
                  csds_re_sig,\
                  csds_im_sig,\
                  csds_ab_sig,\
                  csds_ph_sig,\
                  ns_mer_c,\
                  ns_mer_f):
    check_extension(name_outtxt, EXTENSION_TXT)
    str_header='#Frequency[Hz] Re[CSD_mean] Im[CSD_mean] Abs[CSD_mean][(rms/mean)^2 Hz^-1] Arg[CSD_mean][rad] Re[CSD_sigma] Im[CSD_sigma] Abs[CSD_sigma][(rms/mean)^2 Hz^-1] Arg[CSD_sigma][rad]  N_c[-] N_f[-] \n'
    with open(name_outtxt, 'w') as fout:
        fout.write(str_header)
        for f, csd_re_mea, csd_im_mea, csd_ab_mea, csd_ph_mea, csd_re_sig, csd_im_sig, csd_ab_sig, csd_ph_sig, n_mer_c, n_mer_f\
        in zip(fs, csds_re_mea, csds_im_mea, csds_ab_mea, csds_ph_mea, csds_re_sig, csds_im_sig, csds_ab_sig, csds_ph_sig, ns_mer_c, ns_mer_f):
            str_data=['{:.9g} '.format(f),\
                      '{:.9g} '.format(csd_re_mea),\
                      '{:.9g} '.format(csd_im_mea),\
                      '{:.9g} '.format(csd_ab_mea),\
                      '{:.9g} '.format(csd_ph_mea),\
                      '{:.9g} '.format(csd_re_sig),\
                      '{:.9g} '.format(csd_im_sig),\
                      '{:.9g} '.format(csd_ab_sig),\
                      '{:.9g} '.format(csd_ph_sig),\
                      '{:.0f} '.format(n_mer_c),\
                      '{:.0f}\n'.format(n_mer_f)]
            fout.writelines(str_data)

def write_csd_f(name_outtxt,\
                fs,\
                csds_re_mea,\
                csds_im_mea,\
                csds_ab_mea,\
                csds_ph_mea,\
                csds_ti_mea,\
                csds_re_sig,\
                csds_im_sig,\
                csds_ab_sig,\
                csds_ph_sig,\
                csds_ti_sig,\
                gamma2s,\
                ns_mer_c,\
                ns_mer_f):
    check_extension(name_outtxt, EXTENSION_TXT)
    str_header='#Frequency[Hz] Re[CSD_mean] Im[CSD_mean] Abs[CSD_mean][(rms/mean)^2 Hz^-1] Phase_mean[rad] Time_mean[s] Re[CSD_sigma] Im[CSD_sigma] Abs[CSD_sigma][(rms/mean)^2 Hz^-1] Phase_sigma[rad] Time_sigma[s] Intrinsic coherence(mean)[-] N_c[-] N_f[-] \n'
    with open(name_outtxt, 'w') as fout:
        fout.write(str_header)
        for f,\
            csd_re_mea,\
            csd_im_mea,\
            csd_ab_mea,\
            csd_ph_mea,\
            csd_ti_mea,\
            csd_re_sig,\
            csd_im_sig,\
            csd_ab_sig,\
            csd_ph_sig,\
            csd_ti_sig,\
            gamma2,\
            n_mer_c,\
            n_mer_f\
        in zip(fs,\
               csds_re_mea,\
               csds_im_mea,\
               csds_ab_mea,\
               csds_ph_mea,\
               csds_ti_mea,\
               csds_re_sig,\
               csds_im_sig,\
               csds_ab_sig,\
               csds_ph_sig,\
               csds_ti_sig,\
               gamma2s,\
               ns_mer_c,\
               ns_mer_f):
            str_data=['{:.9g} '.format(f),\
                      '{:.9g} '.format(csd_re_mea),\
                      '{:.9g} '.format(csd_im_mea),\
                      '{:.9g} '.format(csd_ab_mea),\
                      '{:.9g} '.format(csd_ph_mea),\
                      '{:.9g} '.format(csd_ti_mea),\
                      '{:.9g} '.format(csd_re_sig),\
                      '{:.9g} '.format(csd_im_sig),\
                      '{:.9g} '.format(csd_ab_sig),\
                      '{:.9g} '.format(csd_ph_sig),\
                      '{:.9g} '.format(csd_ti_sig),\
                      '{:.9g} '.format(gamma2),\
                      '{:.0f} '.format(n_mer_c),\
                      '{:.0f}\n'.format(n_mer_f)]
            fout.writelines(str_data)

def read_csd_f(name_intxt):
    check_existence(name_intxt)
    check_extension(name_intxt, EXTENSION_TXT)
    f,\
    csd_re_mea,\
    csd_im_mea,\
    csd_ab_mea,\
    csd_ph_mea,\
    csd_ti_mea,\
    csd_re_sig,\
    csd_im_sig,\
    csd_ab_sig,\
    csd_ph_sig,\
    csd_ti_sig,\
    gamma2,\
    n_mer_p,\
    n_mer_f\
    =np.loadtxt(fname=name_intxt, \
                dtype='float', \
                skiprows=1, \
                unpack=True)
    return f,\
           csd_re_mea,\
           csd_im_mea,\
           csd_ab_mea,\
           csd_ph_mea,\
           csd_ti_mea,\
           csd_re_sig,\
           csd_im_sig,\
           csd_ab_sig,\
           csd_ph_sig,\
           csd_ti_sig,\
           gamma2,\
           n_mer_p,\
           n_mer_f

def read_csd_ing(name_intxt):
    check_existence(name_intxt)
    check_extension(name_intxt, EXTENSION_TXT)
    f,\
    csd_re_mea,\
    csd_im_mea,\
    csd_ab_mea,\
    csd_ph_mea,\
    csd_re_sig,\
    csd_im_sig,\
    csd_ab_sig,\
    csd_ph_sig,\
    n_mer_p,\
    n_mer_f\
    =np.loadtxt(fname=name_intxt, \
                dtype='float', \
                skiprows=1, \
                unpack=True)
    return f,\
           csd_re_mea,\
           csd_im_mea,\
           csd_ab_mea,\
           csd_ph_mea,\
           csd_re_sig,\
           csd_im_sig,\
           csd_ab_sig,\
           csd_ph_sig,\
           n_mer_p,\
           n_mer_f
#######################################
########## csd_calc.py (End) ##########
#######################################

#########################################
########## coh_calc.py (Start) ##########
#########################################
def read_psd_ns_ave(name_intxt):
    check_existence(name_intxt)
    check_extension(name_intxt, EXTENSION_TXT)
    f,\
    psd_mean,\
    psd_sigma,\
    n_mer_p,\
    n_mer_f\
    =np.loadtxt(fname=name_intxt, \
                dtype='float', \
                skiprows=1, \
                unpack=True)
    return f,\
           psd_mean,\
           psd_sigma,\
           n_mer_p,\
           n_mer_f

def read_psd_ave_ns(name_intxt):
    check_existence(name_intxt)
    check_extension(name_intxt, EXTENSION_TXT)
    
    with open(name_intxt, 'r') as fin:
        line=fin.readline()
        line=fin.readline()
        line_str=line.split()
        cr_mean=float(line_str[0])
        exposure=float(line_str[1])

    f,\
    psd_mean,\
    psd_sigma,\
    psd_noise,\
    n_mer_p,\
    n_mer_f\
    =np.loadtxt(fname=name_intxt,\
                dtype='float',\
                skiprows=3,\
                unpack=True)
    return f,\
           psd_mean,\
           psd_sigma,\
           psd_noise,\
           n_mer_p,\
           n_mer_f,\
           cr_mean,\
           exposure

def read_csd(name_intxt):
    check_existence(name_intxt)
    check_extension(name_intxt, EXTENSION_TXT)
    f,\
    csd_mean_re,\
    csd_mean_im,\
    csd_sigma_re,\
    csd_sigma_im,\
    n_mer_c,\
    n_mer_f\
    =np.loadtxt(fname=name_intxt,\
                dtype='float',\
                skiprows=1,\
                unpack=True)
    csd_mean=csd_mean_re+1j*csd_mean_im
    csd_sigma=csd_sigma_re+1j*csd_sigma_im
    return f,\
           csd_mean,\
           csd_sigma,\
           n_mer_c,\
           n_mer_f

def bias_calc_coh(psd_raw_1,\
                  psd_noi_1,\
                  psd_raw_2,\
                  psd_noi_2,\
                  n_mer_i,\
                  n_mer_f,\
                  gamma2):
    b2=   ((psd_raw_1*psd_raw_2) - gamma2*(psd_raw_1-psd_noi_1)*(psd_raw_2-psd_noi_2))\
        / (n_mer_i*n_mer_f)
    return b2

def bias_calc(psd_raw_1,\
              psd_noi_1,\
              psd_raw_2,\
              psd_noi_2,\
              n_mer_i,\
              n_mer_f,\
              csd_ab):
    gamma2=1 #intrinsic coherence
    ite=0
    while True:
        b2=bias_calc_coh(psd_raw_1=psd_raw_1,\
                         psd_noi_1=psd_noi_1,\
                         psd_raw_2=psd_raw_2,\
                         psd_noi_2=psd_noi_2,\
                         n_mer_i=n_mer_i,\
                         n_mer_f=n_mer_f,\
                         gamma2=gamma2)
        gamma2_new=intri_coh_calc(psd_raw_1=psd_raw_1,\
                                  psd_noi_1=psd_noi_1,\
                                  psd_raw_2=psd_raw_2,\
                                  psd_noi_2=psd_noi_2,\
                                  csd_ab=csd_ab,\
                                  b2=b2)
        ite=ite+1
        if abs(gamma2_new-gamma2)<10**(-3):
            #print('Iteration: {0}, b2: {1}, gamma2: {2}'.format(ite, b2, gamma2))
            gamma2=gamma2_new
            break
        else:
            gamma2=gamma2_new
        if ite>=20:
            gamma2=gamma2_new
            print('Convergence is not achieved.(b2: {0}, gamma2: {1})'.format(b2, gamma2))
    if n_mer_i*n_mer_f>500 or b2>csd_ab:
        b2=0
    return b2, gamma2

def raw_coh_calc(psd_raw_ref,\
                 psd_raw_coi,\
                 csd_ab,\
                 b2):
    raw_coh=(csd_ab**2-b2)/(psd_raw_ref*psd_raw_coi)
    return raw_coh

def intri_coh_calc(psd_raw_1,\
                   psd_noi_1,\
                   psd_raw_2,\
                   psd_noi_2,\
                   csd_ab,\
                   b2):
    gamma2=(csd_ab**2-b2)/((psd_raw_1-psd_noi_1)*(psd_raw_2-psd_noi_2))
    return gamma2

def coh_calc(psd_raw_1,\
             psd_noi_1,\
             psd_raw_2,\
             psd_noi_2,\
             csd,\
             n_mer_i,\
             n_mer_f):
    csd_ab=abs(csd)
    b2, gamma2=bias_calc(psd_raw_1=psd_raw_1,\
                         psd_noi_1=psd_noi_1,\
                         psd_raw_2=psd_raw_2,\
                         psd_noi_2=psd_noi_2,\
                         n_mer_i=n_mer_i,\
                         n_mer_f=n_mer_f,\
                         csd_ab=csd_ab)
    coh=(csd_ab**2-b2)/(psd_raw_1*psd_raw_2)
    return coh

def write_coh(name_outtxt,\
              fs,\
              cohs,\
              ns_mer_i,\
              ns_mer_f):
    check_extension(name_outtxt, EXTENSION_TXT)
    str_header='#Frequency[Hz] Coherence[-] N_i[-] N_f[-]\n'
    with open(name_outtxt, 'w') as fout:
        fout.write(str_header)
        for f, coh, n_mer_i, n_mer_f\
        in zip(fs, cohs, ns_mer_i, ns_mer_f):
            str_data=['{:.9g} '.format(f),\
                      '{:.9g} '.format(coh),\
                      '{:.0f} '.format(n_mer_i),\
                      '{:.0f}\n'.format(n_mer_f)]
            fout.writelines(str_data)
#######################################
########## coh_calc.py (End) ##########
#######################################

##########################################
########## lagf_calc.py (Start) ##########
##########################################
def read_coh(name_intxt):
    check_existence(name_intxt)
    check_extension(name_intxt, EXTENSION_TXT)
    f,\
    coh_mean,\
    n_mer_i,\
    n_mer_f\
    =np.loadtxt(fname=name_intxt,\
                dtype='float',\
                skiprows=1,\
                unpack=True)
    return f,\
           coh_mean,\
           n_mer_i,\
           n_mer_f

def lagf_calc(f,\
              csd_mean,\
              coh_mean,\
              n_mer_i,\
              n_mer_f):
    n_b=len(f)
    pha_mean=np.zeros(n_b)
    for i_b in range(n_b):
        pha_mean[i_b]=cmath.phase(csd_mean[i_b])

    pha_sigma=np.sqrt((1.-coh_mean)/(2.*coh_mean*n_mer_i*n_mer_f))
    tau_mean=pha_mean/(2.*np.pi*f)
    tau_sigma=pha_sigma/(2.*np.pi*f)
    return pha_mean,\
           pha_sigma,\
           tau_mean,\
           tau_sigma

def write_lagf(name_outtxt,\
               fs,\
               phas_mean,\
               phas_sigma,\
               taus_mean,\
               taus_sigma,\
               ns_mer_i,\
               ns_mer_f):
    check_extension(name_outtxt, EXTENSION_TXT)
    str_header='#Frequency[Hz] Phase_mean[rad] Phase_sigma[rad] Time_mean[s] Time_sigma[s] N_i[-] N_f[-]\n'
    with open(name_outtxt, 'w') as fout:
        fout.write(str_header)
        for f, pha_mean, pha_sigma, tau_mean, tau_sigma, n_mer_i, n_mer_f\
        in zip(fs, phas_mean, phas_sigma, taus_mean, taus_sigma, ns_mer_i, ns_mer_f):
            str_data=['{:.9g} '.format(f),\
                      '{:.9g} '.format(pha_mean),\
                      '{:.9g} '.format(pha_sigma),\
                      '{:.9g} '.format(tau_mean),\
                      '{:.9g} '.format(tau_sigma),\
                      '{:.0f} '.format(n_mer_i),\
                      '{:.0f}\n'.format(n_mer_f)]
            fout.writelines(str_data)
########################################
########## lagf_calc.py (End) ##########
########################################

##########################################
########## lage_calc.py (Start) ##########
##########################################
def pars_cut(x,\
             y,\
             x_min,\
             x_max):
    if not len(x)==len(y):
        print('Error')
        sys.exit()
    x_cut=x[(x_min<x) & (x<x_max)]
    y_cut=y[(x_min<x) & (x<x_max)]
    return x_cut, y_cut

def cut_psd_csd(f_ref,\
                psd_raw_ref,\
                psd_noi_ref,\
                n_mer_i_ref,\
                n_mer_f_ref,\
                f_coi,\
                psd_raw_coi,\
                psd_noi_coi,\
                n_mer_i_coi,\
                n_mer_f_coi,\
                f_csd,\
                csd_raw,\
                n_mer_i_csd,\
                n_mer_f_csd,\
                f_min,\
                f_max):
    if (len(f_ref)!=len(f_coi)) or (len(f_ref)!=len(f_csd)):
        print('Error')
        sys.exit()

    ### Reference band ###
    f_cut, n_mer_i_cut=pars_cut(x=f_ref,\
                                y=n_mer_i_ref,\
                                x_min=f_min,\
                                x_max=f_max)
    f_cut, n_mer_f_cut=pars_cut(x=f_ref,\
                                y=n_mer_f_ref,\
                                x_min=f_min,\
                                x_max=f_max)
    f_cut, psd_raw_ref_cut=pars_cut(x=f_ref,\
                                    y=psd_raw_ref,\
                                    x_min=f_min,\
                                    x_max=f_max)
    f_cut, psd_noi_ref_cut=pars_cut(x=f_ref,\
                                    y=psd_noi_ref,\
                                    x_min=f_min,\
                                    x_max=f_max)

    ### Channel-of-Interest band ###
    f_cut, psd_raw_coi_cut=pars_cut(x=f_coi,\
                                    y=psd_raw_coi,\
                                    x_min=f_min,\
                                    x_max=f_max)
    f_cut, psd_noi_coi_cut=pars_cut(x=f_coi,\
                                    y=psd_noi_coi,\
                                    x_min=f_min,\
                                    x_max=f_max)

    ### Cross spectrum ###
    f_cut, csd_raw_cut=pars_cut(x=f_csd,\
                                y=csd_raw,\
                                x_min=f_min,\
                                x_max=f_max)

    return n_mer_i_cut,\
           n_mer_f_cut,\
           f_cut,\
           psd_raw_ref_cut,\
           psd_noi_ref_cut,\
           psd_raw_coi_cut,\
           psd_noi_coi_cut,\
           csd_raw_cut

def mean_calc_1d_weigh(xs, ws):
    if xs.ndim!=1 or ws.ndim!=1:
        print('Error')
        sys.exit()
    if not len(xs)==len(ws):
        print('Error')
        sys.exit()
    n=sum(ws)
    tot=0
    for x, w in zip(xs, ws):
        tot=tot+x*w
    mea=tot/n
    return mea

def ave_psd_csd(psd_raw_ref,\
                psd_noi_ref,\
                psd_raw_coi,\
                psd_noi_coi,\
                csd_raw,\
                n_mer_i,\
                n_mer_f):
    psd_raw_ref_mea=mean_calc_1d_weigh(xs=psd_raw_ref, ws=n_mer_f)
    psd_noi_ref_mea=mean_calc_1d_weigh(xs=psd_noi_ref, ws=n_mer_f)
    psd_raw_coi_mea=mean_calc_1d_weigh(xs=psd_raw_coi, ws=n_mer_f)
    psd_noi_coi_mea=mean_calc_1d_weigh(xs=psd_noi_coi, ws=n_mer_f)
    csd_raw_mea=mean_calc_1d_weigh(xs=csd_raw, ws=n_mer_f)
    n_mer_i_fin=n_mer_i[0]
    n_mer_f_fin=sum(n_mer_f)

    psd_raw_ref_mea=np.array([psd_raw_ref_mea])
    psd_noi_ref_mea=np.array([psd_noi_ref_mea])
    psd_raw_coi_mea=np.array([psd_raw_coi_mea])
    psd_noi_coi_mea=np.array([psd_noi_coi_mea])
    csd_raw_mea=np.array([csd_raw_mea])
    n_mer_i=np.array([n_mer_i])
    n_mer_f=np.array([n_mer_f])

    return psd_raw_ref_mea,\
           psd_noi_ref_mea,\
           psd_raw_coi_mea,\
           psd_noi_coi_mea,\
           csd_raw_mea,\
           n_mer_i_fin,\
           n_mer_f_fin

def csd_e_process(fs,\
                  psds_raw_ref,\
                  psds_noi_ref,\
                  psds_raw_coi,\
                  psds_noi_coi,\
                  csds,\
                  ns_mer_i,\
                  ns_mer_f):
    ### Mean (Start) ###
    f_mea=mean_calc_1d_weigh(xs=fs, ws=ns_mer_f)
    psd_raw_ref_mea=mean_calc_1d_weigh(xs=psds_raw_ref, ws=ns_mer_f)
    psd_noi_ref_mea=mean_calc_1d_weigh(xs=psds_noi_ref, ws=ns_mer_f)
    psd_raw_coi_mea=mean_calc_1d_weigh(xs=psds_raw_coi, ws=ns_mer_f)
    psd_noi_coi_mea=mean_calc_1d_weigh(xs=psds_noi_coi, ws=ns_mer_f)
    csd_mea=mean_calc_1d_weigh(xs=csds, ws=ns_mer_f)
    n_mer_i=ns_mer_i[0]
    n_mer_f=sum(ns_mer_f)

    csd_re_mea=csd_mea.real
    csd_im_mea=csd_mea.imag
    csd_ab_mea=abs(csd_mea)
    csd_ph_mea=cmath.phase(csd_mea)
    csd_ti_mea=csd_ph_mea/(2.*np.pi*f_mea)
    ### Mean (End) ###

    ### Sigma (Start) ###
    csd_re_sig,\
    csd_im_sig,\
    csd_ab_sig,\
    csd_ph_sig,\
    gamma2\
    =csd_e_sig_calc(csd_ab=csd_ab_mea,\
                    psd_raw_ref=psd_raw_ref_mea,\
                    psd_noi_ref=psd_noi_ref_mea,\
                    psd_raw_coi=psd_raw_coi_mea,\
                    psd_noi_coi=psd_noi_coi_mea,\
                    n_i=n_mer_i,\
                    n_f=n_mer_f)
    csd_ti_sig=csd_ph_sig/(2.*np.pi*f_mea)
    ### Sigma (End) ###

    return f_mea,\
           csd_re_mea,\
           csd_im_mea,\
           csd_ab_mea,\
           csd_ph_mea,\
           csd_ti_mea,\
           csd_re_sig,\
           csd_im_sig,\
           csd_ab_sig,\
           csd_ph_sig,\
           csd_ti_sig,\
           gamma2,\
           n_mer_i,\
           n_mer_f

def write_lage(name_outtxt,\
               ref_min,\
               ref_max,\
               cois_min,\
               cois_max,\
               f_min,\
               f_max,\
               phas_mea,\
               phas_sig,\
               taus_mea,\
               taus_sig,\
               cohs_mea):
    check_extension(name_outtxt, EXTENSION_TXT)
    with open(name_outtxt, 'w') as fout:
        str_header='#Ref_min[-] Ref_max[-] F_min[Hz] F_max[rad]\n'
        str_data=['{:.0f} '.format(ref_min),\
                  '{:.0f} '.format(ref_max),\
                  '{:.9g} '.format(f_min),\
                  '{:.9g}\n'.format(f_max)]
        fout.write(str_header)
        fout.writelines(str_data)
        str_header='#CH_min[-] CH_max[-] Frequency[Hz] Phase_mean[rad] Phase_sigma[rad] Time_mean[s] Time_sigma[s] Coherence[-]\n'
        fout.write(str_header)
        for coi_min, coi_max, pha_mea, pha_sig, tau_mea, tau_sig, coh_mea\
        in zip(cois_min, cois_max, phas_mea, phas_sig, taus_mea, taus_sig, cohs_mea):
            str_data=['{:.0f} '.format(coi_min),\
                      '{:.0f} '.format(coi_max),\
                      '{:.9g} '.format(pha_mea),\
                      '{:.9g} '.format(pha_sig),\
                      '{:.9g} '.format(tau_mea),\
                      '{:.9g} '.format(tau_sig),\
                      '{:.9g}\n'.format(coh_mea)]
            fout.writelines(str_data)

def write_csd_e(name_outtxt,\
                ref_min,\
                ref_max,\
                cois_min,\
                cois_max,\
                f_min,\
                f_max,\
                csds_re_mea,\
                csds_im_mea,\
                csds_ab_mea,\
                csds_ph_mea,\
                csds_ti_mea,\
                csds_re_sig,\
                csds_im_sig,\
                csds_ab_sig,\
                csds_ph_sig,\
                csds_ti_sig,\
                gamma2s,\
                ns_mer_i,\
                ns_mer_f):
    check_extension(name_outtxt, EXTENSION_TXT)
    with open(name_outtxt, 'w') as fout:
        str_header='#Ref_min[-] Ref_max[-] F_min[Hz] F_max[rad]\n'
        str_data=['{:.0f} '.format(ref_min),\
                  '{:.0f} '.format(ref_max),\
                  '{:.9g} '.format(f_min),\
                  '{:.9g}\n'.format(f_max)]
        fout.write(str_header)
        fout.writelines(str_data)
        str_header='#CH_min[-] CH_max[-] Re[CSD_mean] Im[CSD_mean] Abs[CSD_mean][(rms/mean)^2 Hz^-1] Phase_mean[rad] Time_mean[s] Re[CSD_sigma] Im[CSD_sigma] Abs[CSD_sigma][(rms/mean)^2 Hz^-1] Phase_sigma[rad] Time_sigma[s] Intrinsic coherence(mean)[-] N_i[-] N_f[-] \n'
        fout.write(str_header)
        for coi_min,\
            coi_max,\
            csd_re_mea,\
            csd_im_mea,\
            csd_ab_mea,\
            csd_ph_mea,\
            csd_ti_mea,\
            csd_re_sig,\
            csd_im_sig,\
            csd_ab_sig,\
            csd_ph_sig,\
            csd_ti_sig,\
            gamma2,\
            n_mer_i,\
            n_mer_f\
        in zip(cois_min,\
               cois_max,\
               csds_re_mea,\
               csds_im_mea,\
               csds_ab_mea,\
               csds_ph_mea,\
               csds_ti_mea,\
               csds_re_sig,\
               csds_im_sig,\
               csds_ab_sig,\
               csds_ph_sig,\
               csds_ti_sig,\
               gamma2s,\
               ns_mer_i,\
               ns_mer_f):
            str_data=['{:.0f} '.format(coi_min),\
                      '{:.0f} '.format(coi_max),\
                      '{:.9g} '.format(csd_re_mea),\
                      '{:.9g} '.format(csd_im_mea),\
                      '{:.9g} '.format(csd_ab_mea),\
                      '{:.9g} '.format(csd_ph_mea),\
                      '{:.9g} '.format(csd_ti_mea),\
                      '{:.9g} '.format(csd_re_sig),\
                      '{:.9g} '.format(csd_im_sig),\
                      '{:.9g} '.format(csd_ab_sig),\
                      '{:.9g} '.format(csd_ph_sig),\
                      '{:.9g} '.format(csd_ti_sig),\
                      '{:.9g} '.format(gamma2),\
                      '{:.0f} '.format(n_mer_i),\
                      '{:.0f}\n'.format(n_mer_f)]
            fout.writelines(str_data)

def read_csd_e(name_intxt):
    check_existence(name_intxt)
    check_extension(name_intxt, EXTENSION_TXT)

    with open(name_intxt, 'r') as fin:
        line=fin.readline()
        line=fin.readline()
        line_str=line.split()
        ref_min=int(line_str[0])
        ref_max=int(line_str[1])
        f_min=float(line_str[2])
        f_max=float(line_str[3])

    cois_min,\
    cois_max,\
    csds_re_mea,\
    csds_im_mea,\
    csds_ab_mea,\
    csds_ph_mea,\
    csds_ti_mea,\
    csds_re_sig,\
    csds_im_sig,\
    csds_ab_sig,\
    csds_ph_sig,\
    csds_ti_sig,\
    gamma2s,\
    ns_mer_i,\
    ns_mer_f\
    =np.loadtxt(fname=name_intxt, \
                dtype='float', \
                skiprows=3, \
                unpack=True)
    return ref_min,\
           ref_max,\
           f_min,\
           f_max,\
           cois_min,\
           cois_max,\
           csds_re_mea,\
           csds_im_mea,\
           csds_ab_mea,\
           csds_ph_mea,\
           csds_ti_mea,\
           csds_re_sig,\
           csds_im_sig,\
           csds_ab_sig,\
           csds_ph_sig,\
           csds_ti_sig,\
           gamma2s,\
           ns_mer_i,\
           ns_mer_f
########################################
########## lage_calc.py (End) ##########
########################################

#########################################
########## rms_calc.py (Start) ##########
#########################################
def cut_psd(f,\
            psd_raw,\
            psd_sig,\
            psd_noi,\
            n_mer_i,\
            n_mer_f,\
            f_min,\
            f_max):
    f_cut, n_mer_i_cut=pars_cut(x=f,\
                                y=n_mer_i,\
                                x_min=f_min,\
                                x_max=f_max)
    f_cut, n_mer_f_cut=pars_cut(x=f,\
                                y=n_mer_f,\
                                x_min=f_min,\
                                x_max=f_max)
    f_cut, psd_raw_cut=pars_cut(x=f,\
                                y=psd_raw,\
                                x_min=f_min,\
                                x_max=f_max)
    f_cut, psd_sig_cut=pars_cut(x=f,\
                                y=psd_sig,\
                                x_min=f_min,\
                                x_max=f_max)
    f_cut, psd_noi_cut=pars_cut(x=f,\
                                y=psd_noi,\
                                x_min=f_min,\
                                x_max=f_max)
    return n_mer_i_cut,\
           n_mer_f_cut,\
           f_cut,\
           psd_raw_cut,\
           psd_sig_cut,\
           psd_noi_cut

def rms_calc(fs,\
             psds_raw,\
             psds_sig,\
             psds_noi,\
             cr_mea,\
             ns_mer_i,\
             ns_mer_f):
    ### Mean ###
    psds=psds_raw-psds_noi
    var_per_mea=int_tra_1d(xs=fs,\
                           ys=psds)
    var_abs_mea=var_per_mea*(cr_mea**2)
    rms_per_mea=np.sqrt(var_per_mea)
    rms_abs_mea=np.sqrt(var_abs_mea)

    ### Sigma ### Assume coherence~1
    n_mer_i=ns_mer_i[0]
    n_mer_f=sum(ns_mer_f)
    var_per_noi=int_tra_1d(xs=fs,\
                           ys=psds_noi)
    var_abs_noi=var_per_noi*(cr_mea**2)
    rms_per_sig=np.sqrt(  (2.*var_per_mea*var_per_noi+(var_per_noi**2))\
                        / (4.*n_mer_i*n_mer_f*var_per_mea))
    rms_abs_sig=np.sqrt(  (2.*var_abs_mea*var_abs_noi+(var_abs_noi**2))\
                        / (4.*n_mer_i*n_mer_f*var_abs_mea))
    return rms_per_mea,\
           rms_per_sig,\
           rms_abs_mea,\
           rms_abs_sig\

def rms_calc_strict(fs,\
                    psds_raw_ref,\
                    psds_noi_ref,\
                    psds_raw_coi,\
                    psds_noi_coi,\
                    csds,\
                    ns_mer_i,\
                    ns_mer_f,\
                    cr_mea_coi):
    ############
    ### Mean ###
    ############
    psds=psds_raw_coi-psds_noi_coi
    var_per_mea=int_tra_1d(xs=fs,\
                           ys=psds)
    var_abs_mea=var_per_mea*(cr_mea_coi**2)
    rms_per_mea=np.sqrt(var_per_mea)
    rms_abs_mea=np.sqrt(var_abs_mea)

    #############
    ### Sigma ###
    #############
    ### Mean (Start) ###
    f_mea=mean_calc_1d_weigh(xs=fs, ws=ns_mer_f)
    psd_raw_ref_mea=mean_calc_1d_weigh(xs=psds_raw_ref, ws=ns_mer_f)
    psd_noi_ref_mea=mean_calc_1d_weigh(xs=psds_noi_ref, ws=ns_mer_f)
    psd_raw_coi_mea=mean_calc_1d_weigh(xs=psds_raw_coi, ws=ns_mer_f)
    psd_noi_coi_mea=mean_calc_1d_weigh(xs=psds_noi_coi, ws=ns_mer_f)
    csd_mea=mean_calc_1d_weigh(xs=csds, ws=ns_mer_f)
    n_mer_i=ns_mer_i[0]
    n_mer_f=sum(ns_mer_f)

    csd_ab_mea=abs(csd_mea)
    ### Mean (End) ###

    ### Intrinsic coherence (Start) ###
    b2, gamma2=bias_calc(psd_raw_1=psd_raw_coi_mea,\
                         psd_noi_1=psd_noi_coi_mea,\
                         psd_raw_2=psd_raw_ref_mea,\
                         psd_noi_2=psd_noi_ref_mea,\
                         n_mer_i=n_mer_i,\
                         n_mer_f=n_mer_f,\
                         csd_ab=csd_ab_mea)
    ### Intrinsic coherence (End) ###

    var_per_noi=int_tra_1d(xs=fs,\
                           ys=psds_noi_coi)
    var_abs_noi=var_per_noi*(cr_mea_coi**2)
    rms_per_sig=np.sqrt( ( (1.-(gamma2**2))*(var_per_mea**2) + \
                           (var_per_noi)**2 + \
                           (2.*var_per_mea*var_per_noi) \
                         ) / \
                         (4.*n_mer_i*n_mer_f*var_per_mea) \
                       )
    rms_abs_sig=np.sqrt( ( (1.-(gamma2**2))*(var_abs_mea**2) + \
                           (var_abs_noi)**2 + \
                           (2.*var_abs_mea*var_abs_noi) \
                         ) / \
                         (4.*n_mer_i*n_mer_f*var_abs_mea) \
                       )
    return rms_per_mea,\
           rms_per_sig,\
           rms_abs_mea,\
           rms_abs_sig

def write_rms(name_outtxt,\
              cois_min,\
              cois_max,\
              f_min,\
              f_max,\
              rmses_per_mea,\
              rmses_per_sig,\
              rmses_abs_mea,\
              rmses_abs_sig,\
              crs_mea,\
              crs_sig,\
              exposure):
    check_extension(name_outtxt, EXTENSION_TXT)
    with open(name_outtxt, 'w') as fout:
        str_head='#F_min[Hz] F_max[rad] Exposure[s]\n'
        str_data=['{:.9g} '.format(f_min),\
                  '{:.9g} '.format(f_max),\
                  '{:.9g}\n'.format(exposure)]
        fout.write(str_head)
        fout.writelines(str_data)
        str_head='#CH_min[-] CH_max[-] RMS_mean[100%] RMS_sigma[100%] RMS_mean[/s] RMS_sigma[/s] Flux_mean[/s] Flux_sigma[/s]\n'
        fout.write(str_head)
        for coi_min, coi_max, rms_per_mea, rms_per_sig, rms_abs_mea, rms_abs_sig, cr_mea, cr_sig\
        in zip(cois_min, cois_max, rmses_per_mea, rmses_per_sig, rmses_abs_mea, rmses_abs_sig, crs_mea, crs_sig):
            str_data=['{:.0f} '.format(coi_min),\
                      '{:.0f} '.format(coi_max),\
                      '{:.9g} '.format(rms_per_mea),\
                      '{:.9g} '.format(rms_per_sig),\
                      '{:.9g} '.format(rms_abs_mea),\
                      '{:.9g} '.format(rms_abs_sig),\
                      '{:.9g} '.format(cr_mea),\
                      '{:.9g}\n'.format(cr_sig)]
            fout.writelines(str_data)
#######################################
########## rms_calc.py (End) ##########
#######################################

#########################################
########## cov_calc.py (Start) ##########
#########################################
def cov_calc(fs,\
             psds_raw_ref,\
             psds_noi_ref,\
             cr_mea_ref,\
             psds_raw_coi,\
             psds_noi_coi,\
             cr_mea_coi,\
             csds_raw,\
             ns_mer_i,\
             ns_mer_f):
    ### Mean ###
    csds_ab=abs(csds_raw)
    b2=np.zeros(len(fs))
    for i_b,\
        (psd_raw_ref,\
         psd_noi_ref,\
         psd_raw_coi,\
         psd_noi_coi,\
         n_mer_i,\
         n_mer_f,\
         csd_ab)\
    in enumerate(zip(psds_raw_ref,\
                     psds_noi_ref,\
                     psds_raw_coi,\
                     psds_noi_coi,\
                     ns_mer_i,\
                     ns_mer_f,\
                     csds_ab)):
        b2[i_b], gamma2=bias_calc(psd_raw_1=psd_raw_coi,\
                                  psd_noi_1=psd_noi_coi,\
                                  psd_raw_2=psd_raw_ref,\
                                  psd_noi_2=psd_noi_ref,\
                                  n_mer_i=n_mer_i,\
                                  n_mer_f=n_mer_f,\
                                  csd_ab=csd_ab)

    csds=(csds_ab**2-b2)/(psds_raw_ref-psds_noi_ref)
    var_per_mea=int_tra_1d(xs=fs,\
                           ys=csds)
    var_abs_mea=var_per_mea*(cr_mea_coi**2)
    cov_per_mea=np.sqrt(var_per_mea)
    cov_abs_mea=np.sqrt(var_abs_mea)

    ### Sigma ### Assume coherence~1
    n_mer_i=ns_mer_i[0]
    n_mer_f=sum(ns_mer_f)
    var_per_raw_ref=int_tra_1d(xs=fs,\
                               ys=psds_raw_ref)
    var_per_noi_ref=int_tra_1d(xs=fs,\
                               ys=psds_noi_ref)
    var_per_noi_coi=int_tra_1d(xs=fs,\
                               ys=psds_noi_coi)
    var_abs_raw_ref=var_per_raw_ref*(cr_mea_ref**2)
    var_abs_noi_ref=var_per_noi_ref*(cr_mea_ref**2)
    var_abs_noi_coi=var_per_noi_coi*(cr_mea_coi**2)

    cov_abs_sig=np.sqrt(  (   var_abs_mea*var_abs_noi_ref\
                            + var_abs_raw_ref*var_abs_noi_coi\
                            + var_abs_noi_coi*var_abs_noi_ref\
                          )\
                        / (2.*n_mer_i*n_mer_f*var_abs_raw_ref)\
                       )
    return cov_abs_mea,\
           cov_abs_sig

def cov_calc_strict(fs,\
                    psds_raw_ref,\
                    psds_noi_ref,\
                    cr_mea_ref,\
                    psds_raw_coi,\
                    psds_noi_coi,\
                    cr_mea_coi,\
                    csds,\
                    ns_mer_i,\
                    ns_mer_f):
    ############
    ### Mean ###
    ############
    csds_ab=abs(csds)
    b2=np.zeros(len(fs))
    for i_b,\
        (psd_raw_ref,\
         psd_noi_ref,\
         psd_raw_coi,\
         psd_noi_coi,\
         n_mer_i,\
         n_mer_f,\
         csd_ab)\
    in enumerate(zip(psds_raw_ref,\
                     psds_noi_ref,\
                     psds_raw_coi,\
                     psds_noi_coi,\
                     ns_mer_i,\
                     ns_mer_f,\
                     csds_ab)):
        b2[i_b], gamma2=bias_calc(psd_raw_1=psd_raw_coi,\
                                  psd_noi_1=psd_noi_coi,\
                                  psd_raw_2=psd_raw_ref,\
                                  psd_noi_2=psd_noi_ref,\
                                  n_mer_i=n_mer_i,\
                                  n_mer_f=n_mer_f,\
                                  csd_ab=csd_ab)

    rands=(csds_ab**2-b2)/(psds_raw_ref-psds_noi_ref)
    cov2_per_mea=int_tra_1d(xs=fs,\
                            ys=rands)
    cov2_abs_mea=cov2_per_mea*(cr_mea_coi**2)
    cov_per_mea=np.sqrt(cov2_per_mea)
    cov_abs_mea=np.sqrt(cov2_abs_mea)

    #############
    ### Sigma ###
    #############
    ### Mean (Start) ###
    f_mea=mean_calc_1d_weigh(xs=fs, ws=ns_mer_f)
    psd_raw_ref_mea=mean_calc_1d_weigh(xs=psds_raw_ref, ws=ns_mer_f)
    psd_noi_ref_mea=mean_calc_1d_weigh(xs=psds_noi_ref, ws=ns_mer_f)
    psd_raw_coi_mea=mean_calc_1d_weigh(xs=psds_raw_coi, ws=ns_mer_f)
    psd_noi_coi_mea=mean_calc_1d_weigh(xs=psds_noi_coi, ws=ns_mer_f)
    csd_mea=mean_calc_1d_weigh(xs=csds, ws=ns_mer_f)
    n_mer_i=ns_mer_i[0]
    n_mer_f=sum(ns_mer_f)

    csd_ab_mea=abs(csd_mea)
    ### Mean (End) ###

    ### Intrinsic coherence (Start) ###
    bb2, gamma2=bias_calc(psd_raw_1=psd_raw_coi_mea,\
                          psd_noi_1=psd_noi_coi_mea,\
                          psd_raw_2=psd_raw_ref_mea,\
                          psd_noi_2=psd_noi_ref_mea,\
                          n_mer_i=n_mer_i,\
                          n_mer_f=n_mer_f,\
                          csd_ab=csd_ab_mea)
    ### Intrinsic coherence (End) ###

    psds_ref=psds_raw_ref-psds_noi_ref
    psds_coi=psds_raw_coi-psds_noi_coi

    rms2_per_mea_ref=int_tra_1d(xs=fs,\
                                ys=psds_ref)
    rms2_per_mea_noi_ref=int_tra_1d(xs=fs,\
                                    ys=psds_noi_ref)
    rms2_per_mea_coi=int_tra_1d(xs=fs,\
                                ys=psds_coi)
    rms2_per_mea_noi_coi=int_tra_1d(xs=fs,\
                                    ys=psds_noi_coi)

    #cr_mea_ref is not needed actually since it is canceled in the calculation of error of covariance.
    rms2_abs_mea_ref=rms2_per_mea_ref*(cr_mea_ref**2)
    rms2_abs_mea_noi_ref=rms2_per_mea_noi_ref*(cr_mea_ref**2)
    rms2_abs_mea_coi=rms2_per_mea_coi*(cr_mea_coi**2)
    rms2_abs_mea_noi_coi=rms2_per_mea_noi_coi*(cr_mea_coi**2)

    cov_per_sig=np.sqrt( (rms2_per_mea_ref+rms2_per_mea_noi_ref)*\
                         ( (1.-gamma2)*rms2_per_mea_coi + rms2_per_mea_noi_coi )/\
                         (2.*n_mer_i*n_mer_f*rms2_per_mea_ref)\
                       )
    cov_abs_sig=np.sqrt( (rms2_abs_mea_ref+rms2_abs_mea_noi_ref)*\
                         ( (1.-gamma2)*rms2_abs_mea_coi + rms2_abs_mea_noi_coi )/\
                         (2.*n_mer_i*n_mer_f*rms2_abs_mea_ref)\
                       )

    return cov_per_mea,\
           cov_per_sig,\
           cov_abs_mea,\
           cov_abs_sig

def write_cov(name_outtxt,\
              cois_min,\
              cois_max,\
              f_min,\
              f_max,\
              covs_per_mea,\
              covs_per_sig,\
              covs_abs_mea,\
              covs_abs_sig,\
              crs_mea,\
              crs_sig,\
              exposure):
    check_extension(name_outtxt, EXTENSION_TXT)
    with open(name_outtxt, 'w') as fout:
        str_head='#F_min[Hz] F_max[rad] Exposure[s]\n'
        str_data=['{:.9g} '.format(f_min),\
                  '{:.9g} '.format(f_max),\
                  '{:.9g}\n'.format(exposure)]
        fout.write(str_head)
        fout.writelines(str_data)
        str_head='#CH_min[-] CH_max[-] Covariance_mean[100%] Covariance_sigma[100%] Covariance_mean[/s] Covariance_sigma[/s] Flux_mean[/s] Flux_sigma[/s]\n'
        fout.write(str_head)
        for coi_min, coi_max, cov_per_mea, cov_per_sig, cov_abs_mea, cov_abs_sig, cr_mea, cr_sig\
        in zip(cois_min, cois_max, covs_per_mea, covs_per_sig, covs_abs_mea, covs_abs_sig, crs_mea, crs_sig):
            str_data=['{:.0f} '.format(coi_min),\
                      '{:.0f} '.format(coi_max),\
                      '{:.9g} '.format(cov_per_mea),\
                      '{:.9g} '.format(cov_per_sig),\
                      '{:.9g} '.format(cov_abs_mea),\
                      '{:.9g} '.format(cov_abs_sig),\
                      '{:.9g} '.format(cr_mea),\
                      '{:.9g}\n'.format(cr_sig)]
            fout.writelines(str_data)
#######################################
########## cov_calc.py (End) ##########
#######################################

########################################
########## rms_pha.py (Start) ##########
########################################
def read_rms(name_intxt):
    check_existence(name_intxt)
    check_extension(name_intxt, EXTENSION_TXT)

    with open(name_intxt, 'r') as fin:
        line=fin.readline()
        line=fin.readline()
        line_str=line.split()
        f_min=float(line_str[0])
        f_max=float(line_str[1])
        expo=float(line_str[2])

    cois_min,\
    cois_max,\
    rmses_per_mea,\
    rmses_per_sig,\
    rmses_abs_mea,\
    rmses_abs_sig,\
    crs_mea,\
    crs_sig\
    =np.loadtxt(fname=name_intxt, \
                dtype='float', \
                skiprows=3, \
                unpack=True)
    return expo,\
           f_min,\
           f_max,\
           cois_min,\
           cois_max,\
           rmses_per_mea,\
           rmses_per_sig,\
           rmses_abs_mea,\
           rmses_abs_sig,\
           crs_mea,\
           crs_sig
######################################
########## rms_pha.py (End) ##########
######################################

########################################
########## cov_pha.py (Start) ##########
########################################
def read_cov(name_intxt):
    check_existence(name_intxt)
    check_extension(name_intxt, EXTENSION_TXT)

    with open(name_intxt, 'r') as fin:
        line=fin.readline()
        line=fin.readline()
        line_str=line.split()
        f_min=float(line_str[0])
        f_max=float(line_str[1])
        expo=float(line_str[2])

    cois_min,\
    cois_max,\
    covs_per_mea,\
    covs_per_sig,\
    covs_abs_mea,\
    covs_abs_sig,\
    crs_mea,\
    crs_sig\
    =np.loadtxt(fname=name_intxt, \
                dtype='float', \
                skiprows=3, \
                unpack=True)

    return expo,\
           f_min,\
           f_max,\
           cois_min,\
           cois_max,\
           covs_per_mea,\
           covs_per_sig,\
           covs_abs_mea,\
           covs_abs_sig,\
           crs_mea,\
           crs_sig
######################################
########## cov_pha.py (End) ##########
######################################
