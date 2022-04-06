import numpy as np
import sys
from sys import argv
from lightcurve import *
# ---------------------------------------------------------- #
# Extract data which can be compared between LE, ME and HE   #
# due to common sampling times.                              #
# ---------------------------------------------------------- #
def main():
    # --------------------------- #
    # ----- Setting (begin) ----- #
    # --------------------------- #
    if not len(sys.argv)==7:
        print('Error: Number of inputs is wrong(0).(must be {1})'\
            .format(len(sys.argv), 7))
        sys.exit()
    # Filename of input light curve (LE)
    name_inlc_le=argv[1]
    # Filename of input light curve (ME)
    name_inlc_me=argv[2]
    # Filename of input light curve (HE)
    name_inlc_he=argv[3]
    # Filename of output common index (LE)
    name_outtxt_le=argv[4]
    # Filename of output common index (ME)
    name_outtxt_me=argv[5]
    # Filename of output common index (HE)
    name_outtxt_he=argv[6]
    # --------------------------- #
    # ----- Setting (end)   ----- #
    # --------------------------- #

    # ------------------------------------ #
    # ----- Read light curve (begin) ----- #
    # ------------------------------------ #
    lc_le=LightCurve()
    lc_le.read_data(name_file=name_inlc_le)

    lc_me=LightCurve()
    lc_me.read_data(name_file=name_inlc_me)

    lc_he=LightCurve()
    lc_he.read_data(name_file=name_inlc_he)

    if lc_le.dt!=lc_me.dt or lc_le.dt!=lc_he.dt:
        print('Error: Time resolution is different.')
        sys.exit()
    dt=lc_le.dt
    # ---------------------------------- #
    # ----- Read light curve (end) ----- #
    # ---------------------------------- #

    # ------------------------------------- #
    # ----- Print information (begin) ----- #
    # ------------------------------------- #
    lc_le.print_obs_info()
    lc_me.print_obs_info()
    lc_he.print_obs_info()
    # ------------------------------------- #
    # ----- Print information (end)   ----- #
    # ------------------------------------- #

    # ---------------------------------------------- #
    # ----- Extract common light curve (begin) ----- #
    # ---------------------------------------------- #
    lc_le.data['TIME']
    first=True
    n_le=len(lc_le.data['TIME'])
    is_le=np.empty(0)
    is_me=np.empty(0)
    is_he=np.empty(0)
    for i_le, time_le in enumerate(lc_le.data['TIME']):
        if i_le%10000==0:
            print('{0:.0f}/{1:.0f}'.format(i_le, n_le))
        i_me=np.where(np.abs(time_le-lc_me.data['TIME'])<dt/2.)[0]
        i_he=np.where(np.abs(time_le-lc_he.data['TIME'])<dt/2.)[0]
        if (len(i_me)==1) & (len(i_he)==1):
            is_le=np.append(is_le, i_le)
            is_me=np.append(is_me, i_me)
            is_he=np.append(is_he, i_he)
    # ---------------------------------------------- #
    # ----- Extract common light curve (end)   ----- #
    # ---------------------------------------------- #

    # --------------------------------- #
    # ----- Write results (begin) ----- #
    # --------------------------------- #
    with open(name_outtxt_le, 'w') as fout:
        str_head='#Index(min)[-] Index(max)[-]\n'
        fout.write(str_head)

        first=True
        n=len(is_le)
        for i, i_com in enumerate(is_le):
            if first==True:
                first=False
                i_start=i_com
                i_pre=i_com
                continue
            else:
                if (i_com-i_pre==1) & (i!=n-1):
                    i_pre=i_com
                    continue
                else:
                    if i==n-1:
                        i_end=i_com
                    else:
                        i_end=i_pre
                    str_data=[\
                        '{:.0f} ' .format(i_start),\
                        '{:.0f}\n' .format(i_end)]
                    fout.writelines(str_data)
                    i_start=i_com
                    i_pre=i_com

    with open(name_outtxt_me, 'w') as fout:
        str_head='#Index(min)[-] Index(max)[-]\n'
        fout.write(str_head)

        first=True
        n=len(is_me)
        for i, i_com in enumerate(is_me):
            if first==True:
                first=False
                i_start=i_com
                i_pre=i_com
                continue
            else:
                if (i_com-i_pre==1) & (i!=n-1):
                    i_pre=i_com
                    continue
                else:
                    if i==n-1:
                        i_end=i_com
                    else:
                        i_end=i_pre
                    str_data=[\
                        '{:.0f} ' .format(i_start),\
                        '{:.0f}\n' .format(i_end)]
                    fout.writelines(str_data)
                    i_start=i_com
                    i_pre=i_com

    with open(name_outtxt_he, 'w') as fout:
        str_head='#Index(min)[-] Index(max)[-]\n'
        fout.write(str_head)

        first=True
        n=len(is_he)
        for i, i_com in enumerate(is_he):
            if first==True:
                first=False
                i_start=i_com
                i_pre=i_com
                continue
            else:
                if (i_com-i_pre==1) & (i!=n-1):
                    i_pre=i_com
                    continue
                else:
                    if i==n-1:
                        i_end=i_com
                    else:
                        i_end=i_pre
                    str_data=[\
                        '{:.0f} ' .format(i_start),\
                        '{:.0f}\n' .format(i_end)]
                    fout.writelines(str_data)
                    i_start=i_com
                    i_pre=i_com

    # --------------------------------- #
    # ----- Write results (end)   ----- #
    # --------------------------------- #

    print('Common light curves were investigated successfully.')
    #print('Results are stored in {0}.'.format(name_outfits))

if __name__=='__main__':
    main()
