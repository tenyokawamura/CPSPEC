import xspec

names_outtxt=[\
    'ni1200120130_0mpu7_cl_fpm_0101_0260_0261_0480_csd_phase.txt'\
    ]

xspec.AllData('\
    1:1 ni1200120130_0mpu7_cl_fpm_0101_0260_0261_0480_csd_phase.pha \
    ')
xspec.AllData.ignore('1:1 52')
m1=xspec.Model('powerlaw')
xspec.Plot('uf')
fs=xspec.Plot.x()
dfs=xspec.Plot.xErr()
for i, name_outtxt in enumerate(names_outtxt):
    data=xspec.Plot.y(i+1)
    ddata=xspec.Plot.yErr(i+1)
    with open(name_outtxt, 'w') as fout:
        str_header='#f[Hz] df[Hz] Phase[rad] dPhase[rad]\n'
        fout.write(str_header)
        for f, df, datum, ddatum, in zip(fs, dfs, data, ddata):
            str_data=['{:.9f} ' .format(f),\
                      '{:.9f} ' .format(df),\
                      '{:.9f} ' .format(datum),\
                      '{:.9f}\n'.format(ddatum)\
                     ]
            fout.writelines(str_data)

