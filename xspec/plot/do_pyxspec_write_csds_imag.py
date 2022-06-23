import xspec
from include import inputs_xspec_plot

for i_n, (ch_min, ch_max) in enumerate(zip(inputs_xspec_plot.chs_min, inputs_xspec_plot.chs_max)):
    name_inpha =inputs_xspec_plot.name_inevt.replace('.evt', '_{0:04}_{1:04}_{2:04}_{3:04}_csd_imag.pha'.format(int(ch_min), int(ch_max), int(inputs_xspec_plot.ch_ref_min), int(inputs_xspec_plot.ch_ref_max)))
    name_outtxt=inputs_xspec_plot.name_inevt.replace('.evt', '_{0:04}_{1:04}_{2:04}_{3:04}_csd_imag_xspec.txt'.format(int(ch_min), int(ch_max), int(inputs_xspec_plot.ch_ref_min), int(inputs_xspec_plot.ch_ref_max)))

    s1=xspec.Spectrum(name_inpha)
    chan_min=s1.noticed[0]
    chan_max=s1.noticed[-1]
    xspec.AllData.ignore('1:{0} {1}'.format(int(chan_min), int(chan_max)))
    m1=xspec.Model('powerlaw')
    xspec.Plot('uf')
    fs=xspec.Plot.x()
    dfs=xspec.Plot.xErr()

    data=xspec.Plot.y(1)
    ddata=xspec.Plot.yErr(1)
    with open(name_outtxt, 'w') as fout:
        str_header='#f[Hz] df[Hz] Im[C(f)] dIm[C(f)}\n'
        fout.write(str_header)
        for f, df, datum, ddatum, in zip(fs, dfs, data, ddata):
            str_data=['{:.9f} ' .format(f),\
                      '{:.9f} ' .format(df),\
                      '{:.9f} ' .format(datum),\
                      '{:.9f}\n'.format(ddatum)\
                     ]
            fout.writelines(str_data)

    xspec.AllData.clear()

