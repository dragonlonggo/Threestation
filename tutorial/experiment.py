#加载库
from os.path import join

import matplotlib.pyplot as plt
import numpy as np
import obspy
from obspy import read, Stream, Trace, UTCDateTime
from obspy.signal.cross_correlation import correlate
from obspy.geodetics.base import gps2dist_azimuth
import pandas as pd


plt.rcParams.update({
    'figure.dpi': 200,
})


DATADIR = './data'
#读位置数据
meta = pd.read_csv(join(DATADIR, 'meta.csv'))
# meta
#经纬度坐标绘图
rec_sta_1 = 'M07A'
rec_sta_2 = 'M15A'

lon_r1 = meta[meta['sta'] == rec_sta_1]['lon'].iloc[0]
lat_r1 = meta[meta['sta'] == rec_sta_1]['lat'].iloc[0]
lon_r2 = meta[meta['sta'] == rec_sta_2]['lon'].iloc[0]
lat_r2 = meta[meta['sta'] == rec_sta_2]['lat'].iloc[0]


fig, ax = plt.subplots()
ax.scatter(meta['lon'], meta['lat'], s=5, marker='o', c='r')
for s in [rec_sta_1, rec_sta_2]:
    m = meta[meta.sta == s]
    ax.scatter(m['lon'], m['lat'], s=20, marker='^', c='b')
#读地震数据，绘水平截面图
I2_r1 = obspy.read(join(DATADIR, 'I2', rec_sta_1, '*.SAC'))
I2_r2 = obspy.read(join(DATADIR, 'I2', rec_sta_2, '*.SAC'))

t0 = UTCDateTime('2021-03-24T00:00:00')

for I2 in [I2_r1, I2_r2]:
    for tr in I2:
        tr.stats.distance = tr.stats.sac.dist * 1e3
        tr.stats.starttime = t0

I2_r1 = sort_stream(I2_r1, 'dist')
I2_r2 = sort_stream(I2_r2, 'dist')

for I2 in [I2_r1, I2_r2]:
    
#     I2.filter('bandpass', freqmax=1/20, freqmin=1/50, zerophase=True)

    fig = I2[::40].plot(
        type='section', orientation='horizontal', dpi=200,
        reftime=t0+3000, recordstart=-1000, recordlength=2000,
        linewidth=.5, alpha=1, grid_linewidth=0, scale=.5,
#         offset_min=0, offset_max=1500e3,
    )
