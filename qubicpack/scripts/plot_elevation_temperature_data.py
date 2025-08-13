#!/usr/bin/env python3
'''
$Id: plot_elevation_temperature_data.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Thu 11 Jul 16:11:40 CEST 2019
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.


plot the Pulse Tube inclination angle and temperature
'''
import sys,os,pickle
import numpy as np
from matplotlib import pyplot as plt
import datetime as dt
from qubicpack.utilities import figure_window_title

if len(sys.argv)==1:
    whichdata = ''
else:
    whichdata = '_%s' % sys.argv[1]

temperature_file = 'temperatures%s.dat' % whichdata
if not os.path.isfile(temperature_file):
    print('temperature data not found: %s' % temperature_file)
    quit()

elevation_file = 'elevation%s.dat' % whichdata
if not os.path.isfile(elevation_file):
    print('elevation data not found: %s' % elevation_file)
    quit()

          



# elevation data
h = open(elevation_file,'rb')
el_rec = pickle.load(h)
h.close()
el_npts = len(el_rec)
el = np.zeros(el_npts)
el_date = np.zeros(el_npts)
for idx,rec in enumerate(el_rec):
    el[idx] = rec.elevation
    el_date[idx] = rec.timestamp - 7200 # QubicStudio on localtime
el_t0 = min(el_date)

# temperature data
h = open(temperature_file,'rb')
t_rec = pickle.load(h)
h.close()
t_npts = len(t_rec)

rec_names = {}
rec_names['pt1s1'] = 'Temp_5'
rec_names['pt2s1'] = 'Temp_4'
rec_names['pt1s2'] = 'Temp_11'
rec_names['pt2s2'] = 'Temp_12'
rec_names['bath'] = 'AVS47_1_CH2'

t_pts = {}
t_date = np.zeros(t_npts)
for key in rec_names.keys():
    t_pts[key] = np.zeros(t_npts)

for idx,rec in enumerate(t_rec):
    t_date[idx] = rec.RaspberryDate
    for key in rec_names.keys():
        cmd = 't_pts[key][idx] = rec.%s' % rec_names[key]
        exec(cmd)
t_t0 = min(t_date)

t0 = min([el_t0,t_t0])
t0_date = dt.datetime.fromtimestamp(t0)

# make plots
plt.ion()
for key in rec_names.keys():
    idx_good = np.where(t_pts[key] != -1)[0]
    pngname = '%s_vertical-offset_temperature_%s.png' % (key.upper(),t0_date.strftime('%Y%m%d'))
    ttl = 'Temperture during Synthetic Beam Mapping: %s' % t0_date.strftime('%Y-%m-%d')
    t_label = '%s Temperature' % key.upper()
    
    fig = plt.figure(figsize=(16,8))
    figure_window_title(fig,'Temperature vs Inclination')
    plt.title(ttl)  

    date_pts = (t_date[idx_good] - t0)/60
    temp_pts = t_pts[key][idx_good]
                                                         
    plt.plot(date_pts,temp_pts,ls='none',marker='*',color='green',label=t_label)

    t_ax = fig.axes[0]
    el_ax = t_ax.twinx()

    date_pts = (el_date-t0)/60
    el_pts = np.abs(el - 50)

    el_ax.plot(date_pts,el_pts,ls='none',marker='D',color='red',label='Tilt from Vertical')
    el_ax.set_ylabel('Tilt angle from vertical / degrees')           
    el_ax.legend(loc=(0.05,0.075))        

    t_ax.set_xlabel('Time / minutes since %s' % t0_date.strftime('%Y-%m-%d %H:%M:%S UT'))
    t_ax.set_ylabel('%s / K' % t_label)
    t_ax.legend(loc=(0.05,0.025))

    plt.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    



ans = input('press return to exit ')


