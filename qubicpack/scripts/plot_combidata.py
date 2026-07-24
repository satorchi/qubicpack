#!/usr/bin/env python3
'''
$Id: plot_combidata.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Thu 02 Sep 2021 16:37:16 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

plot compressor pressure together with PT s2 temperatures
'''
import sys,os,subprocess
import datetime as dt
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
from satorchipy.plotfunctions import plot_flags, labelprops
from satorchipy.datefunctions import utcnow, utcfromtimestamp, str2dt
from qubicpack.utilities import hostname
from qubicpack.housekeeping.utilities import read_hk_file, read_compressor_log, read_hk_flags, read_hk_labels, download_hk, qc_hk_dir
from qubicpack.housekeeping.plot_options import boxprops, plot_options

flag = read_hk_flags(plot_options['events'])

this_year = utcnow().strftime('%Y')
start_date = list(flag.keys())[0]
start_tstamp = float(start_date.strftime('%s.%f'))


basenames = ['AVS47_2_ch0',
             'AVS47_2_ch1']

hk_dir = plot_options['hk_dir']
if plot_options['download']:
    download_hk(basenames,hk_dir,remote_machine=plot_options['qubic-central'])
    
clog = {}
Pmin_list = []
Pmax_list = []
for cnum in [1,2]:
    clog[cnum] = read_compressor_log(hk_dir+'/compressor%i_log.txt' % cnum)
    Pmin_list.append(min(clog[cnum]['Pin']))
    Pmax_list.append(max(clog[cnum]['Pin']))

P_min = min(Pmin_list)
P_max = max(Pmax_list)
if 'Pin_min' in plot_options.keys():
    P_min = plot_options['Pin_min']
    print('option Pin_min=%.1f' % plot_options['Pin_min'])
if 'Pin_max' in plot_options.keys():
    P_max = plot_options['Pin_max']
    
all_labels = read_hk_labels()
label = all_labels['housekeeping']
sensor = {}
date = {}
timestamp = {}
for sensorname in label.keys():
    sensor[sensorname] = None
    date[sensorname] = []
    timestamp[sensorname] = []

files = []
for name in basenames:
    files.append(hk_dir + os.sep + name + '.txt')
minima = []
maxima = []
start_date_list = []
end_date_list = []
start_temp = {}

for cnum in [1,2]:
    start_date_list.append(clog[cnum]['date'][0])
    end_date_list.append(clog[cnum]['date'][-1])
    print('PT%i: end date = %s' % (cnum,clog[cnum]['date'][-1]))

    
for F in files:
    key = os.path.basename(F).replace('.txt','')
    dat = read_hk_file(F)
    if len(dat[0])==0: continue
    idxrange = np.where(dat[0] > start_tstamp)[0]
    if len(idxrange)==0: continue
    start_idx = idxrange[0]
    sensor[key] = dat[1][start_idx:]
    minima.append(min(sensor[key]))
    maxima.append(max(sensor[key]))
    for tstamp in dat[0][start_idx:]:
        date[key].append(utcfromtimestamp(tstamp))
    start_date_list.append(date[key][0])
    end_date_list.append(date[key][-1])
    start_temp[key] = dat[1][start_idx]
    print('%s: end date = %s' % (key,date[key][-1]))


if 'cryomin' in plot_options.keys():
    cryomin = plot_options['cryomin']
else:
    cryomin = min(minima)

if 'cryomax' in plot_options.keys():
    cryomax = plot_options['cryomax']
else:
    cryomax = max(maxima)

date_margin = dt.timedelta(minutes=5)
if 'tstart' in plot_options.keys() and plot_options['tstart'] is not None:
    start_date = plot_options['tstart']
else:
    start_date = min(start_date_list) - date_margin

if 'tend' in plot_options.keys() and plot_options['tend'] is not None:
    end_date = plot_options['tend']
else:
    end_date = max(end_date_list) + date_margin
    
fig = plt.figure()
fig.canvas.manager.set_window_title('plt: combiplot_%s' % max(end_date_list).strftime('%Y%m%d'))
fig.suptitle('4K cold head temperature and compressor return pressure',fontsize=20)
ax = fig.add_axes((0.05,0.1,0.88,0.82))
curves = []

colour = {}
colour[1] = 'blue'
colour[2] = 'green'
for key in sensor.keys():
    if sensor[key] is None: continue
    if label[key].find('PT1')>=0:
        curves += ax.plot(date[key],sensor[key],ls='none',marker='d',label=label[key],color=colour[1])
    elif label[key].find('PT2')>=0:
        curves += ax.plot(date[key],sensor[key],ls='none',marker='d',label=label[key],color=colour[2])
    else:
        curves += ax.plot(date[key],sensor[key],ls='none',marker='d',label=label[key])
ax.set_ylabel('temperature / K')
ax.set_xlabel('date / UT')
ax.tick_params(axis='both',labelsize=20)
ax.set_xlim(start_date,end_date)
ax.set_ylim(cryomin,cryomax)

axcomp = ax.twinx()
for cnum in [1,2]:
    curves += axcomp.plot(clog[cnum]['date'],clog[cnum]['Pin'],ls='none',marker='.',label='Compressor %i' % cnum,color=colour[cnum])
    
axcomp.set_ylabel('Compressor Return Pressure / bar',rotation=270,va='bottom')
axcomp.tick_params(axis='both',labelsize=20)
axcomp.set_ylim(P_min,P_max)


plot_flags(ax,flag)


labels = [l.get_label() for l in curves]
ax.legend(curves, labels, loc='lower left',facecolor='wheat',framealpha=0.5)

pngname = 'combiplot_%s.png' % max(end_date_list).strftime('%Y%m%d')
fig.savefig(pngname,format='png',dpi=300,bbox_inches='tight')
if hostname.find('qubic-central')<0: plt.show()
