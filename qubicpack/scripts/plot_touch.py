#!/usr/bin/env python3
'''
$Id: plot_touch.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Wed 09 Mar 2022 13:23:03 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

plot the mechanical heat switch positions and the touch resistor together
'''
import os,sys,re
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
from satorchipy.datefunctions import str2dt, utcfromtimestamp
from satorchipy.plotfunctions import labelprops,plot_flags,mouse_click_date,get_colour
from qubicpack.housekeeping.plot_options import boxprops,plot_options
from qubicpack.housekeeping.utilities import read_hk_flags, read_hk_labels, read_hk_file, find_pt_start, download_hk, qc_hk_dir
flag = read_hk_flags(plot_options['events'])

touchkey = 'AVS47_1_ch0'
#mechkeys = ['MHS1','MHS2']
mechkeys = ['MHS1']
tempkeys = ['AVS47_2_ch0','AVS47_1_ch1','AVS47_2_ch4']
'''
tempkeys = ['AVS47_2_ch1',
            'AVS47_2_ch0',
            'TEMPERATURE12',
            'TEMPERATURE10',
            'TEMPERATURE04',
            'TEMPERATURE05',
            'AVS47_2_ch2']
'''
#tempkeys = []
tempkeys = ['AVS47_2_ch0','TEMPERATURE16']
upskey = 'ups_log'
pressurekey = 'PRESSURE1'
heaterkeys = ['HEATER5_Amp','HEATER3_Amp']
#basenames = mechkeys + tempkeys + [touchkey] + [pressurekey] + [upskey]
basenames = mechkeys + [touchkey] + [pressurekey] + tempkeys + [upskey] + heaterkeys

hk_dir = plot_options['hk_dir']
if plot_options['download']:
    download_hk(basenames,hk_dir,remote_machine=plot_options['qubic-central'])

all_labels = read_hk_labels()
label = all_labels['housekeeping']
sensor = {}
date = {}
for key in label.keys():
    sensor[key] = None
    date[key] = []

# start date
start_date = None
date_margin = dt.timedelta(seconds=10)

# find most recent pulse tube start
for labeldate in flag.keys():
    match = re.search('(PT|pt|(P|p)ulse (T|t)ube).* (S|s)tart',flag[labeldate])
    if match: start_date = labeldate

if 'tstart' in plot_options.keys():
    start_date = plot_options['tstart']

end_date = None
if 'tend' in plot_options.keys():
    end_date = plot_options['tend']
    
if start_date is None: start_date = list(flag.keys())[0]
print('plotting points beginning %s' % start_date.strftime('%Y-%m-%d %H:%M:%S'))
start_tstamp = start_date.timestamp()
    
files = []
for name in basenames:
    files.append(hk_dir + os.sep + name + '.txt')
minima = []
maxima = []
start_date_list = []
end_date_list = []
for F in files:
    key = os.path.basename(F).replace('.txt','')
    print('reading file: %s' % F)
    dat = read_hk_file(F)
    if len(dat[0])==0: continue
    idxrange = np.where(dat[0] > start_tstamp)[0]
    if len(idxrange)==0: continue
    start_idx = idxrange[0]
    sensor[key] = dat[1][start_idx:]
    if key.find('MHS')>=0:
        sensor[key] /= 1000
    minima.append(min(sensor[key]))
    maxima.append(max(sensor[key]))
    for tstamp in dat[0][start_idx:]:
        date[key].append(utcfromtimestamp(tstamp))
    start_date_list.append(date[key][0])
    end_date_list.append(date[key][-1])

# read UPS input voltage
key = 'ups'
upslog = hk_dir + os.sep + 'ups_log.txt'
sensor[key] = None
Vin = []
Vdate = []
if os.path.exists(upslog):
    print('reading file: %s' % upslog)
    h = open(upslog,'r')
    lines = h.read().split('\n')
    h.close()
    del(lines[-1])
    for line in lines:
        col = line.split()
        if len(col)<3:continue
        Vdate.append(str2dt(col[0].strip()).replace(tzinfo=dt.UTC))
        val_str = col[1].split('=')[-1]
        Vin.append(eval(val_str))
    if len(Vin)>0:
        sensor[key] = np.array(Vin)
        date[key] = np.array(Vdate)
        start_date_list.append(Vdate[0])
        end_date_list.append(Vdate[-1])
        label[key] = 'V$_\\mathrm{in}$'

if 'mechmin' in plot_options.keys():
    mechmin = plot_options['mechmin']
else:
    mechmin = min(minima)

if 'mechmax' in plot_options.keys():
    mechmax = plot_options['mechmax']
else:
    mechmax = max(maxima)

if 'Tmin' in plot_options.keys():
    Tmin = plot_options['Tmin']
else:
    Tmin = min(minima)
if 'cryomin' in plot_options.keys():
    Tmin = plot_options['cryomin']

if 'Tmax' in plot_options.keys():
    Tmax = plot_options['Tmax']
else:
    Tmax = max(maxima)
if 'cryomax' in plot_options.keys():
    Tmax = plot_options['cryomax']

if touchkey not in sensor.keys() or sensor[touchkey] is None:
    plot_touch = False
    touchmin = 1e-4
    touchmax = 2e5
else:
    plot_touch = True
    touchmin = min(sensor[touchkey])
    touchmax = max(sensor[touchkey])
    
if 'touchmin' in plot_options.keys():
    touchmin = plot_options['touchmin']    

if 'touchmax' in plot_options.keys():
    touchmax = plot_options['touchmax']

if 'pmin' in plot_options.keys():
    pmin = plot_options['pmin']
else:
    pmin = 1e-7

if 'pmax' in plot_options.keys():
    pmax = plot_options['pmax']
else:
    pmax = 1e3

if 'Vmin' in plot_options.keys():
    Vmin = plot_options['Vmin']
else:
    Vmin = 0

if 'Vmax' in plot_options.keys():
    Vmax = plot_options['Vmax']
else:
    Vmax = 240

plot_heater = False
heater_mins = []
heater_maxs = []
for key in heaterkeys:
    if key not in sensor.keys() or sensor[key] is None: continue
    plot_heater = True
    heater_mins.append(sensor[key].min())
    heater_maxs.append(sensor[key].max())
if plot_heater:
    heater_min = min(heater_mins)
    heater_max = max(heater_maxs)
if 'heater_min' in plot_options.keys():
    heater_min = plot_options['heater_min']
if 'heater_max' in plot_options.keys():
    heater_max = plot_options['heater_max']
    
if end_date is None:
    end_date = max(end_date_list) + date_margin
    
    
print('using start date: %s' % start_date.strftime('%Y-%m-%d %H:%M:%S'))
print('using end date: %s' % end_date.strftime('%Y-%m-%d %H:%M:%S'))

padding = 0

curves = []
fig = plt.figure()
figname = 'MHS-touch_%s' % end_date.strftime('%Y%m%d')
figttl = 'MHS - Touch %s' % end_date.strftime('%Y-%m-%d')
fig.canvas.manager.set_window_title('plt: '+figname)
ax = fig.add_axes((0.06,0.14,0.63,0.80))
ax.text(0.5,1.02,figttl,va='bottom',ha='center',transform=ax.transAxes,fontsize=25)
ax.set_xlim(start_date,end_date)
ax.set_ylim(touchmin,touchmax)
ax.set_xlabel('Date / UT')

if plot_touch:
    key = touchkey
    curves += ax.plot(date[key],sensor[key],ls='none',color='red',marker='+',label=label[key])
    ax.set_ylabel('Touch / $\\Omega$',color='red')
    ax.set_yscale('log')
    ax.tick_params(axis='y', labelcolor='red')
    ax.tick_params(axis='x', rotation=315, left=True, labelleft=True,labelsize=12)
    #ax.set_xticklabels(ax.get_xticklabels(),ha='left')

# plot the MHS positions
axmech = ax.twinx()
marker = 'd'
colours = ['blue','green',]
for idx,key in enumerate(mechkeys):
    if sensor[key] is None: continue
    labeltxt = '%s' % (label[key])
    curves += axmech.plot(date[key],sensor[key],ls='none',color=colours[idx],marker=marker,label=labeltxt)
axmech.set_ylabel('MHS position / kilosteps',rotation=270,color=colours[0],va='bottom',ha='left')
axmech.set_ylim(mechmin,mechmax)
axmech.tick_params(axis='y', labelcolor=colours[0],pad=padding)
padding += 100

# plot the temperatures
axtemp = ax.twinx()
marker = 'v'
colours = ['green','magenta','black','#1f77b4ff','#a20cffff','olive','red','blue','purple','cyan']
ntemperatures = 0
for idx,key in enumerate(tempkeys):
    if sensor[key] is None: continue
    labeltxt = '%s' % (label[key])
    print('plotting: %s' % labeltxt)
    curves += axtemp.plot(date[key],sensor[key],ls='none',color=colours[idx],marker=marker,label=labeltxt)
    ntemperatures += 1
if ntemperatures>0:
    axtemp.set_ylabel('Temperature / K',rotation=270,color=colours[0],va='bottom',ha='center')
    axtemp.tick_params(axis='y', labelcolor=colours[0],pad=padding)
    padding += 100
    #axtemp.yaxis.set_label_coords(1.12,0.5)
    axtemp.set_ylim(Tmin,Tmax)
else:
    axtemp.set_visible(False)

if upskey in basenames and sensor['ups'] is not None:
    # plot the Voltage supply
    marker = 'p'
    colour = 'orange'
    key = 'ups'
    axups = ax.twinx()
    labeltxt = '%s' % (label[key])
    curves += axups.plot(date[key],sensor[key],ls='none',color=colour,marker=marker,label=labeltxt)
    axups.set_ylabel('Supply Voltage / VAC',rotation=270,color=colour,va='bottom',ha='right')
    axups.tick_params(axis='y', labelcolor=colour,pad=padding)
    padding += 100
    # axups.yaxis.set_label_coords(1.20,0.5)
    axups.set_ylim(Vmin,Vmax)

if pressurekey in basenames and sensor[pressurekey] is not None:
    # plot the pressure
    key = pressurekey
    axpressure = ax.twinx()
    marker = 'p'
    colour = '#1f77b4ff'
    labeltxt = '%s' % (label[key])
    curves += axpressure.plot(date[key],sensor[key],ls='none',color=colour,marker=marker,label=labeltxt)
    axpressure.set_ylabel('Pressure / mbar',rotation=270,color=colour,va='bottom',ha='right')
    axpressure.tick_params(axis='y', labelcolor=colour,pad=padding)
    padding += 100
    # axpressure.yaxis.set_label_coords(1.20,0.5)
    axpressure.set_ylim(pmin,pmax)
    axpressure.set_yscale('log')

if plot_heater:
    axheater = ax.twinx()
    nheaters = 0
    marker = 'o'
    for idx,key in enumerate(heaterkeys):
        if sensor[key] is None: continue
        labeltxt = '%s' % (label[key])
        print('plotting: %s' % labeltxt)
        heater_colour = colours[idx]
        curves += axheater.plot(date[key],sensor[key],ls='none',color=heater_colour,marker=marker,label=labeltxt)
        nheaters += 1
    if nheaters>0:
        axheater.set_ylabel('Heater / Amp',rotation=270,color=heater_colour,va='bottom',ha='center')
        axheater.tick_params(axis='y', labelcolor=heater_colour,pad=padding)
        padding += 100
        axheater.set_ylim(heater_min,heater_max)
    else:
        axheater.set_visible(False)
    

# the combined legend
labels = [l.get_label() for l in curves]
ax.legend(curves, labels, loc='upper left',facecolor='wheat',framealpha=0.5)

minmax = ax.axis()[2:]
if 'flagpos' not in plot_options.keys():
    flagpos = np.log10(minmax[0]) + 0.3*( np.log10(minmax[1]) - np.log10(minmax[0]))
    flagpos = 10**flagpos
else:
    flagpos = plot_options['flagpos']

plot_flags(ax,flag,flagpos)
pngname = '%s.png' % figname
fig.savefig(pngname,format='png',dpi=300,bbox_inches='tight')
plt.connect('button_press_event',mouse_click_date)
plt.show()
