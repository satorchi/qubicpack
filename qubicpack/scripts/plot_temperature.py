#!/usr/bin/env python3
import sys,os,subprocess,re
import datetime as dt
from glob import glob
from copy import copy

from satorchipy.datefunctions import str2dt, utcnow, utcfromtimestamp, tstamp2dt
from satorchipy.plotfunctions import labelprops, plot_flags, make_legend_label, mouse_click_date, plot_dayboundaries
from qubicpack.utilities import hostname
from qubicpack.housekeeping.utilities import read_hk_file, find_pt_start, read_hk_labels, read_hk_flags, download_hk
from qubicpack.housekeeping.plot_options import boxprops,plot_options

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

hk_dir = plot_options['hk_dir']
flag = read_hk_flags(plot_options['events'])

all_labels = read_hk_labels()
label = all_labels['housekeeping']
basename = {}
for bname in label.keys():
    lbl = label[bname]
    basename[lbl] = bname

datefmt = '%Y-%m-%d %H:%M:%S'

def line_model(x,m,b):
    return m*x + b

T_final = {}
T_final['TEMPERATURE04'] = 40.0
T_final['TEMPERATURE05'] = 40.0
T_final['TEMPERATURE11'] = 4.0
T_final['TEMPERATURE12'] = 4.0
#T_final['AVS47_1_ch2'] = 4.0 # TES stage
T_final['AVS47_2_ch0'] = 4.0
T_final['AVS47_2_ch1'] = 4.0

# when we can start cryo pumps:  see https://elog-qubic.in2p3.fr/demo/1249
T_final['AVS47_2_ch2'] = 10.0 # Fridge plate when we can start the cryo pumps
T_final['AVS47_1_ch1'] = 25.0 # 1K stage, when we can start the cryo pumps


# bad sensors
# 'TEMPERATURE05',
# 'TEMPERATURE10',
# 'AVS47_2_ch3', # a bit wonky

temperatures = plot_options['temperatures']
if temperatures is None:
    basenames = ['TEMPERATURE17',
                 'TEMPERATURE14',
                 'TEMPERATURE16',
                 'TEMPERATURE04',
                 'AVS47_1_ch2',
                 'AVS47_2_ch0',
                 'AVS47_2_ch1',
                 'AVS47_1_ch4',
                 'AVS47_2_ch4',
                 'AVS47_1_ch1',
                 'AVS47_2_ch2',
                 'TEMPERATURE18',
                 ]
    temperatures = []
    for bname in basenames:
        temperatures.append(label[bname])
else:
    basenames = []
    for lbl in temperatures:
        basenames.append(basename[lbl])

print('\n    '.join(['making plot for the following temperatures:']+temperatures))
    
#heaters = ['HEATER3']
heaters = [] # takes too long, needs debugging

if plot_options['pressure']:
    pressurekey = 'PRESSURE1'
    basenames += [pressurekey]
    
for b in heaters:
    for keytype in ['Volt','Amp']:
        basenames.append(b+'_'+keytype)

if plot_options['download']:
    download_hk(basenames,plot_options['hk_dir'],remote_machine=plot_options['qubic-central'])
        
sensor = {}
date = {}
timestamp = {}
for key in label.keys():
    sensor[key] = None
    date[key] = []
    timestamp[key] = []
    

# temperature files
files = []
for name in basenames:
    files.append(hk_dir + os.sep + name + '.txt')
# heater files
for name in heaters:
    files.append(hk_dir + os.sep + name + '_Volt.txt')
    files.append(hk_dir + os.sep + name + '_Amp.txt')

# pressure
if plot_options['pressure']:
    date[pressurekey] = []
    label[pressurekey] = 'pressure'

# find start time and end time
pt_start_date = find_pt_start(flag)
start_date = plot_options['tstart']
if start_date is None:
    if pt_start_date is None:
        start_date = list(flag.keys())[0]
    else:
        start_date = pt_start_date
print('plotting from %s' % start_date.strftime(datefmt))
start_tstamp = start_date.timestamp()
    
minima = []
maxima = []
start_date_list = []
end_date_list = []
start_temp = {}
for F in files:
    print('[%s] processing: %s' % (dt.datetime.now().strftime(datefmt),F))
    key = os.path.basename(F).replace('.txt','')
    print('[%s] reading file' % dt.datetime.now().strftime(datefmt),end='...',flush=True)
    dat = read_hk_file(F)
    print('done')
    if dat[0] is None:
        print('[%s] no timestamp data for %s' % (dt.datetime.now().strftime(datefmt),key))
        continue
    if dat[1] is None:
        print('[%s] no data for %s' % (dt.datetime.now().strftime(datefmt),key))
        continue
    if len(dat[0])==0:
        print('[%s] zero length time data for %s' % (dt.datetime.now().strftime(datefmt),key))
        continue
    print('[%s] range mask' % dt.datetime.now().strftime(datefmt),end='...',flush=True)
    rangemask = dat[0] > start_tstamp
    print('done')
    #idxrange = np.where(dat[0] > start_tstamp)[0]
    if rangemask.sum()==0:
        print('[%s] no data in range for %s' % (dt.datetime.now().strftime(datefmt),key))
        continue
    #start_idx = idxrange[0]
    
    sensor[key] = dat[1][rangemask]
    timestamp[key] = dat[0][rangemask]

    print('[%s] finding min/max' % dt.datetime.now().strftime(datefmt),end='...',flush=True)
    if key.find('HEATER')==0:
        onoff = dat[2][rangemask]
        val = dat[1][rangemask]
        tstamps = dat[0][rangemask]
        sensor[key] = val[onoff]
        timestamp[key] = tstamps[onoff]
    else:
        minima.append(min(sensor[key]))
        maxima.append(max(sensor[key]))
        start_temp[key] = dat[1][rangemask][0]
    print('done')

    print('[%s] converting timestamps to date' % dt.datetime.now().strftime(datefmt),end='...',flush=True)
    date[key] = tstamp2dt(timestamp[key])
    print('done')
    
    if len(date[key])==0:
        print('[%s] timestamps could not be converted to dates for %s' % (dt.datetime.now().strftime(datefmt),key))
        continue
    print('[%s] processed:  %s' % (dt.datetime.now().strftime(datefmt),F))
    start_date_list.append(date[key][0])
    end_date_list.append(date[key][-1])

date_margin = dt.timedelta(minutes=5)
if 'tstart' in plot_options.keys() and plot_options['tstart'] is not None:
    start_date = plot_options['tstart']
else:
    start_date = min(start_date_list) - date_margin

if 'tend' in plot_options.keys() and plot_options['tend'] is not None:
    end_date = plot_options['tend']
else:
    end_date = max(end_date_list) + date_margin

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

if 'pmin' in plot_options.keys():
    pmin = plot_options['pmin']
else:
    pmin = 1e-8

if 'pmax' in plot_options.keys():
    pmax = plot_options['pmax']
else:
    pmax = 1000.0

plot_heaters = plot_options['heaters']
if len(heaters)==0: plot_heaters = False
if plot_heaters:
    # make the power from the Amp and Volt data
    for key in heaters:
        volt_key = key+'_Volt'
        amp_key = key+'_Amp'

        volt_interp = np.interp(timestamp[amp_key], timestamp[volt_key], sensor[volt_key])
        sensor[key] = 0.001*volt_interp*sensor[amp_key]
        date[key] = date[amp_key]
    
    
# fit a straight line to the last fit_npts points and estimate time to 4K
if plot_options['estimate cold date']:
    fit_npts = 500
    estimate_list = ['estimates:']
    estimate_linelen = [len(estimate_list[0])]
    model_slope = {}
    model_offset = {}
    model_date_final = {}
    for fit_sensor in T_final.keys():
        if fit_sensor not in sensor.keys(): continue
        if sensor[fit_sensor] is None: continue

        T_now = sensor[fit_sensor][-1]
        if 'deltafit' in plot_options.keys() and plot_options['deltafit'] is not None:
            tstamp_end = timestamp[fit_sensor][-1]
            tstamp_fitstart = tstamp_end - plot_options['deltafit']
            mask = timestamp[fit_sensor]>=tstamp_fitstart
            fit_npts_tmp = mask.sum()
            if fit_npts_tmp>0:
                fit_npts = fit_npts_tmp
        print('%s: using the last %i points for fitting' % (label[fit_sensor].ljust(20,' '),fit_npts))

        
        ypts = sensor[fit_sensor][-fit_npts:]
        fit_npts = len(ypts)
        xpts = np.empty(fit_npts)
        #print('length xpts,ypts = %i,%i' % (len(xpts),len(ypts)))
        for idx,d in enumerate(date[fit_sensor][-fit_npts:]):
            xpts[idx] = d.timestamp()

        if T_now<T_final[fit_sensor]:
            estimate_list.append('%s is below %.0fK' % (label[fit_sensor],T_final[fit_sensor]))

        else:
            # first guess from first and last points in the set
            slope = (ypts[-1] - ypts[0])/(xpts[-1] - xpts[0])
            offset = -slope*xpts[0] + ypts[0]
            first_guess = np.array((slope,offset))    
            popt,pcov=curve_fit(line_model,xpts,ypts,p0=first_guess)
            slope = popt[0]
            offset = popt[1]
            model_slope[fit_sensor] = slope
            model_offset[fit_sensor] = offset

            
            if np.abs(slope)<1e-8:
                tstamp_final = -1
            else:
                tstamp_final = (T_final[fit_sensor] - offset)/slope

            if slope>0:
                estimate_list.append('%s temperature is rising' % (label[fit_sensor]))
            elif tstamp_final<=0:
                estimate_list.append('%s will never reach %.0fK' % (label[fit_sensor],T_final[fit_sensor]))
            else:
                date_final = utcfromtimestamp(tstamp_final)
                model_date_final[fit_sensor] = date_final
                estimate_list.append('%s will be %.0fK at %s' % (label[fit_sensor],
                                                                 T_final[fit_sensor],
                                                                 date_final.strftime(datefmt)))
        estimate_linelen.append(len(estimate_list[-1]))
        
    linelen = max(estimate_linelen)
    for idx,line in enumerate(estimate_list):
        estimate_list[idx] = line.ljust(linelen)
    estimate_str = '\n'.join(estimate_list)
fig = plt.figure()
fig.canvas.manager.set_window_title('plt: cooldown_%s' % end_date.strftime('%Y%m%d'))
ax = fig.add_axes((0.05,0.1,0.88,0.82))
curves = []
if plot_heaters: ax_heaters = ax.twinx()

temp_diff = {}
temp_diff_txtlist = ['Temperature differences from start of cooling']
heater_colour = ['red','black','green','magenta','cyan','orange']
heater_idx = 0
for idx,key in enumerate(sensor.keys()):
    if sensor[key] is None: continue
    if len(sensor[key])!=len(date[key]):
        print('PROBLEM! %s: %i != %i' % (key,len(sensor[key]),len(date[key])))
        continue

    if key.find('AVS')==0:
        marker = 'd'
    elif key.find('TEMPERATURE')==0:
        marker = 'x'
    else:
        marker = '^'
    # if key.find('HEATER')==0 and key.find('Volt')<0 and key.find('Amp')<0:
    if key.find('HEATER')==0:
        if key.find('Amp')<0 and key.find('Volt')<0:
            colour =heater_colour[heater_idx]
            labeltxt = '%7.03fW: %s' % (sensor[key][-1],label[key])
            curve = ax_heaters.plot(date[key],sensor[key],ls='none',marker='v',color=colour,label=labeltxt)
            curves += curve
            heater_idx += 1
    elif key.find('PRESSURE')<0:
        labeltxt = '%7.03fK: %s' % (sensor[key][-1],label[key])
        curve = ax.plot(date[key],sensor[key],ls='none',marker=marker,label=labeltxt)
        curves += curve
        if plot_options['show differences']:
            colour = curve[0].get_color()
            ax.plot([date[key][0],date[key][-1]],[start_temp[key],start_temp[key]],ls='dotted',color=colour)
            temp_diff[key] = sensor[key][-1] - start_temp[key]
            temp_diff_txtlist.append('%s: %.03fK' % (label[key],temp_diff[key]))
            

if plot_options['show differences']:
    infotxt = '\n'.join(temp_diff_txtlist)
    ax.text(0.25,0.01,infotxt,ha='left',va='bottom',transform=ax.transAxes,bbox=boxprops,fontsize=20)
    
if plot_options['estimate cold date']:
    for fit_sensor in T_final.keys():
        # print('T_final sensor: %s [verifying' % fit_sensor,end=' ...')
        if fit_sensor not in model_date_final.keys():
            # print(' No, not a requested sensor.]')
            continue
        if len(date[fit_sensor])<fit_npts:
            # print(' No, not enough points.]')
            continue
        # print(' OK]')
        tstamp_1 = date[fit_sensor][-fit_npts].timestamp()
        Tmodel_1 = model_slope[fit_sensor]*tstamp_1 + model_offset[fit_sensor]
        tstamp_2 = model_date_final[fit_sensor].timestamp()
        Tmodel_2 = model_slope[fit_sensor]*tstamp_2 + model_offset[fit_sensor]
        ax.plot([date[fit_sensor][-fit_npts],model_date_final[fit_sensor]],
                [Tmodel_1,Tmodel_2],
                ls='dashed',marker='+',markersize=25,color='black')

        ax.plot([date[fit_sensor][-fit_npts],date[fit_sensor][-1]],
                [sensor[fit_sensor][-fit_npts],sensor[fit_sensor][-1]],
                ls='dashed',marker='+',markersize=25,color='green')
        
    ax.text(0.9,0.99,estimate_str,ha='right',va='top',transform=ax.transAxes,bbox=boxprops,fontsize=20)

ax.set_ylim(Tmin,Tmax)
ax.set_xlim(start_date,end_date)
ax.text(0.5,1.01,'QUBIC temperatures up to %s' % max(end_date_list).strftime(datefmt),
        fontsize=24,ha='center',va='bottom',
        transform=ax.transAxes)
ax.set_ylabel('temperature / K',fontsize=24)
ax.set_xlabel('date (%s-DD HH:MM) / UT' % start_date.strftime('%Y-%m'))
ax.tick_params(axis='x', rotation=315, left=True, labelleft=True,labelsize=12)

minmax = ax.axis()[2:]
if 'flagpos' not in plot_options.keys():
    flagpos = minmax[0] + 0.3*(minmax[1] - minmax[0])
else:
    flagpos = plot_options['flagpos']

if 'noflags' not in plot_options.keys():
    plot_flags(ax,flag,flagpos)

plot_dayboundaries(ax)

if plot_heaters:
    ax_heaters.set_ylim(0,2.5)
    ax_heaters.set_ylabel('Heater power / W')



if plot_options['pressure']:
    key = 'PRESSURE1'
    axpressure = ax.twinx()
    marker = 'v'
    colour = 'gold'
    labeltxt = 'pressure'
    curve = axpressure.plot(date[key],sensor[key],ls='none',marker=marker,color=colour,label=labeltxt)
    curves += curve
    axpressure.set_yscale('log')
    axpressure.set_ylabel('Pressure / mbar',rotation=270,va='bottom',ha='center',color=colour)

if 'log' in plot_options.keys() and plot_options['log']:
    ax.set_yscale('log')
    
labels = [l.get_label() for l in curves]
ax.legend(curves, labels, loc='lower left',facecolor='wheat',framealpha=0.8)
pngname = 'cooldown_%s.png' % end_date.strftime('%Y%m%d')
fig.savefig(pngname,format='png',dpi=300,bbox_inches='tight')
plt.connect('button_press_event',mouse_click_date)
if hostname.find('qubic-central')<0: plt.show()
#plt.close(fig)
#ans = input('enter to exit')
