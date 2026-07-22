#!/usr/bin/env python3
import os,sys
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np

from satorchipy.datefunctions import str2dt, utcnow, utcfromtimestamp
from satorchipy.plotfunctions import mouse_click_date, labelprops, make_legend_label, plot_dayboundaries, plot_flags

from qubicpack.utilities import fmt4latex
from qubicpack.housekeeping.utilities import read_hk_file, find_pt_start, read_hk_labels, read_hk_flags, download_hk
from qubicpack.housekeeping.plot_options import boxprops,plot_options
hk_dir = plot_options['hk_dir']
flag = read_hk_flags(plot_options['events'])

basenames = ['PRESSURE1',
             'weather',
             'inside_weather',
             'CRYOSTAT',
             'TEMPERATURE01',
             'TEMPERATURE03',
             'AVS47_1_ch1',
             'AVS47_1_ch4',
             'AVS47_2_ch3',
             'AVS47_2_ch2',
             'AVS47_2_ch4',
             'TEMPERATURE04',
             'AVS47_2_ch0', 
             'AVS47_2_ch1']

basenames = ['PRESSURE1',
             'weather',
             'inside_weather',
             'CRYOSTAT',
             'AVS47_2_ch1',
             'AVS47_2_ch3']

    

def read_weather_file(filename):
    '''
    read the weather data (outside or inside)
    '''
    weatherdat = ([],[],[])
    h = open('%s/%s' % (hk_dir,filename))
    lines = h.read().split('\n')
    h.close()
    del(lines[-1])
    for line in lines:
        col = line.split()
        if len(col)<3: continue

        datpt = np.empty(3)
        ok = False
        for idx in range(3):
            try:
                datpt[idx] =  (eval(col[idx]))
                ok = True
            except:
                break
        
        if not ok: continue
        for idx in range(3):
            weatherdat[idx].append(datpt[idx])
            
    weatherdat = np.array(weatherdat)        
    return weatherdat

if plot_options['download']:
    download_hk(basenames,plot_options['hk_dir'],remote_machine=plot_options['qubic-central'])

print('[%s] reading pressure' % dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
dat = read_hk_file('%s/PRESSURE1.txt' % hk_dir)

# read the temperature labels
all_labels = read_hk_labels()
label = all_labels['housekeeping']

# read the cryo temperatures, if any
print('[%s] reading temperatures' % dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
full_cryodat = {}
for bname in basenames:
    if (bname.find('TEMPERATURE')<0) and (bname.find('AVS47')<0): continue
    temp_dat = read_hk_file('%s/%s.txt' % (hk_dir,bname))
    full_cryodat[bname] = temp_dat
    

# read the weather
print('[%s] reading weather.txt' % dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
weatherdat = read_weather_file('weather.txt')
print('[%s] reading inside_weather.txt' % dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
inside_weatherdat = read_weather_file('inside_weather.txt')

# read the cryostat shell temperature
print('[%s] reading cryoshell temperature' % dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
cryoshell_dat = read_hk_file('%s/CRYOSTAT.txt' % hk_dir)    
print('[%s] done reading' % dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

# apply shift to clock adjustment
clock_adjustment = False
tstamps = dat[0]
if tstamps[0]==1637083004.810148:
    clock_adjustment = True
    offset_idx1 = 9527
    delta1 = tstamps[offset_idx1+1] - tstamps[offset_idx1]
    shift1 = tstamps[offset_idx1] - tstamps[offset_idx1-1] + 240
    offset_idx2 = 9601
    delta2 = tstamps[offset_idx2+1] - tstamps[offset_idx2]
    shift2 = tstamps[offset_idx2] - tstamps[offset_idx1-2] + 240

    tstamps[0:offset_idx1] += shift2+delta2
    tstamps[offset_idx1:offset_idx2] += shift2+delta2 - (shift1-delta1)
    
    

date_margin = dt.timedelta(minutes=5)
    
start_date = list(flag.keys())[0] - date_margin
if 'tstart' in plot_options.keys() and plot_options['tstart'] is not None:
    start_date = plot_options['tstart']
start_tstamp = start_date.timestamp()

start_idx = np.where(dat[0] > start_tstamp)[0][0]

if clock_adjustment:
    offset_idx1 -= start_idx
    offset_idx2 -= start_idx
date = []
for tstamp in dat[0][start_idx:]:
    date.append(utcfromtimestamp(tstamp))

 
end_date = date[-1] + date_margin
if 'tend' in plot_options.keys():
    end_date = plot_options['tend']

p = dat[1][start_idx:]
pmin = p.min() - 0.5*p.min()
pmax = p.max() 
if 'pmin' in plot_options.keys():
    pmin = plot_options['pmin']
if 'pmax' in plot_options.keys():
    pmax = plot_options['pmax']
    

Tmins = []
Tmaxs = []
RHmins = []
RHmaxs = []
plot_weather = False
if len(weatherdat[0])>0 and (weatherdat[0]>start_tstamp).sum()>0:
    plot_weather = True
    weather_start_idx = np.where(weatherdat[0] > start_tstamp)[0][0]
    Toutside = np.array(weatherdat[1][weather_start_idx:]) + 273.15 # convert to K
    RH = np.array(weatherdat[2][weather_start_idx:])
    Tmins.append(Toutside.min())
    Tmaxs.append(Toutside.max())
    RHmins.append(RH.min())
    RHmaxs.append(RH.max())
    weatherdate = []
    for tstamp in weatherdat[0][weather_start_idx:]:
        weatherdate.append(utcfromtimestamp(tstamp))

plot_inside_weather = False
if len(inside_weatherdat[0])>0 and (inside_weatherdat[0] > start_tstamp).sum()>0:
    plot_inside_weather = True
    inside_weather_start_idx = np.where(inside_weatherdat[0] > start_tstamp)[0][0]
    Tinside = np.array(inside_weatherdat[1][inside_weather_start_idx:]) + 273.15 # convert to K
    Tmins.append(Tinside.min())
    Tmaxs.append(Tinside.max())
    insideRH = np.array(inside_weatherdat[2][inside_weather_start_idx:])
    RHmins.append(insideRH.min())
    RHmaxs.append(insideRH.max())
    inside_weatherdate = []
    for tstamp in inside_weatherdat[0][inside_weather_start_idx:]:
        inside_weatherdate.append(utcfromtimestamp(tstamp))
    

plot_cryoshell = False
if cryoshell_dat[0] is not None and (cryoshell_dat[0] > start_tstamp).sum()>0:
    plot_cryoshell = True
    shell_start_idx = np.where(cryoshell_dat[0] > start_tstamp)[0][0]
    Tshell = np.array(cryoshell_dat[1][shell_start_idx:])
    Tmins.append(Tshell.min())
    Tmaxs.append(Tshell.max())
    shelldate = []
    for tstamp in cryoshell_dat[0][shell_start_idx:]:
        shelldate.append(utcfromtimestamp(tstamp))
    
if len(Tmins)>0:
    Tmin = min(Tmins)
    Tmax = max(Tmaxs)
if 'Tmin' in plot_options.keys():
    Tmin = plot_options['Tmin']
if 'Tmax' in plot_options.keys():
    Tmax = plot_options['Tmax']

if len(RHmins)>0:
    rhmin = min(RHmins)
    rhmax = max(RHmaxs)
if 'rhmin' in plot_options.keys():
    rhmin = plot_options['rhmin']
if 'rhmax' in plot_options.keys():
    rhmax = plot_options['rhmax']


cryodat = {}
cryomins = []
cryomaxs = []
for key in full_cryodat.keys():
    if full_cryodat[key] is not None and (full_cryodat[key][0] is not None) and (full_cryodat[key][0] > start_tstamp).sum()>0:
        cryo_start_idx = np.where(full_cryodat[key][0] > start_tstamp)[0][0]
        tstamp_cryo = full_cryodat[key][0][cryo_start_idx:]
        Tcryo = full_cryodat[key][1][cryo_start_idx:]
        date_list = []
        for tstamp in tstamp_cryo:
            date_list.append(utcfromtimestamp(tstamp))

        cryodat[key] = (date_list,Tcryo)
        cryomins.append(min(Tcryo))
        cryomaxs.append(max(Tcryo))

if len(cryomins)>0: cryomin = min(cryomins)
else: cryomin = 0
if len(cryomaxs)>0: cryomax = max(cryomaxs)
else: cryomax = 300
if 'cryomin' in plot_options.keys():
    cryomin = plot_options['cryomin']
if 'cryomax' in plot_options.keys():
    cryomax = plot_options['cryomax']
    

    
fig = plt.figure()
fig.canvas.manager.set_window_title('plt: pumpdown_%s' % date[-1].strftime('%Y%m%d'))
curves = []
ax = fig.add_axes((0.06, 0.1,0.77, 0.85))
p_str = fmt4latex(p[-1],3)+' mbar'
lbl = make_legend_label('pressure',p_str)
if clock_adjustment:
    ax.plot(date[0:offset_idx1],p[0:offset_idx1],ls='none',marker='+',color='green')
    ax.plot(date[offset_idx1:offset_idx2],p[offset_idx1:offset_idx2],ls='none',marker='x',color='red')
    ax.plot(date[offset_idx2:],p[offset_idx2:],ls='none',marker='d',color='blue')
else:
    curves += ax.plot(date,p,ls='none',marker='d',label=lbl)

pmin_idx = np.argmin(p)
pmin_str = fmt4latex(p.min(),3)
minpressure_txt = 'Minimum pressure %s mbar at %s' % (pmin_str,date[pmin_idx].strftime('%Y-%m-%d %H:%M:%S UT'))

ax.set_ylim(pmin,pmax)
ax.set_yscale('log')
ax.set_xlim(start_date,end_date)
ax.text(0.5,1.01,'Pressure at Alto Chorrillos',fontsize=24,ha='center',va='bottom',transform=ax.transAxes)
ax.set_ylabel('pressure / mbar',fontsize=24,color='#1f77b4ff',)
ax.set_xlabel('date (%s-DD HH:MM) / UT' % date[0].strftime('%Y-%m'))
ax.tick_params(axis='y',labelsize=24,labelcolor='#1f77b4ff')
ax.tick_params(axis='x', rotation=315, left=True, labelleft=True,labelsize=12)

    
pressure_txt = 'pressure at %s is %s' % (date[-1].strftime('%Y-%m-%d %H:%M:%S UT'),p_str)
axlegend = ax

axweather = ax.twinx()
if plot_weather:
    val_str = '%.2f C' % (Toutside[-1] - 273.15)
    lbl = make_legend_label('outside temperature',val_str)
    curves += axweather.plot(weatherdate,Toutside,ls='none',color='green',marker='^',label=lbl)
if plot_inside_weather:
    val_str = '%.2f C' % (Tinside[-1] - 273.15)
    lbl = make_legend_label('inside temperature',val_str)
    curves += axweather.plot(inside_weatherdate,Tinside,ls='none',color='olive',marker='v',label=lbl)
axweather.set_ylim(Tmin,Tmax)
axweather.set_ylabel('Ambient Temperature / K',rotation=270,va='bottom',ha='center',color='green')
axweather.tick_params(axis='y',labelcolor='green')
if plot_cryoshell:
    val_str = '%.2f C' % (Tshell[-1]-273.15)
    lbl = make_legend_label('cryostat shell temperature',val_str)
    curves += axweather.plot(shelldate,Tshell,ls='none',marker='o',color='#a20cffff',label=lbl)
axlegend = axweather
    
axrh = ax.twinx()
if plot_weather:
    val_str = '%.1f%%' % RH[-1]
    lbl = make_legend_label('Outside Relative Humidity',val_str)
    curves += axrh.plot(weatherdate,RH,color='red',ls='none',marker='^',label=lbl)
if plot_inside_weather:
    val_str = '%.1f%%' % insideRH[-1]
    lbl = make_legend_label('Inside Relative Humidity',val_str)
    curves += axrh.plot(inside_weatherdate,insideRH,ls='none',marker='v',color='red',label=lbl)
axrh.set_ylim(rhmin,rhmax)
axrh.set_ylabel('Relative Humidity / %',rotation=270,va='bottom',ha='center',color='red')
axrh.tick_params(axis='y',labelcolor='red',pad=90)
axlegend = axrh

if len(cryodat)>0:
    axcryo = ax.twinx()
    for key in cryodat.keys():    
        val_str = '%.2f K' % cryodat[key][1][-1]
        lbl = make_legend_label(label[key],val_str)
        if label[key]=='1K HS' or label[key]=='PT2 S2 CH':
            marker = 'x'
        else:
            marker = '.'
        curves += axcryo.plot(cryodat[key][0],cryodat[key][1],ls='none',marker=marker,label=lbl)
    axcryo.set_ylim(cryomin,cryomax)
    colour = curves[-1].get_color()
    axcryo.set_ylabel('Cryogenic Temperature / K',rotation=270,va='bottom',ha='center',color='blue')
    axcryo.tick_params(axis='y',labelcolor='blue',pad=160)

    # plot 77K and 273K as indicators
    curves += axcryo.plot([cryodat[key][0][0],cryodat[key][0][-1]],[77,77],ls='dashed',color='grey',label='77K')
    curves += axcryo.plot([cryodat[key][0][0],cryodat[key][0][-1]],[273.15,273.15],ls='dashed',color='lightgrey',label='273K')
    axlegend = axcryo
        
labels = [l.get_label() for l in curves]

minmax = axlegend.axis()[2:]
if 'flagpos' not in plot_options.keys():
    if axlegend.get_yscale()=='log':
        flagpos = np.log10(minmax[0]) + 0.3*( np.log10(minmax[1]) - np.log10(minmax[0]))
        flagpos = 10**flagpos
    else:
        flagpos = minmax[0] + 0.3*( minmax[1] - minmax[0] )
        
else:
    flagpos = plot_options['flagpos']

if 'noflags' not in plot_options.keys():
    plot_flags(axlegend,flag,flagpos)
axlegend.legend(curves, labels, loc='lower left',facecolor='wheat',framealpha=0.9)
axlegend.text(0.18,0.10,minpressure_txt,ha='left',va='bottom',transform=ax.transAxes,fontsize=24,color='black',bbox=boxprops)
axlegend.text(0.18,0.02,pressure_txt,   ha='left',va='bottom',transform=ax.transAxes,fontsize=24,color='black',bbox=boxprops)

plot_dayboundaries(axlegend)

pngname = 'pumpdown_%s.png' % date[-1].strftime('%Y%m%d')
fig.savefig(pngname,format='png',dpi=300,bbox_inches='tight')
plt.connect('button_press_event',mouse_click_date)
plt.show()
#ans = input('enter to exit')
