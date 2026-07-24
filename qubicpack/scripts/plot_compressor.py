#!/usr/bin/env python3
'''
$Id: plot_compressor.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Wed 01 Sep 2021 13:39:05 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

plot compressor log
'''
import os
import datetime as dt
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from satorchipy.datefunctions import str2dt, tstamp2dt, utcnow, utcfromtimestamp
from satorchipy.plotfunctions import labelprops, plot_flags, plot_dayboundaries

from qubicpack.utilities import hostname
from qubicpack.housekeeping.utilities import read_compressor_log, read_hk_flags, download_hk
from qubicpack.housekeeping.plot_options import boxprops,plot_options
flag = read_hk_flags(plot_options['events'])

def line_model(x,m,b):
    return m*x + b


def operational_percentage(timestamps,start_date,end_date,offmargin=180):
    '''
    calculate the percentage time that the compressor was operating
    '''
    maskrange = ( (timestamps>=start_date.timestamp()) & (timestamps<=end_date.timestamp()) )
    tstamp_start = timestamps[maskrange][0]
    tstamp_end = timestamps[maskrange][-1]
    tsep = (timestamps[maskrange] - np.roll(timestamps[maskrange],1))[1:]
    maskoff = (tsep>offmargin)
    offtime = tsep[maskoff].sum()
    total_time = tstamp_end - tstamp_start
    ontime = total_time - offtime
    time_range = end_date.timestamp() - start_date.timestamp()
    on_percentage = 100*ontime/time_range
    return on_percentage,tstamp_start,tstamp_end
    

# main program
if __name__=='__main__':
    date_fmt = '%Y-%m-%d %H:%M:%S'
    hk_dir = plot_options['hk_dir']

    if plot_options['download']:
        basenames = ['compressor*']
        download_hk(basenames,hk_dir,remote_machine=plot_options['qubic-central'])

    clog = {}
    clog[1] = read_compressor_log('%s/compressor1_log.txt' % hk_dir)
    clog[2] = read_compressor_log('%s/compressor2_log.txt' % hk_dir)

    end_date = max( [clog[1]['date'][-1], clog[2]['date'][-1]] )
    start_date = min( [clog[1]['date'][0], clog[2]['date'][0]] )
    if 'tstart' in plot_options.keys() and plot_options['tstart'] is not None:
        start_date = plot_options['tstart']
    if 'tend' in plot_options.keys() and plot_options['tend'] is not None:
        end_date = plot_options['tend']        
    
    fig = plt.figure()
    fig.canvas.manager.set_window_title('plt: compressors_%s' % end_date.strftime('%Y%m%d'))
    tstamp_start_list = []
    tstamp_end_list = []
    op_percentage_list = []
    for compressornum in [1,2]:
        timestamps = clog[compressornum]['timestamp']
        op_percentage,tstamp_start,tstamp_end = operational_percentage(timestamps,start_date,end_date)
        tstamp_start_list.append(tstamp_start)
        tstamp_end_list.append(tstamp_end)
        start_date_str = utcfromtimestamp(tstamp_start).strftime(date_fmt)
        end_date_str = utcfromtimestamp(tstamp_end).strftime(date_fmt)
        op_percentage_txt = 'Compressor %i: Operating percentage = %4.1f %%'\
            % (compressornum,op_percentage)
        op_percentage_list.append(op_percentage_txt)

    #start_date = utcfromtimestamp(min(tstamp_start_list))
    #end_date = utcfromtimestamp(max(tstamp_end_list))
    ttl_list = ['Pulse tube compressor data from %s to %s' % (start_date.strftime(date_fmt),end_date.strftime(date_fmt))]\
        + op_percentage_list
    ttl = '\n'.join(ttl_list)
    fig.suptitle(ttl)
    print(ttl)

    key = 'Pin'
    ax = fig.add_axes((0.09, 0.64,0.9, 0.27))
    ax.set_xlim(start_date,end_date)
    min_list = []
    max_list = []
    for compressornum in [1,2]:
        curve = ax.plot(clog[compressornum]['date'],clog[compressornum][key],ls='none',marker='.',label='Compressor %i' % compressornum)
        colour = curve[0].get_color() 
        min_list.append(min(clog[compressornum][key]))
        max_list.append(max(clog[compressornum][key]))

        if 'deltafit' in plot_options.keys() and plot_options['deltafit'] is not None:
            T_now = clog[compressornum]['date'][-1]
            T_fitstart = T_now - dt.timedelta(seconds=plot_options['deltafit'])
            mask = np.array(clog[compressornum]['date'])>=T_fitstart
            fit_npts = mask.sum()
            if fit_npts <=0:
                fit_npts = 60
            print('using the last %i points for fitting' % fit_npts)
            ypts = np.array(clog[compressornum][key][-fit_npts:])
            Pin_mean = ypts.mean()
            fit_npts = len(ypts)
            xpts = np.empty(fit_npts)
            for idx,d in enumerate( clog[compressornum]['date'][-fit_npts:] ):
                xpts[idx] = d.timestamp()
            slope = 0.0 # first guess
            offset = -slope*xpts[0] + ypts[0]
            first_guess = np.array((slope,offset))    
            popt,pcov = curve_fit(line_model,xpts,ypts,p0=first_guess)
            slope = popt[0]
            offset = popt[1]
            print('%s compressor-%i %.3f bar, slope: %e bar/sec' % (key,compressornum,Pin_mean,slope))
            model_ypts = line_model(xpts,slope,offset)
            lbl = 'Compressor %s: %.2f bar, slope = %.1f mbar/hr' % (compressornum,Pin_mean,3600*1000*slope)
            ax.plot(clog[compressornum]['date'][-fit_npts:],model_ypts,color=colour,label=lbl,lw=5)
            
                  

            
        
    ax.set_ylabel('Return Pressure / bar')
    ax.tick_params(axis='both',labelbottom=False,labelsize=24)
    ax.legend(loc='upper left')
    P_min = min(min_list)
    P_max = max(max_list)
    if 'Pin_min' in plot_options.keys():
        P_min = plot_options['Pin_min']
        print('option Pin_min=%.1f' % plot_options['Pin_min'])
    if 'Pin_max' in plot_options.keys():
        P_max = plot_options['Pin_max']
    ax.set_ylim(P_min,P_max)
    plot_dayboundaries(ax)
    

    ax = fig.add_axes((0.09, 0.36, 0.9, 0.27))
    ax.set_xlim(start_date,end_date)
    min_list = []
    max_list = []
    for compressornum in [1,2]:
        for key in ['Tin','Tout']:
            ax.plot(clog[compressornum]['date'],clog[compressornum][key],ls='none',marker='.',label='%s compressor %i' % (key,compressornum))
            min_list.append(min(clog[compressornum][key]))
            max_list.append(max(clog[compressornum][key]))
    ax.tick_params(axis='both',labelbottom=False,labelsize=24)
    ax.set_ylabel('temperature / C')
    ax.legend(loc='upper left')
    T_min = min(min_list)
    T_max = max(max_list)
    if 'Tinout_min' in plot_options.keys():
        T_min = plot_options['Tinout_min']
    if 'Tinout_max' in plot_options.keys():
        T_max = plot_options['Tinout_max']
    ax.set_ylim(T_min,T_max)
    plot_dayboundaries(ax)
    

    key = 'T_He'
    ax = fig.add_axes((0.09, 0.08, 0.9, 0.27))
    ax.set_xlim(start_date,end_date)
    min_list = []
    max_list = []
    for compressornum in [1,2]:
        ax.plot(clog[compressornum]['date'],clog[compressornum][key],ls='none',marker='.',label='%s compressor %i' % (key,compressornum))
        min_list.append(min(clog[compressornum][key]))
        max_list.append(max(clog[compressornum][key]))
    ax.tick_params(axis='both',labelsize=24)
    ax.set_xlabel('date / UT')
    ax.set_ylabel('temperature / C')
    ax.legend(loc='upper left')
    he_min = min(min_list)
    he_max = max(max_list)
    
    if 'T_He_min' in plot_options.keys():
        he_min = plot_options['T_He_min']
    if 'T_He_max' in plot_options.keys():
        he_max = plot_options['T_He_max']
    ax.set_ylim(he_min,he_max)
    plot_dayboundaries(ax)

    if 'noflags' not in plot_options.keys():
        plot_flags(ax,flag)
    else:
        print('not showing flags: noflags=%s' % plot_options['noflags'])
        
    
    pngname = 'compressors_%s.png' % end_date.strftime('%Y%m%d')
    fig.savefig(pngname,format='png',dpi=300,bbox_inches='tight')
    if hostname.find('qubic-central')<0: plt.show()
