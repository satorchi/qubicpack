#!/usr/bin/env python3
'''$Id: hwp_analysis.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Fri 16 Oct 2020 11:59:51 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

analysis of Half Wave Plate data.  
This started as a translation from the Jupyter notebooks by JCH: 
   HWP-2019-12-23.Rmd
   Analyse-HWP-Results.Rmd


# Location of the Data (notes by JCH)
### 2019-12-23: 4 (2019-12-23_19.04.56)
- Modulator: Amplitude = 2V. Offset = 1.5 V
- Nice data although source data does not have the same shape as TES
  data, probably the source measurement was not configured
  correctly. The data can however be exploited using Simulated CalSrc.

### 2019-12-24: 0 (2019-12-24_09.46.19)
- Modulator: Amplitude = 500 mV ; Offest = 250 mV
- Nice data with 180 sec/pos and 3 cycles.
- SrcData not there...
- Can be used with Simulated Cal Src
- This is the data used for Figure 11 of Paper-3

### 2019-12-26: 0
- Modulator: Amplitude: 500mV, Offset 2. V
- Only one cycle but good quality data, The source is ON and seems weell configured

### 2019-12-26: 1
- Modulator: Amplitude: 500mV, Offset 2.5 V
- Long overnight acquisition - to be looked at closely

### 2019-12-27: 2
- Modulator: Amplitude = 500mV ; Offest = 2.5 V
- Excellent data


***Important Remark:
The mount has moved on Dec. 26th 18h30 CET, so this means that the
data: 
  [2019-12-23_Dataset_4, 
   2019-12-24_Dataset_0,
   2019-12-26_Dataset_0] 
have the same pointing while
   [2019-12-26_Dataset_1, 
    2019-12-27_Dataset_2] 
are with another pointing.***

### 2020-10/11 (note by ST)
- We used the SCAN function built into the HWP controller
  it goes from position 2 to 6, and there's a pause, and then back
  we couldn't go to the end positions for fear of breaking the thing

'''
import scipy.ndimage.filters as filters
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt

# need 42 colours (of course)
colour = ['red',
          'blue',
          'green',
          'black',
          'magenta',
          'olive',
          'orange',
          'c',
          'm',
          'y',
          'darkcyan',
          'firebrick',
          'yellowgreen',
          'darkorchid',
          'dimgray',
          'steelblue',
          'purple']


def assign_hwp_chunkinfo(self,TES=None,asic=None):
    '''
    assign the interval and position id's for specific datasets.  
    This was done manually
    '''
    args = self.args_ok(TES,asic)
    if args is None:return
    TES,asic = args

    dset = self.dataset_name
    print(dset)
    retval = {}
    retval['dataset'] = dset

    if dset=='2019-12-23_19.04.56__HWP_Scanning_Vtes_2TimePerPos_60':
        ### intervals for the data from 2019-12-23_19.04.56__HWP_Scanning_Vtes_2TimePerPos_60
        intervals = np.array([(70,10851),
                              (11475,20225),
                              (20865,29601),
                              (30228,38977),
                              (39447,48353),
                              (48979,57729),
                              (58510,67104),
                              (69136,78200),
                              (78979,87728),
                              (88359,96945),
                              (97574,106325),
                              (107107,115701),
                              (116483,125391),
                              (125858,134453),
                              (136952,145859),
                              (146326,155077),
                              (156330,164458),
                              (165077,173986),
                              (174455,183363),
                              (183674,192582),
                              (193207,201958)
                              ])
        positions = np.empty((intervals.shape[0]),dtype=np.int)
        posidx = 1
        for idx in range(intervals.shape[0]):
            positions[idx] = posidx
            posidx += 1
            if posidx>7: posidx = 1
            
        timeaxis = self.timeaxis(datatype='sci',asic=asic)
        timeline = self.timeline(asic=1,TES=TES)
        t1 = timeaxis[intervals[0][0]]
        t2 = timeaxis[intervals[-1][1]]
        chunkinfo = chunkify_hwp_scan(timeaxis,timeline,t1,t2,padding=0,intervals=intervals)
        chunkinfo['dataset'] = self.dataset_name
        chunkinfo['TES'] = TES
        chunkinfo['asic'] = asic
        chunkinfo['positions'] = positions
        return chunkinfo
    
    if dset=='2019-12-24_09.46.19__HWP_Scanning_Vtes_2TimePerPos_180':
        intervals = np.array([(2047,29707),
                              (30483,57993),
                              (58937,86028),
                              (86898,114237),
                              (114867,142373),
                              (142994,170495),
                              (171124,198620),
                              (200652,228313),
                              (229257,256447),
                              (257413,284508),
                              (285498,312846),
                              (313627,340973),
                              (341600,369106),
                              (369725,397225),
                              (399260,427076),
                              (427695,455197),
                              (455918,482541),
                              (484265,511461),
                              (512388,539582),
                              (540358,567708),
                              (568650,582561)
                              ])
        positions = np.empty((intervals.shape[0]),dtype=np.int)
        posidx = 1
        for idx in range(intervals.shape[0]):
            positions[idx] = posidx
            posidx += 1
            if posidx>7: posidx = 1
            
        timeaxis = self.timeaxis(datatype='sci',asic=asic)
        timeline = self.timeline(asic=1,TES=TES)
        t1 = timeaxis[intervals[0][0]]
        t2 = timeaxis[intervals[-1][1]]
        chunkinfo = chunkify_hwp_scan(timeaxis,timeline,t1,t2,padding=0,intervals=intervals)
        chunkinfo['dataset'] = self.dataset_name
        chunkinfo['TES'] = TES
        chunkinfo['asic'] = asic
        chunkinfo['positions'] = positions
        return chunkinfo
        


    if dset=='2019-12-26_16.04.14__HWP_Scanning_Vtes_2TimePerPos_60':
        intervals = np.array([(67,10849),
                              (11474,20221),
                              (21766,29682),
                              (30226,38976),
                              (39446,48354),
                              (48978,57888),
                              (58354,67104)
                              ])
        positions = np.empty((intervals.shape[0]),dtype=np.int)
        posidx = 1
        for idx in range(intervals.shape[0]):
            positions[idx] = posidx
            posidx += 1
            if posidx>7: posidx = 1

        timeaxis = self.timeaxis(datatype='sci',asic=asic)
        timeline = self.timeline(asic=1,TES=TES)
        t1 = timeaxis[intervals[0][0]]
        t2 = timeaxis[intervals[-1][1]]
        chunkinfo = chunkify_hwp_scan(timeaxis,timeline,t1,t2,padding=0,intervals=intervals)
        chunkinfo['dataset'] = self.dataset_name
        chunkinfo['TES'] = TES
        chunkinfo['asic'] = asic
        chunkinfo['positions'] = positions
        return chunkinfo
    
    if dset=='2019-12-27_11.51.37__HWP_Scanning_Vtes_2_TimePerPos_60_Cycle_0_over_5':
        ### intervals for data from 2019-12-27_11.51.37__HWP_Scanning_Vtes_2_TimePerPos_60_Cycle_0_over_5
        ### NOTE:  these are timestamps, but the scripts expect array index.
        ###        Use np.where() to find the appropriate index
        interval_start = np.array([1577447496.9374256,
                                   1577447560.0006204,
                                   1577447620.819865,
                                   1577447680.1653929,
                                   1577447740.180023,
                                   1577447800.083986,
                                   1577447859.9818566])

        interval_end = np.array([1577447556.9966493,
                                 1577447616.1648827,
                                 1577447676.9783556,
                                 1577447737.062811,
                                 1577447795.9415889,
                                 1577447857.9210281,
                                 1577447917.1985645])

        intervals = np.empty((len(interval_start),2),dtype=np.int)
        positions = np.empty((intervals.shape[0]),dtype=np.int)
        posidx = 1
        timeaxis = self.timeaxis(datatype='sci',asic=asic)
        timeline = self.timeline(asic=1,TES=TES)
        for idx,tstamp_start in enumerate(interval_start):
            idxrange = np.where(timeaxis>tstamp_start)[0]
            intervals[idx][0] = idxrange[0]

            tstamp_end = interval_end[idx]
            idxrange = np.where(timeaxis>tstamp_end)[0]
            intervals[idx][1] = idxrange[0]

            positions[idx] = posidx
            posidx += 1
            if posidx>7: posidx = 1
        
        t1 = timeaxis[intervals[0][0]]
        t2 = timeaxis[intervals[-1][1]]
        chunkinfo = chunkify_hwp_scan(timeaxis,timeline,t1,t2,padding=0,intervals=intervals)
        chunkinfo['dataset'] = self.dataset_name
        chunkinfo['TES'] = TES
        chunkinfo['asic'] = asic
        chunkinfo['positions'] = positions
        return chunkinfo



    if dset=='2020-11-05_09.21.43__HWP_Scanning_150GHz_az0.07_el49.99':
        ### positions for the data from 2020-11-05_09.21.43__HWP_Scanning_150GHz_az0.07_el49.99
        ### intervals are derived in chunkify_hwp_scan() with period=10.64,padding=3.0
        ###
        t1 = 1604568170.1192396
        t2 = 1604568731.1094625
        positions = np.array([ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,
                               18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
                               35, 36, 37, 38, 39, 40, 41, 42, 42, 42, 42, 42, 41, 40, 39, 38, 37,
                               36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20,
                               19, 18, 17, 16, 15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,
                               2,  1,  1,  1,  1,  1])
        timeaxis = self.timeaxis(datatype='sci',asic=1)
        timeline = self.timeline(asic=1,TES=95)
        chunkinfo = chunkify_hwp_scan(timeaxis,timeline,t1,t2,period=10.64,padding=3.0)
        chunkinfo['dataset'] = self.dataset_name
        chunkinfo['TES'] = 95
        chunkinfo['asic'] = 1
        chunkinfo['positions'] = positions
        return chunkinfo

    # if it's not a known dataset, then try to chunkify by reading the HWP position
    return self.chunkify_by_position(TES,asic)
                 
def chunkify_by_position(self,TES=None,asic=None,padding=1):
    '''
    find the HWP chunk information using the HWP position
    the intervals are cropped by padding at both ends
    '''

    args = self.args_ok(TES,asic)
    if args is None:return
    TES,asic = args

    chunkinfo = {}

    hwp = self.hwp_position()
    t_hwp = self.timeaxis(datatype='hwp')
    tline = self.timeline(TES=TES,asic=asic)
    t_sci = self.timeaxis(datatype='sci',asic=asic)

    chunkinfo['dataset'] = self.dataset_name
    chunkinfo['TES'] = TES
    chunkinfo['asic'] = asic

    known_position_idxs = np.where(hwp != 0)[0] # indexes where HWP is in a known position
    if len(known_position_idxs)==0:
        print('ERROR!  No known HWP positions.')
        return None
    
    positions = []
    tstart_list = []
    tend_list = []
    npos = 0
    curpos = 0
    for idx,posidx in enumerate(known_position_idxs):
        pos = hwp[posidx]
        if pos!=curpos:
            positions.append(pos)
            tstart_list.append(t_hwp[posidx])
            if npos>0:
                iend = known_position_idxs[idx-1]
                tend_list.append(t_hwp[iend])
            curpos = pos
            npos += 1
    tend_list.append(t_hwp[known_position_idxs[-1]])
    positions = np.array(positions)
    npositions = len(positions)
    
    #hwp_time_intervals = np.empty((npositions,2),dtype=float)
    #intervals = np.empty((npositions,2),dtype=int)
    hwp_time_intervals = []
    intervals = []
    good_positions = []
    for idx in range(npositions):
        tstart = tstart_list[idx]
        tend = tend_list[idx]
        if tend<=tstart+2*padding:
            print('skipping interval with separation %.3f seconds' % (tend-tstart))
            continue

        good_positions.append(positions[idx])
        hwp_time_intervals.append(np.array([tstart,tend],dtype=float))

        tstart_sci = tstart + padding
        tend_sci = tend - padding
        istart = np.where(t_sci>=tstart_sci)[0][0]
        iend = np.where(t_sci>=tend_sci)[0][0]
        intervals.append(np.array([istart,iend],dtype=int))
    
    chunkinfo['positions'] = np.array(good_positions)
    chunkinfo['npositions'] = len(good_positions)
    chunkinfo['hwp time intervals'] = np.array(hwp_time_intervals)
    chunkinfo['intervals'] = np.array(intervals)

    additional_chunkinfo = chunkify_hwp_scan(t_sci,tline,
                                             padding=padding,
                                             intervals=intervals,
                                             positions=positions)

    for key in additional_chunkinfo.keys():
        if key not in chunkinfo.keys():
            chunkinfo[key] = additional_chunkinfo[key]
    
    return chunkinfo

def chunkify_search(timeline,tol=1e-3,jitter_tol=1e-2):
    '''
    chunkify the HWP scan by looking for changes in the min peak-to-peak

    THIS METHOD DOESN'T WORK YET
    '''
    retval = {}
    interval_list = []
    going_up_list = []
    going_down_list = []
    min_list = []
    
    val = timeline[0]
    prev = val
    min_curr = val
    idx_interval_start = -1
    idx_interval_end = 0
    going_down = False
    for idx,val in enumerate(timeline):
        
        if val>prev:
            if prev==0.0:
                norm = 1.0
            else:
                norm = np.abs(prev)
            jitter_compare = (val-prev)/norm
            if jitter_compare>jitter_tol and going_down:
                going_up_list.append(idx)
                min_list.append(idx-1)
                going_down = False
                min_compare = (min_curr-prev)/norm
                if prev<min_curr and min_compare>tol:
                    if idx_interval_start==-1:
                        idx_interval_start = idx - 1
                    else:
                        idx_interval_end = idx - 1
                        interval_list.append((idx_interval_start,idx_interval_end))
                        idx_interval_start = -1
                        min_curr = prev
        else:
            going_down = True
            going_down_list.append(idx)
        prev = val

    retval['intervals'] = np.array(interval_list)
    retval['going up'] = np.array(going_up_list)
    retval['going down'] = np.array(going_down_list)
    retval['minimums'] = np.array(min_list)
    return retval

def chunkify_hwp_scan(timeaxis,
                      timeline,
                      t1=None,
                      t2=None,
                      period=None,
                      padding=1,
                      intervals=None,
                      positions=None):
    '''
    The SCAN function went from position 2 to 6 (starts at time t1)
    and then back (starts at time t2)
    it waits 10 seconds and then moves on (period=10)
    change the period to take into account the time to move (add a fraction of a second)
    crop the interval at the beginning and the end by the amount given (padding in seconds)
    if the positions are already known, use them for the colour index
    '''
    retval = {} # return stuff, including the input arguments
    retval['period'] = period
    retval['padding'] = padding
    retval['t1'] = t1
    retval['t2'] = t2
    retval['timeaxis'] = timeaxis
    
    
    p2p_forward_list = []
    p2p_backward_list = []
    p2p_list = p2p_forward_list
    interval_list = []
    position_list = []

    chunk_list = []

    
    t_end = timeaxis[-1]
    idx_mid = len(timeaxis) // 2
    if intervals is None:
        t_end1 = timeaxis[idx_mid] - 20
    else:
        t_end1 = t_end

    if t1 is None: t1 = timeaxis[0]
    if t2 is None: t2 = timeaxis[-1]
    t = t1
    interval_idx = 0
    minvelope_t = []
    minvelope_a = []
    maxvelope_a = []

    position_increment = 1
    position_counter = 1
    while t < t_end:
        if t>t_end1 and t<t2:
            p2p_list = p2p_backward_list
            t = t2
            position_increment = -1

        if intervals is None:
            t_start = t + padding
            t_stop = t + period - padding
            
            idx_range = np.where(timeaxis>t_start)[0]
            if len(idx_range)==0: break
            idx_start = idx_range[0]
        
            idx_range = np.where(timeaxis>t_stop)[0]
            if len(idx_range)==0: break
            idx_stop = idx_range[0]
        else:
            idx_start = intervals[position_counter-1][0]
            idx_stop = intervals[position_counter-1][1]

        chunk = timeline[idx_start:idx_stop]
        if len(chunk)==0:
            print('skipping zero length interval: idx_start,idx_stop = (%i,%i)' % (idx_start,idx_stop))
            t = max((timeaxis[idx_start],timeaxis[idx_stop]))
            position_counter += position_increment
            continue

        interval_list.append((idx_start,idx_stop))
        if positions is None:
            position_list.append(position_counter)
        else:
            position_list.append(positions[position_counter-1])
        
        p2p = np.nanmax(chunk) - np.nanmin(chunk)
        p2p_list.append(p2p)
        
        offset = 0.5*( np.nanmax(chunk) + np.nanmin(chunk) )
        reduced = chunk - offset
        chunk_t = timeaxis[idx_start:idx_stop]

        minvelope_a.append(reduced.min())
        maxvelope_a.append(reduced.max())
        idx_min = reduced.argmin()
        minvelope_t.append(chunk_t[idx_min])

        chunk_list.append(np.array((chunk_t,reduced)))

        interval_idx += 1
        position_counter += position_increment
        if period is None:
            t = chunk_t[-1]
        else:
            t += period
        if intervals is not None and position_counter>len(intervals): break
        
        

    retval['p2p forward'] = np.array(p2p_forward_list)
    retval['p2p backward'] = np.array(p2p_backward_list)
    retval['intervals'] = np.array(interval_list)
    retval['positions'] = np.array(position_list)
    retval['npositions'] = len(retval['positions'])
    retval['minimums'] = np.array(minvelope_a)
    retval['maximums'] = np.array(maxvelope_a)
    retval['chunks'] = chunk_list
    
    return retval


##################################################################
#### plot methods                                            #####
##################################################################
def plot_interval_borders(timeaxis,intervals,p2p=None,positions=None,ax=None):
    '''
    show the intervals with dotted lines
    '''
    if ax is None: ax = plt.gca()
    if positions is None: positions = np.arange(len(intervals)) + 1
    if p2p is None: p2p = ax.axis()[2:]
    for idx,interval in enumerate(intervals):
        interval_txt = '%i' % positions[idx]
        
        interval_start = timeaxis[interval[0]]
        interval_end = timeaxis[interval[1]]
        
        cidx = positions[idx] - 1
        while cidx>=len(colour):
            cidx -= len(colour)

        ax.plot([interval_start,interval_start],p2p,color=colour[cidx],ls='dotted')
        ax.plot([interval_end,interval_end],p2p,color=colour[cidx],ls='dotted')

        interval_mid = 0.5*(interval_end + interval_start)
        ax.text(interval_mid,p2p[1],interval_txt,va='bottom',ha='center',fontsize=14)

    return

def plot_hwp_scan(chunkinfo,plot_intervals=True,positions=None):
    '''
    plot the HWP scan chunks.
    the chunkinfo is a dictionary generated by chunkify_hwp_scan()
    '''
    retval = {}
    retval['pngname'] = '%s_ASIC%02i_TES%03i_HWP-chunks.png' % (chunkinfo['dataset'],chunkinfo['asic'],chunkinfo['TES'])
    
    t0 = chunkinfo['chunks'][0][0][0]
    t1 = chunkinfo['t1']
    t2 = chunkinfo['t2']
    t_end = chunkinfo['chunks'][-1][0][-1]
    period = chunkinfo['period']
    padding = chunkinfo['padding']

    info_list = [chunkinfo['dataset']]
    info_list.append('TES%3i, asic%2i' % (chunkinfo['TES'],chunkinfo['asic']))
    if t0>1e9:
        timelabel = 'UTC Date / seconds since 1970-01-01'
        if t1 is not None: info_list.append('t1=%s' % dt.datetime.utcfromtimestamp(t1).strftime('%Y-%m-%d %H:%M:%S.%f'))
        if t2 is not None: info_list.append('t2=%s' % dt.datetime.utcfromtimestamp(t2).strftime('%Y-%m-%d %H:%M:%S.%f'))
    else:
        timelabel = 'time / seconds'
        if t1 is not None: info_list.append('t1=%.2fs' % t1)
        if t2 is not None: info_list.append('t2=%.2fs' % t2)

    if period is not None: info_list.append('period=%.3fs' % period)
    if padding is not None: info_list.append('padding=%.3fs' % padding)
    txt = ', '.join(info_list)
    

    fig = plt.figure()
    retval['fig'] = fig
    fig.canvas.manager.set_window_title('plt: HWP scan')
    ax = fig.add_axes((0.05,0.1,0.93,0.83))
    ax.tick_params(axis='x',labelsize=20)
    ax.tick_params(axis='y',labelsize=20)
    ax.set_xlabel(timelabel,fontsize=20)
    ax.set_ylabel('Detector response / Arbitrary units',fontsize=20)

    # show the start
    p2p_all = chunkinfo['maximums'].max() - chunkinfo['minimums'].min()
    p2p_minmax = [chunkinfo['minimums'].min(),chunkinfo['maximums'].max()+0.05*p2p_all]
    if t1 is not None:
        ax.plot([t1,t1],p2p_minmax,color='black',ls='dashed')
        ax.set_xlim(t1,t_end)
    else:
        ax.set_xlim(t0,t_end)
    if t2 is not None:
        ax.plot([t2,t2],p2p_minmax,color='black',ls='dashed')
        
    ax.set_ylim(p2p_minmax)

    # show the intervals
    if plot_intervals:
        if positions is None:
            positions2plot = chunkinfo['positions']
        else:
            positions2plot = positions
        plot_interval_borders(chunkinfo['timeaxis'],chunkinfo['intervals'],p2p=p2p_minmax,positions=positions2plot)
    
    # show the envelope of the minimum
    #ax.plot(minvelope_t,minvelope_a,color='green')

    # some info
    ax.text(0.5,1.05,txt,va='bottom',ha='center',fontsize=18,transform=ax.transAxes)

    # plot the chunks
    if positions is None: positions = chunkinfo['positions']
    for interval_idx,chunk in enumerate(chunkinfo['chunks']):
        chunk_t = chunk[0]
        chunk_a = chunk[1]
        
        cidx = positions[interval_idx] - 1
        while cidx>=len(colour): cidx -= len(colour)
        ax.plot(chunk_t,chunk_a,color=colour[cidx])

    fig.savefig(retval['pngname'],format='png',dpi=100,bbox_inches='tight')
    return retval

