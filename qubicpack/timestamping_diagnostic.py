#!/usr/bin/env python
'''
$Id: timestamping_diagnostic.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Wed 27 Feb 2019 15:39:15 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

make a diagnostic plot of the derived timestamps

'''
from __future__ import division, print_function
from matplotlib import pyplot as plt
import numpy as np


'''
# this script was originally tested on the following data:
a2=qp()
a2.read_qubicstudio_dataset('2019/2019-02-22/2019-02-22_16.28.01__Scan2/',asic=2)
pps = a2.hk['ASIC_SUMS']['PPS']
gps = a2.hk['ASIC_SUMS']['GPSDate']
tstamps    = a2.timeline_timeaxis(axistype='pps')
indextime  = a2.timeline_timeaxis(axistype='index')
compstamps = a2.timeline_timeaxis(axistype='computertime')
'''

def plot_timestamp_diagnostic(self,hk=None,zoomx=None,zoomy=None):
    '''
    make a diagnostic plot of the derived timestamps
    '''
    hk = self.qubicstudio_filetype_truename(hk)
    if hk is None: hk = 'ASIC_SUMS'
    
    if not hk in self.hk.keys():
        print('Please give a valid HK.  Valid names are: %s' % ', '.join(self.hk.keys()))
        return None
              
    if hk=='ASIC_SUMS':
        pps_title = 'PPS Scientific Data'
        tstamps_title = 'Timestamps diagnostic for Scientific Data'
    elif hk=='INTERN_HK':
        pps_title = 'PPS Platform'
        tstamps_title = 'Timestamps diagnostic for Platform Data'
    else:
        pps_title = 'PPS %s' % hk
        tstamps_title = 'Timestamps diagnostic for %s' % hk
        
    
    pps = self.pps(hk=hk)
    if pps is None: return 
    gps = self.gps(hk=hk)
    if gps is None: return
    
    compstamps  = self.hk[hk]['ComputerDate']
    npts = len(pps)
    if hk=='ASIC_SUMS':
        sample_period = self.sample_period()
    else:
        sample_period = float(compstamps.max() - compstamps.min())/len(compstamps)

    datainfo = self.infotext()
    
    ########## find where we have potential problems
    epsilon = 0.1
    separations = []
    separations_idx = []
    pps_high = np.where(pps==1)[0]
    pps_indexes = []
    prev = gps[0]
    for idx in pps_high:
        if (idx>0 and pps[idx-1]==0) or (idx<npts-1 and pps[idx+1]==0):
            sep = gps[idx] - prev
            if sep <> 0:
                pps_indexes.append(idx)
                separations.append(sep)
                separations_idx.append(idx)
                prev = gps[idx]            

    separations = np.array(separations[1:])
    separations_idx = np.array(separations_idx[1:])
    print('number of separations: %i' % len(separations))
    print('number of samples: %i' % len(pps))

    mean_separation = separations.mean()
    print('mean separation between pulses is %.4f second' % mean_separation)
    max_separation = separations.max()
    print('max separation between pulses is %.4f second' % max_separation)
    min_separation = separations.min()
    print('min separation between pulses is %.4f second' % min_separation)

    jump_indexes = np.where(separations>1+epsilon)[0]
    print('there are %i jumps at:  %s' % (len(jump_indexes),separations_idx[jump_indexes]))

    stick_indexes = np.where(separations<1-epsilon)[0]
    print('there are %i sticks at: %s' % (len(jump_indexes),separations_idx[stick_indexes]))

    tstamps = self.pps2date(pps,gps)
    t0 = tstamps[0]

    # do a line fit to get drift of the sampling clock
    xpts = np.arange(len(tstamps))
    linefit = np.polyfit(xpts,tstamps,1,full=True)
    slope = linefit[0][0]
    offset = linefit[0][1]
    derived_sample_period = slope
    sample_period_txt = 'expected sample period: %.6f msec\nderived sample period: %.6f msec' % (sample_period*1000,slope*1000)
    print(sample_period_txt)

    # subtract the progression to view only the residues (horizontal line)
    indextime = slope*xpts + offset

    # mark problems with a vertical line
    yminmax = (compstamps.min(),compstamps.max())

    plt.ion()
    ##### plot the scientific PPS
    ttl = pps_title
    png_rootname = '%s_%s' % (ttl.lower().replace(' ','_'),self.obsdate.strftime('%Y%m%d-%H%M%S'))
    fig0 = plt.figure(figsize=(16,8))
    fig0.canvas.set_window_title('plt: %s' % ttl)
    fig0.text(0.9,0.9,sample_period_txt,ha='right')
    plt.suptitle(ttl)
    plt.title(datainfo)
    plt.plot(pps)
    ax0 = fig0.axes[0]
    ax0.set_ylabel('PPS Level')
    ax0.set_xlabel('sample number')
    pngname = '%s_full.png' % png_rootname
    plt.savefig(pngname,format='png',dpi=100,bbox_inches='tight')

    if zoomx is not None:
        ax0.set_xlim(zoomx)
        pngname = '%s_zoom.png' % png_rootname
        plt.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    
    
    ##### plot diagnostic #1
    ttl = tstamps_title
    png_rootname = '%s_%s' % (ttl.lower().replace(' ','_'),self.obsdate.strftime('%Y%m%d-%H%M%S'))
    fig1 = plt.figure(figsize=(16,8))
    fig1.canvas.set_window_title('plt: %s' % ttl)
    plt.suptitle(ttl)
    plt.title(datainfo)
    fig1.text(0.9,0.9,sample_period_txt,ha='right')
    plt.plot(indextime,                       ls='none',marker='d',label='index time')
    plt.plot(tstamps,                         ls='none',marker='o',label='derived timestamps')
    plt.plot(compstamps,                      ls='none',marker='*',label='computer time')
    plt.plot(gps,                             ls='none',marker='x',label='GPS date')
    plt.plot(pps_high,tstamps[pps_high],      ls='none',marker='^',label='pps high')
    plt.plot(pps_indexes,tstamps[pps_indexes],ls='none',marker='o',label='pps indexes',markerfacecolor='none',markersize=16,color='black')
    for idx in jump_indexes:
        plt.plot((separations_idx[idx],separations_idx[idx]),yminmax,color='red',ls='dashed')
    for idx in stick_indexes:
        plt.plot((separations_idx[idx],separations_idx[idx]),yminmax,color='magenta',ls='dotted')
    
    ax1 = fig1.axes[0]
    ax1.set_ylabel('seconds')
    ax1.set_ylim(yminmax)
    ax1.set_xlabel('sample number')
    plt.legend()
    pngname = '%s_full.png' % png_rootname
    plt.savefig(pngname,format='png',dpi=100,bbox_inches='tight')

    if zoomx is not None:
        ax1.set_xlim(zoomx)
        ax1.set_ylim(tstamps[zoomx[0]],tstamps[zoomx[1]])
        pngname = '%s_zoom.png' % png_rootname
        plt.savefig(pngname,format='png',dpi=100,bbox_inches='tight')

    #### second plot: slope removed
    ttl = '%s horizontal' % tstamps_title
    png_rootname = '%s_%s' % (ttl.lower().replace(' ','_'),self.obsdate.strftime('%Y%m%d-%H%M%S'))
    fig2 = plt.figure(figsize=(16,8))
    fig2.canvas.set_window_title('plt: %s' % ttl)
    plt.suptitle(ttl)
    plt.title(datainfo)
    fig2.text(0.9,0.9,sample_period_txt,ha='right')
    
    plt.plot([0,len(indextime)],[0,0],                       label='index time')
    plt.plot(tstamps-indextime,         ls='none',marker='o',label='derived timestamps')
    plt.plot(compstamps-indextime,      ls='none',marker='*',label='computer time')
    plt.plot(gps-indextime,             ls='none',marker='x',label='GPS date')

    ax2 = fig2.axes[0]
    ax2.set_ylabel('seconds')
    ax2.set_xlabel('sample number')
    horiz_y = tstamps - indextime
    #yminmax = (horiz_y.min(),horiz_y,max())
    yminmax = (-0.5,0.015)
    ax2.set_ylim(yminmax)
    for idx in jump_indexes:
        plt.plot((separations_idx[idx],separations_idx[idx]),yminmax,color='red',ls='dashed')
    for idx in stick_indexes:
        plt.plot((separations_idx[idx],separations_idx[idx]),yminmax,color='magenta',ls='dotted')
    plt.legend()
    pngname = '%s_full.png' % png_rootname
    plt.savefig(pngname,format='png',dpi=100,bbox_inches='tight')

    if zoomx is not None:
        ax2.set_xlim(zoomx)

    if zoomy is not None:
        ax2.set_ylim(zoomy)

    if zoomx is not None or zoomy is not None:
        pngname = '%s_zoom.png' % png_rootname
        plt.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    
    return


