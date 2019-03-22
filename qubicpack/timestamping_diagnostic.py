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
    if hk is None: hk = 'ASIC_SUMS'
    if hk.upper() == 'PLATFORM': hk = 'INTERN_HK'
    if hk.upper() == 'ASIC': hk = 'CONF_ASIC1'
    if hk.upper() == 'EXTERN': hk = 'EXTERN_HK'
    if hk.upper().find('SCI')==0: hk = 'ASIC_SUMS'
    
    if not hk in self.hk.keys():
        print('Please give a valid HK.  Valid names are: %s' % ', '.join(self.hk.keys()))
        return None
              
    
    pps = self.hk[hk]['PPS']
    gps = self.hk[hk]['GPSDate']
    compstamps  = self.hk[hk]['ComputerDate']
    npts = len(pps)
    if hk=='ASIC_SUMS':
        sample_period = self.sample_period()
    else:
        sample_period = float(compstamps.max() - compstamps.min())/len(compstamps)
    indextime = sample_period*np.arange(npts)

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
            pps_indexes.append(idx)
            separations.append(gps[idx] - prev)
            separations_idx.append(idx)
            prev = gps[idx]            

    separations = np.array(separations[1:])
    separations_idx = np.array(separations_idx[1:])
    print('number of separations: %i' % len(separations))
    print('number of samples: %i' % len(pps))

    mean_separation = separations.mean()
    print('mean separation between pulses is %.3f second' % mean_separation)
    max_separation = separations.max()
    print('max separation between pulses is %.3f second' % max_separation)
    min_separation = separations.min()
    print('min separation between pulses is %.3f second' % min_separation)

    jump_indexes = np.where(separations>1+epsilon)[0]
    print('there are %i jumps at:  %s' % (len(jump_indexes),separations_idx[jump_indexes]))

    stick_indexes = np.where(separations<1-epsilon)[0]
    print('there are %i sticks at: %s' % (len(jump_indexes),separations_idx[stick_indexes]))

    tstamps = self.pps2date(pps,gps)
    t0 = tstamps[0]

    if indextime is not None:
        slope = indextime[-1]/(len(tstamps)-1)
    else:
        slope = (tstamps[-1] - tstamps[0])/(len(tstamps)-1)
    offset = t0

    # subtract the progression to view only the residues (horizontal line)
    xpts = np.arange(len(tstamps))
    ypts = slope*xpts + offset

    # mark problems with a vertical line
    yminmax = (compstamps.min(),compstamps.max())

    plt.ion()
    ##### plot the scientific PPS
    ttl = 'PPS Science Data'
    png_rootname = '%s_%s' % (ttl.lower().replace(' ','_'),self.obsdate.strftime('%Y%m%d-%H%M%S'))
    fig0 = plt.figure(figsize=(16,8))
    fig0.canvas.set_window_title('plt: %s' % ttl)
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
    ttl = 'Timestamps diagnostic'
    png_rootname = '%s_%s' % (ttl.lower().replace(' ','_'),self.obsdate.strftime('%Y%m%d-%H%M%S'))
    fig1 = plt.figure(figsize=(16,8))
    fig1.canvas.set_window_title('plt: %s' % ttl)
    plt.suptitle(ttl)
    plt.title(datainfo)
    if indextime is not None:
        plt.plot(indextime+t0,                ls='none',marker='d',label='index time')
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
    ttl = '%s horizontal' % ttl
    png_rootname = '%s_%s' % (ttl.lower().replace(' ','_'),self.obsdate.strftime('%Y%m%d-%H%M%S'))
    fig2 = plt.figure(figsize=(16,8))
    fig2.canvas.set_window_title('plt: %s' % ttl)
    plt.suptitle(ttl)
    plt.title(datainfo)
    if indextime is not None:
        plt.plot(indextime-ypts+t0,                ls='none',marker='d',label='index time')
    plt.plot(tstamps-ypts,                         ls='none',marker='o',label='derived timestamps')
    plt.plot(compstamps-ypts,                      ls='none',marker='*',label='computer time')
    plt.plot(gps-ypts,                             ls='none',marker='x',label='GPS date')

    ax2 = fig2.axes[0]
    ax2.set_ylabel('seconds')
    ax2.set_xlabel('sample number')
    horiz_y = tstamps - ypts
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


