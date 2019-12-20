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

def plot_timestamp_diagnostic(self,hk=None,zoomx=None,zoomy=None,asic=None):
    '''
    make a diagnostic plot of the derived timestamps
    '''
    hk = self.qubicstudio_filetype_truename(hk)
    if hk is None: hk = 'ASIC_SUMS'
    
    if self.__object_type__!='qubicfp':
        asic = self.asic
        HK = self.hk
    else:
        if hk=='ASIC_SUMS':
            if asic is None:
                self.printmsg('Please enter a valid ASIC.')
                return None
            HK = self.asic_list[asic-1].hk
        else:
            HK = self.hk
            asic = None
            
    if not hk in HK.keys():
        self.printmsg('Please give a valid HK.  Valid names are: %s' % ', '.join(HK.keys()))
        return None
              
    asic_txt = ''
    if asic is not None:
        asic_txt = ' (ASIC %i)' % asic

    if hk=='ASIC_SUMS':
        pps_title = 'PPS Scientific Data'+asic_txt
        tstamps_title = 'Timestamps diagnostic for Scientific Data'+asic_txt
        nsamples_title = 'Number of samples between PPS events for Scientific Data'+asic_txt
    elif hk=='INTERN_HK':
        pps_title = 'PPS Platform'
        tstamps_title = 'Timestamps diagnostic for Platform Data'
        nsamples_title = 'Number of samples between PPS events for Platform Data'
    else:
        pps_title = 'PPS %s%s' % (hk,asic_txt)
        tstamps_title = 'Timestamps diagnostic for %s%s' % (hk,asic_txt)
        nsamples_title = 'Number of samples between PPS events for %s%s' % (hk,asic_txt)

        
    
    pps = self.pps(hk=hk,asic=asic)
    if pps is None: return 
    gps = self.gps(hk=hk,asic=asic)
    if gps is None: return
    
    compstamps  = HK[hk]['ComputerDate']
    npts = len(pps)
    if hk=='ASIC_SUMS':
        sample_period = self.sample_period(asic=asic)
    else:
        sample_period = float(compstamps.max() - compstamps.min())/len(compstamps)

    lost_txt = None
    lost_idx = self.lost_packets(hk=hk)
    if lost_idx is not None:
        lost_txt = '%i lost packets' % len(lost_idx)

    datainfo = self.infotext()
    
    ########## find where we have potential problems
    separations = []
    separations_idx = []
    pps_high = np.where(pps==1)[0]
    pps_indexes = []
    prev = gps[0]
    for idx in pps_high:
        if (idx>0 and pps[idx-1]==0) or (idx<npts-1 and pps[idx+1]==0):
            sep = gps[idx] - prev
            if sep != 0:
                pps_indexes.append(idx)
                separations.append(sep)
                separations_idx.append(idx)
                prev = gps[idx]

    separations = np.array(separations[1:])
    separations_idx = np.array(separations_idx[1:])
    self.printmsg('number of separations: %i' % len(separations))
    self.printmsg('number of samples: %i' % len(pps))

    samples_per_pps = []
    idx_prev = pps_indexes[0]
    for idx in pps_indexes[1:]:
        nsamples = idx - idx_prev
        samples_per_pps.append(nsamples)
        idx_prev = idx
    samples_per_pps = np.array(samples_per_pps,dtype=np.float)
    mean_samples_per_pps = samples_per_pps.mean()
    sigma_samples_per_pps = samples_per_pps.std()
    self.printmsg('mean number of samples between pulses: %i' % mean_samples_per_pps)
    self.printmsg('expected number of samples between pulses: %.1f' % (float(len(pps))/len(separations)))
    self.printmsg('sigma samples between pulses: %.1f' % sigma_samples_per_pps)
    weird_event = []
    weird_idx = []
    for idx,val in enumerate(samples_per_pps):
        sigma = np.abs(val-mean_samples_per_pps)/sigma_samples_per_pps
        if sigma > 2.0:
            weird_event.append(idx)
            weird_idx.append(pps_indexes[idx])
    self.printmsg('number of anomolous events: %i' % len(weird_event))

    mean_separation = separations.mean()
    self.printmsg('mean separation between pulses is %.4f second' % mean_separation)
    max_separation = separations.max()
    self.printmsg('max separation between pulses is %.4f second' % max_separation)
    min_separation = separations.min()
    self.printmsg('min separation between pulses is %.4f second' % min_separation)

    tstamps = self.pps2date(pps,gps)
    t0 = tstamps[0]

    # do a line fit to get drift of the sampling clock
    xpts = np.arange(len(tstamps))
    linefit = np.polyfit(xpts,tstamps,1,full=True)
    slope = linefit[0][0]
    offset = linefit[0][1]
    derived_sample_period = slope
    sample_period_txt = 'expected sample period: %.6f msec\nderived sample period: %.6f msec' % (sample_period*1000,slope*1000)
    self.printmsg(sample_period_txt)

    # subtract the progression to view only the residues (horizontal line)
    indextime = slope*xpts + offset

    # mark problems with a vertical line
    yminmax = (compstamps.min(),compstamps.max())

    plt.ion()
    ##### plot the  PPS
    ttl = pps_title
    png_rootname = '%s_%s' % (ttl.lower().replace(' ','_'),self.obsdate.strftime('%Y%m%d-%H%M%S'))
    fig0 = plt.figure(figsize=(16,8))
    fig0.canvas.set_window_title('plt: %s' % ttl)
    fig0.text(0.9,0.9,sample_period_txt,ha='right')
    plt.suptitle(ttl)
    plt.title(datainfo)
    if lost_txt is not None: fig0.text(0.5,0.92,lost_txt,ha='center')
    plt.plot(pps)
    ax0 = fig0.axes[0]
    ax0.set_ylabel('PPS Level')
    ax0.set_xlabel('sample number')
    for idx in weird_idx:
        plt.plot((idx,idx),(0,1),color='red',ls='dashed')    
    pngname = '%s_full.png' % png_rootname
    plt.savefig(pngname,format='png',dpi=100,bbox_inches='tight')

    if zoomx is not None:
        ax0.set_xlim(zoomx)
        pngname = '%s_zoom.png' % png_rootname
        plt.savefig(pngname,format='png',dpi=100,bbox_inches='tight')


    ##### plot number of samples between each PPS event #####
    ttl = nsamples_title
    png_rootname = '%s_%s' % (ttl.lower().replace(' ','_'),self.obsdate.strftime('%Y%m%d-%H%M%S'))
    fig = plt.figure(figsize=(16,8))
    fig.canvas.set_window_title('plt: %s' % ttl)
    fig.text(0.9,0.9,sample_period_txt,ha='right')
    fig.text(0.1,0.9,'avg number of samples between events: %.2f' % samples_per_pps.mean())
    plt.suptitle(ttl)
    plt.title(datainfo)
    if lost_txt is not None: fig.text(0.5,0.92,lost_txt,ha='center')
    plt.plot(samples_per_pps)
    ax = fig.axes[0]
    ax.set_ylabel('Number of samples per PPS event')
    ax.set_xlabel('PPS event number')
    for idx in weird_event:
        ax.plot((idx,idx),(samples_per_pps.min(),samples_per_pps.max()),color='red',ls='dashed')
    pngname = '%s_full.png' % png_rootname
    plt.savefig(pngname,format='png',dpi=100,bbox_inches='tight')

    if zoomx is not None:
        ax.set_xlim(zoomx)
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
    if lost_txt is not None: fig1.text(0.5,0.92,lost_txt,ha='center')
    plt.plot(indextime,                       ls='none',marker='d',label='index time')
    plt.plot(tstamps,                         ls='none',marker='o',label='derived timestamps')
    plt.plot(compstamps,                      ls='none',marker='*',label='computer time')
    plt.plot(gps,                             ls='none',marker='x',label='GPS date')
    plt.plot(pps_high,tstamps[pps_high],      ls='none',marker='^',label='pps high')
    plt.plot(pps_indexes,tstamps[pps_indexes],ls='none',marker='o',label='pps indexes',markerfacecolor='none',markersize=16,color='black')
    for idx in weird_idx:
        plt.plot((idx,idx),yminmax,color='red',ls='dashed')
    
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
    if lost_txt is not None: fig2.text(0.5,0.92,lost_txt,ha='center')

    tstamps_horiz = tstamps - indextime
    #yminmax = (-0.5,0.015)
    yminmax = (tstamps_horiz.min(),tstamps_horiz.max())
    # measure the peak to peak between the first and last PPS
    tstamps_horiz_between_pps = tstamps_horiz[pps_indexes[0]:pps_indexes[-1]]
    peak2peak = tstamps_horiz_between_pps.max() - tstamps_horiz_between_pps.min()
    peak2peak_txt = 'peak to peak variation: %.2f msec' % (1000*peak2peak)
    self.printmsg(peak2peak_txt)
    fig2.text(0.1,0.92,peak2peak_txt,ha='left')
    
    plt.plot([0,len(indextime)],[0,0],                       label='index time')
    plt.plot(tstamps-indextime,         ls='none',marker='o',label='derived timestamps')
    plt.plot(compstamps-indextime,      ls='none',marker='*',label='computer time')
    plt.plot(gps-indextime,             ls='none',marker='x',label='GPS date')

    ax2 = fig2.axes[0]
    ax2.set_ylabel('seconds')
    ax2.set_xlabel('sample number')
    ax2.set_ylim(yminmax)
    for idx in weird_idx:
        plt.plot((idx,idx),yminmax,color='red',ls='dashed')
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


def lost_packets(self,hk='sci',asic=None):
    '''
    check the QubicStudio science data for lost packets
    '''
    datatype = self.qubicstudio_filetype_truename(hk)

    if self.__object_type__!='qubicfp':
        asic = self.asic
        HK = self.hk
    else:
        if asic is None:
            self.printmsg('Please enter a valid ASIC.')
            return None
        if hk=='ASIC_SUMS':
            HK = self.asic_list[asic-1].hk
        else:
            HK = self.hk

    if datatype not in HK.keys():
        self.printmsg('No QubicStudio data of type %s!' % datatype)
        return None

    counter_key = None
    counter_max = None
    if datatype=='ASIC_SUMS':
        counter_key = 'CN'
        counter_max = 2**7
    if datatype=='INTERN_HK':
        counter_key = 'Platform-acqCount'
        counter_max = 2**16
        
    if counter_key is None:
        self.printmsg('No sample counter.')
        return None
        
    if  counter_key not in HK[datatype].keys():
        self.printmsg('Missing sample counter!')
        return None

    cn = HK[datatype][counter_key]
    npts = len(cn)
    counter = cn[0]
    generated_cn = np.zeros(npts)
    for idx in range(npts):
        if counter==counter_max: counter = 0
        generated_cn[idx] = counter
        counter += 1
        
    delta = cn - generated_cn
    idx_lost = np.where(delta!=0)[0]

    if len(idx_lost)==0:
        self.printmsg('No lost packets!')
    else:
        self.printmsg('%i lost packets.' % len(idx_lost))
    
    return idx_lost
