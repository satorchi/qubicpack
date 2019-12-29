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

def timestamp_diagnostic(self,hk=None,asic=None):
    '''
    calculate the diagnostic of the derived timestamps
    '''
    analysis = {}
    
    hk = self.qubicstudio_filetype_truename(hk)
    if hk is None: hk = 'ASIC_SUMS'
    analysis['hk'] = hk
    
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
    analysis['asic'] = asic
            
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
    analysis['pps_title'] = pps_title
    analysis['tstamps_title'] = tstamps_title
    analysis['nsamples_title'] = nsamples_title

    pps = self.pps(hk=hk,asic=asic)
    if pps is None: return None
    gps = self.gps(hk=hk,asic=asic)
    if gps is None: return None
    analysis['pps'] = pps
    analysis['gps'] = gps
    
    compstamps  = HK[hk]['ComputerDate']
    npts = len(pps)
    if hk=='ASIC_SUMS':
        sample_period = self.sample_period(asic=asic)
    else:
        sample_period = float(compstamps.max() - compstamps.min())/len(compstamps)
    analysis['compstamps'] = compstamps
    analysis['npts'] = npts
    analysis['sample_period'] = sample_period
        
    lost_txt = None
    lost_idx = self.lost_packets(hk=hk,asic=asic)
    if lost_idx is not None:
        lost_txt = '%i lost packets' % len(lost_idx)
    analysis['lost_idx'] = lost_idx
    analysis['lost_txt'] = lost_txt
    
    
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
    self.printmsg('number of separations: %i' % len(separations),verbosity=2)
    self.printmsg('number of samples: %i' % len(pps),verbosity=2)
    analysis['pps_high'] = pps_high
    analysis['separations'] = separations
    analysis['pps_indexes'] = pps_indexes

    if len(pps_indexes)<2:
        analysis['samples_per_pps'] = None
        analysis['sample_period_txt'] = 'Insufficient PPS data'
        analysis['avg samples_per_pps txt'] = 'Insufficient PPS data'
        analysis['weird_idx'] = []
        analysis['weird_event'] = []
        self.printmsg('ERROR! insufficient PPS data',verbosity=2)
        return analysis
    
    samples_per_pps = []
    idx_prev = pps_indexes[0]
    for idx in pps_indexes[1:]:
        nsamples = idx - idx_prev
        samples_per_pps.append(nsamples)
        idx_prev = idx
    samples_per_pps = np.array(samples_per_pps,dtype=np.float)
    mean_samples_per_pps = samples_per_pps.mean()
    sigma_samples_per_pps = samples_per_pps.std()
    self.printmsg('mean number of samples between pulses: %i' % mean_samples_per_pps,verbosity=2)
    self.printmsg('expected number of samples between pulses: %.1f' % (float(len(pps))/len(separations)),verbosity=2)
    self.printmsg('sigma samples between pulses: %.1f' % sigma_samples_per_pps,verbosity=2)
    analysis['samples_per_pps'] = samples_per_pps
    analysis['mean_samples_per_pps'] = mean_samples_per_pps
    analysis['avg samples_per_pps txt'] = 'avg number of samples between events: %.2f' % mean_samples_per_pps
    analysis['sigma_samples_per_pps'] = sigma_samples_per_pps

    weird_event = []
    weird_idx = []
    for idx,val in enumerate(samples_per_pps):
        sigma = np.abs(val-mean_samples_per_pps)/sigma_samples_per_pps
        if sigma > 2.0:
            weird_event.append(idx)
            weird_idx.append(pps_indexes[idx])
    self.printmsg('number of anomolous events: %i' % len(weird_event),verbosity=2)
    analysis['weird_idx'] = weird_idx
    analysis['weird_event'] = weird_event

    mean_separation = separations.mean()
    self.printmsg('mean separation between pulses is %.4f second' % mean_separation,verbosity=2)
    max_separation = separations.max()
    self.printmsg('max separation between pulses is %.4f second' % max_separation,verbosity=2)
    min_separation = separations.min()
    self.printmsg('min separation between pulses is %.4f second' % min_separation,verbosity=2)
    
    tstamps = self.pps2date(pps,gps)
    t0 = tstamps[0]
    analysis['tstamps'] = tstamps

    # do a line fit to get drift of the sampling clock
    xpts = np.arange(len(tstamps))
    linefit = np.polyfit(xpts,tstamps,1,full=True)
    slope = linefit[0][0]
    offset = linefit[0][1]
    derived_sample_period = slope
    analysis['linefit'] = linefit
    analysis['slope'] = slope
    analysis['offset'] = offset
    sample_period_txt = 'expected sample period: %.6f msec\nderived sample period: %.6f msec' % (sample_period*1000,slope*1000)
    analysis['sample_period_txt'] = sample_period_txt
    self.printmsg(sample_period_txt,verbosity=2)

    # subtract the progression to view only the residues (horizontal line)
    indextime = slope*xpts + offset
    analysis['indextime'] = indextime
    return analysis


def plot_pps(self,analysis=None,hk=None,zoomx=None,zoomy=None,asic=None,ax=None,fontsize=12):
    '''
    make a plot of the PPS
    '''

    if analysis is None:
        analysis = self.timestamp_diagnostic(hk=hk,asic=asic)
    if analysis is None: return

    ttl = analysis['pps_title']
    png_rootname = '%s_%s' % (ttl.lower().replace(' ','_'),self.obsdate.strftime('%Y%m%d-%H%M%S'))

    newplot = False
    if ax is None:
        fig = plt.figure(figsize=(16,8))
        fig.canvas.set_window_title('plt: %s' % ttl)
        fig.suptitle(self.infotext(),fontsize=fontsize)
        ax = fig.add_axes((0.05,0.08,0.9,0.8))
        newplot = True

    ax.text(0.01,1.00,ttl,ha='left',va='bottom',fontsize=fontsize,transform=ax.transAxes)
    ax.text(0.99,1.01,analysis['sample_period_txt'],ha='right',va='bottom',fontsize=fontsize,transform=ax.transAxes)
    ax.text(0.01,1.05,analysis['avg samples_per_pps txt'],fontsize=fontsize,transform=ax.transAxes)
    if analysis['lost_txt'] is not None:
        ax.text(0.99,1.1,analysis['lost_txt'],ha='right',va='bottom',fontsize=fontsize,transform=ax.transAxes)

    ax.plot(analysis['pps'])
    ax.set_ylabel('PPS Level')
    ax.set_xlabel('sample number')
    ax.tick_params(axis='both',labelsize=fontsize)
    for idx in analysis['weird_idx']:
        ax.plot((idx,idx),(0,1),color='red',ls='dashed')

    if newplot:
        pngname = '%s_full.png' % png_rootname
        fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')

    if zoomx is not None:
        ax.set_xlim(zoomx)
        if newplot:
            pngname = '%s_zoom.png' % png_rootname
            fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    return

def plot_pps_nsamples(self,analysis=None,hk=None,zoomx=None,zoomy=None,asic=None,ax=None,fontsize=12):
    '''
    make a plot of the number of samples per PPS 
    '''

    if analysis is None:
        analysis = self.timestamp_diagnostic(hk=hk,asic=asic)
    if analysis is None: return

    ttl = analysis['nsamples_title']
    png_rootname = '%s_%s' % (ttl.lower().replace(' ','_'),self.obsdate.strftime('%Y%m%d-%H%M%S'))

    newplot = False
    if ax is None:
        newplot = True
        fig = plt.figure(figsize=(16,8))
        fig.canvas.set_window_title('plt: %s' % ttl)
        fig.suptitle(self.infotext(),fontsize=fontsize)
        ax = fig.add_axes((0.05,0.08,0.9,0.8))

    ax.text(0.01,1.00,ttl,ha='left',va='bottom',fontsize=fontsize,transform=ax.transAxes)
    ax.text(0.99,1.01,analysis['sample_period_txt'],ha='right',va='bottom',fontsize=fontsize,transform=ax.transAxes)
    ax.text(0.01,1.05,analysis['avg samples_per_pps txt'],fontsize=fontsize,transform=ax.transAxes)
    if analysis['lost_txt'] is not None:
        ax.text(0.99,1.1,analysis['lost_txt'],ha='right',va='bottom',fontsize=fontsize,transform=ax.transAxes)

    if analysis['samples_per_pps'] is not None:
        ax.plot(analysis['samples_per_pps'])
    else:
        ax.text(0.5,0.5,'Insufficient PPS data',va='center',ha='center',fontsize=2*fontsize,transform=ax.transAxes)
    ax.set_ylabel('Number of samples per PPS event',fontsize=fontsize)
    ax.set_xlabel('PPS event number',fontsize=fontsize)
    ax.tick_params(axis='both',labelsize=fontsize)
    for idx in analysis['weird_event']:
        ax.plot((idx,idx),(analysis['samples_per_pps'].min(),analysis['samples_per_pps'].max()),color='red',ls='dashed')

    if newplot:
        pngname = '%s_full.png' % png_rootname
        fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')

    if zoomx is not None:
        ax.set_xlim(zoomx)
        if newplot:
            pngname = '%s_zoom.png' % png_rootname
            fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    return

def plot_timestamp_diagnostic_fig1(self,analysis=None,hk=None,zoomx=None,zoomy=None,asic=None,ax=None,fontsize=12):
    '''
    make a plot of the timestamp diagnostic
    '''

    if analysis is None:
        analysis = self.timestamp_diagnostic(hk=hk,asic=asic)
    if analysis is None: return

    ttl = analysis['tstamps_title']
    png_rootname = '%s_%s' % (ttl.lower().replace(' ','_'),self.obsdate.strftime('%Y%m%d-%H%M%S'))

    newplot = False
    if ax is None:
        newplot = True
        fig = plt.figure(figsize=(16,8))
        fig.canvas.set_window_title('plt: %s' % ttl)
        fig.suptitle(self.infotext(),fontsize=fontsize)
        ax = fig.add_axes((0.05,0.08,0.9,0.8))
        
    ax.text(0.50,1.00,ttl,ha='center',va='bottom',fontsize=fontsize,transform=ax.transAxes)
    ax.text(0.99,1.01,analysis['sample_period_txt'],ha='right',va='bottom',fontsize=fontsize,transform=ax.transAxes)
    ax.text(0.01,1.05,'avg number of samples between events: %.2f' % analysis['samples_per_pps'].mean(),fontsize=fontsize,transform=ax.transAxes)
    if analysis['lost_txt'] is not None:
        ax.text(0.99,1.1,analysis['lost_txt'],ha='right',va='bottom',fontsize=fontsize,transform=ax.transAxes)

    # mark problems with a vertical line
    yminmax = (analysis['compstamps'].min(),analysis['compstamps'].max())
    pps_high = analysis['pps_high']
    pps_indexes = analysis['pps_indexes']
    ax.plot(analysis['indextime'],                       ls='none',marker='d',label='index time')
    ax.plot(analysis['tstamps'],                         ls='none',marker='o',label='derived timestamps')
    ax.plot(analysis['compstamps'],                      ls='none',marker='*',label='computer time')
    ax.plot(analysis['gps'],                             ls='none',marker='x',label='GPS date')
    ax.plot(pps_high,analysis['tstamps'][pps_high],      ls='none',marker='^',label='pps high')
    ax.plot(pps_indexes,analysis['tstamps'][pps_indexes],ls='none',marker='o',label='pps indexes',markerfacecolor='none',markersize=16,color='black')
    for idx in analysis['weird_idx']:
        ax.plot((idx,idx),yminmax,color='red',ls='dashed')
    
    ax.set_ylabel('seconds',fontsize=fontsize)
    ax.set_ylim(yminmax)
    ax.set_xlabel('sample number',fontsize=fontsize)
    ax.tick_params(axis='both',labelsize=fontsize)
    ax.legend(fontsize=fontsize)
    if newplot:
        pngname = '%s_full.png' % png_rootname
        fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')

    if zoomx is not None:
        ax.set_xlim(zoomx)
        ax.set_ylim(tstamps[zoomx[0]],tstamps[zoomx[1]])
        if newplot:
            pngname = '%s_zoom.png' % png_rootname
            ax.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    return

def plot_timestamp_diagnostic_fig2(self,analysis=None,hk=None,zoomx=None,zoomy=None,asic=None,ax=None,fontsize=12):
    '''
    make a plot of the timestamp diagnostic with the slope removed
    '''

    if analysis is None:
        analysis = self.timestamp_diagnostic(hk=hk,asic=asic)
    if analysis is None: return

    
    ttl = '%s horizontal' % analysis['tstamps_title']
    png_rootname = '%s_%s' % (ttl.lower().replace(' ','_'),self.obsdate.strftime('%Y%m%d-%H%M%S'))

    newplot = False
    if ax is None:
        fig = plt.figure(figsize=(16,8))
        fig.canvas.set_window_title('plt: %s' % ttl)
        fig.suptitle(ttl,fontsize=fontsize)
        ax = fig.add_axes((0.05,0.08,0.9,0.8))

    ax.text(0.50,1.00,ttl,ha='center',va='bottom',fontsize=fontsize,transform=ax.transAxes)
    ax.text(0.99,1.01,analysis['sample_period_txt'],ha='right',va='bottom',fontsize=fontsize,transform=ax.transAxes)
    ax.text(0.01,1.05,'avg number of samples between events: %.2f' % analysis['samples_per_pps'].mean(),fontsize=fontsize,transform=ax.transAxes)
    if analysis['lost_txt'] is not None:
        ax.text(0.99,1.1,analysis['lost_txt'],ha='right',va='bottom',fontsize=fontsize,transform=ax.transAxes)


    tstamps_horiz = analysis['tstamps'] - analysis['indextime']
    compstamps_horiz = analysis['compstamps'] - analysis['indextime']
    gps_horiz = analysis['gps'] - analysis['indextime']
    pps_indexes = analysis['pps_indexes']
    yminmax = (tstamps_horiz.min(),tstamps_horiz.max())
    # measure the peak to peak between the first and last PPS
    tstamps_horiz_between_pps = tstamps_horiz[pps_indexes[0]:pps_indexes[-1]]
    peak2peak = tstamps_horiz_between_pps.max() - tstamps_horiz_between_pps.min()
    peak2peak_txt = 'peak to peak variation: %.2f msec' % (1000*peak2peak)
    self.printmsg(peak2peak_txt,verbosity=2)
    ax.text(0.01,1.0,peak2peak_txt,ha='left',va='bottom',fontsize=fontsize,transform=ax.transAxes)

    
    
    ax.plot([0,len(analysis['indextime'])],[0,0],                       label='index time')
    ax.plot(tstamps_horiz,         ls='none',marker='o',label='derived timestamps')
    ax.plot(compstamps_horiz,      ls='none',marker='*',label='computer time')
    ax.plot(gps_horiz,             ls='none',marker='x',label='GPS date')

    ax.set_ylabel('seconds',fontsize=fontsize)
    ax.set_xlabel('sample number',fontsize=fontsize)
    ax.set_ylim(yminmax)
    for idx in analysis['weird_idx']:
        ax.plot((idx,idx),yminmax,color='red',ls='dashed')
    ax.legend(fontsize=fontsize)

    if newplot:
        pngname = '%s_full.png' % png_rootname
        fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')

    if zoomx is not None:
        ax.set_xlim(zoomx)

    if zoomy is not None:
        ax.set_ylim(zoomy)

    if zoomx is not None or zoomy is not None:
        if newplot:
            pngname = '%s_zoom.png' % png_rootname
            fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')

    return

def plot_timestamp_diagnostic(self,analysis=None,hk=None,zoomx=None,zoomy=None,asic=None,ax=None,fontsize=12):
    '''
    make a diagnostic plot of the derived timestamps
    '''

    if analysis is None:
        analysis = self.timestamp_diagnostic(hk=hk,asic=asic)
    if analysis is None: return

    plt.ion()
    ##### plot the  PPS
    self.plot_pps(analysis=analysis)

    ##### plot number of samples between each PPS event #####
    self.plot_pps_nsamples(analysis=analysis)
            
    ##### plot diagnostic #1
    self.plot_timestamp_diagnostic_fig1(analysis=analysis)

    #### second plot: slope removed
    self.plot_timestamp_diagnostic_fig2(analysis=analysis)
    
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
        if hk=='ASIC_SUMS':
            if asic is None:
                self.printmsg('Please enter a valid ASIC.')
                return None
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
        self.printmsg('No lost packets!',verbosity=2)
    else:
        self.printmsg('%i lost packets.' % len(idx_lost),verbosity=2)
    
    return idx_lost
