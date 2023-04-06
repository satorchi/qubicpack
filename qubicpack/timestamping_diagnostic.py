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
from matplotlib import pyplot as plt
import numpy as np
from qubicpack.utilities import figure_window_title

'''
   the first version of this script was tested on the following data:
   2019/2019-02-22/2019-02-22_16.28.01__Scan2
'''

def find_pps_peaks_scidata(pps):
    '''
    find the pps peaks for scientific data
    for the scientific data, the PPS is a square wave, 1 second on, 1 second off
    '''
    retval = {}
        
    # when the gradient is non-zero we have a pulse rise or fall detected
    ppsgrad = np.gradient(pps)
    idxhigh = np.where(ppsgrad>0)[0]
    idxlow  = np.where(ppsgrad<0)[0]

    # the gradient is non-zero for the sample at the start of the rise, and the next one at the end of the rise
    # we take every other one, which is when the high pulse is detected
    oddidx = 1 + 2*np.arange(len(idxhigh)//2)
    pps_high = idxhigh[oddidx]
    retval['pps_high'] = pps_high
    
    # we do the same for the end of the pulse
    oddidx = 2*np.arange(len(idxlow)//2)
    pps_low_idx = idxlow[oddidx]
    retval['pps_low'] = pps_low_idx

    # the PPS is a square wave 1 second on, 1 second off
    # samples per pps is the number of samples per second
    x1 = np.empty(len(pps_high)+1,dtype=int)
    x2 = np.empty(len(pps_high)+1,dtype=int)
    x1[0:len(pps_high)] = pps_high
    x2[1:len(pps_high)+1] = pps_high
    samples_per_pps = (x1 - x2)[1:-1]
    sample_rate = 0.5*samples_per_pps.mean()
    retval['samples_per_pps'] = samples_per_pps
    retval['mean_samples_per_second'] = sample_rate
    retval['pps_indexes'] = pps_high

    return retval

def find_pps_peaks_hkdata(pps):
    '''
    find the pps peaks for housekeeping data
    for housekeeping data, the PPS is a pulse (single sample) followed by ~6 samples
    '''
    retval = {}

    pps_high = np.where(pps==1)[0]
    retval['pps_high'] = pps_high
    retval['pps_indexes'] = pps_high
    
    pps_low = np.where(pps==0)[0]
    retval['pps_low'] = pps_low

    # samples per pps is the number of samples per second
    x1 = np.empty(len(pps_high)+1,dtype=int)
    x2 = np.empty(len(pps_high)+1,dtype=int)
    x1[0:len(pps_high)] = pps_high
    x2[1:len(pps_high)+1] = pps_high
    samples_per_pps = (x1 - x2)[1:-1]
    sample_rate = 0.5*samples_per_pps.mean()
    retval['samples_per_pps'] = samples_per_pps
    retval['mean_samples_per_second'] = sample_rate    
    
    return retval

def pps2date(pps,gps,gps_sample_offset=0,verbosity=0):
    '''
    convert the gps date to a precise date given the pps
    '''

    if gps is None or gps.min()<1494486000 or gps.max()<1494486000 or len(np.unique(gps))<2:
        if verbosity>1: print('ERROR! Bad GPS data.')
        return None
    
    npts = len(pps)
    pps_separation=1  # exactly one second between pulses
    epsilon = 0.1

    separations = []
    pps_high = np.where(pps==1)[0]
    # select the first/last PPS in each series
    pps_indexes = []
    prev = gps[0]
    for idx in pps_high:
        if (idx>0 and pps[idx-1]==0)\
           or (idx<npts-1 and pps[idx+1]==0):
            sep = gps[idx] - prev
            if sep != 0: # next PPS valid only if we have a non-zero step (modif by MP)
                pps_indexes.append(idx)    
                separations.append(sep)
                prev = gps[idx]

    separations = np.array(separations[1:])
    if separations.size==0:
        if verbosity>0: print('no pps intervals!')
        return None
    
    if verbosity>1: print('mean pps interval is %.4f second' % separations.mean())
    if verbosity>1: print('max pps interval is  %.4f second' % separations.max())
    if verbosity>1: print('min pps interval is  %.4f second' % separations.min())
            
    # find the GPS date corresponding to the PPS
    tstamp = -np.ones(npts)
    prev_gps = gps[0]
    #offset = 50 # delay after PPS for valid GPS (this should be different for scientific and housekeeping) 
    for idx in pps_indexes:
        gps_at_pps = gps[idx]

        ### original algorithm  replaced by the one below (MP & JCH) ##################################
        ## we use the GPS timestamp from a bit later
        # offset_idx = idx + offset
        # if offset_idx>=npts:offset_idx=npts-1
        # next_gps = gps[offset_idx]
        # tstamp[idx] = next_gps
        ###############################################################################################

        ###  assuming the PPS arrives before the corresponding time given by the GPS ##################
        # the pps arrives just before the corresponding gps
        # so we simply add 1 second to current gps value (gps increments in steps of 1 second exactly)
        # (modification by MP & JCH)
        # tstamp[idx] = gps_at_pps + 1
        ###############################################################################################


        ###  assuming the PPS arrives after the corresponding time given by the GPS ###################
        tstamp[idx] = gps_at_pps
        ###############################################################################################

        
    # now finally do the interpolation for the time axis
    first_sample_period = None    
    for idx in range(len(pps_indexes)-1):
        diff_idx = pps_indexes[idx+1] - pps_indexes[idx]
        pps_separation = tstamp[pps_indexes[idx+1]]-tstamp[pps_indexes[idx]]
        sample_period = pps_separation/diff_idx
        if first_sample_period is None:
            first_sample_period = sample_period
        for idx_offset in range(diff_idx):
            tstamp[pps_indexes[idx]+idx_offset] = tstamp[pps_indexes[idx]] + idx_offset*sample_period

    last_sample_period = sample_period

    # do the first bit before the first PPS
    tstamp0 = tstamp[pps_indexes[0]]
    for idx in range(pps_indexes[0]+1):
        tstamp[pps_indexes[0] - idx] = tstamp0 - idx*first_sample_period

    # do the last bit after the last PPS
    tstampF = tstamp[pps_indexes[-1]]
    for idx in range(npts - pps_indexes[-1]):
        tstamp[pps_indexes[-1] + idx] = tstampF + idx*last_sample_period

    return tstamp


def timestamp_diagnostic(self,hk=None,asic=None):
    '''
    calculate the diagnostic of the derived timestamps
    '''
    analysis = {}
    analysis['dataset'] = self.dataset_name
    
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
            if self.asic_list[asic-1] is None:
                self.printmsg('Error! No data for ASIC %i.' % asic)
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
    
    if hk=='ASIC_SUMS':
        ans = find_pps_peaks_scidata(pps)
    else:
        ans = find_pps_peaks_hkdata(pps)
        
    for key in ans.keys():
        analysis[key] = ans[key]
    pps_high = analysis['pps_high']
    pps_indexes = analysis['pps_indexes']
    pps_low = analysis['pps_low']
    samples_per_pps = analysis['samples_per_pps']
    sample_rate = analysis['mean_samples_per_second']
        

    self.printmsg('number of samples: %i' % len(pps),verbosity=2)
    self.printmsg('number of pulses: %i' % len(pps_high),verbosity=2)
    self.printmsg('samples per second: %.3f' % sample_rate,verbosity=2)


    if len(pps_indexes)<2:
        analysis['samples_per_pps'] = None
        analysis['sample_period_txt'] = 'Insufficient PPS data'
        analysis['avg samples_per_pps txt'] = 'Insufficient PPS data'
        analysis['weird_idx'] = []
        analysis['weird_event'] = []
        self.printmsg('ERROR! insufficient PPS data',verbosity=2)
        return analysis

    mean_samples_per_pps = samples_per_pps.mean()
    analysis['mean_samples_per_pps'] = mean_samples_per_pps
    sigma_samples_per_pps = samples_per_pps.std()
    analysis['sigma_samples_per_pps'] = sigma_samples_per_pps
    self.printmsg('mean number of samples between pulses: %i' % mean_samples_per_pps,verbosity=2)
    self.printmsg('expected number of samples between pulses: %.1f' % (float(len(pps))/len(pps_high)),verbosity=2)
    self.printmsg('sigma samples between pulses: %.1f' % sigma_samples_per_pps,verbosity=2)
    analysis['avg samples_per_pps txt'] = 'avg number of samples between events: %.2f (%.3f samples/second)' % (mean_samples_per_pps,sample_rate)

    weird_event = []
    weird_idx = []
    for idx,val in enumerate(samples_per_pps):
        sigma = np.abs(val-mean_samples_per_pps)/sigma_samples_per_pps
        if sigma > 2.0:
            weird_event.append(idx)
            weird_idx.append(pps_indexes[idx])
    self.printmsg('number of anomolous events: %i' % len(weird_event),verbosity=2)
    analysis['weird_idx'] = np.array(weird_idx)
    analysis['weird_event'] = np.array(weird_event)

    mean_separation = mean_samples_per_pps/sample_rate
    self.printmsg('mean separation between pulses is %.4f second' % mean_separation,verbosity=2)
    max_separation = samples_per_pps.max()/sample_rate
    self.printmsg('max separation between pulses is %.4f second' % max_separation,verbosity=2)
    min_separation = samples_per_pps.min()/sample_rate
    self.printmsg('min separation between pulses is %.4f second' % min_separation,verbosity=2)
    
    ########## use the GPS to make the timeaxis
    if gps is None or max(gps)==0 or max(gps)==min(gps):
        refclock = compstamps
        self.printmsg('Bad GPS data.  Using the computer time instead.',verbosity=2)
        analysis['tstamps_title'] += ' (based on computer clock.  Bad GPS data)'
    else:
        refclock = gps

    
    tstamps = pps2date(pps,refclock,verbosity=self.verbosity)
    analysis['tstamps'] = tstamps
    if tstamps is None:
        self.printmsg('Could not assign timestamps from reference clock',verbosity=2)
        return analysis
    t0 = tstamps[0]

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
        figure_window_title(fig,ttl)
        fig.suptitle(self.infotext(),fontsize=fontsize)
        ax = fig.add_axes((0.05,0.08,0.9,0.8))
        newplot = True

    ax.text(0.01,1.00,ttl,ha='left',va='bottom',fontsize=fontsize,transform=ax.transAxes)
    ax.text(0.99,1.01,analysis['sample_period_txt'],ha='right',va='bottom',fontsize=fontsize,transform=ax.transAxes)
    ax.text(0.01,1.05,analysis['avg samples_per_pps txt'],fontsize=fontsize,transform=ax.transAxes)
    if analysis['lost_txt'] is not None:
        ax.text(0.99,1.1,analysis['lost_txt'],ha='right',va='bottom',fontsize=fontsize,transform=ax.transAxes)

    ax.plot(analysis['pps'],label='PPS')
    ax.plot(analysis['pps_high'],np.ones(analysis['pps_high'].size,dtype=int),ls='none',marker='^',color='green',label='pulse start')
    ax.plot(analysis['pps_low'],np.zeros(analysis['pps_low'].size,dtype=int),ls='none',marker='v',color='pink',label='pulse end')
    ax.set_ylabel('PPS Level')
    ax.set_xlabel('sample number')
    ax.tick_params(axis='both',labelsize=fontsize)
    for idx in analysis['weird_idx']:
        ax.plot((idx,idx),(0,1),color='red',ls='dashed',label='anomaly')
    ax.legend(loc='upper right')    

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
    if analysis is None:
        if ax is None:
            return
        ax.text(0.5,0.5,'Insufficient PPS Data',va='center',ha='center',fontsize=2*fontsize,transform=ax.transAxes)
        ax.get_yaxis().set_visible(False)
        ax.get_xaxis().set_visible(False)
        return ax


    ttl = analysis['nsamples_title']
    png_rootname = '%s_%s' % (ttl.lower().replace(' ','_'),self.obsdate.strftime('%Y%m%d-%H%M%S'))

    newplot = False
    if ax is None:
        newplot = True
        fig = plt.figure(figsize=(16,8))
        figure_window_title(fig,ttl)
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
        figure_window_title(fig,ttl)
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
    ax.legend(fontsize=fontsize,loc='upper left')
    if newplot:
        pngname = '%s_full.png' % png_rootname
        fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')

    if zoomx is not None:
        ax.set_xlim(zoomx)
        ax.set_ylim(analysis['tstamps'][zoomx[0]],analysis['tstamps'][zoomx[1]])
        if newplot:
            pngname = '%s_zoom.png' % png_rootname
            fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
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
        figure_window_title(fig,ttl)
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
    ax.legend(fontsize=fontsize,loc='upper right')

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
        delta = cn[idx] - generated_cn[idx]
        if delta!=0:
            counter += delta # lost packet. resync
        counter += 1
        
    delta = cn - generated_cn
    idx_lost = np.where(delta!=0)[0]

    if len(idx_lost)==0:
        self.printmsg('No lost packets!',verbosity=2)
    else:
        self.printmsg('%i lost packets.' % len(idx_lost),verbosity=2)
    
    return idx_lost
