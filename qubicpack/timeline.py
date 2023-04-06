#!/usr/bin/env python
'''
$Id: timeline.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Tue 17 Oct 2017 19:04:04 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

utilities for QUBIC TES timeline data
'''
import numpy as np
import sys,os,time
import datetime as dt
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from qubicpack.utilities import TES_index, figure_window_title
from qubicpack.pix2tes import assign_pix2tes,pix2tes,tes2pix
from qubicpack.plot_fp import plot_fp
from qubicpack.timestamping_diagnostic import pps2date

def exist_timeline_data(self):
    '''
    check if we have timeline data
    '''
    if not isinstance(self.tdata,list):
        return False
    if len(self.tdata)==0:
        return False
    if not 'TIMELINE' in self.tdata[0].keys():
        return False
    if not isinstance(self.tdata[0]['TIMELINE'],np.ndarray):
        return False
    return True

def ntimelines(self):
    '''
    return the number of timelines collected
    '''
    if not self.exist_timeline_data():
        return 0
    return len(self.tdata)

def timeline(self,TES,timeline_index=0):
    '''
    return the timeline for a given TES and timeline index

    timeline_index is leftover from the early days when there were multiple timelines per fits file (2017/18)
    '''
    if not self.exist_timeline_data():return None
    ntimelines=self.ntimelines()
    if timeline_index>=ntimelines:
        self.printmsg('ERROR! timeline index out of range.  Enter an index between 0 and %i' % (ntimelines-1))
        return None
    TES_idx=TES_index(TES)
    if TES_idx is None:return None
    timeline=self.tdata[timeline_index]['TIMELINE'][TES_idx,:]
    return timeline

def timeline_array(self,timeline_index=0):
    '''
    return an array with the timelines for all the TES
    '''
    if not self.exist_timeline_data():return None    
    return self.tdata[timeline_index]['TIMELINE']

def amplitude2DAC(self,amplitude):
    '''
    convert bias voltage amplitude in Volts to DAC to send to QubicStudio
    '''

    '''
    # QubicStudio V2
    if amplitude > 0 and amplitude <= 9:
        DACamplitude = amplitude / 0.001125 - 1 # "-1" is back Mon 30 Oct 2017 19:51:42 CET
        #DACamplitude = amplitude / 0.001125  # "-1" removed Tue 17 Oct 2017 14:09:47 CEST
    else:
        DACamplitude = 65536 + amplitude / 0.001125
    '''
    
    DACamplitude= amplitude / (2*self.DAC2V)      
    DACamplitude = int(np.round(DACamplitude))
    return DACamplitude

def bias_offset2DAC(self,bias):
    '''
    convert bias offset voltage in Volts to DAC to send to QubicStudio
    '''

    '''
    # QubicStudio v2
    # conversion constant DAC <- Volts
    # A = 2.8156e-4
    A = 284.5e-6
    if bias > 0 and bias <= 9:
        DACoffset = bias / A - 1 # "-1" is back Mon 30 Oct 2017 19:51:57 CET
        #DACoffset = bias / A # "-1" removed Tue 17 Oct 2017 15:47:03 CEST
    else:
        DACoffset = 65536 + bias / A
    '''

    '''
    Mon 12 Feb 2018 16:37:49 CET:  max bias offset is 8.6081536 with DAC2V=2.627e-4

    bias offset can go from 0 to 2^15*DAC2V
    DACoffset=0     ->  0V
    DACoffset=32768 ->  8.608V
    DACoffset=32769 ->  0V
    DACoffset=32770 -> -8.608V
                        and decreasing from there

    Tue 13 Feb 2018 13:53:34 CET
    max bias is now 8.837V.  It depends on the power supply of the FPGA card etc
    '''
    max_offset=self.DAC2V * 2**15
    if abs(bias)>max_offset:
        self.printmsg('WARNING! Cannot set bias offset greater than %.3fV.' % max_offset)
        if bias<0.0:
            bias=-max_offset
        else:
            bias= max_offset
        self.printmsg('Setting bias to %.3fV' % bias)


    DACoffset = abs(bias) / self.DAC2V
    DACoffset = int(np.round(DACoffset))
    if bias<0.0: DACoffset = 2**16 - DACoffset

    self.debugmsg("DACoffset=%i" % DACoffset)
    return DACoffset

def sample_period(self,timeline_index=None):
    '''
    the integration time per sample.  This is in seconds.
    2MHz is the sampling rate of the ASIC
    '''
    if not self.exist_timeline_data():return None

    if timeline_index is None: timeline_index = 0

    if 'NSAMPLES' in self.tdata[timeline_index].keys():
        self.nsamples = self.tdata[timeline_index]['NSAMPLES']
        
    if self.nsamples is None:
        time_axis = self.timeline_computertime(timeline_index)
        timeline_npts = len(time_axis)
        if time_axis is None:  return None
            
        measured_sample_period = (time_axis[-1] - time_axis[0])/(timeline_npts-1)
        self.nsamples = measured_sample_period

    if self.nsamples is None:
        self.printmsg('WARNING!  No value for nsamples.  Assuming nsamples = 100',verbosity=1)
        self.nsamples = 100
    
    npixels = self.NPIXELS_sampled
    if npixels is None: npixels = self.NPIXELS
    period = 1.0 / (2e6 / npixels / self.nsamples)
    return period

def timeline_npts(self):
    '''
    the size of the timeline determined from requested parameters.
    This is the number of points in the timeline vector

    *** NOTE: this method is not used anywhere ******************
    '''
    sample_period=self.sample_period()
    if sample_period is None:return None
    timeline_size = int(np.ceil(self.tinteg / sample_period))
    return timeline_size

def timeline_computertime(self,timeline_index=None):
    '''
    return the scientific time axis given by the computer date (not the GPS)
    '''
    if timeline_index is None\
       and self.hk is not None\
       and 'ASIC_SUMS' in self.hk.keys()\
       and 'ComputerDate' in self.hk['ASIC_SUMS'].keys():
        
        timestamps = self.hk['ASIC_SUMS']['ComputerDate']
        if len(timestamps)==0: return None
        return timestamps

    if timeline_index is None: timeline_index = 0    
    if 'DATE-OBS' in self.tdata[timeline_index].keys():
        timeline_date = self.tdata[timeline_index]['DATE-OBS']
        if len(timeline_date)==0: return None
        # fix weirdness
        t = timeline_date[0]
        utcoffset = t.timestamp() - dt.datetime.utcfromtimestamp(t.timestamp()).timestamp()
        timestamps = np.array( [(t.timestamp()+utcoffset) for t in timeline_date] )
        return timestamps

    self.printmsg('ERROR! No Computer Time.',verbosity=2)
    return None

def timeline_timeaxis(self,timeline_index=None,axistype='pps'):
    '''
    the timeline time axis for scientific data (TES)
    This is determined from the sample period and the number of points
    or, possibly given as a list of datetime
    '''
    self.printmsg('DEBUG: call to timeline_timeaxis with axistype=%s' % axistype,verbosity=4)
    
    if not self.exist_timeline_data():return None
    
    computertime = self.timeline_computertime(timeline_index)
    t0 = computertime[0]

    if timeline_index is None:timeline_index = 0

    timeline_npts = self.tdata[timeline_index]['TIMELINE'].shape[1]    
    sample_period = self.sample_period(timeline_index)
    if sample_period is None:
        time_axis_index = None
        time_axis_default = computertime
        default_descr = "computer time"
    else:
        time_axis_index = sample_period*np.arange(timeline_npts)
        time_axis_default = t0 + time_axis_index
        default_descr = "sample rate with start time given by the initial computer time"
    
    if axistype.lower()=='index':
        return time_axis_index

    if axistype.lower()=='pps':
        if 'ASIC_SUMS' in self.hk.keys():
            pps = self.pps(hk='ASIC_SUMS')
            gps = self.gps(hk='ASIC_SUMS')
            time_axis = pps2date(pps,gps,verbosity=self.verbosity)
            if time_axis is None:
                self.printmsg('ERROR! Using default based on %s' % default_descr,verbosity=2)
                return time_axis_default
            return time_axis
        self.printmsg('ERROR! No PPS data.  Using default based on %s' % default_descr,verbosity=2)
        return time_axis_default

    if axistype.lower()=='computertime':
        time_axis = computertime
        return time_axis

    self.printmsg('timeline_timeaxis returning time axis based on %s.' % default_descr,verbosity=2)
    return time_axis_default

def timeaxis(self,datatype=None,axistype='pps',asic=None,TES=None):
    '''
    wrapper to return the time axis for data.
    the datatypes are the various hk or scientific
    '''

    # valid axistypes in order of preference
    valid_axistypes = ['pps','timestamp','index','computertime']
    if axistype.lower() not in valid_axistypes:
        self.printmsg('Invalid axistype request.  Please choose one of: %s' % ', '.join(valid_axistypes))
        return None
    
    datatype = self.qubicstudio_filetype_truename(datatype)
    self.printmsg('timeaxis(): datatype=%s' % datatype,verbosity=3)
    
    if datatype is None: datatype = 'ASIC_SUMS'

    if datatype not in self.hk.keys():
        if self.__object_type__=='qubicfp':
            if asic is None:
                args = self.args_ok(TES,asic,asic_only_required=True)
                if args is None:
                    self.printmsg('Please enter a valid ASIC or TES number')
                    return None
                TES,asic = args
            asic_idx = asic - 1
            return self.asic_list[asic_idx].timeaxis(axistype=axistype,datatype=datatype)

        self.printmsg('No data for %s!' % datatype)
        return None

    if datatype=='ASIC_SUMS':
        self.printmsg('DEBUG: calling timeline_timeaxis from timeaxis() with axistype=%s' % axistype,verbosity=4)
        return self.timeline_timeaxis(axistype=axistype)

    # otherwise, return the time axis from one of the housekeeping sections
    tstart = self.hk[datatype]['ComputerDate'][0]
    tend   = self.hk[datatype]['ComputerDate'][-1]
    span = tend - tstart
    npts = len(self.hk[datatype]['ComputerDate'])
    tindex = (span/npts)*np.arange(npts)
    
    t_default = self.hk[datatype]['ComputerDate']
    t_default_str = 'computer time'

    tstamp = None
    if 'RaspberryDate' in self.hk[datatype].keys():
        tstamp_key = 'RaspberryDate'
    else:
        tstamp_key = 'timestamp'
    
    if tstamp_key in self.hk[datatype].keys():
        tstamp = self.hk[datatype][tstamp_key]
        t_default = tstamp
        t_default_str = 'timestamp'

    if axistype.lower()=='timestamp':
        if 'timestamp' not in self.hk[datatype].keys():
            self.printmsg('No timestamp.  Using %s instead' % t_default_str)
            return t_default            
        return tstamp

    if axistype.lower()=='index':
        return tindex

    if axistype.lower()=='computertime':
        return self.hk[datatype]['ComputerDate']

    # the only remaining option is pps
    pps = self.pps(hk=datatype)
    if pps is None:
        self.printmsg('No PPS.  Using %s instead' % t_default_str)
        return t_default
    if pps.max() == 0:
        self.printmsg('PPS is zero.  Using %s instead' % t_default_str)
        return t_default

    gps = self.gps(hk=datatype)
    if gps is None:
        self.printmsg('No GPS.  Using %s instead' % t_default_str)
        return t_default
    
    if gps.max() == 0.0:
        self.printmsg('GPS is zero.  Using %s instead' % t_default_str)
        return t_default

    pps_time = pps2date(pps,gps,verbosity=self.verbosity)
    if pps_time is None:
        self.printmsg('Could not get timeaxis from PPS/GPS.  Using %s instead' % t_default_str)
        return t_default

    return pps_time
    

def determine_bias_modulation(self,TES,timeline_index=None,timeaxis='pps'):
    '''
    determine the modulation of the bias voltage
    It should be close to the "bias_frequency" which is actually the period.
    '''
    if not self.exist_timeline_data():return None
    if not isinstance(timeline_index,int):timeline_index=0
    ntimelines=self.ntimelines()
    if timeline_index>=ntimelines:
        self.printmsg('Please enter a timeline between 0 and %i' % (ntimelines-1))
        return None

    retval = {}
    retval['TES'] = TES
    retval['timeline_index'] = timeline_index
    retval['timeaxis'] = timeaxis
    
    TES_idx=TES_index(TES)
    timeline=self.timeline(TES,timeline_index)
    timeline_npts=len(timeline)

    sample_period=self.sample_period(timeline_index)
    self.printmsg('DEBUG: calling timeline_timeaxis from determine_bias_modulation() with axistype=%s' % timeaxis, verbosity=4)
    time_axis=self.timeline_timeaxis(timeline_index,axistype=timeaxis)
    measured_sample_period = (time_axis[-1] - time_axis[0])/(timeline_npts-1)
    retval['measured_sample_period'] = measured_sample_period
    if sample_period is None:
        sample_period = measured_sample_period

    # use the bias_phase if it exists
    bias_phase = self.bias_phase()
    if bias_phase is not None:
        self.printmsg('getting bias variation from the saved data',verbosity=2)
        imin = np.argmin(bias_phase)
        imax = np.argmax(bias_phase)
        iperiod = 2*abs(imax-imin)
        if imin==imax:
            imin = 0
            imax = timeline_npts - 1
            iperiod = imax - imin
        ipeak0 = min([imin,imax])
        ipeak1 = ipeak0 + iperiod
        peak0 = time_axis[ipeak0]
        if ipeak1>=timeline_npts:
            peak1 = peak0 + iperiod*measured_sample_period
        else:
            peak1 = time_axis[ipeak1]
        self.bias_period = peak1-peak0
        retval['ipeak0'] = ipeak0
        retval['ipeak1'] = ipeak1
        retval['peak0'] = peak0
        retval['peak1'] = peak1
        return retval

    # the so-called frequency of the bias modulation is, in fact, the period
    # In QubicStudio the selection of bias_frequency=99 changes the significance from frequency to period (confusing)
    if self.bias_frequency is None:
        period_firstguess = 92. # based on experience
    else:
        period_firstguess = self.bias_frequency
        
    bias_period_npts=int(period_firstguess/sample_period) # should we be using measured_sample_period here?
    self.debugmsg('period npts = %i' % bias_period_npts)
        
    # skip the first few seconds which are often noisy
    skip=int(3.0/sample_period)
    self.debugmsg('looking for peaks in I-V timeline.  Skipping the first %i points.' % skip)
    peak0_range=(skip,skip+bias_period_npts)
    peak1_range_end=skip+2*bias_period_npts
    if peak1_range_end>=timeline_npts:
        peak1_range_end=timeline_npts-1
    peak1_range=(skip+bias_period_npts,peak1_range_end)


    # try to find the peaks, otherwise return ipeak0=0, ipeak1=timeline_npts-1
    try:
        ipeak0=np.argmax(timeline[peak0_range[0]:peak0_range[1]])
        ipeak0+=peak0_range[0]
    except:
        ipeak0=0
    peak0=time_axis[ipeak0]

    try:
        ipeak1=np.argmax(timeline[peak1_range[0]:peak1_range[1]])
        ipeak1+=peak1_range[0]
    except:
        ipeak1=timeline_npts-1
    peak1=time_axis[ipeak1]
    self.bias_period = peak1-peak0

    retval['ipeak0'] = ipeak0
    retval['ipeak1'] = ipeak1
    retval['peak0'] = peak0
    retval['peak1'] = peak1
    return retval

def plot_timeline(self,TES,timeline_index=None,fit=False,ipeak0=None,ipeak1=None,plot_bias=True,xwin=True,timeaxis='pps',ax=None,fontsize=12):
    '''
    plot the timeline
    '''
    if not self.exist_timeline_data():
        self.printmsg('ERROR! No timeline data.')
        return None

    if timeline_index is None:
        # by default, plot the first one.  For QubicStudio files, there is only one timeline
        timeline_index=0

    ntimelines=self.ntimelines()
    if timeline_index>=ntimelines:
        self.printmsg('Please enter a timeline between 0 and %i' % (ntimelines-1))
        return None

    tdata = self.tdata[timeline_index]
    keys = tdata.keys()
    
    warning_str = ''
    if 'WARNING' in keys and tdata['WARNING']:
        warning_str = '\n'.join(tdata['WARNING'])

    if 'R_FEEDBK' in keys:
        self.Rfeedback = tdata['R_FEEDBK']
        
    if 'NSAMPLES' in keys:
        self.nsamples = tdata['NSAMPLES']
        
    if 'DATE-OBS' in keys:
        timeline_date=tdata['DATE-OBS']
    else:
        timeline_date=self.obsdate

    if 'BEG-OBS' in keys:
        timeline_start=tdata['BEG-OBS']
    else:
        timeline_start=timeline_date

    if 'BIAS_MIN' in keys:
        self.min_bias=tdata['BIAS_MIN']
    if 'BIAS_MAX' in keys:
        self.max_bias=tdata['BIAS_MAX']

    biasphase = self.bias_phase()
        
    ttl=str('QUBIC Timeline curve for TES#%3i (%s)' % (TES,timeline_start.strftime('%Y-%b-%d %H:%M UTC')))

    if 'TES_TEMP' in keys:
        tempstr='%.0f mK' % (1000*tdata['TES_TEMP'])
    else:
        if self.temperature is None:
            tempstr='unknown'
        else:
            tempstr=str('%.0f mK' % (1000*self.temperature))

    fbstr = ''
    if 'R_FEEDBK' in keys:
        if tdata['R_HEATER']==1:
            onoff = 'ON'
        else:
            onoff = 'OFF'
        fbstr = ', Feedback Relay: %.0fk$\Omega$, Heater %s' % (tdata['R_FEEDBK']*1e-3,onoff)
        
    subttl=str('Array %s, ASIC #%i, Pixel #%i, Temperature %s%s' %
               (self.detector_name,self.asic,tes2pix(TES,self.asic),tempstr,fbstr))
                                
    if xwin: plt.ion()
    else: plt.ioff()

    if ax is None:
        newplot = True
        fig=plt.figure()
        figure_window_title(fig,ttl) 
        ax=plt.gca()
    else:
        newplot = False
        
    ax.set_xlabel('date UT',fontsize=fontsize)
    ax.set_ylabel('Current  /  $\mu$A',fontsize=fontsize)
    ax.tick_params(axis='both',labelsize=fontsize)
    if warning_str:
        boxprops = {}
        boxprops['alpha'] = 0.4
        boxprops['color'] = 'red'
        boxprops['boxstyle'] = 'round'
        ax.text(0.5,0.5,warning_str,ha='center',va='center',fontsize=2*fontsize,transform=ax.transAxes,bbox=boxprops)

    
    TES_idx=TES_index(TES)
    timeline=self.timeline(TES,timeline_index)
    current=self.ADU2I(timeline) # uAmps
    timeline_npts=len(timeline)

    self.printmsg('DEBUG: calling timeline_timeaxis from plot_timeline() with axistype=%s' % timeaxis,verbosity=4)
    time_axis=self.timeline_timeaxis(timeline_index,axistype=timeaxis)

    fitparms=None
    if fit:
        fitparms=self.fit_timeline(TES,timeline_index,ipeak0,ipeak1)
        

    ipeak0=0
    ipeak1=timeline_npts-1
    peak0=time_axis[ipeak0]
    peak1=time_axis[ipeak1]
    ysine = None
    if plot_bias:
        if self.timeline_conversion is None:
            self.timeline2adu(TES=TES,timeline_index=timeline_index,timeaxis=timeaxis)

        if self.timeline_conversion is None or self.min_bias is None or self.max_bias is None:
            plot_bias = False
        else:    
            ipeak0=self.timeline_conversion['ipeak0']
            ipeak1=self.timeline_conversion['ipeak1']
            peak0=self.timeline_conversion['peak0']
            peak1=self.timeline_conversion['peak1']
            shift=self.timeline_conversion['shift']

        if biasphase is not None:
            self.printmsg('DEBUG: taking ysine from QubicStudio FITS file',verbosity=4)
            ysine = self.timeline_vbias
            sinelabel = 'V$_\mathrm{bias}$ from QubicStudio FITS file'
        elif fitparms is None:
            self.printmsg('DEBUG: taking ysine from peak to peak',verbosity=4)
            bias_period=peak1-peak0
            amplitude=0.5*(self.max_bias-self.min_bias)
            offset=self.min_bias+amplitude
            sinelabel='sine curve period=%.2f seconds\npeaks determined from TES %i' % (bias_period,self.timeline_conversion['TES'])
            ysine=offset+amplitude*np.sin((time_axis-peak0)*2*np.pi/bias_period + 0.5*np.pi + shift*2*np.pi)
        else:
            self.printmsg('DEBUG: taking ysine from timeline fit to sine curve',verbosity=4)
            bias_period=fitparms['period']
            amplitude=fitparms['amplitude']
            offset=fitparms['offset']
            shift=fitparms['phaseshift']
            if bias_period is not None and amplitude is not None:
                sinelabel='best fit sine curve: period=%.2f seconds, amplitude=%.2f $\mu$A' % (bias_period,amplitude)
                ysine=self.model_timeline(time_axis,bias_period,shift,offset,amplitude)
    if ysine is None: plot_bias = False
        

    if newplot:
        fig.suptitle(ttl+'\n'+subttl,fontsize=fontsize)
    else:
        ax.text(0.5,1.0,ttl+'\n'+subttl,va='bottom',ha='center',fontsize=fontsize,transform=ax.transAxes)

    # plot date instead of timestamp
    time_axis_date = np.empty(len(time_axis),dtype=dt.datetime)
    for idx,tstamp in enumerate(time_axis):
        time_axis_date[idx] = dt.datetime.utcfromtimestamp(tstamp)
    peak0_date = time_axis_date[ipeak0]
    peak1_date = time_axis_date[ipeak1]
    
    curve1=ax.plot(time_axis_date,current,label='TES current',color='blue')

    #ymax=max([current[ipeak0],current[ipeak1]])
    ymax=np.nanmax(current)
    if np.isnan(ymax):
        ymax = 1.0
    ymin=np.nanmin(current)
    if np.isnan(ymin):
        ymin = 1.0
    yrange=ymax-ymin
    if np.isnan(yrange) or yrange==0:
        yrange = 0.1
    yminmax=(ymin-0.02*yrange,ymax+0.02*yrange)
    ax.set_ylim(yminmax)
    
    if plot_bias:
        ax.plot([peak0_date,peak0_date],yminmax,color='red',label='sine curve first peak')
        ax.plot([peak1_date,peak1_date],yminmax,color='red',label='sine curve second peak')
        if fitparms is None:
            ax_bias = ax.twinx()
            ax_bias.set_ylabel('Bias / V',rotation=270,va='bottom',fontsize=fontsize)
            if self.min_bias==self.max_bias:
                ax_bias.set_ylim([self.min_bias-1,self.max_bias+1])
            else:
                ax_bias.set_ylim([self.min_bias,self.max_bias])
            ax_bias.tick_params(axis='both',labelsize=fontsize)
            curve2_ax = ax_bias
        else:
            curve2_ax = ax
        self.printmsg('DEBUG: plotting sine curve for bias: len(time_axis)=%i, len(ysine)=%i'
                      % (len(time_axis),len(ysine)),verbosity=4)
        curve2 = curve2_ax.plot(time_axis_date,ysine,label=sinelabel,color='green')
        curves = curve1+curve2
    else:
        curves = curve1

    labs = [l.get_label() for l in curves]
    ax.legend(curves, labs, loc='upper right',fontsize=fontsize)

    pngname=str('TES%03i_array-%s_ASIC%i_timeline_%s.png' % (TES,self.detector_name,self.asic,timeline_start.strftime('%Y%m%dT%H%M%SUTC')))
    pngname_fullpath=self.output_filename(pngname)
    if newplot and isinstance(pngname_fullpath,str):
        plt.savefig(pngname_fullpath,format='png',dpi=100,bbox_inches='tight')
    if xwin:plt.show()
    else: plt.close('all')

    retval = {}
    retval['fitparms'] = fitparms
    retval['ax'] = ax
    retval['curves'] = curves
    retval['plotname'] = pngname
    return retval


def plot_timeline_physical_layout(self,
                                  timeline_index=None,
                                  xwin=True,
                                  imin=None,
                                  imax=None,
                                  tmin=None,
                                  tmax=None,
                                  lutmin=None,
                                  lutmax=None):
    '''
    plot the timeline curves in thumbnails mapped to the physical location of each detector
    '''
    TES2PIX = assign_pix2tes(self.obsdate)
    
    if not self.exist_timeline_data():return None
    ntimelines=self.ntimelines()
    
    if timeline_index is None:
        # by default, plot the first one.
        timeline_index=0
    
    if timeline_index>=ntimelines:
        self.printmsg('Please enter a timeline between 0 and %i' % (ntimelines-1))
        return None

    tdata = self.tdata[timeline_index]
    keys = tdata.keys()    
    timeline_npts = tdata['TIMELINE'].shape[1]
    if lutmax is None:
        lutmax = tdata['TIMELINE'].max() - tdata['TIMELINE'].min()
    if lutmin is None:
        lutmin = 0.0
    

    if 'DATE-OBS' in keys:
        timeline_date=tdata['DATE-OBS']
    else:
        timeline_date=self.obsdate

    if 'BEG-OBS' in keys:
        timeline_start=tdata['BEG-OBS']
    else:
        timeline_start=timeline_date

    ttl=str('QUBIC Timeline curves (%s)' % (timeline_start.strftime('%Y-%b-%d %H:%M UTC')))

    if 'TES_TEMP' in keys:
        tempstr='%.0f mK' % (1000*tdata['TES_TEMP'])
    else:
        if self.temperature is None:
            tempstr='unknown'
        else:
            tempstr=str('%.0f mK' % (1000*self.temperature))
    subttl=str('Array %s, ASIC #%i, T$_\mathrm{bath}$=%s' % (self.detector_name,self.asic,tempstr))

    # use the plot_fp algorithm to plot the focal plane
    asic_key = 'ASIC%i' % self.asic
    args = {}
    args['title'] = ttl
    args['subtitle'] = subttl
    
    pngname=str('QUBIC_Array-%s_ASIC%i_timeline_%s.png' % (self.detector_name,self.asic,timeline_start.strftime('%Y%m%dT%H%M%SUTC')))
    pngname_fullpath=self.output_filename(pngname)
    args['pngname'] = pngname_fullpath


    # plot subsection of timeline
    tlim=[0,timeline_npts]
    if tmin is None:
        tlim[0] = 0
    else:
        tlim[0] = tmin
    if tmax is None:
        tlim[1] = timeline_npts
    else:
        tlim[1] = tmax    
    args[asic_key] = self.timeline_array(timeline_index=timeline_index)[:,tlim[0]:tlim[1]]

    plot_fp(args)

    return args


def timeline2adu(self,TES=None,ipeak0=None,ipeak1=None,timeline_index=0,shift=0.0,timeaxis='pps'):
    '''
    transfer timeline data with I-V curves to the ADU matrix 
    this is done so that we can directly use all the I-V methods
    '''
    self.printmsg('timeline2adu()',verbosity=4)
    if not self.exist_timeline_data():return None
    ntimelines=self.ntimelines()
    if timeline_index>=ntimelines:
        self.printmsg('Please enter a timeline between 0 and %i' % (ntimelines-1))
        return None

    if not isinstance(TES,int):
        self.printmsg('Please enter a TES which is the reference for extracting the bias timeline')
        return None

    biasmod = self.determine_bias_modulation(TES,timeline_index,timeaxis=timeaxis)
    ip0 = biasmod['ipeak0']
    ip1 = biasmod['ipeak1']
    peak0 = biasmod['peak0']
    peak1 = biasmod['peak1']
    
    if not isinstance(ipeak0,int):ipeak0=ip0
    if not isinstance(ipeak1,int):ipeak1=ip1
    self.debugmsg('timeline2adu: ipeak0=%i' % ipeak0)
    self.debugmsg('timeline2adu: ipeak1=%i' % ipeak1)
    self.timeline_conversion={}
    self.timeline_conversion['ipeak0']=ipeak0
    self.timeline_conversion['ipeak1']=ipeak1
    self.timeline_conversion['TES']=TES
    self.timeline_conversion['timeline_index']=timeline_index
    self.timeline_conversion['shift']=shift
    
    self.printmsg('DEBUG: calling timeline_timeaxis from timeline2adu() with axistype=%s' % timeaxis,verbosity=4)
    time_axis=self.timeline_timeaxis(timeline_index,axistype=timeaxis)
    peak0=time_axis[ipeak0]
    if ipeak1<len(time_axis):
        peak1=time_axis[ipeak1]
    self.timeline_conversion['peak0']=peak0
    self.timeline_conversion['peak1']=peak1
    bias_period=peak1-peak0
    self.timeline_conversion['bias_period']=bias_period


    # find the number of bias cycles (one cycle is back to the original Vbias)
    # one full bias_period is a return to the same bias voltage, so one period is 1 cycle
    if bias_period==0.0:
        ncycles = 0
    else:
        ncycles = int( (peak1 - peak0)/bias_period )
    if ncycles==0:
        # it's not down/up
        self.nbiascycles=1
        self.cycle_vbias=False
    else:
        self.nbiascycles=ncycles
        self.cycle_vbias=True    


    if (self.max_bias is None) or (self.min_bias is None):
        self.printmsg('No Bias Voltage information.',verbosity=2)
        return False
    
    amplitude=0.5*(self.max_bias-self.min_bias)
    offset=self.min_bias+amplitude

    # if bias phase was saved by QubicStudio then we use the vbias calculated by bias_phase()
    biasphase = self.bias_phase()
    if biasphase is None:
        # the last term is if we're applying a shift in terms of period
        ysine = offset+amplitude*np.sin((time_axis-peak0)*2*np.pi/bias_period + 0.5*np.pi + shift*2*np.pi)
        self.printmsg('DEBUG: timeline2adu() : setting vbias to derived sine curve',verbosity=4)
        self.vbias = ysine[ipeak0:ipeak1]
    else:
        self.printmsg('DEBUG: timeline2adu() : setting vbias to subset between peaks',verbosity=4)
        self.vbias = self.timeline_vbias[ipeak0:ipeak1]
        
        
    self.min_bias = np.nanmin(self.vbias)
    self.max_bias = np.nanmax(self.vbias)

    tdata = self.tdata[timeline_index]
    keys = tdata.keys()
    if 'BEG-OBS' in keys:
        self.obsdate=tdata['BEG-OBS']
    else:
        self.obsdate=tdata['DATE-OBS']
    if 'TES_TEMP' in keys:
        self.temperature=tdata['TES_TEMP']
    else:
        self.temperature=None

    if 'R_FEEDBK' in keys:
        self.Rfeedback = tdata['R_FEEDBK']
        
    if 'NSAMPLES' in keys:
        self.nsamples = tdata['NSAMPLES']
        
        
    npts=len(self.vbias)
    self.adu=np.empty((self.NPIXELS,npts))
    for idx in range(self.NPIXELS):
        TES=idx+1
        self.adu[idx,:]=self.timeline(TES,timeline_index)[ipeak0:ipeak1]
        

    return True


def model_timeline(self,t,period,phaseshift,offset,amplitude):
    '''
    a sine function to fit to the timeline data
    '''
    ysine=offset + amplitude*np.sin( 2*np.pi * (t/period + phaseshift) )
    return ysine
    
def fit_timeline(self,TES,timeline_index=None,ipeak0=None,ipeak1=None,timeaxis='pps'):
    '''
    fit the timeline to a sine curve
    '''
    # return a dictionary
    fit={}
    fit['TES']=TES
    fit['DET_NAME']=self.detector_name
    fit['ASIC']=self.asic
    
    if timeline_index is None:timeline_index=0    
    ntimelines=self.ntimelines()
    if timeline_index>=ntimelines:
        self.printmsg('Please enter a timeline between 0 and %i' % (ntimelines-1))
        return None
    fit['timeline_index']=timeline_index
    fit['date']=self.tdata[timeline_index]['DATE-OBS']
    fit['Tbath']=self.tdata[timeline_index]['TES_TEMP']
    
    TES_idx=TES_index(TES)
    timeline=self.timeline(TES,timeline_index)
    current=self.ADU2I(timeline)
    timeline_npts=len(timeline)

    self.printmsg('DEBUG: calling timeline_timeaxis from fit_timeline() with axistype=%s' % timeaxis,verbosity=4)
    time_axis=self.timeline_timeaxis(timeline_index,axistype=timeaxis)

    # first guess;  use the peak search algorithm
    biasmod = self.determine_bias_modulation(TES,timeline_index,timeaxis=timeaxis)
    if biasmod is None:
        self.printmsg('ERROR! Could not determine bias modulation.',verbosity=2)
        return None
    fit['ipeak0'] = biasmod['ipeak0']
    fit['ipeak1'] = biasmod['ipeak1']
    fit['peak0'] = biasmod['peak0']
    fit['peak1'] = biasmod['peak1']
    self.timeline_conversion['ipeak0']=ipeak0 #MP
    self.timeline_conversion['ipeak1']=ipeak1 #MP
    self.timeline_conversion['peak0']=peak0 #MP
    self.timeline_conversion['peak1']=peak1 #MP

    peak0 = biasmod['peak0']
    peak1 = biasmod['peak1']
    period = peak1-peak0
    amplitude = 0.5*(max(current)-min(current))
    offset = min(current)+amplitude
    phaseshift = peak0/period
    
    p0=[period,phaseshift,offset,amplitude]

    try:
        popt,pcov=curve_fit(self.model_timeline,time_axis,current,p0=p0)
    except:
        fit['period'] = None
        fit['phaseshift'] = None
        fit['offset'] = None
        fit['amplitude'] = None
        fit['R amplitude'] = None
        return fit

        
    period,phaseshift,offset,amplitude=popt
    fit['period']=period
    fit['phaseshift']=phaseshift
    fit['offset']=offset
    fit['amplitude']=amplitude # this is in microAmps

    Vtes=self.Rshunt*( (self.max_bias*self.bias_factor)/self.Rbias - 1e-6*abs(amplitude) )    
    fit['R amplitude']=abs(Vtes/amplitude)
    return fit
