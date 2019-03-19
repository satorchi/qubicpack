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
from __future__ import division, print_function
import numpy as np
import sys,os,time
import datetime as dt
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def exist_timeline_data(self):
    '''
    check if we have timeline data
    '''
    if self.tdata is None:
        return False
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
    '''
    if not self.exist_timeline_data():return None
    ntimelines=self.ntimelines()
    if timeline_index>=ntimelines:
        print('ERROR! timeline index out of range.  Enter an index between 0 and %i' % (ntimelines-1))
        return None
    TES_index=self.TES_index(TES)
    if TES_index is None:return None
    timeline=self.tdata[timeline_index]['TIMELINE'][TES_index,:]
    return timeline

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
        print('WARNING! Cannot set bias offset greater than %.3fV.' % max_offset)
        if bias<0.0:
            bias=-max_offset
        else:
            bias= max_offset
        print('Setting bias to %.3fV' % bias)


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

    if timeline_index is not None and 'NSAMPLES' in self.tdata[timeline_index].keys():
        self.nsamples = self.tdata[timeline_index]['NSAMPLES']        
    if self.nsamples is None:return None
    
    npixels=self.NPIXELS_sampled
    if npixels is None:npixels=self.NPIXELS
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

def timeline_timeaxis(self,timeline_index=None,axistype='pps'):
    '''
    the timeline time axis.  
    This is determined from the sample period and the number of points
    or, possibly given as a list of datetime
    '''
    if not self.exist_timeline_data():return None
    if timeline_index is None:timeline_index = 0

    if 'DATE-OBS' in self.tdata[timeline_index].keys():
        timeline_date=self.tdata[timeline_index]['DATE-OBS']
    else:
        timeline_date = self.obsdate

    timeline_npts = self.tdata[timeline_index]['TIMELINE'].shape[1]    
    sample_period = self.sample_period(timeline_index)
    time_axis_index = sample_period*np.arange(timeline_npts)

    if axistype.lower()=='index':
        return time_axis_index

    if axistype.lower()=='pps':
        if 'ASIC_SUMS' in self.hk.keys():
            pps = self.pps(hk='ASIC_SUMS')
            gps = self.gps(hk='ASIC_SUMS')
            if gps is None or gps.min()<1494486000 or gps.max()<1494486000:
                print('ERROR! Bad GPS data.  Using sample rate to make a time axis.')
                return time_axis_index
            time_axis = self.pps2date(pps,gps)
            return time_axis
        print('ERROR! No PPS data.')
        return time_axis_index

    if axistype.lower()=='computertime':
        if isinstance(timeline_date,list) and isinstance(timeline_date[0],dt.datetime):
            time_axis = np.array([ eval(t.strftime('%s.%f')) for t in timeline_date ])
            t0 = time_axis[0]
            return time_axis
        print('ERROR! No Computer Time.')
        return time_axis_index

    return time_axis_index

def timeaxis(self,datatype=None,axistype='pps'):
    '''
    wrapper to return the time axis for data.
    the two datatypes are: hk or scientific
    '''
    if datatype is None: datatype = 'sci'
    if datatype.lower()[0:3]=='sci':
        return self.timeline_timeaxis(axistype=axistype)

    # otherwise, return the platform time axis
    hktype = 'INTERN_HK'
    pps = self.pps(hk=hktype)
    if pps is None: return None    
    gps = self.gps(hk=hktype)
    if gps is None: return None

    tstamp = pps2date(pps,gps)
    if axistype=='index':
        return tstamp - tstamp[0]
    return tstamp
    

def determine_bias_modulation(self,TES,timeline_index=None):
    '''
    determine the modulation of the bias voltage
    It should be close to the "bias_frequency" which is actually the period.
    '''
    if not self.exist_timeline_data():return None
    if not isinstance(timeline_index,int):timeline_index=0
    ntimelines=self.ntimelines()
    if timeline_index>=ntimelines:
        print('Please enter a timeline between 0 and %i' % (ntimelines-1))
        return None
    
    TES_index=self.TES_index(TES)
    timeline=self.timeline(TES,timeline_index)
    timeline_npts=len(timeline)

    sample_period=self.sample_period(timeline_index)
    time_axis=self.timeline_timeaxis(timeline_index)

    # the so-called frequency of the bias modulation is, in fact, the period
    # In QubicStudio the selection of bias_frequency=99 changes the significance from frequency to period (confusing)
    if self.bias_frequency is None:
        period_firstguess = 92. # based on experience
    else:
        period_firstguess = self.bias_frequency
        
    bias_period_npts=int(period_firstguess/sample_period)
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
    return (ipeak0,ipeak1)

def plot_timeline(self,TES,timeline_index=None,fit=False,ipeak0=None,ipeak1=None,plot_sine=True,xwin=True,timeaxis='pps'):
    '''
    plot the timeline
    '''
    if not self.exist_timeline_data():
        print('ERROR! No timeline data.')
        return None

    if timeline_index is None:
        # by default, plot the first one.  For QubicStudio files, there is only one timeline
        timeline_index=0

    ntimelines=self.ntimelines()
    if timeline_index>=ntimelines:
        print('Please enter a timeline between 0 and %i' % (ntimelines-1))
        return None

    tdata = self.tdata[timeline_index]
    keys = tdata.keys()
    
    warning_str = ''
    if 'WARNING' in keys and tdata['WARNING']:
        warning_str = '\n'.join(tdata['WARNING'])
                                 
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

    
        
    ttl=str('QUBIC Timeline curve for TES#%3i (%s)' % (TES,timeline_start.strftime('%Y-%b-%d %H:%M UTC')))

    if 'TES_TEMP' in keys:
        tempstr='%.0f mK' % (1000*tdata['TES_TEMP'])
    else:
        if self.temperature is None:
            tempstr='unknown'
        else:
            tempstr=str('%.0f mK' % (1000*self.temperature))
    subttl=str('Array %s, ASIC #%i, Pixel #%i, Temperature %s' % (self.detector_name,self.asic,self.tes2pix(TES),tempstr))
    if warning_str: subttl += '\n'+warning_str
                                
    if xwin: plt.ion()
    else: plt.ioff()

    fig=plt.figure(figsize=self.figsize)
    fig.canvas.set_window_title('plt: '+ttl) 
    ax=plt.gca()
    ax.set_xlabel('time  /  seconds')
    ax.set_ylabel('Current  /  $\mu$A')
    
    TES_index=self.TES_index(TES)
    timeline=self.timeline(TES,timeline_index)
    current=self.ADU2I(timeline) # uAmps
    timeline_npts=len(timeline)

    time_axis=self.timeline_timeaxis(timeline_index,axistype=timeaxis)

    fitparms=None
    if fit:
        fitparms=self.fit_timeline(TES,timeline_index,ipeak0,ipeak1)

    ipeak0=0
    ipeak1=timeline_npts-1
    if plot_sine:
        if self.timeline_conversion is None:
            try:
                self.timeline2adu(TES=TES,timeline_index=timeline_index)
                ipeak0=self.timeline_conversion['ipeak0']
                ipeak1=self.timeline_conversion['ipeak1']
                shift=self.timeline_conversion['shift']
            except:
                plot_sine=False
        else:
            ipeak0=self.timeline_conversion['ipeak0']
            ipeak1=self.timeline_conversion['ipeak1']
            shift=self.timeline_conversion['shift']

            
    peak0=time_axis[ipeak0]
    peak1=time_axis[ipeak1]

    if plot_sine:
        if fitparms is None:
            bias_period=peak1-peak0
            amplitude=0.5*(self.max_bias-self.min_bias)
            offset=self.min_bias+amplitude
            sinelabel='sine curve period=%.2f seconds\npeaks determined from TES %i' % (bias_period,self.timeline_conversion['TES'])
            ysine=offset+amplitude*np.sin((time_axis-peak0)*2*np.pi/bias_period + 0.5*np.pi + shift*2*np.pi)
        else:
            bias_period=fitparms['period']
            amplitude=fitparms['amplitude']
            offset=fitparms['offset']
            shift=fitparms['phaseshift']
            sinelabel='best fit sine curve: period=%.2f seconds, amplitude=%.2f $\mu$A' % (bias_period,amplitude)
            ysine=self.model_timeline(time_axis,bias_period,shift,offset,amplitude)
        

    fig.suptitle(ttl+'\n'+subttl,fontsize=16)
    
    curve1=ax.plot(time_axis,current,label='I-V timeline',color='blue')

    #ymax=max([current[ipeak0],current[ipeak1]])
    ymax=max(current)
    ymin=min(current)
    yrange=ymax-ymin
    yminmax=(ymin-0.02*yrange,ymax+0.02*yrange)
    ax.plot([peak0,peak0],yminmax,color='red')
    ax.plot([peak1,peak1],yminmax,color='red')
    ax.set_ylim(yminmax)

    if plot_sine:
        if fitparms is None:
            ax_bias = ax.twinx()
            ax_bias.set_ylabel('Bias / V',rotation=270,va='bottom')
            ax_bias.set_ylim([self.min_bias,self.max_bias])
            curve2_ax=ax_bias
        else:
            curve2_ax=ax        
        curve2=curve2_ax.plot(time_axis,ysine,label=sinelabel,color='green')
        curves = curve1+curve2
    else:
        curves = curve1
        
    labs = [l.get_label() for l in curves]
    ax.legend(curves, labs, loc=0)

    pngname=str('TES%03i_array-%s_ASIC%i_timeline_%s.png' % (TES,self.detector_name,self.asic,timeline_start.strftime('%Y%m%dT%H%M%SUTC')))
    pngname_fullpath=self.output_filename(pngname)
    if isinstance(pngname_fullpath,str): plt.savefig(pngname_fullpath,format='png',dpi=100,bbox_inches='tight')
    if xwin:plt.show()
    else: plt.close('all')
    
    return fitparms


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
    if not self.exist_timeline_data():return None
    ntimelines=self.ntimelines()
    
    if timeline_index is None:
        # by default, plot the first one.
        timeline_index=0
    
    if timeline_index>=ntimelines:
        print('Please enter a timeline between 0 and %i' % (ntimelines-1))
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

    nrows=self.pix_grid.shape[0]
    ncols=self.pix_grid.shape[1]

    if xwin: plt.ion()
    else: plt.ioff()
    # need a square figure for this plot to look right
    figlen=max(self.figsize)
    fig,ax=plt.subplots(nrows,ncols,figsize=[figlen,figlen])
    pngname=str('QUBIC_Array-%s_ASIC%i_timeline_%s.png' % (self.detector_name,self.asic,timeline_start.strftime('%Y%m%dT%H%M%SUTC')))
    pngname_fullpath=self.output_filename(pngname)
    if xwin: fig.canvas.set_window_title('plt:  '+ttl)
    fig.suptitle(ttl+'\n'+subttl,fontsize=16)
    

    # the pixel number is between 1 and 248
    TES_translation_table=self.TES2PIX[self.asic_index()]
    ilim=[None,None]
    tlim=[0,timeline_npts]
    
    text_y=0.0
    text_x=1.0
    for row in range(nrows):
        for col in range(ncols):
            ax[row,col].get_xaxis().set_visible(False)
            ax[row,col].get_yaxis().set_visible(False)

            # the pixel identity associated with its physical location in the array
            physpix=self.pix_grid[row,col]
            pix_index=physpix-1
            self.debugmsg('processing PIX %i' % physpix)

            if physpix==0:
                pix_label = 'EMPTY'
                label_colour = 'black'
                face_colour = 'black'
            elif physpix in TES_translation_table:
                TES = self.pix2tes(physpix)
                pix_label = str('%i' % TES)
                label_colour = 'black'
                face_colour = 'white'
                TES_index = self.TES_index(TES)
                timeline = self.timeline(TES,timeline_index)
                I = self.ADU2I(timeline)
                self.debugmsg('plotting TES %i' % TES)
                ax[row,col].plot(I,color='blue')

                # plot subsection of timeline
                if tmin is None:
                    tlim[0] = 0
                else:
                    tlim[0] = tmin
                if tmax is None:
                    tlim[1] = timeline_npts
                else:
                    tlim[1] = tmax
                ax[row,col].set_xlim(tlim)

                # get min/max from timeline window
                negpeak = min(I[ tlim[0]:tlim[1] ])
                pospeak = max(I[ tlim[0]:tlim[1] ])
                if imin is None:
                    ilim[0] = negpeak
                else:
                    ilim[0] = imin
                if imax is None:
                    ilim[1] = pospeak
                else:
                    ilim[1] = imax
                ax[row,col].set_ylim(ilim)

                # get face colour from peak-to-peak
                peak2peak = pospeak - negpeak
                #face_colour = self.lut(peak2peak,lutmin,lutmax)
                              
                

            else:
                pix_label='other\nASIC'
                label_colour='yellow'
                face_colour='blue'

            ax[row,col].set_facecolor(face_colour)
            ax[row,col].text(text_x,text_y,pix_label,va='bottom',ha='right',
                             color=label_colour,fontsize=8,transform=ax[row,col].transAxes)
            
    if isinstance(pngname_fullpath,str): plt.savefig(pngname_fullpath,format='png',dpi=100,bbox_inches='tight')
    if xwin: plt.show()
    else: plt.close('all')

    return


def timeline2adu(self,TES=None,ipeak0=None,ipeak1=None,timeline_index=0,shift=0.0):
    '''
    transfer timeline data with I-V curves to the ADU matrix 
    this is done so that we can directly use all the I-V methods
    '''
    if not self.exist_timeline_data():return None
    ntimelines=self.ntimelines()
    if timeline_index>=ntimelines:
        print('Please enter a timeline between 0 and %i' % (ntimelines-1))
        return None

    if not isinstance(TES,int):
        print('Please enter a TES which is the reference for extracting the bias timeline')
        return None
    
    ip0,ip1=self.determine_bias_modulation(TES,timeline_index)
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
    
    time_axis=self.timeline_timeaxis(timeline_index)
    peak0=time_axis[ipeak0]
    peak1=time_axis[ipeak1]
    bias_period=peak1-peak0

    # find the number of cycles from the bias modulation "frequency" which is the period
    ### what's this? it doesn't make sense.
    # ncycles=int( np.round((peak1-peak0)/self.bias_frequency) )
    ncycles = int( (time_axis[-1] - time_axis[0])/bias_period )
    if ncycles==0:
        # it's not down/up
        self.nbiascycles=1
        self.cycle_vbias=False
    else:
        self.nbiascycles=ncycles
        self.cycle_vbias=True    
    
    amplitude=0.5*(self.max_bias-self.min_bias)
    offset=self.min_bias+amplitude

    # the last term is if we're applying a shift in terms of period
    ysine=offset+amplitude*np.sin((time_axis-peak0)*2*np.pi/bias_period + 0.5*np.pi + shift*2*np.pi)
    self.vbias=ysine[ipeak0:ipeak1]
    self.min_bias=min(self.vbias)
    self.max_bias=max(self.vbias)
    if 'BEG-OBS' in self.tdata[timeline_index].keys():
        self.obsdate=self.tdata[timeline_index]['BEG-OBS']
    else:
        self.obsdate=self.tdata[timeline_index]['DATE-OBS']
    if 'TES_TEMP' in self.tdata[timeline_index].keys():
        self.temperature=self.tdata[timeline_index]['TES_TEMP']
    else:
        self.temperature=None
    
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
    
def fit_timeline(self,TES,timeline_index=None,ipeak0=None,ipeak1=None):
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
        print('Please enter a timeline between 0 and %i' % (ntimelines-1))
        return None
    fit['timeline_index']=timeline_index
    fit['date']=self.tdata[timeline_index]['DATE-OBS']
    fit['Tbath']=self.tdata[timeline_index]['TES_TEMP']
    
    TES_index=self.TES_index(TES)
    timeline=self.timeline(TES,timeline_index)
    current=self.ADU2I(timeline)
    timeline_npts=len(timeline)

    time_axis=self.timeline_timeaxis(timeline_index)

    # first guess;  use the peak search algorithm
    i0,i1=self.determine_bias_modulation(TES,timeline_index)
    if ipeak0 is None:ipeak0=i0
    if ipeak1 is None:ipeak1=i1
    peak0=time_axis[ipeak0]
    peak1=time_axis[ipeak1]
    period=peak1-peak0
    amplitude=0.5*(max(current)-min(current))
    offset=min(current)+amplitude
    phaseshift=peak0/period
    
    p0=[period,phaseshift,offset,amplitude]

    popt,pcov=curve_fit(self.model_timeline,time_axis,current,p0=p0)
    period,phaseshift,offset,amplitude=popt
    fit['period']=period
    fit['phaseshift']=phaseshift
    fit['offset']=offset
    fit['amplitude']=amplitude # this is in microAmps

    Vtes=self.Rshunt*( (self.max_bias*self.bias_factor)/self.Rbias - 1e-6*abs(amplitude) )    
    fit['R amplitude']=abs(Vtes/amplitude)
    return fit
