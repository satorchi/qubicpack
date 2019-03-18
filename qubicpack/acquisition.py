'''
$Id: acquisition.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Mon 07 Aug 2017 07:35:33 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

acquisition methods.  These require a connection to QubicStudio

     the following methods are originally by Michel Piat:
      set_VoffsetTES()
      set_slowDAC() [removed from QubicStudio v3]
      set_diffDAC() [removed from QubicStudio v3]
      get_amplitude()
      get_mean()
      integrate_scientific_data()


Note on the use of sendSetCalibPolar()
name change in QubicStudio v3: sendSetCalibPolar -> sendSetTESDAC

see comments by Pierre Chanial in dispatcheraccess.pyx
translated here...

sendSetCalibPolar(int asic, int mode, int shape, int frequency, int amplitude, int offset)

this sets the shape of the Bias Voltage.  
The French for "Bias Voltage" is "Polarisation", which is very confusing.
and then the word is further truncated to "Polar"

arguments:
  asic:  
       asic=0xFF : all ASIC's receive the command
       asic<16 : command sent to the specified ASIC
       asic bits 8 to 23 : send to a list of ASICs
                           example: 0x00FF00 corresponds to ASIC 0 to 7
                           example: 0x0AD200 corresponds to ASIC 0,3,6,9

  mode: (this was deleted in QubicStudio v3)
       mode=0 : no signal
       mode=1 : apply Bias voltage

  shape: 
       shape=0 : sinus
       shape=1 : triangle
       shape=2 : continuous

  frequency:
       frequency=99 defaults to lowest possible value.  i.e. f=1Hz


'''
from __future__ import division, print_function
import numpy as np
import pystudio
import sys,os,time
import datetime as dt
import matplotlib.pyplot as plt

def connect_QubicStudio(self,client=None, ip=None, silent=False):
    if ip is None:
        ip=self.QubicStudio_ip
    else:
        self.QubicStudio_ip=ip
    
    if client is None:
        client = pystudio.get_client()

    if client is None:
        if not silent:print("connecting to QubicStudio on host: ",self.QubicStudio_ip)
        client = pystudio.DispatcherAccess(self.QubicStudio_ip, 3002)
        if not silent:print('wait 3 seconds before doing anything')
        time.sleep(3)

    if not client.connected:
        if not silent:print("ERROR: could not connect to QubicStudio")
        return None

    client.waitingForAckMode = True
    # client.sendSetScientificDataTfUsed(1) # data in Volts
    return client

def verify_QS_connection(self):
    '''
    verify that we have a valid connection to QubicStudio
    This is useful for verifying that we're trying to connect to the correct ASIC
    '''
    client = self.connect_QubicStudio()
    if client is None:return False

    parameter = 'QUBIC_PixelScientificDataTimeLine_%i' % self.QS_asic_index
    req = client.request(parameter)

    try:
        timeline=req.next()
    except Exception as e:
        msg= '\n*****************************************************************************'
        msg+='\nERROR! Could not get data from QubicStudio.  Did you choose the correct ASIC?'
        msg+='\n       ASIC=%i' % self.asic
        msg+='\n       %s' % e
        msg+='\n*****************************************************************************\n'
        print(msg)
        return False

    return True

def get_NPIXELS(self):
    '''
    get the number of pixels currently being samples
    '''
    if self.AVOID_HANGUP:
        # HACK: do not request NPIXELS again if we already have it.
        if not self.NPIXELS_sampled is None: return self.NPIXELS_sampled

    client = self.connect_QubicStudio()
    if client is None:return None

    # flush the request queue just in case
    #q=client.abort_requests()
    
    self.debugmsg('getting NPIXELS...')
    NPIXELS_all = client.fetch('ASIC_NbPixelSelected')
    self.debugmsg('got NPIXELS.')
    NPIXELS=NPIXELS_all[self.QS_asic_index]
    self.debugmsg('NPIXELS=%i' % NPIXELS)
    self.NPIXELS_sampled=NPIXELS
    self.NPIXELS_requested=True
    return NPIXELS

def get_PID(self):
    '''
    get the current Flux Lock Loop P,I,D and state
    '''
    client = self.connect_QubicStudio()
    if client is None:return None

    if not self.AVOID_HANGUP or self.FLL_state is None:
        self.debugmsg('getting FLL State...')
        req=client.fetch('QUBIC_FLL_State')
        self.FLL_state=req[self.QS_asic_index]
    
    if not self.AVOID_HANGUP or self.FLL_P is None:
        self.debugmsg('getting FLL P...')
        req=client.fetch('QUBIC_FLL_P')
        self.FLL_P=req[self.QS_asic_index]
    
    if not self.AVOID_HANGUP or self.FLL_I is None:
        self.debugmsg('getting FLL I...')
        req=client.fetch('QUBIC_FLL_I')
        self.FLL_I=req[self.QS_asic_index]
    
    if not self.AVOID_HANGUP or self.FLL_D is None:
        self.debugmsg('getting FLL D...')
        req=client.fetch('QUBIC_FLL_D')
        self.FLL_D=req[self.QS_asic_index]

    return self.FLL_state,self.FLL_P,self.FLL_I,self.FLL_D


def configure_PID(self,P=0,I=20,D=0):
    '''
    configure the FLL (Flux Lock Loop) PID
    '''
    client = self.connect_QubicStudio()
    if client is None:return False

    # first switch off the loop
    client.sendActivatePID(self.QS_asic_index,0)

    # wait a bit
    self.wait_a_bit()

    # next, set the parameters
    client.sendConfigurePID(self.QS_asic_index, P,I,D)

    # wait a bit
    self.wait_a_bit()

    # and reactivate the loop
    client.sendActivatePID(self.QS_asic_index,1)

    # wait a bit before returning
    self.wait_a_bit()

    return True


def set_Rfeedback(self,Rfeedback):
    '''
    set the resistance for the feedback loop
    it can be 10kOhm or 100kOhm.

    for the QubicStudio command:
          val=1 -> 10kOhm
          val=0 -> 100kOhm

    '''
    client = self.connect_QubicStudio()
    if client is None:return None

    # by default, we set it to 10kOhm.
    if Rfeedback<100:
        client.sendSetFeedbackRelay(self.QS_asic_index,1)
        self.Rfeedback=10e3
    else:
        client.sendSetFeedbackRelay(self.QS_asic_index,0)
        self.Rfeedback=100e3

    return self.Rfeedback

def compute_offsets(self,count=10,consigne=0.0):
    '''
    measure the offsets and upload them to the table for future use
    '''
    client = self.connect_QubicStudio()
    if client is None:return False

    # first switch off the loop
    client.sendActivatePID(self.QS_asic_index,0)

    # make sure relay=10kOhm
    self.set_Rfeedback(10)
    
    # set sampling frequency 400Hz
    freq=400.
    # set sampling amplitude 1V
    amplitude=1.0
    # set sampling offset 6V
    bias=6.0
    # set shape to sinus
    shape=0
    if not self.set_VoffsetTES(bias, amplitude, freq, shape):return False

    # to begin, assign zero offset
    offsets = np.zeros(self.NPIXELS)
    client.sendSetOffsetTable(self.QS_asic_index, offsets)
    self.wait_a_bit()

    # to begin, get the current offsets
    #parameter='QUBIC_OffsetDACValues_%i' % self.QS_asic_index
    #offsets=client.fetch(parameter)

    k=1.0 # the first step is big
    for counter in range(count):

        print('count: %i/%i' % (counter+1,count))
        timeline = self.integrate_scientific_data(save=False)
        this_data_avg=timeline.mean(axis=-1)
        prev_offsets=offsets
        offsets=-k*(this_data_avg-consigne)+prev_offsets
        client.sendSetOffsetTable(self.QS_asic_index, offsets)
        self.wait_a_bit()
        k=0.5 # and subsequent steps are smaller
    return True

def feedback_offsets(self,count=10,consigne=0.0,Rfb=10,I=10,k=1.0): #MP
    '''
    measure the feedback offsets and upload them to the table for future use
    '''
    client = self.connect_QubicStudio()
    if client is None:return False

    ## switch off the feedback loop
    client.sendActivatePID(self.QS_asic_index,0)

    # make sure relay=Rfb #MP
    self.set_Rfeedback(Rfb)

    # set sampling frequency 10Hz
    freq=10.
    # set sampling amplitude 0.0V
    amplitude=0.0
    # set sampling offset 6V
    bias=6.0
    # set shape to sinus
    shape=0
    if not self.set_VoffsetTES(bias, amplitude, freq, shape):return False

    # to begin, assign zero offset
    offsets = np.zeros(self.NPIXELS)
    client.sendSetFeedbackTable(self.QS_asic_index, offsets)
    self.wait_a_bit(1.0)

    ## switch on the feedback loop
    self.configure_PID(P=0,I=I,D=0) #MP
    self.wait_a_bit(5.0)
    self.assign_pausetime(0.5)

    # correction direction changes with ASIC
    # Tue 13 Feb 2018 15:16:48 CET:  switchover ASIC 1 & 2
    if self.QS_asic_index==0:
        correction_direction = -1
    else:
        correction_direction =  1
    

    # k=1.0 # the first step is big #MP
    for counter in range(count):

        self.debugmsg('count %i/%i: integrating...' % (counter+1,count))
        timeline = self.integrate_scientific_data(save=False)
        self.debugmsg('count %i/%i: finished integrating' % (counter+1,count))
        this_data_avg=timeline.mean(axis=-1)
        prev_offsets=offsets
        offsets = correction_direction*k*(this_data_avg-consigne)+prev_offsets
        self.debugmsg('count %i/%i: applying feedback offsets...' % (counter+1,count))
        client.sendSetFeedbackTable(self.QS_asic_index, offsets)
        self.debugmsg('count %i/%i: feedback offsets applied.' % (counter+1,count))
        self.wait_a_bit()

        self.debugmsg('count %i/%i: data for TES 37: %.5e' % (counter+1,count,this_data_avg[36]))
        k=k/5.0 # 0.2 # and subsequent steps are smaller #MP
    
    return True

def get_amplitude(self):
    """
    Parameters
    ----------
    integration_time: float
    integration time, in seconds
    asic : int, optional
    ASIC number.
        
    """
    timeline = self.integrate_scientific_data(save=False)
    min_timeline = np.min(timeline, axis=-1)
    max_timeline = np.max(timeline, axis=-1)
    return max_timeline - min_timeline

def get_mean(self):
    """
    Parameters
    ----------
    integration_time: float
    integration time, in seconds
    asic : int, optional
    ASIC number.

    """
    timeline = self.integrate_scientific_data(save=False)
    return timeline.mean(axis=-1)

def get_nsamples(self):
    '''
    get the number of samples from QubicStudio
    '''

    # Thu 15 Feb 2018 11:01:09 CET
    # pyStudio hangs sometimes.  I don't know why.  It's not the network switch.
    # it seems to happen when I request nsamples, but not all the time.

    # is this because of "waitingForAckMode" ?
    # what is the difference between client.fetch and client.request ?

    if self.AVOID_HANGUP:
        # HACK: do not request nsamples again if we already have it.
        if not self.nsamples is None: return self.nsamples

    # Mon 19 Feb 2018 10:59:53 CET
    # I think I got it.  There is a hanging request which must be flushed.
    # Tue 27 Feb 2018 16:13:45 CET
    # the flushing didn't work
    client = self.connect_QubicStudio()
    if client is None:return None

    # flush the request queue just in case
    #q=client.abort_requests()
    
    self.debugmsg('getting nsamples...')
    nsamples_all = client.fetch('QUBIC_Nsamples')
    nsamples=nsamples_all[self.asic_index()]
    self.debugmsg('got nsamples.')
    self.debugmsg('nsamples=%i' % nsamples)
    self.nsamples=nsamples
    return nsamples

def get_chunksize(self):
    '''
    get the timeline chunk size
    '''
    if self.AVOID_HANGUP:
        # HACK: don't ask for it if we've already got it (see get_nsamples())
        if not self.chunk_size is None: return self.chunk_size
    
    client = self.connect_QubicStudio()
    if client is None:return None

    # flush the request queue just in case
    #q=client.abort_requests()
    
    self.debugmsg('getting chunk size...')
    chunk_size = client.fetch('QUBIC_PixelScientificDataTimeLineSize')
    # QubicStudio returns an array of integer of length 1.
    # convert this to a simple integer
    chunk_size = int(chunk_size)
    self.debugmsg('chunk size=%i' % chunk_size)
    self.chunk_size=chunk_size
    return chunk_size

def get_RawMask(self):
    '''
    get the mask which identifies which samples are filtered out
    This is 125 values, each an 8-bit bitmask: 1 -> masked
    for example: 255,0,0,0,....  means that the first 8 samples are masked
    '''
    if self.AVOID_HANGUP:
        # HACK: don't ask for it if we've already got it (see get_nsamples())
        if not self.rawmask is None: return self.rawmask
    
    client = self.connect_QubicStudio()
    if client is None:return None

    # flush the request queue just in case
    #q=client.abort_requests()

    rawmasks=client.fetch('QUBIC_RawMasks')
    rawmask=rawmasks[self.asic_index(),:]
    self.rawmask=rawmask
    return rawmask
    
def get_bias(self):
    '''
    get the bias configuration value from QubicStudio
    '''
    client = self.connect_QubicStudio()
    if client is None:return None

    # flush the request queue just in case
    #q=client.abort_requests()

    if not self.AVOID_HANGUP or self.min_bias is None:
        # HACK! if we already have these values, don't ask again
        # still searching for the cause of the hangups
        self.debugmsg('getting bias offset')
        DACoffset_all=client.fetch('QUBIC_TESDACOffset')
        DACoffset=DACoffset_all[self.QS_asic_index]
        bias_offset=DACoffset*self.DAC2V

        self.debugmsg('getting bias amplitude')
        DACamplitude_all=client.fetch('QUBIC_TESDACAmplitude')
        DACamplitude=DACamplitude_all[self.QS_asic_index]
        bias_amplitude=DACamplitude*self.DAC2V

        self.min_bias=bias_offset-bias_amplitude
        self.max_bias=bias_offset+bias_amplitude
    else:
        bias_amplitude=0.5*(self.max_bias-self.min_bias)
        bias_offset=self.min_bias+bias_amplitude

    if not self.AVOID_HANGUP or self.bias_mode is None:
        self.debugmsg('getting bias mode')
        mode_all=client.fetch('QUBIC_TESDACShape')
        self.bias_mode=mode_all[self.QS_asic_index]

    if not self.AVOID_HANGUP or self.bias_frequency is None:
        self.debugmsg('getting bias frequency')
        modulation_all=client.fetch('QUBIC_TESDACFreq')
        self.bias_frequency=modulation_all[self.QS_asic_index]

    self.debugmsg('returning bias info')
    return self.bias_mode,bias_offset,bias_amplitude,self.bias_frequency
    
def integrate_scientific_data(self,save=True):
    '''
    get a data timeline
    '''
    client = self.connect_QubicStudio()
    if client is None:return None

    if not self.exist_timeline_data():self.tdata=[]
    tdata={}
    
    self.debugmsg('calling integrate_scientific_data for ASIC %i' % self.asic)
    self.debugmsg ('integration_time=%.2f' % self.tinteg)

    npixels=self.get_NPIXELS()
    tdata['NPIXSAMP']=npixels
    
    nsamples = self.get_nsamples()
    if nsamples is None: return None
    tdata['NSAMPLES']=nsamples
    
    period = self.sample_period()
    self.debugmsg('period=%.3f msec' % (1000*period))
    
    timeline_size = int(np.ceil(self.tinteg / period))
    self.debugmsg('timeline size=%i' % timeline_size)

    chunk_size=self.get_chunksize()
    if chunk_size is None: return None
    tdata['CHUNK']=chunk_size

    bias_config=self.get_bias()
    tdata['BIAS_MIN']=self.min_bias
    tdata['BIAS_MAX']=self.max_bias
    tdata['BIAS_MOD']=self.bias_frequency
    tdata['BIASMODE']=self.bias_mode
    
    timeline = np.empty((self.NPIXELS, timeline_size))

    # date of the observation
    self.assign_obsdate()
    tdata['DATE-OBS']=self.obsdate
    
    # bath temperature
    self.oxford_read_bath_temperature()
    tdata['TES_TEMP']=self.temperature

    # feedback loop resistance (relay)
    tdata['R_FEEDBK']=self.Rfeedback

    FLL_state,FLL_P,FLL_I,FLL_D=self.get_PID()
    tdata['FLL_STAT']=FLL_state
    tdata['FLL_P']=FLL_P
    tdata['FLL_I']=FLL_I
    tdata['FLL_D']=FLL_D
        
    # integration time
    tdata['INT-TIME']=self.tinteg

    # get the rawmask from which we calculate n_masked
    rawmask=self.get_RawMask()
    
    self.debugmsg('requesting scientific data timeline...')
    parameter = 'QUBIC_PixelScientificDataTimeLine_%i' % self.QS_asic_index
    req = client.request(parameter)
    self.debugmsg('scientific data requested.')
    istart = 0
    for i in range(int(np.ceil(timeline_size / chunk_size))):
        delta = min(chunk_size, timeline_size - istart)
        self.debugmsg('getting next data chunk...',level=2)
        timeline[:, istart:istart+delta] = req.next()[:, :delta]
        self.debugmsg('got data chunk.',level=2)
        istart += chunk_size
    req.abort()
    tdata['TIMELINE']=timeline

    if not self.AVOID_HANGUP:
        for count in range(10):
            req.abort() # maybe we need this more than once?
        del(req) # maybe we need to obliterate this?
        
    if save:self.tdata.append(tdata)
    return timeline


def set_VoffsetTES(self, bias, amplitude, frequency=99, shape=0):
    '''
    command the bias voltage for the TES array
    integration time, asic, should be selected previously with the appropriate assign_() method
    ''' 
    client = self.connect_QubicStudio()
    if client is None:return None

    # check that we don't surpass the max bias permitted by the ADC
    max_offset=self.DAC2V * 2**15
    if bias+amplitude>max_offset:
        print('ERROR!  This combination of offset and amplitude will surpass the maximum bias available: max=%.3fV, offset+amplitude=%.3fV' % (max_offset,bias+amplitude))
        return None

    DACoffset=self.bias_offset2DAC(bias)
    DACamplitude=self.amplitude2DAC(amplitude)

    self.min_bias=bias-amplitude
    self.max_bias=bias+amplitude
    
    # arguments (see comments at top of file):
    #                                  asic, shape, frequency, amplitude,   offset
    self.bias_frequency=frequency
    client.sendSetTESDAC(self.QS_asic_index, shape, frequency, DACamplitude, DACoffset)
    # wait and send the command again to make sure
    self.wait_a_bit()
    client.sendSetTESDAC(self.QS_asic_index, shape, frequency, DACamplitude, DACoffset)
    return True

def get_iv_data(self,replay=False,TES=None,monitor=False):
    '''
    get IV data and make a running plot
    optionally, replay saved data.

    you can monitor the progress of a given TES by the keyword TES=<number>

    setting monitor=True will monitor *all* the TES, but this slows everything down
    enormously!  Not recommended!!

    '''
    client = self.connect_QubicStudio()
    if client is None:return None

    monitor_iv=False
    if isinstance(TES,int):
        monitor_TES_index=self.TES_index(TES)
        monitor_iv=True

    if replay:
        if not isinstance(self.adu,np.ndarray):
            print('Please read an I-V data file, or run a new measurement!')
            return None
        if not isinstance(self.vbias,np.ndarray):
            print('There appears to be I-V data, but no Vbias info.')
            print('Please run make_Vbias() with the correct max and min values')
            return None
        adu=self.adu
    else:
        client = self.connect_QubicStudio()
        if client is None: return None
        self.assign_obsdate(dt.datetime.utcnow())
        if not isinstance(self.vbias,np.ndarray):
            vbias=make_Vbias()
        nbias=len(self.vbias)
        adu = np.empty((self.NPIXELS,nbias))
        self.oxford_read_bath_temperature()
        
    vbias=self.vbias
    nbias=len(self.vbias)

    # figavg=self.setup_plot_Vavg()
    if monitor_iv:figiv,axiv=self.setup_plot_iv(TES)
    if monitor:
        nrows=16
        ncols=8
        figmulti,axmulti=self.setup_plot_iv_multi()

    for j in range(nbias) :
        self.debugmsg("Vbias=%gV " % vbias[j])
        if not replay:
            self.set_VoffsetTES(vbias[j],0.0)
            self.wait_a_bit()
            Vavg= self.get_mean()
            adu[:,j]=Vavg
            self.oxford_read_bath_temperature()
        else:
            Vavg=adu[:,j]

        # print ("a sample of V averages :  %g %g %g " %(Vavg[0], Vavg[43], Vavg[73]) )
        # plt.figure(figavg.number)
        # self.plot_Vavg(Vavg,vbias[j])
        if monitor_iv:
            plt.figure(figiv.number)
            I_tes=adu[monitor_TES_index,0:j+1]
            Iadjusted=self.ADU2I(I_tes)
            self.draw_iv(Iadjusted,axis=axiv)

        if monitor:
            # monitor all the I-V curves:  Warning!  Extremely slow!!!
            TES_index=0
            for row in range(nrows):
                for col in range(ncols):
                    axmulti[row,col].get_xaxis().set_visible(False)
                    axmulti[row,col].get_yaxis().set_visible(False)

                    Iadjusted=self.ADU2I(self.adu[TES_index,0:j+1])
                    self.draw_iv(Iadjusted,colour='blue',axis=axmulti[row,col])
                    text_y=min(Iadjusted)
                    axmulti[row,col].text(max(self.vbias),text_y,str('%i' % (TES_index+1)),va='bottom',ha='right',color='black')
            
                    TES_index+=1
        


    # plt.show()
    self.endobs=dt.datetime.utcnow()
    self.assign_ADU(adu)
    if not replay:
        self.write_fits()
    
    return adu


def get_iv_timeline(self,vmin=None,vmax=None,frequency=None,shape=0):
    '''
    get timeline data with Bias set to sinusoid shape and then extract the I-V data

    integration time should be set previously with assign_integration_time()
    if vmin,vmax are not given, try to get them from self.vbias
    '''
    client = self.connect_QubicStudio()
    if client is None:return None

    if vmin is None:
        if not isinstance(self.vbias,np.ndarray):
            vbias=self.make_Vbias()
        vmin=min(self.vbias)
    if vmax is None:
        if not isinstance(self.vbias,np.ndarray):
            vbias=self.make_Vbias()
        vmax=max(self.vbias)

    self.min_bias=vmin
    self.max_bias=vmax
    
    amplitude=0.5*(vmax-vmin)
    offset=vmin+amplitude
    
    if frequency is None:frequency=99
    self.bias_frequency=frequency
    #amplitude=2*amplitude # BUG CHECK: is this peak-to-peak or amplitude?
    self.debugmsg('amplitude=%.2f, offset=%.2f, frequency=%.2f' % (amplitude,offset,frequency))
    if self.set_VoffsetTES(offset, amplitude, frequency=frequency, shape=shape) is None:return None

    timeline=self.integrate_scientific_data()
    if not isinstance(timeline,np.ndarray):
        print('ERROR! could not acquire timeline data')
        return None
        
    npts_timeline=timeline.shape[1]
    self.debugmsg('number of points in timeline: %i' % npts_timeline)

    ntimelines=self.ntimelines()
    timeline_index=ntimelines-1
    return timeline

def get_ASD(self,TES=1,tinteg=None,ntimelines=10,nbins=1):
    '''
    get timeline data and plot the Amplitude Spectral Density
    timeline data is saved in FITS file, unless in monitor mode.
    in monitor mode, the plots will refresh indefinitely.  exit with Ctrl-C
    '''
    client = self.connect_QubicStudio()
    if client is None:return None

    TES_index=self.TES_index(TES)
    monitor_mode=False
    if not isinstance(ntimelines,int) or ntimelines<=0:
        monitor_mode=True
    save=not monitor_mode
        
    self.assign_integration_time(tinteg)
    self.assign_obsdate()

    # for noise measurements, we set the feedback resistance to 100kOhm
    self.set_Rfeedback(100)

    idx=0
    ax_timeline=None
    ax_asd=None
    while monitor_mode or idx<ntimelines:
        self.debugmsg('ASD monitoring loop count: %i' % idx)

        # read the bath temperature at each loop
        Tbath=self.oxford_read_bath_temperature()
        
        timeline = self.integrate_scientific_data(save=True) # have to save in memory for plotting afterwards
        self.debugmsg('ASD monitoring: ntimelines=%i' % self.ntimelines())
        timeline_index=self.ntimelines()-1
        result=self.plot_ASD(TES,timeline_index,ax_timeline=ax_timeline,ax_asd=ax_asd,save=save,nbins=nbins)
        ax_asd=result['ax_asd']
        ax_timeline=result['ax_timeline']

        # if only monitoring, get rid of the one just plotted
        if monitor_mode:del(self.tdata[-1])
        idx+=1

    if not monitor_mode: self.write_fits()

    return 

