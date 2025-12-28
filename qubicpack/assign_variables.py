"""
$Id: assign_variables.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Thu 13 Jul 2017 14:11:07 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

default values for various parameters in qubicpack

"""
import numpy as np
import os,subprocess,pickle
import datetime as dt
import matplotlib

from qubicpack.utilities import asic_reversal_date, NPIXELS, TES_index, ASIC_index
from qubicpack.pix2tes import assign_pix_grid, assign_pix2tes, tes2pix, pix2tes, TES2PIX
from qubicpack import __file__
from satorchipy.datefunctions import utcnow

### the rest of the defs are methods of the qubicasic object

def assign_defaults(self):
    self.assign_constants()
    self.logfile=None
    self.NPIXELS = NPIXELS
    self.assign_obsdate()
    #self.assign_datadir() # already called from assign_obsdate() above
    self.endobs=None
    self.NPIXELS_sampled=None
    self.detector_name='undefined'
    self.FLL_state=None
    self.FLL_P=None
    self.FLL_I=None
    self.FLL_D=None
    self.Rfeedback=10e3 # 10kOhm, this is selectable between 10kOhm and 100kOhm (also called "relay" resistance)
    self.asic=None
    self.assign_integration_time()
    self.adu=None
    self.vbias=None
    self.timeline_vbias=None
    self.timeline_Ites = None
    self.timeline_Vtes = None
    self.timeline_Ptes = None
    self.cycle_vbias=True
    self.nbiascycles=None
    self.bias_frequency=None
    self.bias_period=None
    self.max_permitted_bias=3.0
    self.max_bias=None
    self.min_bias=None
    self.max_bias_position=None
    self.bias_factor=1.0
    self.bias_mode=None
    self.pausetime=0.3
    self.observer='APC LaboMM'
    self.nsamples=None
    self.chunk_size=None
    self.rawmask=None
    self.timeline_conversion=None
    self.tdata = [{}]
    self.tdata[0]['WARNING'] = []
    self.pix_grid = assign_pix_grid()
    self.filtersummary=[]
    for idx in range(self.NPIXELS): self.filtersummary.append(None)
    self.assign_lookup_table()
    self.temperature=None
    self.calsource_LF=None
    self.calsource_HF=None
    self.modulator=None
    self.datafiletype=None
    self.dataset_name=None
    self.hornswitch_files = None
    self.hk = {}
    self.assign_fitsblurbs()
    self.temperature_labels = None
    return

def assign_fitsblurbs(self):
    '''
    keynames and descriptions for FITS files
    '''
    
    self.fitsblurbs={}
    self.fitsblurbs['TELESCOP']='Telescope used for the observation'
    self.fitsblurbs['OBSERVER']='name of the observer'
    self.fitsblurbs['DATE-OBS']='date of the observation in UTC'
    self.fitsblurbs['END-OBS'] ='end time of the observation in UTC'
    self.fitsblurbs['NSAMPLES']='number of samples per integration time'
    self.fitsblurbs['INT-TIME']='integration time in seconds'
    self.fitsblurbs['NPIXELS'] ='number of TES detectors in the array'
    self.fitsblurbs['NPIXSAMP']='number of pixels sampled'
    self.fitsblurbs['ASIC']    ='ASIC id (one quarter of the full QUBIC array)'
    self.fitsblurbs['QUBIC-IP']='address of the QUBIC Local Control Computer'
    self.fitsblurbs['NCYCLES'] ='number of cycles of the Bias voltage'
    self.fitsblurbs['CYCBIAS'] ='ramp return Bias,  yes or no'
    self.fitsblurbs['TES_TEMP']='TES physical temperature in K'
    self.fitsblurbs['BIAS_MOD']='bias modulation frequency'
    self.fitsblurbs['BIAS_MIN']='minimum bias in V'
    self.fitsblurbs['BIAS_MAX']='maximum bias in V'
    self.fitsblurbs['BIAS_FAC']='multiplicative factor for bias'
    self.fitsblurbs['BIASMODE']='sinusoidal, triangle, continuous'
    self.fitsblurbs['FLL_STAT']='Flux Lock Loop state (on/off)'
    self.fitsblurbs['FLL_P']   ='Flux Lock Loop P level'
    self.fitsblurbs['FLL_I']   ='Flux Lock Loop I level'
    self.fitsblurbs['FLL_D']   ='Flux Lock Loop D level'
    self.fitsblurbs['DET_NAME']='ID of the detector array'
    self.fitsblurbs['R_FEEDBK']='Feedback resistance in Flux Lock Loop'
    self.fitsblurbs['R_HEATER']='Heater in feedback loop:  On/Off'
    self.fitsblurbs['CHUNK']='data chunk size delivered by QubicStudio'
    return

def assign_constants(self):
    '''
    assign some constant values used in calculations and measurements
    '''
    self.asic_reversal_date = asic_reversal_date
    self.zero=1e-9
    #self.QubicStudio_ip='134.158.186.233'
    #self.QubicStudio_ip='134.158.187.21'
    self.QubicStudio_ip='192.168.2.8'
    #self.DAC2V=2.627e-4    # email from Michel Piat 2018/02/09 17:14 CET
    self.DAC2V=9.404/2**15 # measured Tue 13 Feb 2018 15:25:11 CET
    self.kBoltzmann=1.3806485279e-23
    self.Rshunt=10.e-3  # 10mOhm, mail from M.Piat to M.Salatino 2017-08-10
    self.Rbias =10.e3   # 10kOhm, mail from M.Piat to M.Salatino 2017-08-10
    self.colours=['blue','green','red','cyan','magenta','yellow','black']
    self.Vinfinity=9.0 # used to calculate the I-V offset (force the line through I,Vinfinity to get R=1 at infinity)
    return

def assign_observer(self,observer='APC LaboMM'):
    if not isinstance(observer,str):
        observer=str(observer)
    self.observer=observer
    return

def assign_asic(self,asic=1):
    if asic is None:asic=self.asic
    if not isinstance(asic,int) or asic<1 or asic>2:
        self.printmsg('asic should have an integer value: 1 or 2.  assigning default asic=1')
        self.asic=1
    else:
        self.asic=asic

    if self.asic==1 or self.asic==2:
        self.detector_name='P87'


    ######### Is this correct?
    #### QubicStudio has a reverse nomenclature for the ASIC index
    #### compared to the translation tables (eg. Correspondance.xlsx)
    # so define here a specific QubicStudio ASIC index which should be used in the acquisition methods
    # see, for example, integrate_scientific_data() in acquisition.py
    asic_index=self.asic_index()
    self.QS_asic_index=asic_index
    if self.obsdate < self.asic_reversal_date:
        if asic_index==0:
            self.QS_asic_index=1
        else:
            self.QS_asic_index=0

    # Wed 02 Aug 2017 15:48:15 CEST
    # during lab tests, QubicStudio is always using asic_index=0
    # we change the asic by physically switching a cable to another connector
    #self.QS_asic_index=0
    # Fri 04 Aug 2017 13:38:10 CEST
    # in fact, we should change a jumper on the FPGA board to change the ID of the ASIC
    
    return

def asic_index(self):
    return ASIC_index(self.asic)

def assign_integration_time(self,tinteg=0.1):
    if tinteg is None:tinteg=self.tinteg
    if tinteg is None or tinteg<0.0:
        self.printmsg('integration time should be a positive number of seconds.  Assigning default Tinteg=0.1')
        self.tinteg=0.1
    else:
        self.tinteg=tinteg
    return

def assign_ADU(self,adu):
    if (not isinstance(adu,np.ndarray)):
        self.printmsg('Please enter a 2 dimensional numpy array with the first dimension=%i' % self.NPIXELS)
        return None
    self.adu=adu
    return

def assign_pausetime(self,pausetime):
    if (not isinstance(pausetime,int)) and (not isinstance(pausetime,float)):
        self.printmsg('pause time should be a number of seconds.  Assigning default pausetime=%.3f seconds' % self.pausetime)
    else:
        self.pausetime=pausetime
    return

    
def assign_ip(self,ip):
    if (not isinstance(ip,str)):
        self.printmsg('please give an IP address for QubicStudio in the form xxx.xxx.xxx.xxx:')
        self.printmsg('assigning default IP address: %s' % self.QubicStudio_ip)
        return None
    self.QubicStudio_ip=ip
    return


def assign_temperature(self,temp):
    '''
    assign the bath temperature to the asic object
    '''
    try:
        temperature = float(temp)
    except:
        self.printmsg('ERROR! Temperature should be a number in Kelvin (not milliKelvin)')
        temperature = None
    
    self.temperature = temperature
    if self.tdata is None: self.tdata = [{}]
    self.tdata[-1]['TES_TEMP'] = temperature
    return self.temperature

def assign_temperature_labels(self):
    '''
    read temperature labels
    '''
    pkg_dir = os.path.dirname(__file__)
    label_file = os.sep.join([pkg_dir,'data','TD_TEMPERATURE_LABELS.txt'])
    if not os.path.isfile(label_file):
        self.temperature_labels = None
        self.printmsg('could not find temperature labels: %s' % label_file,verbosity=3)
        return
    h = open(label_file)
    lines = h.read().split('\n')
    h.close()
    self.temperature_labels = {}
    for line in lines:
        if line=='': continue
        keyval = line.split('=')
        if len(keyval)<2: continue
        key = keyval[0].strip()
        val = keyval[1].strip()
        #self.temperature_labels[key] = val

        # QubicStudio assigned different keynames to these
        if key.find('AVS')==0:
            qskey = key.upper()
            self.temperature_labels[qskey] = val
        
        if key.find('HEATER')==0:
            heater_num = int(key.replace('HEATER',''))
            qskey = 'Heaters_Amp_%i' % (heater_num - 1)
            self.temperature_labels[qskey] = val
            qskey = 'Heaters_Volt_%i' % (heater_num - 1)
            self.temperature_labels[qskey] = val
        
        if key.find('TEMPERATURE')==0:
            temp_num = int(key.replace('TEMPERATURE',''))
            qskey = 'Temp_%i' % (temp_num - 1)
            self.temperature_labels[qskey] = val
            
    # assign the temperature labels to the asic objects
    for asicobj in self.asic_list:
        if asicobj is None: continue
        asicobj.temperature_labels = self.temperature_labels
        
    return

def assign_obsdate(self,d=None):
    '''
    assign the observation date
    '''
    global TES2PIX
    if not isinstance(d,dt.datetime):
        self.obsdate=utcnow()
    else:
        self.obsdate=d

    self.assign_datadir()
    TES2PIX = assign_pix2tes(self.obsdate)
    return self.obsdate
    

def assign_datadir(self,d=None):
    '''
    find a place to write data
    '''

    # valid possibilities, in order of preference
    cwd=os.getcwd()
    datadirs=['/data/qubic',
              '/home/qubic/data',
              '/home/qubic/Bureau/PyStudio Work/data',
              '/home/work/qubic/data',
              cwd]
    if 'HOME' in os.environ.keys():
        home=os.environ['HOME']
        datadirs.append(home+'/data')
    datadirs.append(cwd)
    fallback_dir='/tmp/qubic'

    # if a directory is given, make this the priority possibility
    if isinstance(d,str): datadirs=[d]+datadirs

    # make sure we have write permission
    tmpfile=utcnow().strftime('tmpFile_%Y%m%dT%M%H%S.%f.UTC.tmp')
    tmpdir_ok=False
    for datadir in datadirs:
        try:
            tmpfile_full=datadir+os.sep+tmpfile
            tmpfile_full.replace('/',os.sep)
            h=open(tmpfile_full,'w')
            h.write('check if we have write permission\n')
            h.close()
            os.remove(tmpfile_full)
            tmpdir_ok=True
            self.datadir=datadir
            break
        except:
            tmpdir_ok=False
            

    if not tmpdir_ok:
        # try the fall back directory
        datadir=fallback_dir.replace('/',os.sep)
        try:
            if not os.path.exists(fallback_dir): os.system('mkdir -p %s' % fallback_dir)
            tmpfile_full=datadir+os.sep+tmpfile
            h=open(tmpfile_full,'w')
            h.write('check if we have write permission\n')
            h.close()
            os.remove(tmpfile_full)
            tmpdir_ok=True
            self.datadir=datadir
        except:
            self.printmsg('ERROR! Could not find a suitable data directory!')
            return None

    msg='Data will be written to directory: %s' % self.datadir
    self.debugmsg(msg)

    # check how much space is available
    try:
        cmd='df %s' % self.datadir
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out,err=proc.communicate()
        gigs_available=eval(str(out.decode('UTF-8')).split('\n')[1].split()[3])/float(1024**2)
        if gigs_available<1:
            self.printmsg('WARNING! running out of disk space.  Only %.1f GiB space left on disk' % gigs_available)
    except:
        self.debugmsg('WARNING! Could not determine disk space available.')
        
    
    return self.datadir


def assign_bias_factor(self,factor):
    '''
    assign the multiplicative factor for the bias voltage
    '''
    if not (isinstance(factor,float) or isinstance(factor,int)):
        self.printmsg('ERROR! Bias factor should be a number.  Assigning default: 1.0')
        self.bias_factor=1.0
        return
    self.bias_factor=factor
    return

def assign_detector_name(self,det_name):
    '''
    assign the name of the detector array
    examples: P73 (tested July-October 2017)
              P82 (tested October-November 2017)
    '''
    if not isinstance(det_name,str):
        self.printmsg('ERROR! Please enter a valid name of the detector array (eg. P73, P82)')
        det_name = 'undefined'

        
    self.detector_name = det_name
    if self.__object_type__=='qubicfp':
        for asicobj in self.asic_list:
            if asicobj is not None:
                asicobj.assign_detector_name(det_name)
    return

              
def guess_detector_name(self):
    '''
    guess which detector matrix it is for this data
    the DET_NAME keyword was added Tue 21 Nov 2017
    We tested P73 from July to 4 Nov 2017 and then P82
    '''
    if not self.detector_name=='undefined':
        return self.detector_name

    if self.obsdate is None:
        self.printmsg('no observation date!')
        return self.detector_name

    P73_lastdate = dt.datetime.strptime('2017-11-05','%Y-%m-%d').replace(tzinfo=dt.timezone.utc)
    if self.obsdate<P73_lastdate:
        self.assign_detector_name('P73')
        return self.detector_name

    P82_lastdate = dt.datetime.strptime('2017-11-30','%Y-%m-%d').replace(tzinfo=dt.timezone.utc)
    if self.obsdate<P82_lastdate:
        self.assign_detector_name('P82')
        return self.detector_name

    QS_firstdate=dt.datetime.strptime('2018-11-19','%Y-%m-%d').replace(tzinfo=dt.timezone.utc)
    if self.datafiletype!='QP_FITS':
        if self.obsdate>QS_firstdate:            
            if self.asic==1 or self.asic==2:
                self.assign_detector_name('P87')
                return self.detector_name

    self.printmsg('Detector array is: %s' % self.detector_name,verbosity=2)
    return self.detector_name

def assign_logfile(self,rootname=None):
    '''
    assign a filename for the log file
    '''
    if rootname is None:rootname='logfile'
    
    logfile='QUBIC_%s_%s.txt' % (rootname,utcnow().strftime('%Y%m%dT%H%M%SUTC'))
    logfile_fullpath=self.output_filename(logfile)
    self.logfile=logfile_fullpath
    return self.logfile


### lookup table stuff was moved from pix2tes Fri 24 May 2019 12:22:52 CEST
def assign_lookup_table(self):
    '''
    make the lookup table with results for comparison

    this was only done for array P73, so don't ask for this table otherwise
    '''
    filename='TES_translation_table.pickle'
    cwd=os.getcwd()
    dirs=[cwd]
    if not isinstance(self.datadir,str):
        self.assign_datadir()
    if isinstance(self.datadir,str):
        dirs.append(self.datadir)
        
    gotit=False
    for d in dirs:
        filename_fullpath='%s/%s' % (d,filename)
        if os.path.exists(filename_fullpath):
            gotit=True
            break
    if not gotit:
        if self.detector_name=='P73':
            self.printmsg('WARNING! Cannot find translation table file: %s' % filename_fullpath,verbosity=2)
            self.printmsg('Open loop and Room Temperature tests will not be noted in plots etc.',verbosity=2)
        self.transdic=None
        return None

    h=open(filename_fullpath,'rb')
    self.transdic=pickle.load(h)
    h.close()
    return self.transdic

def lookup_TEStable(self,key='PIX',value=100):
    if self.transdic is None:
        self.printmsg('No translation table.  Please load.',verbosity=2)
        return None
    
    if not key in self.transdic[0].keys():
        self.printmsg('Please enter a valid key.  Choose from:',verbosity=0)
        for k in transdic[0].keys():
            print('    %s' % k)
        return None
    
    ret=[entry for entry in self.transdic if entry[key]==value]
    if len(ret)==1:
        ret=ret[0]
    return ret
