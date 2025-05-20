'''
$Id: fpmethods.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Tue 28 May 2019 12:43:24 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

methods for the QUBIC focal plane class qubicfp
'''
import numpy as np
import datetime as dt
import sys,os,time,re
from copy import copy

from qubicpack.qubicasic import qubicasic
from qubicpack.plot_fp import plot_fp
from qubicpack.utilities import NPIXELS,TES_index
from qubicpack.pixel_translation import make_id_focalplane

def assign_defaults(self):
    '''default values for object variables
    '''
    self.assign_constants()

    self.obsdate = None
    self.endobs = None
    self.logfile = None
    self.observer = 'APC LaboMM'
    self.tdata = None
    self.temperature = None
    self.hk  =  {}
    self.hornswitch_files = None
    self.detector_name = 'undefined'
    self.dataset_name = None
    self.datadir = None
    self.assign_fitsblurbs()
    self.mapinfo_list = None
    self.assign_temperature_labels()
    return

def assign_verbosity(self,verbosity):
    '''
    assign the verbosity level for messages on the screen
    '''
    self.verbosity = verbosity
    for asic_obj in self.asic_list:
        if asic_obj is not None:
            asic_obj.verbosity = verbosity
    return

def assign_temperature(self,temp):
    '''
    assign the bath temperature to the asic objects
    '''
    for asic_obj in self.asic_list:
        if asic_obj is not None:
            asic_obj.assign_temperature(temp)
    return

def infotext(self):
    '''
    some basic info to put on plots
    '''

    txt = self.dataset_name
    if self.temperature is not None:
        txt += ' T$_\\mathrm{bath}$=%.1fmK' % (1000*self.temperature)

    return txt

def calsource_oldinfo(self):
    '''
    return calsource info for data before the implementation of MsgStr
    (see calsource_info() below)
    '''
    if 'CALSOURCE-CONF' not in self.hk.keys():
        return None

    conf = self.hk['CALSOURCE-CONF']

    info = {}

    if 'timestamp' in conf.keys():
        info_tstamp = conf['timestamp'][0]
        info_date = dt.datetime.utcfromtimestamp(info_tstamp)
        info['date'] = info_date
    else:
        info['date'] = None

    keytranslation = {}

            
    keytranslation['calsource'] = {'status'    :'CalSource',
                                   'frequency' :'Cal_freq',
                                   'synth_freq':'Syn_freq'}

    keytranslation['modulator'] = {'status'    :'Modulator',
                                   'frequency' :'Mod_freq',
                                   'amplitude' :'Mod_ampl',
                                   'duty_cycle':'Mod_duty',
                                   'offset'    :'Mod_offs',
                                   'shape'     :'Mod_shap'}

    keytranslation['amplifier'] = {'status'               :'Amplifier',
                                   'mode'                 :'Amp_mode',
                                   'filter low frequency' :'Amp_lfreq',
                                   'filter high frequency':'Amp_hfreq',
                                   'coupling'             :'Amp_coup',
                                   'dynamic range'        :'Amp_rang',
                                   'gain'                 :'Amp_gain'}
    

    idx_translation = {}
    idx_translation['shape'] = ['square','sine','DC']
    idx_translation['mode'] = ['bypass',
                                '6db_low_pass',
                                '12db_low_pass',
                                '6db_high_pass',
                                '12db_high_pass',
                                'bandpass']
    idx_translation['status'] = ['OFF','ON']
    idx_translation['coupling'] = ['GROUND','DC','AC']

    for dev in keytranslation.keys():
        
        info[dev] = {}

        for parm in keytranslation[dev].keys():
            confkey = keytranslation[dev][parm]
            if confkey in conf.keys():
                confval = conf[confkey][0]
                try:
                    confidx = int(confval)
                except:
                    confidx = None

                if confidx is not None and parm in idx_translation.keys():
                    info[dev][parm] = idx_translation[parm][confidx]
                else:
                    info[dev][parm] = confval

    # for some reason, on 2019-11-12 and 2019-11-14, the modulation amplitude and frequency are given in mV and mHz
    if info_tstamp>=1573572852 and info_tstamp<=1573750907:
        if 'frequency' in info['modulator'].keys():
            modfreq = info['modulator']['frequency']
            if modfreq>100:
                info['modulator']['frequency'] = modfreq/1000
        if 'amplitude' in info['modulator'].keys():
            modamp = info['modulator']['amplitude']
            if modamp>100:
                info['modulator']['amplitude'] = modamp/1000
                


    return info

def calsource_info(self):
    '''
    return a dictionary of calibration source configuration information
    '''
    if 'CALSOURCE-CONF' not in self.hk.keys():
        return None

    # go back to "oldinfo" because of a bug with MsgStr which gets cutoff at 512 characters
    #return self.calsource_oldinfo()
    
    if 'MsgStr' not in self.hk['CALSOURCE-CONF'].keys():
        return self.calsource_oldinfo()
    
    info_txt = self.hk['CALSOURCE-CONF']['MsgStr'][0]
    if type(info_txt)!=str: return self.calsource_oldinfo()

    # device_list = ['amplifier','modulator','calsource','cf']
    
    info_rawlist = info_txt.split()
    device_list = []
    for item in info_rawlist:
        if item.find(':')<0: continue
        dev = item.split(':')[0]
        if dev not in device_list: device_list.append(dev)
            
    info = {}

    info_tstamp = float(info_rawlist[0])
    info_date = dt.datetime.utcfromtimestamp(info_tstamp)
    info['date'] = info_date

    for dev in device_list:
        info[dev] = {}

    if info_txt.find('busy')>0:
        for dev in device_list:
            info[dev]['status'] = 'busy'
        return info

    munits = ['mHz','mVpp','mVdc']
    units = ['GHz','HZ','Hz','hz','Vpp','Vdc','V','%']
    for item in info_rawlist[2:]:
        if item.find(':')<0:
            info[item] = 'incomplete'
            continue
        
        cols = item.split(':')

        if len(cols)==1:
            info['calsource']['status'] = info_txt
            continue
        
        dev = cols[0]
        if dev=='lamp' or dev=='arduino' or dev=='synthesiser' or dev=='synthesizer':
            continue

        val_list = cols[1].split('=')
        if len(val_list)==1:
            status_str = val_list[0].upper()
            if status_str=='ON' or status_str=='OFF':
                info[dev]['status'] = val_list[0]
            else:
                info[dev]['incomplete'] = val_list[0]
            continue

        parm = val_list[0].lower()
        val = val_list[1]
        
        if val=='--' or val.upper()=='UNKNOWN':
            info[dev][parm.lower()] = -1
            continue

        goto_next_item = False
        for munit in munits:
            if val.find(munit)>0:
                info[dev][parm.lower()] = 0.001*float(val.replace(munit,''))
                goto_next_item = True
                break
        if goto_next_item: continue
        
        if parm=='gain':
            info[dev][parm] = int(val)
            continue

        if parm=='duty_cycle':
            val_stripped = val.strip().lower().replace('+','').replace('-','')
            if len(val_stripped)==0 or re.search('[a-z]',val_stripped):
                info[dev][parm] = 'none'
                continue
            info[dev][parm] = eval(val_stripped)
            continue

        for unit in units:
            if val.find(unit)>0:
                info[dev][parm.lower()] = float(val.replace(unit,''))
                goto_next_item = True
                break
        if goto_next_item: continue

        if val.upper()=='ON' or val.upper()=='OFF':
            info[dev][parm.lower()] = val.upper()
            continue

        info[dev][parm.lower()] = val.lower()


    return info

def calsource_infotext(self):
    '''
    return a calsource info in a string suitable for plot subtitle
    '''
    info = self.calsource_info()
    if info is None:
        return 'Calsource: No information'

    if info['calsource']['status'] == 'OFF':
        return 'Calsource %s' % info['calsource']['status']

    
    calsrc_txt = 'Calsource: '
    if info['calsource']['status'] == 'UNKNOWN':
        calsrc_txt = 'Calsource:UNKNOWN '

    if info['calsource']['frequency'] > 0:
        calsrc_txt += 'frequency=%.2fGHz' % info['calsource']['frequency']
    else:
        calsrc_txt += 'frequency=UNKNOWN'

        

    modulator_units = {'frequency':'%.3fHz', 'shape':'%s','amplitude':'%.3fVpp','offset':'%.3fVdc','duty_cycle':'%.1f%%'}
    if info['modulator']['status'] == 'OFF':
        calsrc_txt += ' modulator OFF'
    elif info['modulator']['shape'] == 'DC':
        calsrc_txt += ' No modulation, offset=%.2fVdc' % info['modulator']['offset']
    else:
        txt_list = []
        for key in modulator_units.keys():
            if key in info['modulator'].keys():
                if info['modulator'][key]=='none':
                    txt = key+'=N/A'
                else:
                    txt = key+'='+modulator_units[key] % info['modulator'][key]
                txt_list.append(txt)
            else:
                txt_list.append(key+'=unknown')
        calsrc_txt += '\nmodulator: '+' '.join(txt_list)

    if info['amplifier']['status'] == 'OFF':
        calsrc_txt += 'amplifier OFF'
    else:
        calsrc_txt += '\namplifier:'
        if info['amplifier']['status'] == 'UNKNOWN':
            calsrc_txt += 'UNKNOWN'
        for parm in info['amplifier'].keys():
            if parm!='status':
                calsrc_txt += ' %s=%s' % (parm,info['amplifier'][parm])

    return calsrc_txt

def read_qubicstudio_science_fits(self,hdu):
    '''
    read the science data for an ASIC
    The HDU passed here as the argument should already have been identified as the Science HDU
    '''
    self.printmsg('DEBUG: read_qubicstudio_science_fits object type is %s' % self.__object_type__,verbosity=3)
        
    asic_no = hdu.header['ASIC_NUM']
    asic_ctr = asic_no - 1
    qubicasic.verbosity = self.verbosity
    if self.asic_list[asic_ctr] is None:
        self.asic_list[asic_ctr] = qubicasic()
    self.asic_list[asic_ctr].read_qubicstudio_science_fits(hdu)

    obsdate = self.asic_list[asic_ctr].obsdate
    if obsdate is None:
        obsdate_str = 'None'
    else:
        obsdate_str = obsdate.strftime('%Y-%m-%d %H:%M:%S.%f')
    self.printmsg('ASIC%i     Observation date: %s' % (asic_no,obsdate_str))
    if self.obsdate is None:
        self.printmsg('DEBUG: setting obsdate which was None: %s' % obsdate_str,verbosity=3)
        self.obsdate = obsdate
    if obsdate is not None and self.obsdate>obsdate:
        self.obsdate = obsdate
        self.printmsg('Observation date is not the same between ASIC data.  Using the earlier date',verbosity=3)

    endobs = self.asic_list[asic_ctr].endobs
    if endobs is None:
        endobs_str = 'None'
    else:
        endobs_str = endobs.strftime('%Y-%m-%d %H:%M:%S.%f')
    self.printmsg('ASIC%i Observation end date: %s' % (asic_no,endobs_str))
    if self.endobs is None:
        self.printmsg('DEBUG: setting endobs which was None: %s' % endobs_str,verbosity=3)
        self.endobs = endobs
    if endobs is not None and self.endobs<endobs:
        self.endobs = endobs
        self.printmsg('End Observation is not the same between ASIC data.  Using the later date',verbosity=3)
        

    return

def read_qubicstudio_raw_fits(self,hdu):
    '''
    read the RAW data for a given ASIC
    The HDU passed here as the argument should already have been identified as the RAW HDU
    '''
    asic_num = hdu.header['ASIC_NUM']
    asic_idx = asic_num - 1
    qubicasic.verbosity = self.verbosity
    if self.asic_list[asic_idx] is None:
        self.asic_list[asic_idx] = qubicasic()
        
    return self.asic_list[asic_idx].read_qubicstudio_raw_fits(hdu)

def read_qubicstudio_asic_fits(self,hdulist):
    '''
    read the data giving the ASIC configuration
    The HDU passed here as the argument should already have been identified as the ASIC HDU

    we should read the science data first, and then read the corresponding ASIC table
    so we read this file for each of the defined asic objects
    '''
    for asic_obj in self.asic_list:
        if asic_obj is not None:
            asic_obj.read_qubicstudio_asic_fits(hdulist)
    return

def read_qubicpack_fits(self,hdulist):
    '''
    read a FITS file that was written by QubicPack
    argument is an hdulist after opening a FITS file
    and confirming that it really is a QubicPack fits file
    '''
    
    self.datafiletype = 'QP_FITS'
    self.observer = hdulist[0].header['OBSERVER']
    self.assign_obsdate(dt.datetime.strptime(hdulist[0].header['DATE-OBS'],'%Y-%m-%d %H:%M:%S UTC'))
            
    asic_idx = hdulist[0].header['ASIC'] - 1
    self.QubicStudio_ip=hdulist[0].header['QUBIC-IP']

    if 'TES_TEMP' in hdulist[0].header.keys():
        self.temperature=hdulist[0].header['TES_TEMP']
    else:
        self.temperature=None

    if 'END-OBS' in hdulist[0].header.keys():
        self.endobs=dt.datetime.strptime(hdulist[0].header['END-OBS'],'%Y-%m-%d %H:%M:%S UTC')
    else:
        self.endobs=None

    if 'DET_NAME' in hdulist[0].header.keys():
        self.assign_detector_name(hdulist[0].header['DET_NAME'])
    # in case detector name is undefined...
    self.guess_detector_name()

    self.asic_list[asic_idx] = qubicasic()
    self.asic_list[asic_idx].read_qubicpack_fits(hdulist)

    return

#### check if there is data
def exist_data(self):
    '''
    check if there's any data
    '''
    if self.hk: return True
    for asicobj in self.asic_list:
        if asicobj is None:continue
        if asicobj.exist_timeline_data(): return True
        if asicobj.exist_iv_data(): return True
        if asicobj.hk: return  True

    self.printmsg('No data!')
    return False

        

#### wrappers to return values
def nasics(self):
    '''
    count the number of asics that have data
    '''
    nasics = 0
    for asicobj in self.asic_list:
        if asicobj is None: continue
        if asicobj.exist_timeline_data(): nasics += 1
    self.printmsg('Number of ASICs with data: %i' % nasics,verbosity=3)
    return nasics
    
def args_ok(self,TES=None,asic=None,allow_multiple_TES=False,asic_only_required=False):
    '''
    check if arguments are okay for the wrapper

    if TES or asic is None, then it is requested
    if TES is greater than 128, then return the ASIC number and the TES number for that ASIC
    return None if not valid
    return tuple (TES,asic)
    '''

    self.printmsg('checking arguments:  TES=%s, asic=%s' % (TES,asic),verbosity=4)

    # first, count the number of asics that have data
    nasics = self.nasics()

    if nasics==0:
        self.printmsg('There is no data!')
        return None

    # at least one of asic or TES has to be specified
    if asic is None and TES is None:
        if asic_only_required:
            self.printmsg('Please give an asic number')
        else:
            self.printmsg('Please give a TES number')
        return None            


    # unless we only want the asic number, we must have something for TES
    if TES is None and not asic_only_required:
        self.printmsg('Please give a TES number')
        return None

    # if TES is given as a list or tuple, change it to a numpy array
    if isinstance(TES,list) or isinstance(TES,tuple):
        TES = np.array(TES)
        
    if allow_multiple_TES:

        if isinstance(TES,str) and TES=='all':
            if asic is not None: # asking for all TES from a particular ASIC
                TES = np.ones((NPIXELS),dtype=bool)
                return (TES,asic)
            # asking for all TES from all ASICs
            TES = np.ones((NPIXELS*nasics),dtype=bool)
            return (TES,None)
        
        errmsg = 'The TES to be used should be defined by an array of Boolean with a size that is a multiple of %i' % NPIXELS
        if isinstance(TES,np.ndarray):
            if not TES.dtype==bool\
               or not len(TES.shape)==1:
                self.printmsg(errmsg)
                return None            
       
            nTES = TES.size
            if nTES%NPIXELS != 0:
                self.printmsg(errmsg)
                return None
       
            nasics_requested = nTES//NPIXELS 
            self.printmsg('The number of ASICs required for %i TES is %i' % (nTES,nasics_requested),verbosity=3)
            if nasics_requested==1:
                if asic is None:
                    self.printmsg('Please give an asic number')
                    return None
                return (TES,asic)
            
            if nTES!=nasics*NPIXELS:
                errmsg = 'The size of the requested set of TES is %i and it does not match the total number of TES=%i' % (nTES,nasics*NPIXELS)
                self.printmsg(errmsg)
                return None
            return (TES,None)
    
    # from here on, TES should not be an array
    if not allow_multiple_TES and isinstance(TES,np.ndarray):
        self.printmsg('Multiple TES are not permitted for this application')
        return None

    if TES is not None:
        TESidx = TES-1
    
    if asic is None:
        asic_idx = TESidx // NPIXELS
        asic = asic_idx + 1
    else:
        asic_idx = asic - 1

    self.printmsg('DEBUG: type(asic) should be int: %s' % (type(asic)),verbosity=5)
    if asic<=0:
        self.printmsg('Please give a valid asic number')
        return None

    if asic>len(self.asic_list) or self.asic_list[asic_idx] is None or not self.asic_list[asic_idx].exist_timeline_data():
        self.printmsg('No data for ASIC %i' % asic)
        return None

    if asic_only_required:
        return (None,asic)
    
    if TES < 1:
        self.printmsg('Please give a valid TES number')
        return None

    TESidx = TESidx % NPIXELS
    TES = TESidx + 1

    return (TES,asic)

def asic(self,asic=None,TES=None):
    '''
    return the requested asic object
    '''
    args =self.args_ok(None,asic,asic_only_required=True)
    if args is None:return
    TES,asic = args
    asic_idx = asic - 1
    self.printmsg('return object for asic=%i' % asic,verbosity=3)
    return self.asic_list[asic_idx]
    
def Rfeedback(self,asic=None):
    '''
    return the feedback resistance for a given asic
    if they are all the same, it doesn't matter which asic
    '''
    Rfb_list = []
    for asicobj in self.asic_list:
        if asicobj is not None:
            Rfb_list.append(asicobj.Rfeedback)
    Rfb_values = np.unique(Rfb_list)
    if len(Rfb_values)==1:
        return Rfb_values[0]
    
    if asic is None:
        self.printmsg('Please enter an asic number')
        return None

    asicobj = self.asic(asic)
    if asicobj is None: return None
    return asicobj.Rfeedback

def relay_heater(self,asic=None,timeline_index=0):
    '''
    return the state of the relay heater:  On or Off
    '''
    if asic is None:
        self.printmsg('Please enter an asic number')
        return None

    asic_idx = asic-1
    asicobj = self.asic_list[asic_idx]
    tdata = asicobj.tdata[timeline_index]
    if 'R_HEATER' in tdata.keys():
        if tdata['R_HEATER']==1:
            return 'ON'
        return 'OFF'
    return 'UNKNOWN'

def FLL_State(self,asic=None):
    '''
    return the Flux Lock Loop state
    '''
    if asic is None:
        self.printmsg('Please enter an asic number')
        return None

    asic_idx = asic-1
    if asic_idx>len(self.asic_list):
        self.printmsg('ASIC%i does not exist' % asic)
        return None
    
    asicobj = self.asic_list[asic_idx]
    if asicobj is None:return None
    
    return asicobj.FLL_State()

#### timeline methods
def bias_phase(self,asic=0,TES=None):
    '''
    return the bias on the detectors
    '''
    args = self.args_ok(TES,asic,asic_only_required=True)
    if args is not None:
        TES,asic = args
        asic_idx = asic-1
        self.printmsg('Returning bias phase for ASIC %i' % asic,verbosity=2)
        return self.asic_list[asic_idx].bias_phase()

    # otherwise, compare to see if Vbias is different between ASICs
    # check all the values against the first non-None ASIC object
    bp0 = None
    warning_msg = 'WARNING! Bias phase value not equal to that of ASIC %i at index %i: %f != %f'
    warn = False
    for asic_idx,asicobj in enumerate(self.asic_list):
        if asicobj is None: continue
        
        bp = asicobj.bias_phase()
        if bp is None:
            warn = True
            self.printmsg('WARNING! no bias phase info.',verbosity=3)
            continue
            
        if bp0 is None:
            bp0 = bp
            compare_idx = asic_idx
            continue

        if len(bp)!=len(bp0):
            warn = True
            self.printmsg('WARNING! not the same number of samples between ASICs',verbosity=3)
            continue
            
        for idx,chk in enumerate(bp==bp0):
            if not chk and bp[idx]!=-bp0[idx]:
                warn = True
                self.printmsg(warning_msg % (compare_idx+1,idx,bp[idx],bp0[idx]),verbosity=5)
                    

    if warn:
        self.printmsg('WARNING! The bias phase is different between the ASICs.  To see where, please set verbosity>2 and rerun bias_phase()')
    return bp0

def timeline_vbias(self,asic=None,TES=None):
    '''
    wrapper to get the bias voltage for the TES timeline
    '''
    args =self.args_ok(TES,asic,asic_only_required=True)
    if args is not None:
        TES,asic = args
        asic_idx = asic-1
        self.printmsg('Returning bias voltage for ASIC %i' % asic,verbosity=2)
        return self.asic_list[asic_idx].timeline_vbias

    # otherwise, compare to see if Vbias is different between ASICs
    vb0 = None
    warning_msg = 'WARNING! Bias voltage value not equal to that of ASIC %i at index %i: %f != %f'
    warn = False
    for asic_idx,asic_obj in enumerate(self.asic_list):
        if asic_obj is not None:
            vb = asic_obj.timeline_vbias
            if vb0 is None:
                vb0 = vb
                compare_idx = asic_idx
                continue
            
            for idx,chk in enumerate(vb==vb0): 
                if not chk:
                    warn = True
                    self.printmsg(warning_msg % (compare_idx+1,idx,vb[idx],vb0[idx]),verbosity=3)
    if warn:
        self.printmsg('WARNING! The Vbias is different between the ASICs.  To see where, please set verbosity>2 and rerun timeline_vbias()')
    return vb0

def sample_period(self,TES=None,asic=None):
    '''
    wrapper to get the sample period for an asic
    '''
    args =self.args_ok(TES,asic,asic_only_required=True)
    if args is None:return
    TES,asic = args

    asic_idx = asic-1
    return self.asic_list[asic_idx].sample_period()

def timeline_array(self,asic=None,TES=None):
    '''
    wrapper to get the timeline array for an asic
    '''
    args =self.args_ok(TES,asic,asic_only_required=True)
    if args is None:return
    TES,asic = args

    asic_idx = asic-1
    return self.asic_list[asic_idx].timeline_array()

def nsamples4interpolation(nsamples,epsilon=0.01,verbosity=0):
    '''
    find a good number of samples to use for the interpolated results
    it should be bigger than nsamples, but not too much
    try to get a power of 2, or nearly
    '''
    solution = {}
    if nsamples%2==0:
        fallback_nsamples = nsamples
    else:
        fallback_nsamples = nsamples + 1
    solution['fallback'] = fallback_nsamples

    max_delta = int(epsilon*nsamples)
    bitpower = int(np.log(fallback_nsamples)/np.log(2)) + 1
    new_nsamples = fallback_nsamples
    
    primes = np.array((3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
                       101,103,107,109,113,127,131,137,139,149,151,157,163,167,173))

    attempted_nsamples = []
    delta = max_delta + 1
    counter = 1
    gotit = False
    while not gotit and counter<10000:
        for prime in primes:
            for primepower in range(4):
                new_nsamples = (prime**primepower) * (2**bitpower)
                attempted_nsamples.append(new_nsamples)
                delta = new_nsamples - nsamples
                if delta<0: continue
                if delta<max_delta:
                    gotit = True
                    solution['bitpower'] = bitpower
                    solution['prime'] = prime
                    solution['primepower'] = primepower
                    break
                
            if gotit: break            
        if gotit: break

        for prime1 in primes:
            for prime2 in primes:
                new_nsamples = (prime1*prime2) * (2**bitpower)
                attempted_nsamples.append(new_nsamples)
                delta = new_nsamples - nsamples
                if delta<0: continue
                if delta<max_delta:
                    gotit = True
                    solution['bitpower'] = bitpower
                    solution['prime1'] = prime1
                    solution['prime2'] = prime2
                    break
            if gotit: break            
        if gotit: break
                                
        counter += 1
        bitpower -=1
            
    new_nsamples = nsamples + delta
    if not gotit: new_samples = fallback_nsamples
    
    if verbosity>0:
        if gotit:
            print('suggest nsamples=%i which is ' % new_nsamples,end=' ')
            if 'primepower' in solution.keys():
                print('2^%i*%i' % (bitpower,prime))
            else:
                print('2^%i*%i*%i' % (bitpower,prime1,prime2))
        else:
            print("didn't find anything... returning default nsamples=%i" % new_nsamples)

    
    solution['new_nsamples'] = new_nsamples
    solution['attempted'] = np.array(attempted_nsamples)
    return solution

def tod(self,axistype='pps',indextype='TES'):
    '''
    return a tuple containing the time axis, and the array of all TES timelines
    this is the timeaxis for all ASIC interpolated to the first ASIC

    with indextype=='TES' (default), the TOD array indexes are in the order of TES
    with indextype=='QS', the TOD array indexes are in the order compatible with qubicsoft simulations
    see for more info: http://qubic.in2p3.fr/wiki/pmwiki.php/TD/TEStranslation
    '''
    FPidentity = make_id_focalplane()
    
    timeaxis_list = []
    nsamples_list = []
    t0_list = []
    tfinal_list = []
    asic_ctr = 0
    qsidx_minlist = []
    for asic_idx,asicobj in enumerate(self.asic_list):
        if asicobj is None: continue
        asic_ctr += 1
        asic = asicobj.asic
        
        tstamps = asicobj.timeaxis(datatype='sci',axistype=axistype)
        timeaxis_list.append(tstamps)
        nsamples_list.append(tstamps.size)
        t0_list.append(tstamps[0])
        tfinal_list.append(tstamps[-1])

        # if we are using the indextype='QS' option
        # we have to adjust the QSindex so it starts at 0 regardless of the location of the detector array in the focal plane
        # for example, for the TD, it's in Quadrant-3, and the QSindex begins at 496
        mask = (FPidentity.ASIC==asic) & (FPidentity.TES>0)
        qsidx_minlist.append(FPidentity.QSindex[mask].min())

    # this won't work for scattered ASICs.
    # Hopefully, we only ever have to deal with the TD (2 asics in quadrant-3) or the Full Instrument
    QSoffset = min(qsidx_minlist)
    self.printmsg('tod(): QSoffset=%i' % QSoffset,verbosity=3)
        

    # choose number of samples for interpolation
    nsample_finder = nsamples4interpolation(max(nsamples_list),verbosity=self.verbosity)
    nsamples = nsample_finder['new_nsamples']

    # make the new timestamps axis
    t_tod = np.empty(nsamples,dtype=float)
    t0 = min(t0_list)
    tfinal = max(tfinal_list)
    t_tod = t0 + (tfinal-t0)*np.arange(nsamples)/(nsamples-1)

    # prepare the numpy array
    n_asics = copy(asic_ctr)
    if indextype.upper().find('QS')==0:
        ndets = n_asics*124
    else:
        ndets = n_asics*NPIXELS
    todarray = np.empty((ndets,nsamples),dtype=float)
    todarray[:] = np.nan

    self.printmsg('number of ASICs with data: %i' % asic_ctr,verbosity=3)
            
    # interpolate the timelines
    asic_ctr = 0
    for asic_idx,asicobj in enumerate(self.asic_list):
        if asicobj is None: continue
        asic_ctr += 1
        
        asic = asicobj.asic
        self.printmsg('asic index = %i, asic counter = %i, asic = %i' % (asic_idx,asic_ctr,asic),verbosity=3)
        tline_array = self.asic(asic).timeline_array()
        
        tstamps = timeaxis_list[asic_ctr-1]
        for TESidx in range(NPIXELS):
            TESnum = TESidx + 1
            tline_interp = np.interp(t_tod, tstamps, tline_array[TESidx,:])
            if indextype.upper().find('QS')==0:
                mask = (FPidentity.ASIC==asic) & (FPidentity.TES==TESnum)
                if mask.sum()!=1:
                    self.printmsg('ERROR! Dark pixel? Could not find the index for ASIC%02i TES%03i' % (asic,TESnum),verbosity=3)
                    continue
                QSindex = FPidentity.QSindex[mask][0]
                self.printmsg('tod(): QSindex=%i' % QSindex,verbosity=3)
                tod_index = QSindex - QSoffset
            else:
                tod_index = (asic_ctr-1)*NPIXELS+TESidx
            self.printmsg('tod(): tod_index = %03i' % tod_index,verbosity=3)
            todarray[tod_index,:] = tline_interp            

    return (t_tod,todarray)

def timeline(self,TES=None,asic=None):
    '''
    wrapper to get a timeline for a TES from an asic object
    '''
    args =self.args_ok(TES,asic)
    if args is None:return
    TES,asic = args
    return self.asic(asic).timeline(TES=TES)


def plot_timeline(self,TES=None,asic=None,timeaxis='pps',ax=None,fontsize=12,
                  plot_bias=False,
                  plot_calsource=False,
                  plot_azel=False,
                  plot_Tbath=False,
                  plot_hwp=False,
                  plot_raw=False):
    '''
    wrapper to plot timeline of the asic object 
    and plotted together with various possible housekeeping information
    '''
    args =self.args_ok(TES,asic)
    if args is None:return
    TES,asic = args
    self.printmsg('plotting timeline for asic=%i, TES=%i' % (asic,TES),verbosity=2)


    if ax is None:
        newplot = True
    else:
        newplot = False

    if plot_raw:
        ret = self.asic(asic).plot_raw(TES=TES,timeaxis=timeaxis,ax=ax)
        curves = [ret['curves'][0],ret['curves'][-1]]
    else:
        ret = self.asic(asic).plot_timeline(TES=TES,plot_bias=plot_bias,timeaxis=timeaxis,ax=ax,fontsize=fontsize)
        curves = ret['curves']
        
    ax = ret['ax']
    fig = ax.get_figure()
    if plot_calsource:
        t_src,v_src = self.calsource()
        d_src = np.empty(len(t_src),dtype=dt.datetime)
        for idx,tstamp in enumerate(t_src):
            d_src[idx] = dt.datetime.utcfromtimestamp(tstamp)            
        axsrc = ax.twinx()
        curvesrc = axsrc.plot(d_src,-v_src,color='red',label='calibration source')
        axsrc.text(0.5,1.01,self.calsource_infotext(),va='bottom',ha='center',fontsize=fontsize,transform=axsrc.transAxes)
        curves += curvesrc

    if plot_azel:
        az = self.azimuth()
        el = self.elevation()
        t_hk = self.timeaxis(datatype='platform')
        d_hk = np.empty(len(t_hk),dtype=dt.datetime)
        for idx,tstamp in enumerate(t_hk):
            d_hk[idx] = dt.datetime.utcfromtimestamp(tstamp)            
        
        axaz = ax.twinx()
        curves += axaz.plot(d_hk,az,color='red',label='azimuth')
        axaz.tick_params(axis='y',labelcolor='red')
        axaz.set_ylim(az.min(),az.max())
        axaz.set_ylabel('azimuth',rotation=270,ha='right',va='bottom',color='red')

        axel = ax.twinx()
        curves += axel.plot(d_hk,el,color='green',label='elevation')
        axel.tick_params(axis='y',labelcolor='green',pad=80)
        axel.set_ylim(el.min(),el.max())
        axel.set_ylabel('elevation',rotation=270,ha='left',va='bottom',color='green')

    if plot_Tbath:
        Tbath = self.Tbath[1]
        t_Tbath = self.Tbath[0]
        d_Tbath = np.empty(len(t_Tbath),dtype=dt.datetime)
        for idx,tstamp in enumerate(t_Tbath):
            d_Tbath[idx] = dt.datetime.utcfromtimestamp(tstamp)            
        
        axbath = ax.twinx()
        curves += axbath.plot(d_Tbath,Tbath,color='magenta',label='T$_\\mathrm{bath}$')
        axbath.tick_params(axis='y',labelcolor='magenta')
        axbath.set_ylim(Tbath.min(),Tbath.max())
        axbath.set_ylabel('T$_\\mathrm{bath}$ / K',rotation=270,ha='right',va='bottom',color='magenta')
    
    if plot_hwp:
        hwp = self.hwp_position()
        t_hwp = self.timeaxis(datatype='hwp')
        d_hwp = np.empty(len(t_hwp),dtype=dt.datetime)
        for idx,tstamp in enumerate(t_hwp):
            d_hwp[idx] = dt.datetime.utcfromtimestamp(tstamp)            
        
        axhwp = ax.twinx()
        curves += axhwp.plot(d_hwp,hwp,color='purple',label='HWP')
        axhwp.tick_params(axis='y',labelcolor='purple')
        axhwp.set_ylim(hwp.min(),hwp.max())
        axhwp.set_ylabel('HWP position',rotation=270,ha='right',va='bottom',color='purple')

    

    labs = [l.get_label() for l in curves]
    ax.legend(curves, labs, loc='upper right',facecolor='white',framealpha=0.7)
    pngname = '%s_TES%03i_ASIC%02i_timeline.png' % (self.dataset_name,TES,asic)
    ret['pngname'] = pngname
    if newplot: fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
            
    return ret
    
def plot_timeline_focalplane(self,xwin=True):
    '''
    plot all the timelines in the focal plane
    '''

    args= {}
    args['title'] = 'QUBIC Focal Plane: %s' % self.dataset_name
    args['xwin'] = xwin
    subttl_list = []
    obsdates = []
    for idx,asic_obj in enumerate(self.asic_list):
        if asic_obj is None: continue
        if not asic_obj.exist_timeline_data(): continue

        key = 'ASIC%i' % (idx+1)
            
        subttl_list.append(asic_obj.infotext())

        args[key] = asic_obj.timeline_array()
        obsdates.append(asic_obj.obsdate)
            
    args['subtitle'] = '\n'.join(subttl_list)
    args['obsdate'] = min(obsdates)
    args['pngname'] = 'QUBIC_focal_plane_timeline_%s.png' % args['obsdate'].strftime('%Y%m%dT%H%M%S')

    plot_fp(args)
    return args


#### I-V methods 
def plot_iv(self,TES=None,asic=None,multi=False,xwin=True,best=True):
    '''
    wrapper to plot I-V of the asic object
    '''
    args =self.args_ok(TES,asic)
    if args is None:return
    TES,asic = args
    return self.asic(asic).plot_iv(TES=TES,multi=multi,xwin=xwin,best=best)

def plot_responsivity(self,TES=None,asic=None,xwin=True,npts_region=500,window_size=51,filter_sigma=10,
                      plot_model=True,rmax=None,rmin=None):
    '''
    wrapper to plot responsivity of a TES
    '''
    args =self.args_ok(TES,asic)
    if args is None:return
    TES,asic = args
    return self.asic(asic).plot_responsivity(TES,xwin,npts_region,window_size,filter_sigma,plot_model,rmax,rmin)


        
def plot_iv_focalplane(self,labels=True,xwin=True):
    '''
    plot all the I-V curves in the focal plane
    '''

    args= {}
    args['title'] = 'QUBIC Focal Plane I-V curves: %s' % self.dataset_name
    if not labels: args['nolabels'] = True
    args['xwin'] = xwin
    
    subttl_list = []
    obsdates = []
    ngood = []
    tot_ngood = 0
    tot_npixels = 0
    for idx,asic_obj in enumerate(self.asic_list):
        if asic_obj is None: continue
        if not asic_obj.exist_iv_data(): continue
        obsdates.append(asic_obj.obsdate)

        key = 'ASIC%i' % (idx+1)
            
        subttl_list.append(asic_obj.infotext())

        bias,adu = asic_obj.best_iv_curve()
        args[key] = adu
            
        keyx = '%s x-axis' % key
        args[keyx] = bias

        keygood = '%s good' % key
        args[keygood] = asic_obj.is_good_iv()

        keybg = '%s bg' % key
        args[keybg] = asic_obj.turnover()
        filtersummary = asic_obj.filterinfo()
        for idx,f in enumerate(filtersummary):
            if f['ignore_turnover']:
                args[keybg][idx] = None

        ngood = asic_obj.ngood()
        if ngood is not None:
            tot_ngood += ngood
            subttl_list.append('%i flagged as bad pixels : yield = %.1f%%' %
                               (asic_obj.NPIXELS-ngood,100.0*ngood/asic_obj.NPIXELS))
        tot_npixels += asic_obj.NPIXELS

    if tot_npixels>0:
        subttl_list.append('overall yield %i/%i = %.1f%%' % (tot_ngood,tot_npixels,100.0*tot_ngood/tot_npixels))
    args['subtitle'] = '\n'.join(subttl_list)

    if len(obsdates)>0:
        args['obsdate'] = min(obsdates)
    args['pngname'] = 'QUBIC_focal_plane_IV_%s.png' % args['obsdate'].strftime('%Y%m%dT%H%M%S')

    plot_fp(args)
    return args

def ngood(self):
    '''
    return the number of good TES
    '''
    tot_ngood = 0
    for asic_obj in self.asic_list:
        if asic_obj is not None:
            tot_ngood += asic_obj.ngood()
    return tot_ngood

def is_good(self):
    '''
    return a boolean array of evaluation for each TES
    '''
    is_good_list = []
    for asic_obj in self.asic_list:
        if asic_obj is not None:
            asic_is_good_list = asic_obj.is_good_iv()
            if asic_is_good_list is None: continue
            is_good_list.append(asic_is_good_list)
    return np.concatenate(is_good_list)
    
    
def filter_iv_all(self,
                  R1adjust=1.0,
                  R1_limit=10,
                  residual_limit=3.0,
                  abs_amplitude_limit=0.01,
                  rel_amplitude_limit=0.1,
                  bias_margin=-3,
                  ignore_turnover=False,
                  jumplimit=None,
                  fitfunction='COMBINED',
                  Vsuper=None,
                  Vnormal=None,
                  istart=None,
                  iend=None,
                  flip=False):
    '''
    find which TES are good
    '''

    for asic_obj in self.asic_list:
        if asic_obj is not None:
            asic_obj.filter_iv_all(R1adjust=R1adjust,
                                   R1_limit=R1_limit,
                                   residual_limit=residual_limit,
                                   abs_amplitude_limit=abs_amplitude_limit,
                                   rel_amplitude_limit=rel_amplitude_limit,
                                   ignore_turnover=ignore_turnover,
                                   bias_margin=bias_margin,
                                   jumplimit=jumplimit,
                                   fitfunction=fitfunction,
                                   Vsuper=Vsuper,
                                   Vnormal=Vnormal,
                                   istart=istart,
                                   iend=iend,
                                   flip=flip)

    return
    

def make_iv_tex_report(self,tableonly=False):
    '''
    produce a latex report on I-V curves
    '''
    for asic_obj in self.asic_list:
        if asic_obj is not None:
            asic_obj.make_iv_tex_report(tableonly=tableonly)
    return


def make_iv_report(self):
    '''
    do all the business to generate the I-V report document
    '''
    for asic_obj in self.asic_list:
        if asic_obj is not None:
            asic_obj.make_iv_report()
    return
