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
import sys,os,time

from qubicpack.qubicasic import qubicasic
from qubicpack.plot_fp import plot_fp

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
        txt += ' T$_\mathrm{bath}$=%.1fmK' % (1000*self.temperature)

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

    info_tstamp = conf['timestamp'][0]
    info_date = dt.datetime.utcfromtimestamp(info_tstamp)
    info['date'] = info_date

    devlist = ['calsource','modulator','amplifier']
    statuskeys = ['CalSource','Modulator','Amplifier']
    for idx,dev in enumerate(devlist):
        info[dev] = {}
        statuskey = statuskeys[idx]
        status = int(conf[statuskey][0])
        if status==1:
            info[dev]['status'] = 'ON'
        elif status==0:
            info[dev]['status'] = 'OFF'
        else:
            info[dev]['status'] = 'UNKNOWN'

    info['calsource']['frequency'] = conf['Cal_freq'][0]
    info['calsource']['synth_freq'] = conf['Syn_freq'][0]
    
    info['modulator']['frequency'] = conf['Mod_freq'][0]
    info['modulator']['amplitude'] = conf['Mod_ampl'][0]
    info['modulator']['duty_cycle'] = conf['Mod_duty'][0]
    info['modulator']['offset'] = conf['Mod_offs'][0]
    shape_idx = int(conf['Mod_shap'][0])
    shapes = ['square','sine','DC']
    info['modulator']['shape'] = shapes[shape_idx]

    ampkeys =  ['Amp_mode', 'Amp_lfreq', 'Amp_hfreq', 'Amp_coup', 'Amp_rang', 'Amp_gain']
    amp_infokeys = ['mode','filter low frequency','filter high frequency','coupling','dynamic range','gain']
    for idx,key in enumerate(ampkeys):
        if key in conf.keys():
            localkey = amp_infokeys[idx]
            info['amplifier'][localkey] = conf[key][0]

    # convert mode number to human understandable
    amp_modes = ["bypass",
                 "6db_low_pass",
                 "12db_low_pass",
                 "6db_high_pass",
                 "12db_high_pass",
                 "bandpass"]
    if 'mode' in info['amplifier'].keys():
        mode_idx = int(info['amplifier']['mode'])
        info['amplifier']['mode'] = amp_modes[mode_idx]

    # convert coupling to human understandable
    coupling_modes = ['GROUND','DC','AC']
    if 'coupling' in info['amplifier'].keys():
        mode_idx = int(info['amplifier']['coupling'])
        info['amplifier']['coupling'] = coupling_modes[mode_idx]

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
    info_rawlist = info_txt.split()
    info = {}

    info_tstamp = float(info_rawlist[0])
    info_date = dt.datetime.utcfromtimestamp(info_tstamp)
    info['date'] = info_date

    info['amplifier'] = {}
    info['modulator'] = {}
    info['calsource'] = {}

    if info_txt.find('busy')>0:
        for dev in ['amplifier','modulator','calsource']:
            info[dev]['status'] = 'busy'
        return info

    munits = ['mHz','mVpp','mVdc']
    units = ['GHz','Hz','Vpp','Vdc','%']
    for item in info_rawlist[2:]:
        cols = item.split(':')

        if len(cols)==1:
            info['calsource']['status'] = info_txt
            continue
        
        dev = cols[0]
        if dev=='lamp' or dev=='arduino' or dev=='synthesiser':
            continue

        val_list = cols[1].split('=')
        if len(val_list)==1:
            status_str = val_list[0].upper()
            if status_str=='ON' or status_str=='OFF':
                info[dev]['status'] = val_list[0]
            else:
                info[dev]['incomplete'] = val_list[0]
            continue

        parm = val_list[0]
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

        for unit in units:
            if val.find(unit)>1:
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
                txt_list.append(key+'='+modulator_units[key] % info['modulator'][key])
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
    self.printmsg('ASIC%i     Observation date: %s' % (asic_no,obsdate.strftime('%Y-%m-%d %H:%M:%S.%f')))
    if self.obsdate is None:
        self.printmsg('DEBUG: setting obsdate which was None: %s' % obsdate.strftime('%Y-%m-%d %H:%M:%S.%f'),verbosity=3)
        self.obsdate = obsdate
    if self.obsdate>obsdate:
        self.obsdate = obsdate
        self.printmsg('Observation date is not the same between ASIC data.  Using the earlier date',verbosity=3)

    endobs = self.asic_list[asic_ctr].endobs
    self.printmsg('ASIC%i Observation end date: %s' % (asic_no,endobs.strftime('%Y-%m-%d %H:%M:%S.%f')))
    if self.endobs is None:
        self.printmsg('DEBUG: setting endobs which was None: %s' % endobs.strftime('%Y-%m-%d %H:%M:%S.%f'),verbosity=3)
        self.endobs = endobs
    if self.endobs<endobs:
        self.endobs = endobs
        self.printmsg('End Observation is not the same between ASIC data.  Using the later date',verbosity=3)
        

    return

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
def args_ok(self,TES=None,asic=None):
    '''
    check if arguments are okay for the wrapper

    if TES or asic is None, then it is requested
    if TES or asic is 0, then it is ignored (argument not required)
    '''
    if asic is None:
        self.printmsg('Please give an asic number')
        return False

    asic_idx = asic - 1
    if asic>len(self.asic_list) or self.asic_list[asic_idx] is None:
        self.printmsg('No data for ASIC %i' % asic)
        return False

    if TES is None:
        self.printmsg('Please give a TES number')
        return False

    if TES=='no TES number required':
        return True
    
    if TES < 1:
        self.printmsg('Please give a valid TES number')
        return False

    return True

def asic(self,asic=None):
    '''
    return the requested asic object
    '''
    if not self.args_ok('no TES number required',asic):return
    asic_idx = asic - 1
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


#### timeline methods
def bias_phase(self,asic=0):
    '''
    return the bias on the detectors
    '''
    if not self.args_ok(TES='no TES number required',asic=asic):return None
    asic_idx = asic-1
    if asic_idx >= 0:
        self.printmsg('Returning bias phase for ASIC %i' % (asic_idx+1),verbosity=2)
        return self.asic_list[asic_idx].bias_phase()

    # check all the values against the first non-None ASIC object
    bp0 = None
    warning_msg = 'WARNING! Bias phase value not equal to that of ASIC %i at index %i: %f != %f'
    warn = False
    for asic_idx,asicobj in enumerate(self.asic_list):
        if asicobj is not None:
            bp = asicobj.bias_phase()
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

def timeline_vbias(self,asic=0):
    '''
    wrapper to get the bias voltage for the TES timeline
    '''
    if not self.args_ok(TES='no TES number required',asic=asic):return None
    asic_idx = asic-1
    if asic_idx >= 0:
        self.printmsg('Returning bias voltage for ASIC %i' % (asic_idx+1),verbosity=2)
        return self.asic_list[asic_idx].timeline_vbias
    
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

def sample_period(self,asic=None):
    '''
    wrapper to get the sample period for an asic
    '''
    if asic is None:
        self.printmsg('Please enter an asic number')
        return None

    asic_idx = asic-1
    return self.asic_list[asic_idx].sample_period()

def timeline_array(self,asic=None):
    '''
    wrapper to get the timeline array for an asic
    '''
    if asic is None:
        self.printmsg('Please enter an asic number')
        return None

    asic_idx = asic-1
    return self.asic_list[asic_idx].timeline_array()
    

def timeline(self,TES=None,asic=None):
    '''
    wrapper to get a timeline for a TES from an asic object
    '''
    if not self.args_ok(TES,asic):return None
    asic_idx = asic-1
    return self.asic_list[asic_idx].timeline(TES=TES)


def plot_timeline(self,TES=None,asic=None,plot_bias=True,timeaxis='pps',ax=None,fontsize=12):
    '''
    wrapper to plot timeline of the asic object
    '''
    if not self.args_ok(TES,asic):return None
    asic_idx = asic-1
    return self.asic_list[asic_idx].plot_timeline(TES=TES,plot_bias=plot_bias,timeaxis=timeaxis,ax=ax,fontsize=fontsize)
    
def plot_timeline_focalplane(self):
    '''
    plot all the timelines in the focal plane
    '''

    args= {}
    args['title'] = 'QUBIC Focal Plane: %s' % self.dataset_name
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

    plot_fp(args)
    return args


#### I-V methods 
def plot_iv(self,TES=None,asic=None,multi=False,xwin=True,best=True):
    '''
    wrapper to plot I-V of the asic object
    '''

    if not self.args_ok(TES,asic):return    
    asic_idx = asic - 1
    return self.asic_list[asic_idx].plot_iv(TES=TES,multi=multi,xwin=xwin,best=best)

def plot_responsivity(self,TES=None,asic=None,xwin=True,npts_region=500,window_size=51,filter_sigma=10,
                      plot_model=True,rmax=None,rmin=None):
    if not self.args_ok(TES,asic):return None
    asic_idx = asic - 1
    return self.asic_list[asic_idx].plot_responsivity(TES,xwin,npts_region,window_size,filter_sigma,plot_model,rmax,rmin)


        
def plot_iv_focalplane(self,labels=True):
    '''
    plot all the I-V curves in the focal plane
    '''

    args= {}
    args['title'] = 'QUBIC Focal Plane I-V curves: %s' % self.dataset_name
    if not labels: args['nolabels'] = True
    
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
    
def filter_iv_all(self,
                  R1adjust=1.0,
                  R1_limit=10,
                  residual_limit=3.0,
                  abs_amplitude_limit=0.01,
                  rel_amplitude_limit=0.1,
                  bias_margin=0.2,
                  ignore_turnover=False,
                  jumplimit=None,
                  fitfunction='COMBINED',
                  Vsuper=None,
                  Vnormal=None,
                  istart=None,
                  iend=None):
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
                                   iend=iend)

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
