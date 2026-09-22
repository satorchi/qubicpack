'''
$Id: calsource.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Thu 03 Sep 2026 14:01:36 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

methods for reading/interpreting the calibration source and the carbon fibre data

most of these were originally in tools.py and fpmethods.py
'''
import os
import datetime as dt
from glob import glob
from .utilities import TZUTC
from satorchipy.datefunctions import utcfromtimestamp

def read_calsource_fits(self,hdu):
    '''
    read the calibration source data from the given HDU of a fits file
    '''
    
    self.hk['CALSOURCE'] = {}
    self.hk['CALSOURCE']['timestamp'] = hdu.data.field(0)
    self.hk['CALSOURCE']['Value'] = hdu.data.field(1)
    
    return

def read_calsource_infofile(self,datadir):
    '''
    read the calsource information from CALINFO.txt
    we set this up in a way that is compatible with earlier versions
    '''
    calinfo_file = os.sep.join([datadir,'Hks','CALINFO.txt'])
    if not os.path.isfile(calinfo_file):
        self.printmsg('WARNING! Did not find calsource information file: %s' % calinfo_file,verbosity=1)
        return False

    h = open(calinfo_file,'r')
    calinfo_rawtxt = h.read()
    h.close()

    # replace newlines
    calinfo_txt = calinfo_rawtxt.replace('\n',' ')

    if 'CALSOURCE-CONF' not in self.hk.keys():
        self.hk['CALSOURCE-CONF'] = {}
    if 'MsgStr' not in self.hk['CALSOURCE-CONF'].keys():
        self.hk['CALSOURCE-CONF']['MsgStr'] = ['INIT MSG STR']

    self.hk['CALSOURCE-CONF']['MsgStr'][0] = calinfo_txt
    return True

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
        info_date = utcfromtimestamp(info_tstamp)
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

    # get the various timestamps:  info received, info sent, command received
    # date is considered to be "info sent"
    date_list = []
    for tstamp_str in info_rawlist:
        if tstamp_str.find(':')>=0: continue
        try:
            tstamp = float(tstamp_str)
        except:
            info[tstamp_str] = 'incomplete'
            continue
        date = utcfromtimestamp(tstamp)
        date_list.append(date)

    start_idx = len(date_list)
    if len(date_list)==2:
        info['date'] = date_list[0]
        info['command received'] = date_list[1]
    elif len(date_list)>=3:
        info['date'] = date_list[1]
        info['command received'] = date_list[2]
        info['info received'] = date_list[0]
    elif len(date_list)==0:
        info['date'] = self.obsdate
    else:
        info['date'] = date_list[0]

    # create dictionaries for each device
    for dev in device_list:
        info[dev] = {}

    if info_txt.find('busy')>0:
        for dev in device_list:
            info[dev]['status'] = 'busy'
        return info

    munits = ['mHz','mVpp','mVdc']
    units = ['GHz','HZ','Hz','hz','Vpp','Vdc','V','%']
    for item in info_rawlist:
        if item.find(':')<0: continue
        
        cols = item.split(':')

        if len(cols)==1:
            info['calsource']['status'] = info_txt
            continue
        
        dev = cols[0]

        #### 2026-09-22 14:58:46 I don't remember why I skipped all these.
        # if dev=='lamp' or dev=='arduino' or dev=='synthesiser' or dev=='synthesizer':
        #    continue

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
            try:
                info[dev][parm] = int(val)
            except:
                info[dev][parm] = val
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

def calsource(self):
    '''
    return the calibration source data
    '''    
    if 'CALSOURCE' not in self.hk.keys():
        self.printmsg('No calibration source data',verbosity=2)
        return None, None
    
    t_src = self.hk['CALSOURCE']['timestamp']
    data_src = self.hk['CALSOURCE']['Value']
    return t_src,data_src

def find_calsource(self,datadir):
    '''
    try to find, and then read the calsource file corresponding to the dataset
    '''
    # look for files within the last hour, and then take the closest one to the start time
    # the files are in FITS format as of Wed 10 Apr 2019 10:21:35 CEST
    self.printmsg('trying to find calsource data corresponding to %s' % self.dataset_name,verbosity=2)

    if self.obsdate is None:
        self.printmsg('No date for observation!',verbosity=1)
        return

    # calsource directory is normally two up
    calsource_dir = '%s/calsource' % os.path.dirname(os.path.dirname(datadir))
    filetype = 'calsource'
    datadir = calsource_dir
    search_start = self.obsdate - dt.timedelta(minutes=30)
    pattern = []
    pattern.append('%s/calsource_%s*.fits' % (calsource_dir,search_start.strftime('%Y%m%dT%H')))
    pattern.append('%s/calsource_%s*.fits' % (calsource_dir,self.obsdate.strftime('%Y%m%dT%H')))
    files = []
    for p in pattern:
        files += glob(p)
    if len(files)==0:
        self.printmsg('No %s data found in directory: %s' % (filetype,datadir),verbosity=1)
        return
    files.sort()

    # find the file which starts before and nearest to obsdate
    filename = None
    file_delta = 1e6
    for f in files:
        basename = os.path.basename(f)
        file_date = dt.datetime.strptime(basename,'calsource_%Y%m%dT%H%M%S.fits').replace(tzinfo=TZUTC)
        delta = (self.obsdate - file_date).total_seconds()
        if np.abs(delta)<file_delta:
            file_delta = np.abs(delta)
            filename = f

    if file_delta>30:
        self.printmsg('Did not find a corresponding calsource file.')
        return
    
    self.printmsg('found calsource file which started %.1f seconds before the data acquisition' % file_delta)
    self.printmsg('reading calsource file: %s' % filename)
    hdulist=pyfits.open(filename)
    nhdu=len(hdulist)
    if nhdu!=2:
        self.printmsg("This doesn't look like a calsource file!")
        hdulist.close()
        return
    hdu = hdulist[1]
    if 'EXTNAME' not in hdu.header.keys()\
       and hdu.header['EXTNAME']!='CALSOURCE':
        self.printmsg("This is not a calsource FITS file!")
        hdulist.close()
        return
    
    self.read_calsource_fits(hdu)
    hdulist.close()
    return

