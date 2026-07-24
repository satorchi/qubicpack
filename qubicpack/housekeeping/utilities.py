'''
$Id: utilities.py<housekeeping>
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Wed 22 Jul 2026 10:20:51 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

utilities for dealing with QUBIC housekeeping data
'''
import os,re
import datetime as dt
TZUTC = dt.timezone.utc

import numpy as np
from satorchipy.datefunctions import str2dt
from .. import __file__
from ..utilities import hostname
pkg_dir = os.path.dirname(__file__)

qc_hk_dir = '/home/qubic/data/temperature/broadcast'

def read_hk_labels():
    '''
    read the housekeeping labels associated with each HK data file
    assign these for use in QubicStudio datasets and also directly from HK files
    '''
    label_file = os.sep.join([pkg_dir,'data','TD_TEMPERATURE_LABELS.txt'])
    if not os.path.isfile(label_file):
        print('could not find temperature labels: %s' % label_file)
        return None
    
    h = open(label_file)
    lines = h.read().split('\n')
    h.close()
    qs_labels = {} # compatible with QubicStudio datasets
    hk_labels = {} # compatible with housekeeping raw files
    for line in lines:
        if line=='': continue
        keyval = line.split('=')
        if len(keyval)<2: continue
        key = keyval[0].strip()
        val = keyval[1].strip()

        # QubicStudio assigned different keynames to these
        if key.find('AVS')==0:
            qskey = key.upper()
            hkkey = key
            qs_labels[qskey] = val
            hk_labels[hkkey] = val
            continue
        
        if key.find('HEATER')==0:
            heater_num = int(key.replace('HEATER',''))
            for meas in ['Amp','Volt']:
                qskey = 'Heaters_%s_%i' % (meas,heater_num - 1)
                hkkey = 'HEATER%i_%s' % (heater_num,meas)
                qs_labels[qskey] = val
                hk_labels[hkkey] = val
            continue
        
        if key.find('TEMPERATURE')==0:
            temp_num = int(key.replace('TEMPERATURE',''))
            qskey = 'Temp_%i' % (temp_num - 1)
            hkkey = 'TEMPERATURE%02i' % temp_num
            qs_labels[qskey] = val
            hk_labels[hkkey] = val
            continue

        if key.find('PRESSURE')==0:
            temp_num = int(key.replace('PRESSURE',''))
            qskey = 'Pressure_%i' % (temp_num - 1)
            hkkey = 'PRESSURE%i' % temp_num
            qs_labels[qskey] = val
            hk_labels[hkkey] = val
            continue
        

        qs_labels[key] = val
        hk_labels[key] = val

    labels = {}
    labels['housekeeping'] = hk_labels
    labels['QubicStudio'] = qs_labels
    return labels


def read_hk_file(filename):
    '''
    return the date,data from the Housekeeping broadcast file
    '''
    if not os.path.isfile(filename):
        print('ERROR! File not found: %s' % filename)
        return None,None

    h = open(filename,'r')
    txt = h.read()
    h.close()
    txt_clean = txt.replace('\x00','').replace('inf','64218') # infinity = 0xfada
    txt_cleanclean = re.sub('\n.*\\.[0-9]*\\..*\n','\n',txt_clean)
    lines = txt_cleanclean.split('\n')
    del(lines[-1])

    # try to use numpy loadtxt which is fastest
    try:
        dat = np.loadtxt(lines)
    except:
        dat = None

    if dat is not None:
        ncols = dat.shape[-1]
        npts = dat.shape[0]
        t = dat[:,0]
        v = dat[:,1]
        if ncols==3:
            onoff = dat[:,2]
        else:
            onoff = np.zeros(npts,dtype=bool)
        return t,v,onoff
        
    
    npts = len(lines)
    t = np.zeros(npts)
    v = np.zeros(npts)
    onoff = np.zeros(npts,dtype=bool)
    idx=0
    badpattern = re.compile('[a-zA-Z]')
    for line_idx,line in enumerate(lines):
        cols = line.strip().split()
        if len(cols)<2: continue
        tstamp_str = cols[0]
        val_str = cols[1]
        if badpattern.match(tstamp_str): continue
        if badpattern.match(val_str): continue
        try:
            tstamp = float(tstamp_str)
            reading = eval(val_str)
        except:
            print("ERROR! Couldn't read line: %i) %s" % (line_idx+1,line))
            continue

        if tstamp>4e9: continue
        if tstamp<1.4e9: continue
        v[idx] = reading
        t[idx] = tstamp
        

        if len(cols)>2:
            if cols[2]=='ON': onoff[idx] = True
        idx+=1
        
    if idx<npts:
        t = t[0:idx]
        v = v[0:idx]
        onoff = onoff[0:idx]
        print('%s: idx,npts = %i,%i' % (filename,idx,npts))
    return t,v,onoff

def find_pt_start(flag):
    '''
    find the most recent start date of the pulse tubes

    argument:
       flag : this is the dictionary of events
    '''
    for labeldate in flag.keys():
        match = re.search('(PT(1|2|s)|PT|pt|(P|p)ulse (T|t)ube).* ((S|s)tart|(ON|on))',flag[labeldate])
        if match: pt_start_date = labeldate

    print('most recent pulse tube start: %s' % pt_start_date.strftime('%Y-%m-%d %H:%M:%S'))
    return pt_start_date

def read_compressor_log(filename):
    '''
    read the pulse tube compressor data
    '''

    if not os.path.isfile(filename):
        print('file not found: %s' % filename)
        return None

    h = open(filename,'r')
    txt = h.read()
    h.close()
    lines = txt.replace('\0','').split('\n')
    del(lines[-1])
    
    compressorlog = {}
    compressorlog['date'] = []
    compressorlog['timestamp'] = []

    for line in lines:
        if line.find('OFFLINE')>0: continue
        sample_line = line
        break
    
    col = sample_line.split(' ')
    for item in col:
        if item.find('=')<0: continue
        key = item.split('=')[0]
        compressorlog[key] = []
        
    for line in lines:
        col = line.split(' ')
        if col[0] is None: continue
        date = str2dt(col[0]).replace(tzinfo=TZUTC)
        if date is None:continue
        tstamp = date.timestamp()
        keyval = None
        for item in col[1:]:
            keyval = None
            if item.find('=')<0:continue
            keyval = item.split('=')
            key = keyval[0]
            val = float(keyval[1])
            compressorlog[key].append(val)
        if keyval is not None:
            compressorlog['date'].append(date)
            compressorlog['timestamp'].append(tstamp)


    compressorlog['timestamp'] = np.array(compressorlog['timestamp'])
    return compressorlog

def read_hk_flags(flagfile=None):
    '''
    read the list of events
    '''
    flagfiles = []
    
    default_flagfile = os.sep.join([pkg_dir,'data','events.txt'])
    if os.path.isfile(default_flagfile):
        flagfiles.append(default_flagfile)
    else:
        print('could not find the default event list: %s' % default_flagfile)

    if flagfile is not None and os.path.isfile(flagfile):
        flagfiles.append(flagfile)

    flag = {}
    for flagfile in flagfiles:
        h = open(flagfile,'r')
        lines = h.read().split('\n')
        h.close()

        for line in lines:
            if line.find('=')<0: continue
            col = line.split('=')
            datekey = str2dt(col[0])
            if datekey is None: continue
            eventmsg = col[1]
            flag[datekey] = eventmsg

    return flag

def download_hk(basenames,hk_dir,remote_machine='qubic'):
    '''
    download using rsync from qubic-central
    '''
    if hostname.find('qubic-central')==0:
        print('NOT DOWNLOADING FROM qubic-central to qubic-central!')
        return

    if len(basenames)==0:
        print('No files to download!')
        return
    
    filenames = []
    for b in basenames:
        filenames.append(b+'.txt')
        
    print('download file: %s' % '\ndownload file: '.join(filenames))
    print('[%s] downloading...' % dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    if len(filenames)==1:
        cmd = 'rsync -Pavtz %s:%s/%s %s' % (remote_machine,qc_hk_dir,filenames[0], hk_dir)
    else:
        cmd = 'rsync -Pavtz %s:%s/{%s} %s' % (remote_machine,qc_hk_dir,','.join(filenames), hk_dir)
    os.system(cmd)
    print('[%s] files downloaded' % dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    return


    
