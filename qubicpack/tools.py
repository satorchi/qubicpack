#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
$Id: tools.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

tools for reading/writing and organizing data

"""
import numpy as np
import sys,os,time,subprocess,struct,re
import datetime as dt
from glob import glob
import pickle
from collections import OrderedDict
from astropy.io import fits as pyfits
from qubicpack.utitlities import obsmount_implemented

qubicasic_hk_keys = ['Apol',
                     'CN',
                     'ColumnEnd',
                     'ColumnStart',
                     'FLL_D',
                     'FLL_I',
                     'FLL_P',
                     'FLL_State',
                     'Feedback-DAC-values',
                     'FrequencyASICSerialLink',
                     'NETQUIC synchro error',
                     'NbSamplesPerSum',
                     'Offset-DAC-values',
                     'Raw-mask',
                     'Relays state',
                     'RowEnd',
                     'RowStart',
                     'Spol',
                     'TES Sinus phase',
                     'TESAmplitude',
                     'TESFreq',
                     'TESOffset',
                     'TESShapeMode',
                     'Undersampling raw mode',
                     'Vicm',
                     'Vocm',
                     'nsample',
                     'testPatternMode']


def debugmsg(self,msg,verbosity=3):
    if verbosity<=self.verbosity:
        if self.logfile is None:
            print('DEBUG %s : %s' % (dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC'),msg))
        else:
            self.writelog('DEBUG: %s' % msg)
    return

def printmsg(self,msg,verbosity=1):
    '''
    print a message to screen
    '''
    if verbosity<=self.verbosity:
        print(msg)
    return

def read_date_from_filename(self,filename):
    '''
    read the date from the filename. 
    this hasn't been used in a long time
    '''
    try:
        datestr=filename.split('_')[-1].split('.')[0]
        date=dt.datetime.strptime(datestr,'%Y%m%dT%H%M%SUTC')
    except:
        date=None
    return date

def data_subdir(self):
    '''
    make a subdirectory for output files based on the date of the data acquisition
    '''
    if not isinstance(self.obsdate,dt.datetime):
        self.printmsg('ERROR! No date for this data.')
        return None

    if not isinstance(self.datadir,str):
        datadir=self.assign_datadir()
        if datadir is None:return None

    subdir=self.obsdate.strftime('%Y/%Y%m%d')
    fullpath='%s%s%s' % (self.datadir,os.sep,subdir)
    # make the subdirectory if necessary
    if not os.path.exists(fullpath):
        cmd='mkdir -p %s' % fullpath
        proc=subprocess.Popen(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out,err=proc.communicate()

    if not os.path.exists(fullpath):
        self.printmsg('ERROR! Could not create subdirectory: %s' % fullpath)
        return None
    
    return subdir

def output_filename(self,rootname):
    '''
    assign a filename for plots, data, etc with the full path
    '''
    if not isinstance(rootname,str):
        return None

    if not isinstance(self.datadir,str):
        self.assign_datadir()

    if not isinstance(self.datadir,str):
        self.printmsg('ERROR! no appropriate location for file save.')
        return None

    subdir=self.data_subdir()
    if isinstance(subdir,str):
        filename='%s%s%s%s%s' % (self.datadir,os.sep,subdir,os.sep,rootname)
    else:
        filename='%s%s%s' % (self.datadir,os.sep,rootname)

    return filename


def keyvals(self):
    '''
    assign the FITS keyword values for the primary header
    the keyword descriptions are done in assign_variables.py
    '''
    keyvals={}
    keyvals['TELESCOP']='QUBIC'
    keyvals['OBSERVER']=self.observer
    keyvals['DATE-OBS']=self.obsdate.strftime('%Y-%m-%d %H:%M:%S UTC')
    keyvals['END-OBS']=self.endobs.strftime('%Y-%m-%d %H:%M:%S UTC')
    keyvals['NSAMPLES']=self.nsamples
    keyvals['INT-TIME']=self.tinteg
    keyvals['NPIXELS'] =self.NPIXELS
    keyvals['NPIXSAMP']=self.NPIXELS_sampled
    keyvals['ASIC']    =self.asic
    keyvals['QUBIC-IP']=self.QubicStudio_ip
    keyvals['NCYCLES'] =self.nbiascycles
    keyvals['CYCBIAS'] =self.cycle_vbias
    keyvals['TES_TEMP']=self.temperature
    keyvals['BIAS_MOD']=self.bias_frequency
    keyvals['BIAS_MIN']=self.min_bias
    keyvals['BIAS_MAX']=self.max_bias
    keyvals['BIAS_FAC']=self.bias_factor
    keyvals['BIASMODE']=self.bias_mode
    keyvals['FLL_STAT']=self.FLL_state
    keyvals['FLL_P']   =self.FLL_P
    keyvals['FLL_I']   =self.FLL_I
    keyvals['FLL_D']   =self.FLL_D
    keyvals['DET_NAME']=self.detector_name
    keyvals['R_FEEDBK']=self.Rfeedback
    keyvals['CHUNK']=self.chunk_size
    return keyvals
    

def write_fits(self):
    '''
    write data to file
    it could be timeline data or I-V data, or both
    '''
    datefmt='%Y%m%dT%H%M%SUTC'
    if self.obsdate is None: self.assign_obsdate()
    datestr=self.obsdate.strftime(datefmt)

    if self.endobs is None:
        self.endobs=self.obsdate

    keyvals=self.keyvals()
    prihdr = pyfits.Header()
    for key in self.fitsblurbs.keys():
        prihdr[key]=(keyvals[key],self.fitsblurbs[key])
    prihdu = pyfits.PrimaryHDU(header=prihdr)

    tbhdu0=None
    if not self.rawmask is None:
        col  = pyfits.Column(name='RawMask',format='I', unit='bitmask', array=self.rawmask)
        cols  = pyfits.ColDefs([col])
        tbhdu0 = pyfits.BinTableHDU.from_columns(cols)
    
    if isinstance(self.adu,np.ndarray):
        fitsfile=str('QUBIC_TES_%s.fits' % datestr)
        fitsfile_fullpath=self.output_filename(fitsfile)
        if os.path.exists(fitsfile_fullpath):
            self.printmsg('file already exists! %s' % fitsfile_fullpath)
            fitsfile=dt.datetime.utcnow().strftime('resaved-%Y%m%dT%H%M%SUTC__')+fitsfile
            fitsfile_fullpath=self.output_filename(fitsfile)
            self.printmsg('instead, saving to file: %s' % fitsfile_fullpath)

        fmtstr=str('%iD' % self.adu.shape[1])
        dimstr=str('%i' % self.adu.shape[0])
        #self.printmsg('format=',fmtstr)
        #self.printmsg('dim=',dimstr)
        col1  = pyfits.Column(name='V_tes', format=fmtstr, dim=dimstr, unit='ADU', array=self.adu)
        cols  = pyfits.ColDefs([col1])
        tbhdu1 = pyfits.BinTableHDU.from_columns(cols)

        col2  = pyfits.Column(name='V_bias',format='D', unit='V', array=self.vbias)
        cols  = pyfits.ColDefs([col2])
        tbhdu2 = pyfits.BinTableHDU.from_columns(cols)

        if tbhdu0 is None:
            thdulist = pyfits.HDUList([prihdu, tbhdu1, tbhdu2])
        else:
            thdulist = pyfits.HDUList([prihdu, tbhdu0, tbhdu1, tbhdu2])
        thdulist.writeto(fitsfile_fullpath)
        self.printmsg('FITS file written: %s' % fitsfile_fullpath)

    if self.exist_timeline_data():
        fitsfile=str('QUBIC_timeline_%s.fits' % datestr)
        fitsfile_fullpath=self.output_filename(fitsfile)
        if os.path.exists(fitsfile_fullpath):
            self.printmsg('file already exists! %s' % fitsfile_fullpath)
            fitsfile=dt.datetime.utcnow().strftime('resaved-%Y%m%dT%H%M%SUTC__')+fitsfile
            fitsfile_fullpath=self.output_filename(fitsfile)
            self.printmsg('instead, saving to file: %s' % fitsfile_fullpath)

        ntimelines=self.ntimelines()

        hdulist=[prihdu]
        if not tbhdu0 is None:hdulist.append(tbhdu0)
                    
        for timeline_index in range(ntimelines):
            timeline_array=self.tdata[timeline_index]['TIMELINE']
            fmtstr=str('%iD' % timeline_array.shape[1])
            dimstr=str('%i' % timeline_array.shape[0])
            col1  = pyfits.Column(name='timelines', format=fmtstr, dim=dimstr, unit='ADU', array=timeline_array)
            cols  = pyfits.ColDefs([col1])
            tbhdu = pyfits.BinTableHDU.from_columns(cols)
            for keyword in self.fitsblurbs.keys():
                if keyword in self.tdata[timeline_index].keys():
                    val=self.tdata[timeline_index][keyword]
                    if isinstance(val,dt.datetime):
                        tbhdu.header[keyword]=(val.strftime('%Y-%m-%d %H:%M:%S UTC'),self.fitsblurbs[keyword])
                    else:
                        tbhdu.header[keyword]=(val,self.fitsblurbs[keyword])
                
            hdulist.append(tbhdu)
            
        thdulist = pyfits.HDUList(hdulist)
        thdulist.writeto(fitsfile_fullpath)
        self.printmsg('FITS file written: %s' % fitsfile_fullpath)

    return

def read_fits(self,filename):
    '''
    open a FITS file and determine whether it is QubicStudio or QubicPack
    '''
    if not isinstance(filename,str):
        self.printmsg('ERROR! please enter a valid filename.')
        return None
    
    if not os.path.exists(filename):
        self.printmsg('ERROR! file not found: %s' % filename)
        return None

    try:
        hdulist = pyfits.open(filename)
        nhdu = len(hdulist)
    except:
        self.printmsg('FITS Error! %s' % filename,verbosity=1)
        return False

    # check if it's a QUBIC file
    if nhdu>1\
       and ('TELESCOP' in hdulist[0].header.keys())\
       and (hdulist[0].header['TELESCOP'].strip()=='QUBIC'):
        self.printmsg('Reading QubicPack file: %s' % filename)
        self.dataset_name = os.path.basename(filename).replace('.fits','')
        self.read_qubicpack_fits(hdulist)
        hdulist.close()
        return True

    # check if it's a QubicStudio file
    # QubicStudio FITS files always have at least 2 HDUs, with nothing in the primary header
    nogood_msg = 'Unrecognized FITS file!'
    if 'INSTRUME' not in hdulist[1].header.keys():
        self.printmsg("'INSTRUME' keyword not found\n%s" % nogood_msg)
        hdulist.close()
        return False
    if hdulist[1].header['INSTRUME'].strip() !=  'QUBIC':
        self.printmsg('Instrument is not QUBIC\n%s' % nogood_msg)
        hdulist.close()
        return False
    if 'EXTNAME' not in hdulist[1].header.keys():
        self.printmsg("'EXTNAME' keyword not found\n%s" % nogood_msg)
        hdulist.close()
        return False
    
    self.printmsg('Reading QubicStudio FITS file: %s' % filename,verbosity=2)
    chk = self.read_qubicstudio_fits(hdulist)
    hdulist.close()
    return chk 

def read_fits_field(self,hdu,fieldname):
    '''
    check if a field exists, and if so read the data
    '''
    nfields = hdu.header['TFIELDS']
    fieldnames = []
    for field_idx in range(nfields):
        fieldno = field_idx+1
        fieldnames.append(hdu.header['TTYPE%i' % fieldno].strip())

    if fieldname in fieldnames:
        field_idx = fieldnames.index(fieldname)
        units = hdu.header['TUNIT%i' % (field_idx+1)]
        if units=='ms since 1970-01-01T00:00:00': # convert timestamp to seconds instead of milliseconds
            return 0.001*hdu.data.field(field_idx)
        else:
            return hdu.data.field(field_idx)

    return None


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
        file_date = dt.datetime.strptime(basename,'calsource_%Y%m%dT%H%M%S.fits')
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

def find_hornswitch(self,datadir):
    '''
    try to find hornswitch files corresponding to the observation date
    '''
    self.printmsg('trying to find hornswitch data corresponding to %s' % self.dataset_name,verbosity=2)
    if self.obsdate is None:
        self.printmsg('No observation date')
        return None

    if self.endobs is None:
        self.printmsg('No end observation date')
        return None

    # hornswitch directory is normally two up
    hornswitch_dir = '%s/hornswitch' % os.path.dirname(os.path.dirname(datadir))
    filetype = 'hornswitch'
    datadir = hornswitch_dir
    search_start = self.obsdate
    search_end = self.endobs
    datefmt = '%Y-%m-%d %H:%M:%S'
    self.printmsg('looking for files in time period: %s to %s'
                  % (search_start.strftime(datefmt),search_end.strftime(datefmt)),verbosity=3)
    one_minute = dt.timedelta(minutes=1)
    pattern_time = search_start
    files = []
    while pattern_time<=search_end:
        pattern = '%s/hornswitch_?_%s*.fits' % (hornswitch_dir,pattern_time.strftime('%Y%m%dT%H%M'))
        self.printmsg('looking for files with pattern: %s' % pattern,verbosity=3)
        files += glob(pattern)
        pattern_time += one_minute
        
    if len(files)==0:
        self.printmsg('No %s data found in directory: %s' % (filetype,datadir),verbosity=1)
        return None
    files.sort()
    
    self.printmsg('found %i hornswitch files' % len(files),verbosity=2)
    
    return files

def assign_bath_temperature(self):
    '''
    assign the TES bath temperature from the housekeeping data
    '''
    # the bath temperature probes, in order of preference
    bath_probes = ['TES stage','MGC3_PID_0_Mes','MMR3_CH2_R']
    for probe in bath_probes:
        testemp = None
        hkdat = self.get_hk(probe)
        if hkdat is None: continue

        if probe=='MMR3_CH2_R':
            # transfer function recalibrated 2020-06-18
            # R=Ro*exp(sqrt(Tg/T): Ro=8.44521 Ohm Tg=31.18887 K
            # see elog: https://elog-qubic.in2p3.fr/demo/443
            Ro=8.44521
            Tg=31.18887
            R = hkdat
            try:
                testemp = Tg/(np.log(R/Ro))**2
            except:
                testemp = None
        else:
            testemp = hkdat

        if testemp is not None: break

    # we need the time axis to make the clean Tbath array
    self.printmsg('assign_bath_temperature using probe: %s' % probe,verbosity=2)
    hktstamps = self.timeaxis(datatype=probe)

    # filter spike artefacts in the HK data
    if testemp is not None:
        idxok = (testemp<0.5) & (testemp>0) & (hktstamps<dt.datetime.utcnow().timestamp())
        if idxok.sum()==0:
            idxok = (testemp<300) & (testemp>0)

    if testemp is None or idxok.sum()==0:
        # in desperation, temperature given in the dataset title
        if self.dataset_name is not None:
            match = re.search('[0-9][0-9][0-9]mK',self.dataset_name)
            if match:
                tempstr = match[0].replace('mK','')
                temperature = 1e-3*float(tempstr)
                testemp = np.array([temperature])
                idxok = np.array([True])
                utcoffset = self.obsdate.timestamp() - dt.datetime.utcfromtimestamp(self.obsdate.timestamp()).timestamp()
                hktstamps = np.array([self.obsdate.timestamp()+utcoffset])
                self.temperature = temperature                
                self.printmsg('Assigning TES temperature from the dataset name: %.1fmK' % (1000*self.temperature),verbosity=1)

    if testemp is None or idxok.sum()==0:
        # that's it.  give up.  the temperature is unknown.
        self.printmsg('WARNING!  Bath temperature is unknown!',verbosity=1)
        return
        
    min_temp = testemp[idxok].min()
    max_temp = testemp[idxok].max()
    temperature = testemp[idxok].mean()
        
    self.printmsg('TES temperature varies between %.1fmK and %.1fmK during the measurement' % (1000*min_temp,1000*max_temp))
    self.printmsg('Using TES temperature %.1fmK' % (1000*temperature),verbosity=2)
    self.temperature = temperature

    # assign the bath temperature to the asic objects
    self.assign_temperature(temperature)

    # assign a clean Tbath array with its timestamps
    self.Tbath = (hktstamps[idxok],testemp[idxok])
    
    return

def read_qubicstudio_dataset(self,datadir,asic=None):
    '''
    read a QubicStudio data set which consists of a number of FITS files in a directory
    '''
    if not os.path.isdir(datadir):
        self.printmsg('Please enter a valid directory.')
        return None

    if self.__object_type__=='qubicfp':
        asic = 'ALL'

    if asic is None:
        asic = 1
        self.printmsg('If you would like to read data for a specific ASIC, please include the keyword asic=<N>')

    # if directory name ends with a slash, remove it
    if datadir[-1]==os.sep: datadir = datadir[:-1]
    self.dataset_name = os.path.basename(datadir)
        
    subdir = OrderedDict()
    subdir['science'] = 'Sums'
    subdir['asic'] = 'Hks'
    subdir['hkextern'] = 'Hks'
    subdir['raw'] = 'Raws'
    subdir['MMR'] = 'Hks'
    subdir['hkintern'] = 'Hks'
    subdir['MGC'] = 'Hks'
    subdir['calsource'] = 'Hks'
    subdir['calsource-conf'] = 'Hks'
    
    pattern = OrderedDict()
    if asic=='ALL':
        pattern['science'] = '%s/%s/science-asic*-*.fits' % (datadir,subdir['science'])
        pattern['raw'] = '%s/%s/raw-asic*-*.fits' % (datadir,subdir['raw'])
    else:
        pattern['science'] = '%s/%s/science-asic%i-*.fits' % (datadir,subdir['science'],asic)
        pattern['raw'] = '%s/%s/raw-asic%i-*.fits' % (datadir,subdir['raw'],asic)


    pattern['asic'] = '%s/%s/conf-asics-*.fits' % (datadir,subdir['asic'])
    pattern['hkextern'] = '%s/%s/hk-extern-*.fits' % (datadir,subdir['hkextern'])
    pattern['hkintern'] = '%s/%s/hk-intern-*.fits' % (datadir,subdir['hkintern'])
    pattern['MMR'] = '%s/%s/hk-MMR-*.fits' % (datadir,subdir['MMR'])
    pattern['MGC'] = '%s/%s/hk-MGC-*.fits' % (datadir,subdir['MGC'])
    pattern['calsource-conf'] = '%s/%s/calibConf-*.fits' % (datadir,subdir['calsource'])
    pattern['calsource'] = '%s/%s/calibData-*.fits' % (datadir,subdir['calsource'])

    # check for files, and read if found
    for filetype in pattern.keys():
        files = glob(pattern[filetype])
        if len(files)==0:
            self.printmsg('No %s data found in directory: %s/%s' % (filetype,datadir,subdir[filetype]),verbosity=1)
            continue
        files.sort()

        # we expect only one file of each type (per ASIC) unless we're reading all ASIC
        if asic!='ALL' and len(files)>1:
            self.printmsg('WARNING! There are %i %s data files!' % (len(files),filetype))
            self.printmsg('         There should only be 1 file.')

        for filename in files:
            chk = self.read_fits(filename)

        
    # check if there's an observation date.
    if self.obsdate is None:
    
        # If not, get it from the first available housekeeping
        if self.hk is not None:
            for hktype in self.hk.keys():
                if 'ComputerDate' in self.hk[hktype].keys():
                    self.obsdate = dt.datetime.fromtimestamp(self.hk[hktype]['ComputerDate'][0])
                    self.endobs = dt.datetime.fromtimestamp(self.hk[hktype]['ComputerDate'][-1])
        else:    
            self.printmsg('Error! No HK Data!',verbosity=1)


    if self.obsdate is None:
        # still unsuccessful, so try to get the date from the dataset name
        try:
            self.obsdate = dt.datetime.strptime(self.dataset_name.split('__')[0],'%Y-%m-%d_%H.%M.%S')
            self.endobs = self.obsdate
            self.printmsg('Assigning observation date from dataset name',verbosity=2)
        except:
            self.printmsg('Error! Could not find an observation date.',verbosity=2)
            self.obsdate = dt.datetime.strptime('2017-05-11T09:00:00','%Y-%m-%dT%H:%M:%S')
            self.endobs = self.obsdate

    for asicobj in self.asic_list:
        if asicobj is None: continue
        if asicobj.obsdate is None:
            asicobj.obsdate = self.obsdate
        if asicobj.endobs is None:
            asicobj.endobs = self.endobs
            


    # assign temperature labels
    self.assign_temperature_labels()

    # now try to find the corresponding calsource file, if necessary
    if 'CALSOURCE' not in self.hk.keys():
        self.find_calsource(datadir)

    # and try to find corresponding hornswitch files
    self.hornswitch_files = self.find_hornswitch(datadir)

    # run biasphase() in order to assign timeline_vbias
    biasphase = self.bias_phase()

    # read IV evaluation if it exists, and also assign the dataset name to the qubicasic objects
    # and assign the default gps sample offset for the timestamp algorithm
    if self.__object_type__=='qubicfp':
        for asic_obj in self.asic_list:
            if asic_obj is None: continue
            asic_obj.read_filter()
            asic_obj.dataset_name = self.dataset_name
            asic_obj.assign_default_gps_sample_offset()
    else:
        self.read_filter()

    # assign bath temperature
    self.assign_bath_temperature()
    
    return True

def read_calsource_fits(self,hdu):
    '''
    read the calibration source data from the given HDU of a fits file
    '''
    
    self.hk['CALSOURCE'] = {}
    self.hk['CALSOURCE']['timestamp'] = hdu.data.field(0)
    self.hk['CALSOURCE']['Value'] = hdu.data.field(1)
    
    return

def read_qubicstudio_fits(self,hdulist):
    '''
    read a QubicStudio FITS file
    the argument h is an hdulist after opening the FITS file
    and having confirmed that it really is a QubicStudio FITS file
    '''

    nhdu = len(hdulist)
    hdu = hdulist[1] # the primary header has nothing in it
    keys = hdu.header.keys()

    # what kind of QubicStudio file?
    QS_filetypes = ['ASIC_SUMS',
                    'CONF_ASIC1',
                    'EXTERN_HK',
                    'ASIC_RAW',
                    'INTERN_HK',
                    'MMR_HK',
                    'MGC_HK',
                    'CALSOURCE',
                    'CALSOURCE-CONF']
    extname = hdu.header['EXTNAME'].strip()
    if extname not in QS_filetypes:
        self.printmsg('ERROR! Unrecognized QubicStudio FITS file: %s' % extname)
        return None

    self.datafiletype = extname

    # read the timeline data from each detector
    if extname=='ASIC_SUMS':
        chk = self.read_qubicstudio_science_fits(hdu)
        hdulist.close()
        return chk

    # the EXTNAME of the hdulist[1] is ASIC1, but the hdulist has the config for all the ASICS
    if extname=='CONF_ASIC1': 
        chk = self.read_qubicstudio_asic_fits(hdulist)
        hdulist.close()
        return chk

    if extname=='ASIC_RAW':
        chk = self.read_qubicstudio_raw_fits(hdu)
        hdulist.close()
        return chk
        

    if extname=='EXTERN_HK':
        chk = self.read_qubicstudio_hkextern_fits(hdu)
        hdulist.close()
        return chk

    chk = self.read_qubicstudio_hkfits(hdu)
    hdulist.close()
    return chk

def read_qubicstudio_science_fits(self,hdu):
    '''
    read the science data
    The HDU passed here as the argument should already have been identified as the Science HDU
    '''
    self.printmsg('DEBUG: read_qubicstudio_science_fits object type is %s' % self.__object_type__,verbosity=3)

    extname = hdu.header['EXTNAME'].strip()
    # save PPS/GPS etc as we do for HK files
    if extname not in self.hk.keys():
        self.hk[extname] = {}

    self.printmsg("existing keys for hk['%s']: %s" % (extname,self.hk[extname].keys()),verbosity=3)

    # dictionary for timeline data.
    # tdata is a list.  This is historical bagage from early development 2017/18
    tdata = self.tdata[-1]
    
    # check which ASIC
    asic = hdu.header['ASIC_NUM']
    if self.asic is None:
        self.asic = asic
    elif self.asic != asic:
        msg='ERROR! ASIC Id does not match: previous=%i, current=%i' % (self.asic, asic)
        tdata['WARNING'].append(msg)
        self.printmsg(msg)
        self.asic = asic
    self.printmsg('Reading science data for ASIC %i' % asic,verbosity=2)



    # read the science data
    npts = hdu.header['NAXIS2']
    adu = np.zeros((self.NPIXELS,npts))
    for pix_idx in range(self.NPIXELS):
        pix_no = pix_idx+1
        fieldname = 'pixel%i' % pix_no
        adu[pix_idx,:] = self.read_fits_field(hdu,fieldname)

    if 'TIMELINE' not in tdata.keys():
        tdata['TIMELINE'] = adu
        self.printmsg('storing new timeline data for ASIC %i' % asic,verbosity=2)
    else: # multi file data set
        tstamp_start = 0.001*hdu.data.field(0)[0]
        start_str = dt.datetime.utcfromtimestamp(tstamp_start).strftime('%Y-%m-%d %H:%M:%S')
        tstamp_end = 0.001*hdu.data.field(0)[-1]
        self.printmsg('concatenating detector timeline data to pre-existing timeline: starting at %s' % start_str,verbosity=2)
        tdata['TIMELINE'] = np.concatenate((tdata['TIMELINE'],adu),axis=1)
        
        

    # get number of samples per sum
    #################################################################
    # mail from Michel 20181221:
    ## nsample=100 est le nombre total de points par TES. Cela
    ## permet de remonter à la fréquence d'échantillonnage:
    ## fs=2E6/128/nsample
    ##
    ## NbSamplesPerSum=64 est le nombre de points
    ## utilisé pour obtenir le signal scientifique, après le
    ## masquage de certains des 100 points par TES. En fait, le
    ## signal scientifique est la somme des 64 points non masqué sur
    ## les 100 échantillons pris sur chaque TES.
    #################################################################
    nbsamplespersum_list  =  self.read_fits_field(hdu,'NbSamplesPerSum')
    if 'NSAMSUM_LST' not in tdata.keys():
        tdata['NSAMSUM_LST'] = nbsamplespersum_list
        self.printmsg('storing new NsamplesPerSum data',verbosity=2)
    else:
        tdata['NSAMSUM_LST'] = np.concatenate((tdata['NSAMSUM_LST'],nbsamplespersum_list))
        self.printmsg('concatenating NsamplesPerSum data to pre-existing',verbosity=2)

    
    self.hk[extname]['NbSamplesPerSum'] = tdata['NSAMSUM_LST']
    if len(nbsamplespersum_list)==0:
        NbSamplesPerSum = None
    else:
        NbSamplesPerSum = nbsamplespersum_list[-1]
    tdata['NbSamplesPerSum'] = NbSamplesPerSum
        
    ## check if they're all the same
    difflist = np.unique(nbsamplespersum_list)
    if len(difflist)==0:
        msg = 'WARNING! nsamples per sum is not recorded!'
        self.printmsg(msg)
        tdata['WARNING'].append(msg)
        
    if len(difflist)>1:
        msg = 'WARNING! nsamples per sum changed during the measurement!'
        self.printmsg(msg)
        tdata['WARNING'].append(msg)


    # get the time axis
    computertime_idx = 0
    gpstime_idx = 1
    gpstime = 1e-3*hdu.data.field(gpstime_idx)
    if len(np.unique(gpstime))<2:
        msg="ERROR!  Bad GPS data!"
        tdata['WARNING'].append(msg)

    if 'GPSDate' not in self.hk[extname].keys():
        self.hk[extname]['GPSDate'] = gpstime
        self.printmsg('storing new GPS data',verbosity=2)
    else:
        self.hk[extname]['GPSDate'] = np.concatenate((self.hk[extname]['GPSDate'],gpstime))
        self.printmsg('concatenating GPSDate to pre-existing',verbosity=2)

    ppstime_idx = 2
    pps = hdu.data.field(ppstime_idx)
    if 'PPS' not in self.hk[extname].keys():
        self.hk[extname]['PPS'] = pps
        self.printmsg('storing new PPS data',verbosity=2)
    else:
        self.printmsg('len(PPS) = %i' % len(self.hk[extname]['PPS']),verbosity=3)
        self.printmsg('concatenating PPS to pre-existing',verbosity=2)
        self.hk[extname]['PPS'] = np.concatenate((self.hk[extname]['PPS'],pps))
                                                 
    dateobs = []
    timestamp = 1e-3*hdu.data.field(computertime_idx)
    if 'ComputerDate' not in self.hk[extname].keys():
        self.hk[extname]['ComputerDate'] = timestamp
        self.printmsg('storing new ComputerDate data',verbosity=2)
    else:
        self.hk[extname]['ComputerDate'] = np.concatenate((self.hk[extname]['ComputerDate'],timestamp))
        self.printmsg('concatenating ComputerDate to pre-existing',verbosity=2)
                                                 
        
    for tstamp in timestamp:
        dateobs.append(dt.datetime.utcfromtimestamp(tstamp))

    if 'DATE-OBS' not in tdata.keys():
        tdata['DATE-OBS'] = dateobs
        if len(dateobs)==0:
            tdata['BEG-OBS'] = None
        else:
            tdata['BEG-OBS'] = dateobs[0]
    else:
        tdata['DATE-OBS'] += dateobs
        
    self.obsdate = tdata['BEG-OBS']
    if len(dateobs)==0:
        tdata['END-OBS'] = None
    else:
        tdata['END-OBS'] = dateobs[-1]
    self.endobs = tdata['END-OBS']

    # other info in the science file
    # CN is the "chronological number" counter tag for each sample so we can see if we lost packets
    # Sine phase is the phase of the bias voltage
    keys = ['CN','TES Sinus phase']
    for key in keys:
        if key not in self.hk[extname].keys():
            self.hk[extname][key] = self.read_fits_field(hdu,key)
            self.printmsg('storing new %s' % key,verbosity=2)
        else:
            self.hk[extname][key] = np.concatenate((self.hk[extname][key],self.read_fits_field(hdu,key)))
            self.printmsg('concatenating %s to pre-existing' % key,verbosity=2)

    # try to guess the name of the detector array (P87, or whatever)
    self.guess_detector_name()
        
    return True

def read_qubicstudio_asic_fits(self,hdulist):
    '''
    read the data giving the ASIC configuration
    The HDU passed here as the argument should already have been identified as the ASIC HDU

    we should read the science data first, and then read the corresponding ASIC table

    '''
    self.printmsg('DEBUG: call to read_qubicstudio_asic.  asic=%s' % self.asic, verbosity=3)
    
    if self.asic is None:
        self.printmsg('ERROR! Please read the science data first (asic is unknown)')
        return None
    
    # dictionary for timeline data.
    # tdata is a list.  This is historical bagage from early development 2017/18
    tdata = self.tdata[-1]

    # check which ASIC
    hdu = hdulist[self.asic]
    asic = hdu.header['ASIC_NUM']
    if self.asic != asic:
        msg="ERROR! I'm reading the wrong ASIC table!  want %i, found %i" % (self.asic, asic)
        tdata['WARNING'].append(msg)
        self.printmsg(msg)
            
    # get the time axis
    computertime_idx = 0
    gpstime_idx = 1
    dateobs = []
    timestamp = 1e-3*hdu.data.field(computertime_idx)
    gpstime = 1e-3*hdu.data.field(gpstime_idx)
    if len(np.unique(gpstime))<2:
        msg="ERROR!  Bad GPS data!"
        tdata['WARNING'].append(msg)
    npts = len(timestamp)
    for tstamp in timestamp:
        dateobs.append(dt.datetime.utcfromtimestamp(tstamp))
    if npts==0:
        tdata['ASICDATE'] = None
        tdata['BEGASIC%i' % asic] = None
        tdata['ENDASIC%i' % asic] = None
        msg = 'ASIC%i: There are no housekeeping measurements!' % asic
        self.printmsg(msg)
        tdata['WARNING'].append(msg)
        return self.read_qubicstudio_hkfits(hdu)

    tdata['ASICDATE'] = dateobs
    tdata['BEGASIC%i' % asic] = dateobs[0]
    tdata['ENDASIC%i' % asic] = dateobs[-1]

    # print some info
    datefmt = '%Y-%m-%d %H:%M:%S'
    self.printmsg('ASIC%i: There are %i housekeeping measurements in the period %s to %s'\
                  % (asic,npts,dateobs[0].strftime(datefmt),dateobs[-1].strftime(datefmt)),verbosity=2)
    
    # get the Raw Mask
    rawmask_lst = self.read_fits_field(hdu,'Raw-mask')
    tdata['RAW-MASK'] = rawmask_lst[0]
    self.rawmask = tdata['RAW-MASK']
    for idx in range(rawmask_lst.shape[0]):
        if not np.array_equal(self.rawmask,rawmask_lst[idx]):
            msg = 'WARNING! Raw-mask varies during the measurement!'
            self.printmsg(msg)
            tdata['WARNING'].append(msg)
            break
            

    # get bias level (this is given directly in Volts.  No need to translate from DAC values)
    # TESAmplitude is in fact, peak-to-peak, so multiply by 0.5
    amplitude = 0.5*self.read_fits_field(hdu,'TESAmplitude')
    offset = self.read_fits_field(hdu,'TESOffset')

    ## check if they're all the same
    bias_min = offset-amplitude
    difflist = np.unique(bias_min)
    if len(difflist)!=1:
        msg = 'WARNING! Minimum Bias changed during the measurement!'
        self.printmsg(msg)
        tdata['WARNING'].append(msg)
    tdata['BIAS_MIN'] = np.nanmin(bias_min)
    self.min_bias = tdata['BIAS_MIN']
    
    bias_max = offset+amplitude
    difflist = np.unique(bias_max)
    if len(difflist)!=1:
        msg = 'WARNING! Maximum Bias changed during the measurement!'
        self.printmsg(msg)
        tdata['WARNING'].append(msg)
    tdata['BIAS_MAX'] = np.nanmax(bias_max)
    self.max_bias = tdata['BIAS_MAX']

    # get the number of samples
    nsamples_list = self.read_fits_field(hdu,'nsample')
    tdata['NSAM_LST'] = nsamples_list
    tdata['NSAMPLES'] = nsamples_list[-1]
    self.nsamples = tdata['NSAMPLES']
    difflist = np.unique(nsamples_list)
    if len(difflist)!=1:
        msg = 'WARNING! nsample changed during the measurement!'
        tdata['WARNING'].append(msg)

    # Relay feedback resistance (bit0: heater on/off and bit1: 10kOhm/100kOhm)
    onoff_list = []
    rfb_list = []
    relay_list = self.read_fits_field(hdu,'Relays state')
    for val in relay_list:
        onoff_list.append(val & 1)
        rfb_bit = val >> 1 & 1
        if rfb_bit==1:
            rfb_list.append(100e3)
        else:
            rfb_list.append(10e3)
            
    tdata['RFB_LST'] = rfb_list
    tdata['R_FEEDBK'] = rfb_list[-1]
    self.Rfeedback = tdata['R_FEEDBK']
    tdata['R_HEATER_LST'] =  onoff_list
    tdata['R_HEATER'] = onoff_list[-1]
    difflist1 = np.unique(rfb_list)
    difflist2 = np.unique(onoff_list)
    if len(difflist1)!=1 or len(difflist2)!=1:
        msg = 'WARNING! Feedback Relay parameters changed during the measurement!'
        tdata['WARNING'].append(msg)

    # check for synchro error
    sync_dat = self.read_fits_field(hdu,'NETQUIC synchro error')
    sync_err = (sync_dat==1).sum()
    if sync_err>0:
        msg = 'WARNING! Lost synch during the measurement!'
        tdata['WARNING'].append(msg)

    # read all the stuff in the asic file as an HK file
    return self.read_qubicstudio_hkfits(hdu)

def read_qubicstudio_hkfits(self,hdu):
    '''
    read a QubicStudio housekeeping FITS file
    '''
    hkname = hdu.header['EXTNAME'].strip()
    self.printmsg('Reading HK fits: %s' % hkname,verbosity=3)
    self.hk[hkname] = {}
    nfields = hdu.header['TFIELDS']
    for idx in range(nfields):
        fieldno = idx + 1
        fieldname = hdu.header['TTYPE%i' % fieldno]
        if fieldname=='Platform-PPS': fieldname = 'PPS'
        self.hk[hkname][fieldname] = hdu.data.field(idx)

    # convert QubicStudio timestamps in ms to s
    for key in ['ComputerDate','GPSDate']:
        if key in self.hk[hkname].keys():
            self.hk[hkname][key] = 1e-3*self.hk[hkname][key]


    # correct misinterpreted RaspberryDate
    key = 'RaspberryDate'
    if key in self.hk[hkname].keys():
        if self.hk[hkname][key][0] > 1e10:
            self.printmsg('correcting misinterpreted RaspberryDate',verbosity=2)
            # QubicStudio misinterpreted the socket broadcast
            raspdate_raw = np.int64(self.hk[hkname][key])
            npts = len(raspdate_raw)
            fmt = '<%iq' % npts
            raspdate_packed = struct.pack(fmt,*raspdate_raw)
            fmt = '<%id' % npts
            raspdate = struct.unpack(fmt,raspdate_packed)
            self.hk[hkname][key] = np.array(raspdate)
    return

def read_qubicstudio_hkextern_fits(self,hdu):
    '''
    read the housekeeping data
    the main item of interest here is the TES bath temperature
    which is read in assign_bath_temperature()
    this wrapper is kept for historical reasons
    '''
    return self.read_qubicstudio_hkfits(hdu)

def read_qubicpack_fits(self,h):
    '''
    read a FITS file that was written by QubicPack
    argument h is an hdulist after opening a FITS file
    and confirming that it really is a QubicPack fits file
    '''

    self.datafiletype='QP_FITS'
    self.observer=h[0].header['OBSERVER']
    self.assign_obsdate(dt.datetime.strptime(h[0].header['DATE-OBS'],'%Y-%m-%d %H:%M:%S UTC'))
            
    self.nsamples=h[0].header['NSAMPLES']
    if self.nsamples=='': self.nsamples=100 # data from 12/13 July 2017
    self.tinteg=h[0].header['INT-TIME']
    self.NPIXELS=h[0].header['NPIXELS']
    self.asic=h[0].header['ASIC']
    self.QubicStudio_ip=h[0].header['QUBIC-IP']

    if 'NCYCLES' in h[0].header.keys():
        self.debugmsg('reading ncycles')
        self.nbiascycles=h[0].header['NCYCLES']

    if 'CYCBIAS' in h[0].header.keys():
        self.debugmsg('reading cycle_vbias')
        self.debugmsg(str(h[0].header.cards[13]))
        self.cycle_vbias=h[0].header['CYCBIAS']
        self.debugmsg('cycle_vbias is %s' % str(self.cycle_vbias))

    if 'TES_TEMP' in h[0].header.keys():
        self.temperature=h[0].header['TES_TEMP']
    else:
        self.temperature=None

    if 'END-OBS' in h[0].header.keys():
        self.endobs=dt.datetime.strptime(h[0].header['END-OBS'],'%Y-%m-%d %H:%M:%S UTC')
    else:
        self.endobs=None

    if 'BIAS_MOD' in h[0].header.keys():
        self.bias_frequency=h[0].header['BIAS_MOD']
    if 'BIASMODE' in h[0].header.keys():
        self.bias_mode=h[0].header['BIASMODE']
    if 'BIAS_MIN' in h[0].header.keys():
        self.min_bias=h[0].header['BIAS_MIN']
    if 'BIAS_MAX' in h[0].header.keys():
        self.max_bias=h[0].header['BIAS_MAX']
    if 'BIAS_FAC' in h[0].header.keys():
        self.bias_factor=h[0].header['BIAS_FAC']

    if 'DET_NAME' in h[0].header.keys():
        self.detector_name=h[0].header['DET_NAME']
    # in case detector name is undefined...
    self.guess_detector_name()

    if 'R_FEEDBK' in h[0].header.keys():
        self.Rfeedback=h[0].header['R_FEEDBK']

    if 'FLL_STAT' in h[0].header.keys():
        self.FLL_state=h[0].header['FLL_STAT']
    if 'FLL_P' in h[0].header.keys():
        self.FLL_P=h[0].header['FLL_P']
    if 'FLL_I' in h[0].header.keys():
        self.FLL_I=h[0].header['FLL_I']
    if 'FLL_D' in h[0].header.keys():
        self.FLL_D=h[0].header['FLL_D']
                
    if 'NPIXSAMP' in h[0].header.keys():
        self.NPIXELS_sampled=h[0].header['NPIXSAMP']
    
    self.debugmsg('Finished reading the primary header.')
    
    self.tdata=[]
    for hdu in h[1:]:
        hdrtype=hdu.header['TTYPE1']

        if hdrtype=='RawMask':
            '''
            this is the mask of the samples (filtered out samples)
            '''
            self.debugmsg('Reading RawMask HDU')

            # number of values should be 125
            nvals=hdu.header['NAXIS2']
            self.debugmsg('RawMask: nvals=%i' % nvals)
            if nvals!=125:
                self.printmsg('WARNING! RawMask has the wrong number of mask values: %i' % nvals)

            # read the raw mask
            self.rawmask=np.zeros(nvals,dtype=int)
            data=hdu.data
            for idx in range(nvals):
                self.rawmask[idx]=data[idx][0]
            self.debugmsg('Finished reading RawMask HDU')

                        
        if hdrtype=='V_tes':
            '''
            this is I-V data
            '''
        
            # number of bias points
            nbias=eval(hdu.header['TFORM1'].strip()[:-1])

            # read the adu matrix
            data=hdu.data
            self.adu=np.empty((self.NPIXELS,nbias))
            for n in range(self.NPIXELS):
                self.adu[n,:]=data[n][0]

        if hdrtype=='V_bias':
            '''
            this is the bias points
            '''

            data=hdu.data
            self.vbias=np.empty(nbias)
            for n in range(nbias):
                self.vbias[n]=data[n][0]
            self.max_bias=max(self.vbias)
            self.min_bias=min(self.vbias)
            self.max_bias_position=np.argmax(self.vbias)
            

        if hdrtype=='timelines':
            '''
            this is the timeline data
            '''
            self.printmsg('reading timeline data')
            tdata={}
            data=hdu.data
            npts=eval(hdu.header['TFORM1'].strip()[:-1])
            timeline=np.empty((self.NPIXELS,npts))
            for n in range(self.NPIXELS):
                timeline[n,:]=data[n][0]
            tdata['TIMELINE']=timeline
            for keyword in self.fitsblurbs.keys():
                if keyword in hdu.header.keys():
                    if keyword=='DATE-OBS' or keyword=='END-OBS':
                        val=dt.datetime.strptime(hdu.header[keyword],'%Y-%m-%d %H:%M:%S UTC')
                    else:
                        val=hdu.header[keyword]
                    if val=='':val=None
                    tdata[keyword]=val
            if 'DATE-OBS' not in tdata.keys():
                tdata['DATE-OBS'] = self.obsdate

            self.tdata.append(tdata)

    h.close()


    if self.exist_iv_data():
        f=self.read_filter()
        #if f is None:f=self.filter_iv_all()

    return

        
def get_from_keyboard(self,msg,default=None):
    ''''
    get interactive input from the keyboard
    '''
    pythonmajor = sys.version_info[0]
    prompt='%s (default: %s) ' % (msg,str(default))
    if pythonmajor==2:
        ans = raw_input(prompt)
    else:
        ans = input(prompt)
    if ans=='':return default
    if type(default)==str:
        return ans
    
    try:
        x=eval(ans)
    except:
        self.printmsg('invalid response.')
        return None
    return x
    

def writelog(self,msg):
    '''
    write some output with a timestamp to a log file
    and also write it on the screen
    '''

    if self.logfile is None: self.assign_logfile()
    
    handle=open(self.logfile,'a')
    timestamp=dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC -- ')
    handle.write(timestamp+msg+'\n')
    handle.close()
    self.printmsg(timestamp+msg)

    # send it to the QubicStudio logbook, if possible
    client=self.connect_QubicStudio(silent=True)
    if not client is None:client.sendAddToLogbook('pyStudio',msg)

    return
         


########################################################################
###### convenient wrappers for returning data ##########################
########################################################################

def RawMask(self,asic=None):
    '''
    return the Raw Mask from the Housekeeping data
    '''
    if self.__object_type__!='qubicfp':
        asic = self.asic

    if asic is None:
        self.printmsg('Please enter a valid ASIC.')
        return None
    
    hktype = 'CONF_ASIC%i' % asic

    if self.__object_type__=='qubicfp':
        HK = self.asic_list[asic-1].hk
    else:
        HK = self.hk
        
    if hktype not in HK.keys():
        self.printmsg('No ASIC housekeeping data!')
        return None

    keyname = 'Raw-mask'
    if keyname not in HK[hktype].keys():
        self.printmsg('No Raw Mask data!')
        return None
    
    rawmask = HK[hktype][keyname]
    return rawmask

def get_hk(self,data=None,hk=None,asic=None):
    '''
    return the data for a requested housekeeping
    '''
    self.printmsg('DEBUG: get_hk() called with data=%s hk=%s asic=%s' % (data,hk,asic),verbosity=4)
    
    data = self.qubicstudio_hk_truename(data)
    if data is None:
        self.printmsg('Please enter a valid housekeeping type')
        return None

    if hk is None:
        hk = self.qubicstudio_filetype_truename(data,asic=asic)
        
    if hk is None: # try to find which Housekeeping
        for key in self.hk.keys():
            if data in self.hk[key].keys():
                hk = key
                break
    
    hk = self.qubicstudio_filetype_truename(hk)
    if hk is None:
        self.printmsg('Please enter a valid hk type, for example "platform" or "sci"')
        return None
    
    if hk not in self.hk.keys():
        if hk=='ASIC_SUMS' and self.__object_type__=='qubicfp':
            if asic is None:
                self.printmsg('Please enter a valid ASIC number')
                return None
            asic_idx = asic - 1
            HK = self.asic_list[asic_idx].hk[hk]
        else:
            self.printmsg('Not a valid housekeeping ID: %s' % hk)
            return None
    else:
        HK = self.hk[hk]

    if data not in HK.keys():
        self.printmsg('No %s in %s' % (data,hk))
        return None

    val = HK[data]
    if data.upper().find('ERROR')<0 and val.max() == 0:
        self.printmsg('get_hk() : Bad %s Data!' % data,verbosity=3)
    return val

def pps(self,hk=None,asic=None):
    '''
    return PPS data for a given HK
    '''
    return self.get_hk('PPS',hk,asic)

def gps(self,hk=None,asic=None):
    '''
    return GPS data for a given HK
    '''
    gps = self.get_hk('GPSDate',hk,asic)
    if gps is None: return None

    ### do not correct for localtime.  this is done elsewhere with utcfromtimestamp
    # t0 = float(self.obsdate.strftime('%s.%f'))
    # delta = np.abs(gps[0] - t0)
    # if delta > 1799: return gps-delta
    return gps

def azimuth_redmount(self):
    '''
    return the Azimuth data timeline from the Red Mount (a.k.a. Calibration Mount)
    '''
    hktype = 'INTERN_HK'
    if hktype not in self.hk.keys():
        self.printmsg('No platform data!')
        return None

    azkey = 'Platform-Azimut'
    if azkey not in self.hk[hktype].keys():
        self.printmsg('No Azimuth data!')
        return None

    azRaw = self.hk[hktype][azkey]
    az = (azRaw.astype(np.int) - 2**15) * 360.0/2**16
    return az

def azimuth(self):
    '''
    return the Azimuth data
    '''
    if self.obsdate < obsmount_implemented:
        return self.azimuth_redmount()

    hktype = 'EXTERN_HK'
    if hktype not in self.hk.keys():
        self.printmsg('No platform data!')
        return None

    azkey = 'Pressure_6' # a bit of trickery here while waiting for QS to be upgraded
    if azkey not in self.hk[hktype].keys():
        self.printmsg('No Azimuth data!')
        return None

    az = self.hk[hktype][azkey]
    return az

def elevation_redmount(self):
    '''
    return the Elevation data timeline from the Red Mount (a.k.a. Calibration Mount)
    '''
    hktype = 'INTERN_HK'
    if hktype not in self.hk.keys():
        self.printmsg('No platform data!')
        return None

    elkey = 'Platform-Elevation'
    if elkey not in self.hk[hktype].keys():
        self.printmsg('No Elevation data!')
        return None

    elRaw = self.hk[hktype][elkey]
    # offset is deduced from beam synthesis mapping on 2019-04-06
    offset = 10131.591
    el = (elRaw.astype(np.int) - offset) * 360.0/2**16
    return el

def elevation(self):
    '''
    return the Elevation data
    '''

    if self.obsdate < obsmount_implemented:
        return self.elevation_redmount()
    
    hktype = 'EXTERN_HK'
    if hktype not in self.hk.keys():
        self.printmsg('No platform data!')
        return None

    azkey = 'Pressure_7' # a bit of trickery here while waiting for QS to be upgraded
    if azkey not in self.hk[hktype].keys():
        self.printmsg('No Elevation data!')
        return None

    el = self.hk[hktype][azkey]
    return el
    

def hwp_position(self):
    '''
    return the HWP position timeline
    '''
    hktype = 'INTERN_HK'
    if hktype not in self.hk.keys():
        self.printmsg('No internal housekeeping data!')
        return None

    key = 'HWP-Position'
    if key not in self.hk[hktype].keys():
        self.printmsg('No HWP data!')
        return None

    hwp_pos = self.hk[hktype][key]
    return hwp_pos

def bias_phase(self):
    '''
    return the bias voltage phase when in Sine mode
    '''
    self.printmsg('DEBUG: call to bias_phase()',verbosity=4)

    hktype = 'ASIC_SUMS'
    
    if hktype not in self.hk.keys():
        self.printmsg('No scientific housekeeping data!')
        return None

    sinekey = 'TES Sinus phase'
    if sinekey not in self.hk[hktype].keys():
        self.printmsg('No bias sine data!')
        return None

    # convert uint to +/- float, and normalize to amplitude +/- 1
    sineraw = self.hk[hktype][sinekey]
    if len(sineraw)==0:
        self.printmsg('bias sine data is empty!')
        return None
        
    if max(sineraw)==0: return None
    sinephase = np.array(sineraw,dtype=np.float64)
    idx_neg = np.where(sineraw > 32767)
    sinephase[idx_neg] -= 65536.0
    norm = max(abs(sinephase.max()),abs(sinephase.min()))
    sinephase = sinephase/norm

    # assign the vbias for the timeline (vbias for the I-V curve will be a subset)
    if self.max_bias is not None and self.min_bias is not None:
        amplitude = 0.5*(self.max_bias - self.min_bias)
        offset = self.min_bias + amplitude
        self.printmsg('DEBUG: bias_phase() : scaling vbias to max/min bias',verbosity=4)
        self.timeline_vbias = amplitude*sinephase + offset
    
    return sinephase

def FLL_State(self):
    '''
    return the FLL state with its timestamps.  
    NOTE:  there is a variable called FLL_state which is still here for historical reasons
           It was used in the acquisition method when we ran the Oxford dilution fridge.
    '''
    hktype = 'CONF_ASIC%i' % self.asic
    fllkey = 'FLL_State'
    if hktype not in self.hk.keys():
        self.printmsg('No ASIC configuration data!')
        return None
    fll_state = self.hk[hktype][fllkey]
    fll_timestamps = self.hk[hktype]['GPSDate']
    return fll_timestamps,fll_state

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

def infotext(self,TES=None):
    '''
    information to put on plots as a subtitle
    '''
    
    ttl = self.obsdate.strftime('%Y-%m-%d %H:%M:%S')
    if self.detector_name!='undefined':
        ttl += ' Array %s' % self.detector_name

    if self.__object_type__!='qubicfp':
        ttl += ' ASIC#%i' % self.asic

    if TES is not None:
        ttl += ' TES#%03i' % TES

    ttl += ' T$_\mathrm{bath}$='
    if self.temperature is None or self.temperature<0:
        ttl += 'unknown'
    else:
        ttl += '%.1fmK' % (1000*self.temperature)

    return ttl

        
def qubicstudio_filetype_truename(self,ftype,asic=None):
    '''
    return the valid key name given a nickname for the QubicStudio
    filetype within the dataset
    '''

    self.printmsg('DEBUG: calling filetype_truename with ftype=%s' % ftype,verbosity=4)
    if ftype is None: return None
    if ftype.upper() == 'PLATFORM': return 'INTERN_HK'
    if ftype.upper() == 'HK': return 'INTERN_HK'
    if ftype.upper().find('AZ')==0: return 'INTERN_HK'
    if ftype.upper().find('EL')==0: return 'INTERN_HK'
    if ftype.upper().find('HWP')==0: return 'INTERN_HK'
    if ftype.upper() == 'EXTERN': return 'EXTERN_HK'
    if ftype.upper() == 'TEMPERATURE': return 'EXTERN_HK'
    if ftype.upper() == 'CALSOURCE': return 'CALSOURCE'
    if ftype.upper().find('SCI')==0: return 'ASIC_SUMS'
    if ftype.upper() == 'TES': return 'ASIC_SUMS'
    if ftype.upper().find('DET')==0: return 'ASIC_SUMS'
    if ftype.upper().find('MMR')==0: return 'MMR_HK'
    if ftype.upper().find('MGC')==0: return 'MGC_HK'
    if ftype.upper().find('TBATH')==0: return 'EXTERN_HK'
    if ftype.upper().find('SWITCH')==0: return 'INTERN_HK'
    if ftype.upper().find('RAW')>=0: return 'ASIC_RAW'

    # if we give a particular HK data name, return the associated filetype
    hktruename = self.qubicstudio_hk_truename(ftype)
    for key in self.hk.keys():
        for datkey in self.hk[key]:
            if hktruename == datkey: return key

    # maybe we gave a temperature label
    if self.temperature_labels is None: return ftype.upper()        
    for key in self.temperature_labels.keys():
        if ftype.upper() == self.temperature_labels[key].upper(): return 'EXTERN_HK'


    # if we are in the qubicasic object, and unsuccessful up to here, just return the ftype
    if self.__object_type__=='qubicasic':
        return ftype.upper()
        
    # is it an ASIC specific request?
    if asic is None:
        if hktruename in qubicasic_hk_keys:
            print('Please enter an ASIC number.')
        return ftype.upper()

    asic_idx = asic-1
    if self.asic_list[asic_idx] is not None:
        return self.asic(asic).qubicstudio_filetype_truename(ftype)

    self.printmsg('DEBUG: failed filetype_truename.  Returning ftype=%s' % ftype.upper(),verbosity=4)
    return ftype.upper()

def qubicstudio_hk_truename(self,hktype):
    '''
    return the valid key name for a given housekeeping nickname
    '''
    self.printmsg('DEBUG: calling hk_truename with hktype=%s' % hktype,verbosity=4)
    if hktype is None: return None
    
    hktype_upper = hktype.upper()
    if hktype_upper == 'SWITCH1': return 'RFSwitch 1 closed'
    if hktype_upper == 'SWITCH2': return 'RFSwitch 2 closed'
    if hktype_upper == 'AZ': return 'Platform-Azimut'
    if hktype_upper == 'EL': return 'Platform-Elevation'
    if hktype_upper.find('AV')==0: return hktype_upper
    if hktype_upper == 'CALSOURCE': return 'CALSOURCE'
    if hktype_upper == 'SOURCE': return 'CALSOURCE'
    if hktype_upper == 'TBATH': return 'AVS47_1_CH2'
    if hktype_upper == 'HWP': return 'HWP-Position'
    if hktype_upper.find('SYNC')>=0: return 'NETQUIC synchro error'

    hktype_spaces = hktype_upper.replace(' ','_')
    if hktype_spaces.find('FLL_P')>=0: return 'FLL_P'
    if hktype_spaces.find('FLL_I')>=0: return 'FLL_I'
    if hktype_spaces.find('FLL_D')>=0: return 'FLL_D'

    if self.temperature_labels is None: return hktype

    for key in self.temperature_labels.keys():
        if hktype_upper == self.temperature_labels[key].upper(): return key
    
    return hktype

def azel_etc(self,TES=None):
    '''
    return Az/El and data in a dictionary for the convenience of JCH
    '''

    msg =    '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    msg += '\nDEPRECATION WARNING! Do you really need to use this method?'
    msg += "\nThe elements in this object (called 'a' in this example) can be retrieved with the following methods:"
    if TES is None:
        msg += "\n'data <n>' : a.timeline_array(asic=<n>)"
    else:
        msg += "\n'data <n>' : a.timeline(asic=<n>,TES=%i" % TES
    msg += "\n't_data <n>' : a.timeaxis(datatype='science',asic=<n>)"
    msg += "\n'az' : a.azimuth()"
    msg += "\n'el' : a.elevation()"
    msg += "\n't_azel' : a.timeaxis(datatype='hk')"
    msg += "\n'hwp_position' : a.hwp_position()"
    msg += "\n'data_src' : a.calsource()[1]"
    msg += "\n't_src' : a.calsource()[0]"
    msg += '\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print(msg)
    
    retval = {}

    az = self.azimuth()
    retval['az'] = az

    el = self.elevation()
    retval['el'] = el

    t_azel = self.timeaxis(datatype='hk',axistype='pps')
    retval['t_azel'] = t_azel
    retval['t_hk'] = t_azel
    retval['hwp_position'] = self.hwp_position()

    calsource = self.calsource()
    t_src = calsource[0]
    data_src = calsource[1]
    retval['t_src'] = t_src
    retval['data_src'] = data_src

    if self.__object_type__!='qubicfp':
        if TES is None:
            data = self.timeline_array()
        else:
            data = self.timeline(TES=TES)
        self.printmsg('DEBUG: calling timeline_timeaxis from azel_etc()',verbosity=4)
        t_data = self.timeline_timeaxis(axistype='pps')
        retval['t_data'] = t_data
        retval['data'] = data
        print(msg) # print message again to be annoying about it
        return retval


    for idx,asic_obj in enumerate(self.asic_list):
        asic_no = idx + 1
        if TES is None:
            data = self.timeline_array(asic=asic_no)
        else:
            data = self.timeline(TES=TES,asic=asic_no)
        self.printmsg('DEBUG: calling timeaxis from azel_etc()',verbosity=4)
        t_data = self.timeaxis(datatype='sci',asic=asic_no)
        retval['t_data %i' % asic_no] = t_data
        retval['data %i' % asic_no] = data

    print(msg) # print message again to be annoying about it
    return retval
    
