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
from __future__ import division, print_function
import numpy as np
import sys,os,time,subprocess,struct
import datetime as dt
from glob import glob
import pickle
from collections import OrderedDict
from astropy.io import fits as pyfits

from satorchipy.datefunctions import tot_seconds

def debugmsg(self,msg,verbosity=2):
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

    hdulist=pyfits.open(filename)
    nhdu=len(hdulist)

    # check if it's a QUBIC file
    if nhdu>1\
       and ('TELESCOP' in hdulist[0].header.keys())\
       and (hdulist[0].header['TELESCOP'].strip()=='QUBIC'):
        self.printmsg('Reading QubicPack file: %s' % filename)
        return self.read_qubicpack_fits(hdulist)

    # check if it's a QubicStudio file
    # QubicStudio FITS files always have at least 2 HDUs, with nothing in the primary header
    nogood_msg = 'Unrecognized FITS file!'
    if 'INSTRUME' not in hdulist[1].header.keys():
        self.printmsg("'INSTRUME' keyword not found\n%s" % nogood_msg)
        return False
    if hdulist[1].header['INSTRUME'].strip() !=  'QUBIC':
        self.printmsg('Instrument is not QUBIC\n%s' % nogood_msg)
        return False
    if 'EXTNAME' not in hdulist[1].header.keys():
        self.printmsg("'EXTNAME' keyword not found\n%s" % nogood_msg)
        return False
    
    self.printmsg('Reading QubicStudio FITS file: %s' % filename)
    return self.read_qubicstudio_fits(hdulist)

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
        return hdu.data.field(field_idx)

    return None


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

    # calsource directory is normally two up
    calsource_dir = '%s/calsource' % os.path.dirname(os.path.dirname(datadir))
        
    scidir = '%s/Sums' % datadir
    hkdir = '%s/Hks' % datadir
    subdir = OrderedDict()
    subdir['science'] = 'Sums'
    subdir['asic'] = 'Hks'
    subdir['hkextern'] = 'Hks'
    subdir['raw'] = 'Raws'
    subdir['MMR'] = 'Hks'
    subdir['hkintern'] = 'Hks'
    subdir['MGC'] = 'Hks'
    
    pattern = OrderedDict()
    if asic=='ALL':
        pattern['science'] = '%s/%s/science-asic*-*.fits' % (datadir,subdir['science'])
    else:
        pattern['science'] = '%s/%s/science-asic%i-*.fits' % (datadir,subdir['science'],asic)
    
    pattern['asic'] = '%s/%s/conf-asics-*.fits' % (datadir,subdir['asic'])
    pattern['hkextern'] = '%s/%s/hk-extern-*.fits' % (datadir,subdir['hkextern'])
    pattern['hkintern'] = '%s/%s/hk-intern-*.fits' % (datadir,subdir['hkintern'])
    pattern['MMR'] = '%s/%s/hk-MMR-*.fits' % (datadir,subdir['MMR'])
    pattern['MGC'] = '%s/%s/hk-MGC-*.fits' % (datadir,subdir['MGC'])

    # check for files, and read if found
    for filetype in pattern.keys():
        files = glob(pattern[filetype])
        if len(files)==0:
            self.printmsg('No %s data found in directory: %s/%s' % (filetype,datadir,subdir[filetype]))
            continue

        # we expect only one file of each type (per ASIC) unless we're reading all ASIC
        if asic<>'ALL' and len(files)>1:
            self.printmsg('WARNING! There are %i %s data files!' % (len(files),filetype))
            self.printmsg('         There should only be 1 file.')

        for filename in files:
            chk = self.read_fits(filename)


    # assign bath temperature to asic objects.  This is useful for infotext()
    if self.__object_type__=='qubicfp':
        for asic_obj in self.asic_list:
            self.printmsg('assigning bath temperature of %.3fK to asic %i' % (self.temperature,asic_obj.asic),verbosity=2)
            asic_obj.tdata[-1]['TES_TEMP'] = self.temperature
            asic_obj.temperature = self.temperature
            

    # now try to find the corresponding calsource file
    # look for files within the last hour, and then take the closest one to the start time
    # the files are in FITS format as of Wed 10 Apr 2019 10:21:35 CEST
    self.printmsg('trying to find calsource data')
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
        self.printmsg('No %s data found in directory: %s' % (filetype,datadir))
        return

    # find the file which starts before and nearest to obsdate
    filename = None
    file_delta = 1e6
    for f in files:
        basename = os.path.basename(f)
        file_date = dt.datetime.strptime(basename,'calsource_%Y%m%dT%H%M%S.fits')
        delta = tot_seconds(self.obsdate - file_date)
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
    if nhdu<>2:
        self.printmsg("This doesn't look like a calsource file!")
        hdulist.close()
        return
    hdu = hdulist[1]
    if 'EXTNAME' not in hdu.header.keys()\
       and hdu.header['EXTNAME']<>'CALSOURCE':
        self.printmsg("This is not a calsource FITS file!")
        hdulist.close()
        return

    self.read_calsource_fits(hdu)
    hdulist.close()

    return 

def read_calsource_fits(self,hdu):
    '''
    read the calibration source data from the given HDU of a fits file
    '''
    
    self.hk['CALSOURCE'] = {}
    self.hk['CALSOURCE']['ComputerDate'] = hdu.data.field(0)
    self.hk['CALSOURCE']['amplitude'] = hdu.data.field(1)
    
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
    QS_filetypes = ['ASIC_SUMS','CONF_ASIC1','EXTERN_HK','ASIC_RAW','INTERN_HK','MMR_HK','MGC_HK','CALSOURCE']
    extname = hdu.header['EXTNAME'].strip()
    if extname not in QS_filetypes:
        self.printmsg('ERROR! Unrecognized QubicStudio FITS file: %s' % extname)
        return None

    self.datafiletype = extname

    # dictionary for the return of timeline data
    if self.tdata is None:self.tdata = [{}]
    tdata = self.tdata[-1]
    if 'WARNING' not in tdata.keys():
        tdata['WARNING'] = []
    
    # get the timeline data from each detector
    if extname=='ASIC_SUMS':
        self.read_qubicstudio_science_fits(hdu)
        hdulist.close()
        return
        
    if extname=='CONF_ASIC1':
        self.read_qubicstudio_asic_fits(hdulist)
        hdulist.close()
        return

    if extname=='EXTERN_HK':
        self.read_qubicstudio_hkextern_fits(hdu)
        hdulist.close()
        return

    if extname=='CALSOURCE':
        self.read_calsource_fits(hdu)
        hdulist.close()
        return

    self.read_qubicstudio_hkfits(hdu)
    hdulist.close()
    return

def read_qubicstudio_science_fits(self,hdu):
    '''
    read the science data
    The HDU passed here as the argument should already have been identified as the Science HDU
    '''
    self.printmsg('DEBUG: read_qubicstudio_science_fits object type is %s' % self.__object_type__,verbosity=3)
    if self.tdata is None:self.tdata = [{}]
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
    self.printmsg('Reading data for ASIC %i' % asic)

    # save PPS/GPS etc as we do for HK files
    extname = hdu.header['EXTNAME'].strip()
    self.hk[extname] = {}

    # read the science data
    npts = hdu.header['NAXIS2']
    adu = np.zeros((self.NPIXELS,npts))
    for pix_idx in range(self.NPIXELS):
        pix_no = pix_idx+1
        fieldname = 'pixel%i' % pix_no
        adu[pix_idx,:] = self.read_fits_field(hdu,fieldname)
    tdata['TIMELINE'] = adu

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
    tdata['NSAMSUM_LST'] = nbsamplespersum_list
    self.hk[extname]['NbSamplesPerSum'] = nbsamplespersum_list
    NbSamplesPerSum = nbsamplespersum_list[-1]
    tdata['NbSamplesPerSum'] = NbSamplesPerSum
    ## check if they're all the same
    difflist = np.unique(nbsamplespersum_list)
    if len(difflist)!=1:
        msg = 'WARNING! nsamples per sum changed during the measurement!'
        self.printmsg(msg)
        tdata['WARNING'].append(msg)


    # get the time axis
    computertime_idx = 0
    gpstime_idx = 1
    self.hk[extname]['GPSDate'] = 1e-3*hdu.data.field(gpstime_idx)

    ppstime_idx = 2
    self.hk[extname]['PPS'] = hdu.data.field(ppstime_idx)

    dateobs = []
    timestamp = 1e-3*hdu.data.field(computertime_idx)
    self.hk[extname]['ComputerDate'] = timestamp

    for tstamp in timestamp:
        dateobs.append(dt.datetime.utcfromtimestamp(tstamp))
    tdata['DATE-OBS'] = dateobs
    tdata['BEG-OBS'] = dateobs[0]
    self.obsdate = tdata['BEG-OBS']
    tdata['END-OBS'] = dateobs[-1]

    # other info in the science file
    # CN is the "chronological number" counter tag for each sample so we can see if we lost packets
    # Sine phase is the phase of the bias voltage
    keys = ['CN','TES Sinus phase']
    for key in keys:
        self.hk[extname][key] = self.read_fits_field(hdu,key) 

    # try to guess the name of the detector array (P87, or whatever)
    self.guess_detector_name()
        
    return

def read_qubicstudio_asic_fits(self,hdulist):
    '''
    read the data giving the ASIC configuration
    The HDU passed here as the argument should already have been identified as the ASIC HDU

    we should read the science data first, and then read the corresponding ASIC table
    '''
    if self.asic is None:
        self.printmsg('ERROR! Please read the science data first (asic is unknown)')
        return None
    
    if self.tdata is None:self.tdata = [{}]
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
    npts = len(timestamp)
    for tstamp in timestamp:
        dateobs.append(dt.datetime.utcfromtimestamp(tstamp))
    tdata['ASICDATE'] = dateobs
    tdata['BEGASIC%i' % asic] = dateobs[0]
    tdata['ENDASIC%i' % asic] = dateobs[-1]

    # print some info
    datefmt = '%Y-%m-%d %H:%M:%S'
    self.printmsg('There are %i housekeeping measurements in the period %s to %s'\
          % (npts,dateobs[0].strftime(datefmt),dateobs[-1].strftime(datefmt)))
    
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
    tdata['BIAS_MIN'] = min(bias_min)
    self.min_bias = tdata['BIAS_MIN']
    
    bias_max = offset+amplitude
    difflist = np.unique(bias_max)
    if len(difflist)!=1:
        msg = 'WARNING! Maximum Bias changed during the measurement!'
        self.printmsg(msg)
        tdata['WARNING'].append(msg)
    tdata['BIAS_MAX'] = max(bias_max)
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

    # read all the stuff in the asic file as an HK file
    return self.read_qubicstudio_hkfits(hdu)

def read_qubicstudio_hkfits(self,hdu):
    '''
    read a QubicStudio housekeeping FITS file
    '''
    hkname = hdu.header['EXTNAME'].strip()
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
            
    return

def read_qubicstudio_hkextern_fits(self,hdu):
    '''
    read the housekeeping data
    '''
    if self.tdata is None:self.tdata = [{}]
    tdata = self.tdata[-1]
    testemp = hdu.data.field(5)
    min_temp = testemp.min()
    max_temp = testemp.max()
    temperature = testemp.mean()
    self.printmsg('TES temperatures varies between %.1fmK and %.1fmK during the measurement' % (1000*min_temp,1000*max_temp))
    self.printmsg('Using TES temperature %.1fmK' % (1000*temperature))
    self.tdata[0]['TES_TEMP'] = temperature
    self.temperature = temperature

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

            self.tdata.append(tdata)

    h.close()


    if self.exist_iv_data():
        f=self.read_filter()
        if f is None:f=self.filter_iv_all()

    return


def read_bins(self,filename):
    if not isinstance(filename,str):
        self.printmsg('ERROR! please enter a valid filename.')
        return None
            
    if not os.path.exists(filename):
        self.printmsg('ERROR! file not found: %s' % filename)
        return None
            
    self.printmsg('reading binary file: %s, I suppose this is a timeline' % filename)

    data=[]
    with open(filename, "rb") as f:
        b = f.read(14)
        data.append(struct.unpack('128i', f.read(4*128)))
        while f.read(14) != "": data.append(struct.unpack('128i', f.read(4*128)))
        
    data=np.asarray(zip(*data))
    self.NPIXELS=128
    npts=int(data.size/128.)

    timeline_array=[]
    timeline_tes=np.empty((npts))
    for n in range(self.NPIXELS):
        timeline_tes=data[n]        
        timeline_array.append(timeline_tes)

    timeline_array=np.array(timeline_array)
    self.tdata=[]
    tdata={}
    tdata['TIMELINE']=timeline_array
    self.tdata.append(tdata)
    f.close()
    return
        
def get_from_keyboard(self,msg,default=None):
    ''''
    get interactive input from the keyboard
    '''
    prompt='%s (default: %s) ' % (msg,str(default))
    ans=raw_input(prompt)
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
         

def pps2date(self,pps,gps):
    '''
    convert the gps date to a precise date given the pps
    '''
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
            if sep <> 0: # next PPS valid only if we have a non-zero step (modif by MP)
                pps_indexes.append(idx)    
                separations.append(sep)
                prev = gps[idx]

    separations = np.array(separations[1:])

    self.printmsg('mean pps interval is %.4f second' % separations.mean())
    self.printmsg('max pps interval is  %.4f second' % separations.max())
    self.printmsg('min pps interval is  %.4f second' % separations.min())
            
    # find the GPS date corresponding to the PPS
    tstamp = -np.ones(npts)
    prev_gps = gps[0]
    #offset = 50 # delay after PPS for valid GPS (this should be different for scientific and housekeeping) 
    for idx in pps_indexes:
        gps_at_pps = gps[idx]

        # the pps arrives just before the corresponding gps
        # so we simply add 1 second to current gps value (gps increments in steps of 1 second exactly)
        # (modification by MP & JCH)
        tstamp[idx] = gps_at_pps + 1

        ''' this is  replaced by the above (MP & JCH)
        # we use the GPS timestamp from a bit later
        offset_idx = idx + offset
        if offset_idx>=npts:offset_idx=npts-1
        next_gps = gps[offset_idx]
        tstamp[idx] = next_gps
        '''
        
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

########################################################################
###### convenient wrappers for returning data ##########################
########################################################################

def RawMask(self,asic=None):
    '''
    return the Raw Mask from the Housekeeping data
    '''
    if self.__object_type__<>'qubicfp':
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

def get_hk(self,data='PPS',hk=None,asic=None):
    '''
    return the data for a requested housekeeping
    '''
    
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
    if val.max() == 0:
        self.printmsg('Bad %s Data!' % data)
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
    return self.get_hk('GPSDate',hk,asic)

def azimuth(self):
    '''
    return the Azimuth data timeline
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

def elevation(self):
    '''
    return the Elevation data timeline
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
    el = (elRaw.astype(np.int) - 2**15) * 360.0/2**16
    return el

def bias_phase(self):
    '''
    return the bias voltage phase when in Sine mode
    '''
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
    sinephase = np.array(sineraw,dtype=np.float)
    idx_neg = np.where(sineraw > 32767)
    sinephase[idx_neg] -= 65536.0
    norm = max(abs(sinephase.max()),abs(sinephase.min()))
    sinephase = sinephase/norm

    # assign the vbias
    if self.max_bias is not None and self.min_bias is not None:
        amplitude = 0.5*(self.max_bias - self.min_bias)
        offset = self.min_bias + amplitude
        self.vbias = amplitude*sinephase + offset
    
    return sinephase

def calsource(self):
    '''
    return the calibration source data
    '''
    if 'CALSOURCE' in self.hk.keys():
        t_src = self.hk['CALSOURCE']['ComputerDate']
        data_src = self.hk['CALSOURCE']['amplitude']
        return t_src,data_src

    return None

def infotext(self,TES=None):
    '''
    information to put on plots as a subtitle
    '''
    ttl = self.obsdate.strftime('%Y-%m-%d %H:%M:%S')
    if self.detector_name!='undefined':
        ttl += ' Array %s' % self.detector_name

    ttl += ' ASIC#%i' % self.asic

    if TES is not None:
        ttl += ' TES#%03i' % TES

    ttl += ' T$_\mathrm{bath}$='
    if self.temperature is None or self.temperature<0:
        ttl += 'unknown'
    else:
        ttl += '%.1fmK' % (1000*self.temperature)

    return ttl

        
def qubicstudio_filetype_truename(self,ftype):
    '''
    return the valid key name given a nickname for the QubicStudio
    filetype within the dataset
    '''
    if ftype is None: return None
    if ftype.upper() == 'PLATFORM': return 'INTERN_HK'
    if ftype.upper() == 'HK': return 'INTERN_HK'
    if ftype.upper().find('AZ')==0: return 'INTERN_HK'
    if ftype.upper().find('EL')==0: return 'INTERN_HK'
    if ftype.upper() == 'ASIC': return 'CONF_ASIC1'
    if ftype.upper() == 'EXTERN': return 'EXTERN_HK'
    if ftype.upper().find('SCI')==0: return 'ASIC_SUMS'
    return ftype.upper()


def azel_etc(self,TES=None):
    '''
    return Az/El and data in a dictionary for the convenience of JCH
    '''
    if TES:
        data = self.timeline(TES=TES)
    else:
        data = self.timeline_array()

    t_data = self.timeline_timeaxis(axistype='pps')
    az = self.azimuth()
    el = self.elevation()
    t_azel = self.timeaxis(datatype='hk',axistype='pps')

    calsource = self.calsource()
    if calsource is None:
        t_src = -1
        data_src = -1
    else:
        t_src = calsource[0]
        data_src = calsource[1]

    retval = {}
    retval['t_data'] = t_data
    retval['data'] = data
    retval['t_azel'] = t_azel
    retval['az'] = az
    retval['el'] = el
    retval['t_src'] = t_src
    retval['data_src'] = data_src
    return retval
    
