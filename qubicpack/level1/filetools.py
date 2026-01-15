'''
$Id: filetools.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Mon 23 Jan 2023 14:37:03 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

utilities to write and read QUBIC Level-1 data files
'''
import os
import numpy as np
import h5py

from satorchipy.datefunctions import utcnow, tstamp2dt
from .flags import nflags, flag_definition

# explanatory notes for the primary header
hdr_comment = {}
hdr_comment['TELESCOP'] ='Telescope used for the observation'
hdr_comment['FILETYPE'] = 'File identification'
hdr_comment['DATASET']  = 'QubicStudio dataset name'
hdr_comment['UNITS']    = 'physical units:  ADU, Watt, or Amp'
hdr_comment['FILENAME'] = 'name of this file'
hdr_comment['FILEDATE'] = 'UT date this file was created'
hdr_comment['OBSDATE']  = 'UT date beginning this dataset'
hdr_comment['ENDOBS']   = 'UT date end of this dataset'
hdr_comment['INDXTYPE'] = 'index ordering: QUBICSOFT or TES'
hdr_comment['INFOINDX'] = 'link to more details about index ordering'
hdr_comment['INFO_FMT'] = 'link to more details about QUBIC Level1 data format'
hdr_comment['INFODSET'] = 'link to more details about this dataset'
hdr_keys = hdr_comment.keys()

# explanatory notes for the housekeeping header
hk_comment = {}
hk_comment['COMMENT'] = "all housekeeping data has been interpolated to the same time axis"
hk_comment['TBATH'] = "The bolometer bath temperature should be around 320mK"
hk_comment['1K Stage'] = "The 1K stage should be around 1.1K"
hk_comment['HWP'] = "Half Wave Plate position"
hk_comment['AZIMUTH'] = "telescope azimuth pointing in degrees"
hk_comment['ELEVATION'] = "telescope elevation pointing in degrees"
hk_comment['ROTATION'] = "telescope bore-sight rotation angle in degrees"
hk_keys = hk_comment.keys()


datefmt = '%Y-%m-%dT%H:%M:%S.%fZ'

def write_level1_header(self,handle,
                        filename='UNNAMED.hdf5',
                        indextype='QUBICSOFT',
                        units='Watt',
                        infolink='https://qubic.in2p3.fr/wiki/TD/DatasetDetails'):
    '''
    write the QUBIC Level-1 primary header to file
    '''
    
    hdr = {}
    for key in hdr_keys:
        hdr[key] = None
    hdr['TELESCOP'] = 'QUBIC'
    hdr['FILETYPE'] = 'QUBIC LEVEL1'
    hdr['FILEDATE'] = utcnow().strftime(datefmt)
    hdr['DATASET']  = self.dataset_name
    hdr['OBSDATE']  = self.obsdate.strftime(datefmt)
    hdr['ENDOBS']   = self.endobs.strftime(datefmt)
    hdr['UNITS']    = units
    hdr['FILENAME'] = filename
    if indextype.upper().find('Q')>=0:
        hdr['INDXTYPE'] = 'QUBICSOFT'
    else:
        hdr['INDXTYPE'] = 'TES'
    hdr['INFOINDX'] = 'https://qubic.in2p3.fr/wiki/TD/TEStranslation'
    hdr['INFO_FMT'] = 'https://qubic.in2p3.fr/wiki/TD/ProcessedDataFormat'
    hdr['INFODSET'] = infolink

    # create a group for supplementary information
    infogrp = handle.create_group('EXPLANATORY NOTES')
    infogrp.attrs['COMMENT'] = 'definitions of header keywords'
    
    # Primary header
    for key in hdr.keys():
        handle.attrs[key] = hdr[key]
        if key in hdr_comment.keys():
            infogrp.attrs[key] = hdr_comment[key]

    # add a comment about flag definitions
    handle.attrs['FLAG DEFINITIONS'] = 'Flag definitions can be found in the attributes of the flag data: /SCIENCE DATA/FLAGS'

    # add a comment about calibration source and carbon fibre info
    handle.attrs['CALIBRATION INFO'] = 'Information about calibration source and carbon fibre settings are in the attributes of CALIBRATION INFO group'
    
    return

def write_level1_housekeeping(self,handle):
    '''
    write some explanatory information about the housekeeping data
    interpolate all housekeeping data to the same time axis
    and save the housekeeping data
    '''

    # make a group for the housekeeping data
    hkgrp = handle.create_group('HOUSEKEEPING')
    for key in hk_comment.keys():
        hkgrp.attrs[key] = hk_comment[key]


    # make a list of all the time axes, and interpolate all to the one with the largest number of samples
    hk = {}
    nsamples_list = []
    tstamp_start_list = []
    tstamp_end_list = []

    # bath temperature
    timeaxis = self.Tbath[0]
    if timeaxis is not None:
        hk['TBATH'] = (timeaxis,self.Tbath[1])
        nsamples_list.append(len(timeaxis))
        tstamp_start_list.append(timeaxis[0])
        tstamp_end_list.append(timeaxis[-1])
        
    for key in hk_keys:
        if key=='COMMENT': continue
        if key=='TBATH': continue # special case because there are multiple sensors. already done above.
        val = self.get_hk(key)
        if val is not None:
            timeaxis = self.timeaxis(datatype=key)
            hk[key] = (timeaxis,val)
            nsamples_list.append(len(timeaxis))
            tstamp_start_list.append(timeaxis[0])
            tstamp_end_list.append(timeaxis[-1])            

    tstamp_start = min(tstamp_start_list)
    tstamp_end = max(tstamp_end_list)
    nsamples = max(nsamples_list)

    # interpolate housekeeping to an evenly sampled time axis
    # the number of samples is derived from the maximum sampling of all housekeeping
    t_interp = tstamp_start + (tstamp_end-tstamp_start)*np.arange(nsamples)/(nsamples-1)
    hkgrp.create_dataset('TIMESTAMP',data=t_interp, compression='gzip', compression_opts=9)
    for key in hk.keys():
        val_interp = np.interp(t_interp, hk[key][0], hk[key][1])
        self.printmsg('write_level1: hk val_interp shape=%s' % (str(val_interp.shape)),verbosity=3)
        hkgrp.create_dataset(key,data=val_interp, compression='gzip', compression_opts=9)

    # add a comment if data is missing
    for key in hk_keys:
        if key=='COMMENT': continue
        if key not in hk.keys():
            error_key = 'ERROR! %s'
            hkgrp.attrs[error_key] = 'ERROR! No data for %s' % key
            
    return hkgrp

def write_level1_calinfo(self,handle):
    '''
    write the calibration source and carbon fibre information
    '''

    calgrp = handle.create_group('CALIBRATION INFO')

    # check for calibration information
    calinfo = self.calsource_info()
    if calinfo is None:
        calgrp.attrs['NO CALIBRATION INFORMATION'] = 'There is no calibration source nor carbon fibre information in this dataset'
        return calgrp

    # an explanatory preamble
    comment_list = ['Calibration source power monitor is found in the SCIENCE DATA.',
                    'The source configuration parameters are found in the attributes of the subsystem groups.',
                    'The DATE is the time that the calibration information was read',
                    'which is usually a bit before the data acquisition start time.']
                    
    calgrp.attrs['COMMENT'] = '\n'.join(comment_list)

    if 'date' in calinfo.keys():
        calgrp.attrs['DATE'] = calinfo['date'].strftime(datefmt)
        
    for key in calinfo.keys():
        if key=='date': continue # already done above
        if key=='lamp': continue # legacy stuff
        if len(calinfo[key])==0: continue # ignore empty 
        if key=='cf':
            grpkey = 'CARBON FIBRE'
        else:
            grpkey = key.upper()
        subgrp = calgrp.create_group(grpkey)
        subinfo = calinfo[key]
        for subkey in subinfo.keys():
            subgrp.attrs[subkey.upper()] = subinfo[subkey]
    
    return calgrp

def write_level1_horninfo(self,handle):
    '''
    write the horn shutter info
    '''

    # !!! TO BE IMPLEMENTED !!!
    
    return None
    

def write_level1(self,indextype='QUBICSOFT',units='Watt',infolink=None,savepath=None):
    '''
    write a HDF5 file with Level-1 data
    '''
    if infolink is None:
        infolink = 'https://qubic.in2p3.fr/wiki/TD/DatasetDetails'
    
    if indextype.upper().find('QUBICSOFT')>=0 or indextype.upper()=='QS':
        is_QSindex = True
        indextype = 'QUBICSOFT'
    else:
        is_QSindex = False
        indextype = 'TES'
    
    # initialize
    filename_suffix = '_level1.hdf5'
    filename = self.dataset_name+filename_suffix
    basename = os.path.basename(filename)
    if savepath is None: # save in current directory
        savepath = '.'

    filename_fullpath = os.sep.join([savepath,filename])
    try:
        handle = h5py.File(filename_fullpath,'w')
    except:
        print('ERROR! Could not create file: %s' % filename_fullpath)
        return None

    # write the primary header
    self.write_level1_header(handle,basename,indextype,units,infolink)

    # create a group for the science data
    scigrp = handle.create_group('SCIENCE DATA')
    
    # scigrp contains the level-1 array and the flagarray
    # with indextype=QS, the dark pixels are in a separate array
    # otherwise, they are in the TOD array in the TES order
    tod_tuple = self.tod(indextype=indextype,units=units)
    t_tod = tod_tuple[0]
    todarray = tod_tuple[1]
    # tod_tuple might be only length 2 for backwards compatibility
    if len(tod_tuple)==2:
        dark_pixels = None
    else:
        dark_pixels = tod_tuple[2]


    scigrp.create_dataset('TIMESTAMP', data=t_tod, compression='gzip', compression_opts=9)
    scigrp['TIMESTAMP'].attrs['COMMENT'] = 'timestamps are the number of seconds since 1970-01-01T00:00:00Z'
    scigrp.create_dataset('TOD', data=todarray, compression='gzip', compression_opts=9)
    scigrp['TOD'].attrs['COMMENT'] = 'index order is %s' % indextype
    if dark_pixels is not None:
        scigrp.create_dataset('DARK PIXELS', data=dark_pixels, compression='gzip', compression_opts=9)
    
    # assign the flags for the TOD
    flag_array = self.assign_flags(indextype,t_tod)

    # write flags
    scigrp.create_dataset('FLAGS',data=flag_array, compression='gzip', compression_opts=9)

    # add flag definitions in attributes, in reverse order
    for flagbit in range(nflags-1,-1,-1):
        flagdef = flag_definition[flagbit]
        if flagdef=='available': continue
        flag_key = 'FLAG_%03i' % flagbit
        scigrp['FLAGS'].attrs[flag_key] = flagdef

    # write the calibration source power monitor data, interpolated to the science data time axis
    t_calsrc,calsrc = self.calsource()
    if calsrc is not None:
        calsrc_interp = np.interp(t_tod,t_calsrc,calsrc)
        scigrp.create_dataset('CALSOURCE',data=calsrc_interp, compression='gzip', compression_opts=9)
        scigrp['CALSOURCE'].attrs['COMMENT'] = 'calsource data has been interpolated to the science data time sampling'
        scigrp['CALSOURCE'].attrs['NSAMPLES'] = len(t_calsrc)
        scigrp['CALSOURCE'].attrs['START DATE'] = tstamp2dt(t_calsrc[0]).strftime(datefmt)
        scigrp['CALSOURCE'].attrs['END DATE'] = tstamp2dt(t_calsrc[-1]).strftime(datefmt)
    
    # write the housekeeping data
    hkgrp = self.write_level1_housekeeping(handle)

    # make a group for the calibration source and carbon fibre info
    calgrp = self.write_level1_calinfo(handle)
    
    # horn status
    horngrp = self.write_level1_horninfo(handle)
    
    
    # write the file
    handle.close()
    print('Level-1 data written to file: %s' % filename_fullpath)
    return


def read_level1(filename):
    '''
    read a QUBIC Level-1 HDF5 file
    '''

    if not os.path.exists(filename):
        print('ERROR! File does not exist: %s' % filename)
        return None
              
    if not os.path.isfile(filename):
        print('ERROR! Not a file: %s' % filename)
        return None

    handle = h5py.open(filename,'r')
    prihdr = handle.attrs
    
    if 'TELESCOP' not in prihdr.keys() or prihdr['TELESCOP']!='QUBIC':
        print('ERROR! Not a QUBIC file: %s' % filename)
        handle.close()
        return None

    if 'FILETYPE' not in prihdr.keys() or prihdr['FILETYPE']!='QUBIC LEVEL1':
        print('ERROR! Not a QUBIC Level-1 file: %s' % filename)
        handle.close()
        return None
    
    return handle
