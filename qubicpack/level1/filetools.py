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

from satorchipy.datefunctions import utcnow
from .flags import nflags, flag_definition

hdr_comment = {}
hdr_comment['TELESCOP'] ='Telescope used for the observation'
hdr_comment['FILETYPE'] = 'File identification'
hdr_comment['DATASET']  = 'QubicStudio dataset name'
hdr_comment['UNITS']    = 'physical units:  ADU, Watt, or Amp'
hdr_comment['FILENAME'] = 'name of this file'
hdr_comment['FILEDATE'] = 'UT date this file was created'
hdr_comment['INDXTYPE'] = 'index ordering: QUBICSOFT or TES'
hdr_comment['INFOINDX'] = 'link to more details about index ordering'
hdr_comment['INFO_FMT'] = 'link to more details about QUBIC Level1 data format'
hdr_comment['INFODSET'] = 'link to more details about this dataset'

hdr_keys = hdr_comment.keys()

def write_level1_header(self,handle,
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
    hdr['FILEDATE'] = utcnow().strftime('%Y-%m-%dT%H:%M:%S')
    hdr['DATASET'] = self.dataset_name
    hdr['UNITS'] = units
    hdr['FILENAME'] = filename+filename_suffix
    if indextype.upper().find('Q')>=0:
        hdr['INDXTYPE'] = 'QUBICSOFT'
    else:
        hdr['INDXTYPE'] = 'TES'
    hdr['INFOINDX'] = 'https://qubic.in2p3.fr/wiki/TD/TEStranslation'
    hdr['INFO_FMT'] = 'https://qubic.in2p3.fr/wiki/TD/ProcessedDataFormat'
    hdr['INFODSET'] = infolink
    
    # Primary header
    for key in hdr.keys():
        handle.attrs[key] = hdr[key]
        if key in hdr_comment.keys():
            comment_key = '%s COMMENT' % key
            handle.attrs[comment_key] = hdr_comment[key]

    # add flag definitions in header, in reverse order
    for flagbit in range(nflags-1,-1,-1):
        flagdef = flag_definition[flagbit]
        if flagdef=='available': continue
        flag_key = 'FLAG_%03i' % flagbit
        handle.attrs[flag_key] = flagdef
    
    return

def write_level1(self,indextype='QUBICSOFT',units='Watt',infolink='https://qubic.in2p3.fr/wiki/TD/DatasetDetails',savepath=None):
    '''
    write a HDF5 file with Level-1 data
    '''

    if indextype.upper().find('QUBICSOFT')>=0 or indextype.upper()=='QS':
        is_QSindex = True
        indextype = 'QUBICSOFT'
    else:
        is_QSindex = False
        indextype = 'TES'
    
    # initialize
    filename_suffix = '_level1.hdf5'
    filename = self.dataset_name
    if savepath is None: # save in current directory
        savepath = '.'

    filename_fullpath = os.sep.join([savepath,filename+filename_suffix])
    try:
        handle = h5py.open(filename_fullpath,'w')
    except:
        print('ERROR! Could not create file: %s' % filename_fullpath)
        return None

    # write the primary header
    self.write_level1_header(handle,indextype,units,infolink)

    # create a group for the science data
    scigrp = h.create_group('SCIENCE DATA')
    
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


    scigrp.create_dataset('TIMESTAMP',data=t_tod)
    scigrp.create_dataset('TOD',data=todarray)
    if dark_pixels is not None:
        scigrp.create_dataset('DARK PIXELS',data=dark_pixels)
    
    # assign the flags for the TOD
    flag_array = self.assign_flags(indextype,t_tod)

    # write flags
    scigrp.create_dataset('FLAGS',data=flagarray)

    
    # make another group with the housekeeping data
    # 300mK stage Temperature
    # 1K stage Temperature
    # HWP position
    # azimuth
    # elevation
    # bore sight rotation
    # calsource info and value interpolated to data time axis
    # carbon fibre parameters
    # horn status
    
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
