'''
$Id: pointing.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Mon 22 Dec 2025 11:12:57 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

utilities for reading the binary data saved from the observation mount PLC
these utilities are also used by class obsmount in package qubichw
'''
import os
import numpy as np
from qubicpack.utilities import fmt_translation

axis_fullname = {'AZ': 'azimuth', 'EL': 'elevation', 'RO': 'boresight rotation', 'TR': 'Little Train'}
axis_names = list(axis_fullname.keys())

header_keys = ['TIMESTAMP',
               'IS_ETHERCAT',
               'IS_SYNC',
               'IS_MAINT',
               'AXES_ASYNC_COUNT']
n_header_keys = len(header_keys)
rec_header_names = ','.join(header_keys)
rec_header_format_list = ['float64','uint8','uint8','uint8','int16']
rec_header_format = ','.join(rec_header_format_list)

data_keys = ['AXIS',
             'ACT_VEL_RES',
             'ACT_VEL_ENC',
             'ACT_POS_RES',
             'ACT_POS_ENC',
             'ACT_TORQUE',
             'IS_ENABLED',
             'IS_HOMING',
             'IS_HOMINGSKIP',
             'IS_OPERATIVE',
             'IS_MOVING',
             'IS_OUTOFRANGE',
             'FAULT']
position_key = {'AZ':'ACT_POS_RES', 'EL':'ACT_POS_ENC', 'RO':'ACT_POS_ENC', 'TR':'ACT_POS_ENC'}
rec_data_names = ','.join(data_keys[1:])
rec_data_format_list = ['float64','float64','float64','float64','float64','uint8','uint8','uint8','uint8','uint8','uint8','uint8']
rec_data_format = ','.join(rec_data_format_list)

n_data_keys = len(data_keys)

delimiter = ':'
STX = bytearray([0xaa,0xaa])

def interpret_pointing_chunk(dat):
    '''
    interpret a bytes data chunk which is returned from the observation mount PLC
    '''
    packet = {}
    packet['ok'] = False
    packet['error'] = 'NONE'
    
    dat_str = dat.decode()
    if len(dat_str)==0:
        packet['error'] = 'POINTING: empty data chunk'
        return packet

    dat_list = dat_str.split('\n')
    if len(dat_list)<5:
        packet['error'] = 'partial data: %s' % dat_str
        return packet

    axis = None
    for line in dat_list:
        if len(line)==0: continue
        col = line.split(delimiter)
        ncols = len(col)

        # data for each axis
        if col[0] in axis_names:
            axis = col[0]
            axis_data = {}
            for subidx,val_str in enumerate(col[1:]):
                idx = subidx + 1
                if idx<n_data_keys:
                    key = data_keys[idx]
                else:
                    key = 'UNKNOWN%02i' % idx
                    
                try:
                    val = eval(val_str)
                except:
                    val = val_str
                    
                axis_data[key] = val
            packet[axis] = axis_data
            continue

        # header data
        for idx,val_str in enumerate(col):
            key = header_keys[idx]
            try:
                packet[key] = eval(val_str)
            except:
                packet[key] = val_str

    packet['TIMESTAMP'] *= 0.001
    packet['ok'] = True
    return packet

def read_pointing_bindat(filename):
    '''
    read the binary data saved by the observation mount PLC
    '''
    dat = {}
    dat['ok'] = False
    dat['header'] = None
    dat['data'] = None
    if not os.path.isfile(filename):
        print('ERROR!  File not found: %s' % filename)
        return dat

    h = open(filename,'rb')
    dat_bytes = h.read()
    h.close()

    chunk_list = dat_bytes.split(STX)
    npts = len(chunk_list) - 1

    headerdat = np.recarray(names=rec_header_names,formats=rec_header_format,shape=(npts))
    axdat = {}
    for axname in axis_names:
        axdat[axname] = np.recarray(names=rec_data_names,formats=rec_data_format,shape=(npts))

    idx = 0
    for chunk in chunk_list:
        # print('[%04i]' % idx)
        packet = interpret_pointing_chunk(chunk)
        if not packet['ok']:
            print('packet not okay: %s' % packet['error'])
            continue

        for header in header_keys:
            cmd = 'headerdat[idx].%s = packet[header]' % header
            exec(cmd)

        for axname in axis_names:
            for datname in data_keys[1:]:
                if datname not in packet[axname].keys():
                    # print('ERROR! Invalid key for packet: %s' % datname)
                    continue
                cmd = 'axdat[axname][idx].%s = packet[axname][datname]' % datname
                exec(cmd)
        idx += 1
                
    dat['header'] = headerdat[0:idx]
    final_axdat = {}
    for axname in axis_names:
        final_axdat[axname] = axdat[axname][0:idx]
        
    dat['data'] = final_axdat
    dat['ok'] = True
    return dat
