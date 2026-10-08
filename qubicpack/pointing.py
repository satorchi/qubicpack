'''
$Id: pointing.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Mon 22 Dec 2025 11:12:57 CET
$auth: Mattia Bianchessi <mbianchess@apc.in2p3.fr>
$date: Wed 07 Oct 2026 14:48:42 CEST starting from commit 

$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

utilities for reading the binary data saved from the observation mount PLC
these utilities are also used by class obsmount in package qubichw
'''
import os
import numpy as np

# position offsets measured
# 'EL': 49.315, # see elog: https://elog-qubic.in2p3.fr/demo/1296
# 'EL': 49.935, # see elog: https://elog-qubic.in2p3.fr/demo/1321
# 'AZ':  9.0    # see elog: https://elog-qubic.in2p3.fr/demo/1322 
position_offset = {'AZ': 9.0, 
                   'EL': 49.935,
                   'RO': 0.0,
                   'TR': 0.0
                   }

axis_fullname = {'AZ': 'azimuth', 'EL': 'elevation', 'RO': 'boresight rotation', 'TR': 'Little Train'}
axis_names = list(axis_fullname.keys())

v1_header_keys = ['TIMESTAMP',
                  'IS_ETHERCAT',
                  'IS_SYNC',
                  'IS_MAINT',
                  'AXES_ASYNC_COUNT']

rec_header_names = ','.join(v1_header_keys)
rec_header_format_list = ['float64','uint8','uint8','uint8','int16']
rec_header_format = ','.join(rec_header_format_list)

v2_header_keys = v1_header_keys
v2_rec_header_names = ','.join(['RX_TIMESTAMP']+v1_header_keys)
v2_rec_header_format_list = ['float64','float64','uint8','uint8','uint8','int16']
v2_rec_header_format = ','.join(v2_rec_header_format_list)

## 2026-05-13
v3_header_keys = ['TIMESTAMP1',
                  'TIMESTAMP2',
                  'IS_ETHERCAT',
                  'IS_SYNC',
                  'IS_MAINT',
                  'AXES_ASYNC_COUNT']
v3_rec_header_names = ','.join(['RX_TIMESTAMP']+v3_header_keys)
v3_rec_header_format_list = ['float64','float64','float64','uint8','uint8','uint8','int16']
v3_rec_header_format = ','.join(v3_rec_header_format_list)

### 2026-07-30 17:59:18 CEST
v4_header_keys = ['TIMESTAMP1',
                  'TIMESTAMP2',
                  'IS_ETHERCAT',
                  'IS_SYNC',
                  'IS_MAINT',
                  'AXES_ASYNC_COUNT',
                  'NTP_RESULT',
                  'SET_RTC_RESULT']
v4_rec_header_names = ','.join(['RX_TIMESTAMP']+v4_header_keys)
v4_rec_header_format_list = ['float64','float64','float64','uint8','uint8','uint8','int16','int16','int16']
v4_rec_header_format = ','.join(v4_rec_header_format_list)


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
#position_key = {'AZ':'ACT_POS_RES', 'EL':'ACT_POS_ENC', 'RO':'ACT_POS_ENC', 'TR':'ACT_POS_ENC'}
### 2026-09-04 19:14:00 discussion with Luciano
### 2026-10-08 18:56:34 discussion with Luciano:  use the resolver until we figure out what's wrong with the encoder
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
    packet['version'] = 1
    
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
                # data_keys[0] is the axis name, so values start from index 1
                idx = subidx + 1
                if idx<n_data_keys:
                    key = data_keys[idx]
                else:
                    key = 'UNKNOWN%02i' % idx

                try:
                    val = float(val_str)
                except ValueError:
                    val = val_str
                
                axis_data[key] = val
            packet[axis] = axis_data
            continue

        # header data
        n_headers = len(col)
        if n_headers==len(v4_header_keys):
            packet['version'] = 4
            header_keys = v4_header_keys
        elif n_headers==len(v3_header_keys):
            packet['version'] = 3
            header_keys = v3_header_keys
        else:
            header_keys = v1_header_keys
        for idx,val_str in enumerate(col):
            key = header_keys[idx]
            try:
                packet[key] = float(val_str)
            except ValueError:
                packet[key] = val_str

    # PLC data packet has timestamp in milliseconds
    for key in packet.keys():
        if key.find('TIMESTAMP')==0:
            packet[key] *= 0.001
            
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

    # check if it's the version 1 or version 2 which includes the timestamp of reception
    STXv2 = STX + 'RX'.encode()
    v2_separator = 'XR'.encode()
    chunk_list = dat_bytes.split(STXv2)
    npts = len(chunk_list) - 1
    if npts>1:
        # the format version is read once, from the first packet
        first_chunk = chunk_list[1].split(v2_separator)[-1]
        packet = interpret_pointing_chunk(first_chunk)
        print('PLC data first chunk packet is version: %i' % packet['version'])
        if packet['version']==4:
            pointing_file_version = 4
            rechdr_names = v4_rec_header_names
            rechdr_fmts = v4_rec_header_format
            header_keys = v4_header_keys
        elif packet['version']==3:
            pointing_file_version = 3
            rechdr_names = v3_rec_header_names
            rechdr_fmts = v3_rec_header_format
            header_keys = v3_header_keys
        else:        
            pointing_file_version = 2
            rechdr_names = v2_rec_header_names
            rechdr_fmts = v2_rec_header_format
            header_keys = v1_header_keys
    else:
        chunk_list = dat_bytes.split(STX)
        pointing_file_version = 1        
        npts = len(chunk_list) - 1
        rechdr_names = rec_header_names
        rechdr_fmts = rec_header_format
        header_keys = v1_header_keys

    print('PLC data format version: %i' % pointing_file_version)

    # column names of the header recarray (RX_TIMESTAMP exists only in v2+)
    header_cols = rechdr_names.split(',')
    data_cols = data_keys[1:]

    # one list of values per valid packet; the recarrays are built after the loop
    header_rows = []
    axis_rows = {axname: [] for axname in axis_names}

    for chunk in chunk_list:
        rx_timestamp = None
        if pointing_file_version>1:
            # split the reception timestamp from the PLC payload
            chunks = chunk.split(v2_separator)
            if len(chunks)<2: continue
            chunk = chunks[-1]
            try:
                rx_timestamp = float(chunks[0].decode())
            except ValueError:
                continue
        
        # from here, both v1 and v2 are the same
        packet = interpret_pointing_chunk(chunk)
        if not packet['ok']:
            print('packet not okay: %s' % packet['error'])
            continue

        # build all rows first, so that an incomplete packet is discarded as a whole
        try:
            hrow = [packet[key] for key in header_keys]
            arows = {axname: [packet[axname][key] for key in data_cols] for axname in axis_names}
        except KeyError as err:
            print('packet not okay: missing field %s' % err)
            continue

        header_rows.append(([rx_timestamp] if rx_timestamp is not None else []) + hrow)
        for axname in axis_names:
            axis_rows[axname].append(arows[axname])

    npts_valid = len(header_rows)
    if npts_valid==0:
        print('ERROR!  No valid packets in %s' % filename)
        return dat

    # allocate the recarrays at their final size and fill them column by column;
    hmat = np.array(header_rows,dtype='float64')
    headerdat = np.recarray(names=rechdr_names,formats=rechdr_fmts,shape=(npts_valid,))
    for icol,name in enumerate(header_cols):
        headerdat[name] = hmat[:,icol]

    final_axdat = {}
    for axname in axis_names:
        amat = np.array(axis_rows[axname],dtype='float64')
        final_axdat[axname] = np.recarray(names=rec_data_names,formats=rec_data_format,shape=(npts_valid,))
        for icol,name in enumerate(data_cols):
            final_axdat[axname][name] = amat[:,icol]

    dat['header'] = headerdat
    dat['data'] = final_axdat
    dat['ok'] = True
    return dat


def compare_pointing_dat(dat_old, dat_new):
    '''
    check that two outputs of read_pointing_bindat* are identical in dtype, shape and values
    '''
    pairs = [('header', dat_old['header'], dat_new['header'])]
    pairs += [(ax, dat_old['data'][ax], dat_new['data'][ax]) for ax in axis_names]
    for label,a,b in pairs:
        assert a.dtype==b.dtype, '%s: dtype mismatch' % label
        assert a.shape==b.shape, '%s: shape mismatch %s vs %s' % (label,a.shape,b.shape)
        for name in a.dtype.names:
            assert np.array_equal(a[name],b[name],equal_nan=True), '%s.%s: values differ' % (label,name)
    print('identical')
