#!/usr/bin/env python3
'''
$Id: collect_elevation_temperature_data.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Thu 11 Jul 16:40:33 CEST 2019
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.


gather the temperature data from a series of QubicStudio datasets given as an argument
'''
import os,sys
import numpy as np
import pickle
from qubicpack.qubicfp import qubicfp
qubicfp.verbosity = 1

if len(sys.argv)==1:
    print('usage: %s <list of data directories>' % sys.argv[0])
    quit()
          
datadirs = sys.argv[1:]

t_keys = ['RaspberryDate',
          'AVS47_1_CH0',
          'AVS47_1_CH1',
          'AVS47_1_CH2',
          'AVS47_1_CH3',
          'AVS47_1_CH4',
          'AVS47_1_CH5',
          'AVS47_1_CH6',
          'AVS47_1_CH7',
          'AVS47_2_CH0',
          'AVS47_2_CH1',
          'AVS47_2_CH2',
          'AVS47_2_CH3',
          'AVS47_2_CH4',
          'AVS47_2_CH5',
          'AVS47_2_CH6',
          'AVS47_2_CH7',
          'Temp_1',
          'Temp_2',
          'Temp_3',
          'Temp_4',
          'Temp_5',
          'Temp_6',
          'Temp_7',
          'Temp_8',
          'Temp_9',
          'Temp_10',
          'Temp_11',
          'Temp_12',
          'Temp_13',
          'Temp_14',
          'Temp_15',
          'Temp_16',
          'Temp_17',
          'Temp_18',
          'Temp_19',
          'Temp_20']



name_str = ','.join(t_keys)
fmt_list = []
for idx,key in enumerate(t_keys):
    fmt_list.append('float64')
fmt_str = ','.join(fmt_list)

temperature_rec = np.recarray(names=name_str,formats=fmt_str,shape=(0))
elevation_rec = np.recarray(names='timestamp,elevation',formats='int64,int64',shape=(0))



hkname = 'EXTERN_HK'
key = t_keys[0]
a = qubicfp()
for d in datadirs:
    try:
        chk = a.read_qubicstudio_dataset(d)
    except:
        print('****** ERROR! Could not read dataset: %s' % d)
        continue
    
    if chk is None:
        print('****** bad directory: %s' % d)
        continue

    if a.temperature == -1:
        print('****** no temperature data: %s' % d)
        continue

    el = a.elevation()
    el_timestamp = a.timeaxis(axistype='pps',datatype='platform')
    if el_timestamp is None:
        el_timestamp = a.timeaxis(axistype='computertime',datatype='platform')
    if el_timestamp is None:
        print('****** no timestamp for elevation: %s' % d)
        continue
    
    el_npts = el.shape[0]
    el_rec = np.recarray(names='timestamp,elevation',formats='int64,int64',shape=(el_npts))
    el_rec.timestamp = el_timestamp
    el_rec.elevation = el
    elevation_rec = np.concatenate([elevation_rec,el_rec])

    
    t_npts = a.hk[hkname][key].shape[0]    
    t_rec = np.recarray(names=name_str,formats=fmt_str,shape=(t_npts))
    for idx,key in enumerate(t_keys):
        cmd = 't_rec.%s = a.hk[hkname][key]' % key
        exec(cmd)
    temperature_rec = np.concatenate([temperature_rec,t_rec])

        

whichdata = a.obsdate.strftime('%Y%m%d')
temperature_file = 'temperatures_%s.dat' % whichdata
h = open(temperature_file,'wb')
pickle.dump(temperature_rec,h)
h.close()

elevation_file = 'elevation_%s.dat' % whichdata
h = open(elevation_file,'wb')
pickle.dump(elevation_rec,h)
h.close()

