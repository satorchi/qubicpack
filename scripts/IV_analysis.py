#!/usr/bin/env python3
'''
$Id: IV_analysis.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Wed 06 Apr 2022 15:16:46 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

run the I-V analysis on a QUBIC dataset
'''
import os,sys
import datetime as dt
import numpy as np
from glob import glob
from satorchipy.datefunctions import str2dt

from qubicpack.utilities import Qubic_DataDir
from qubicpack.qubicfp import qubicfp
qubicfp.verbosity = 0

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
plt.rcParams['figure.figsize'] = [24,12]

n_TES = 128

def help():
    '''
    a help message
    '''
    msg = 'Analysis of I-V curves for QUBIC' 
    msg += '\nusage: %s [options]' % sys.argv[0]
    msg += '\nOPTIONS:'
    msg += '\n --start-time=YYYY-MM-DDTHH:MM:SS'
    msg += '\n       give the start time of the dataset you want to analyse'
    msg += '\n --datadir=/path/to/toplevel/data/directory'
    msg += '\n --TES-selection=<comma separated list of TES numbers for plotting, or "ALL">'
    msg += '\n                  WARNING:  If you select "ALL" you will have 258 plots'
    msg += '\n --help'
    msg += '\n        show this help message'
    msg += '\n'
    print(msg)
    return

def parseargs(argv):
    '''
    parse the arguments, and return a dictionary of parameters
    '''
    parms = {}
    parms['parameters_ok'] = False
    
    parms['do_plot'] = False
    parms['TES_selection'] = None

    parms['data_dir'] = None
    parms['start_time'] = None
    parms['day_str'] = None
    for arg in argv:
        
        if arg.find('help')>=0:
            help()
            return parms

        if arg.find('--start-time=')==0:
            date_str = arg.split('=')[-1]
            parms['start_time'] = str2dt(date_str)
            continue

        if arg.find('--datadir=')==0:
            parms['data_dir'] = arg.split('=')[-1]
            continue

        if arg.find('--TES_selection=')==0 or arg.find('--TES-selection=')==0:
            parms['do_plot'] = True
            TES_selection_strlist = arg.split('=')[-1].split(',')
            if len(TES_selection_strlist)>1:
                TES_selection = []
                for val in TES_selection_strlist:
                    TES_selection.append(int(val))
                TES_selection = np.array(TES_selection)
            else:
                keyword = TES_selection_strlist[0]
                if keyword.lower()=='all':
                    TES_selection = 1 + np.arange(n_TES)
                else:
                    TES_selection = [int(keyword)]

            parms['TES_selection'] = TES_selection
            continue


    if parms['start_time'] is not None:
        parms['day_str'] = parms['start_time'].strftime('%Y-%m-%d')

    if parms['day_str'] is None:
        print('Enter a valid date to find the data.')
        print('   Use the option --start-time=YYYY-MM-DDTHH:MM:SS\n')
        help()
        return parms

    parms['parameters_ok'] = True
    return parms
    
parms = parseargs(sys.argv)
for parm in parms.keys():
    if type(parms[parm])==str:
        cmd = "%s = '%s'" % (parm,parms[parm])
    elif type(parms[parm])==dt.datetime:
        cmd = "%s = str2dt('%s')" % (parm,parms[parm].strftime('%Y-%m-%dT%H:%M:%S'))
    elif type(parms[parm])==np.ndarray:
        cmd = '%s = %s' % (parm,str(parms[parm]).replace(' ',','))
    else:
        cmd = '%s = %s' % (parm,parms[parm])
    exec(cmd)
if not parameters_ok: quit()

# start processing
year_str = day_str.split('-')[0]

# use the command line to specify the data directory that is appropriate on your system
# Try to find the data
if data_dir is None:
    data_dir = os.path.join(Qubic_DataDir(),day_str)
    if not os.path.isdir(data_dir):
        data_dir = os.path.join(Qubic_DataDir(),year_str,day_str)
    if not os.path.isdir(data_dir):
        print('Could not find the data directory! %s' % data_dir)
        data_dir = None
else:
    try_data_dir = os.path.join(data_dir,day_str)
    if not os.path.isdir(try_data_dir):
        print('Could not find data directory: %s' % try_data_dir)
        try_data_dir = os.path.join(data_dir,year_str,day_str)
        if not os.path.isdir(try_data_dir):
            print('Could not find data directory: %s' % try_data_dir)
            try_data_dir = None
    data_dir = try_data_dir

if data_dir is None:  quit()
print('searching for data in: %s' % data_dir)

# search for datasets
pattern = '%s/%s__*' % (data_dir,start_time.strftime('%Y-%m-%d_%H.%M.??'))
dsets = glob(pattern)
if len(dsets)<1:
    print('No datasets found for pattern: %s' % pattern)
    quit()
    
dsets.sort()
dset = dsets[0]
print('Processing dataset: %s' % dset)

# Process data for I-V curves
a = qubicfp()
a.read_qubicstudio_dataset(dset)

# first of all, have a look at all timelines for all the TES
res = a.plot_timeline_focalplane(xwin=False)

# We run the filter which fits a model to each TES I-V curve
# This will take a few minutes
for asicobj in a.asic_list:
    if asicobj is None: continue
    read_saved_filter = asicobj.read_filter()
    if not read_saved_filter:
        f = asicobj.filter_iv_all(bias_margin=-3)

# now look at selected TES: plot timeline and I-V
if TES_selection is not None:
    for asicobj in a.asic_list:
        for TESnum in TES_selection:
            res = asicobj.plot_timeline(TES=TESnum,xwin=False)
            fig = asicobj.plot_iv(TES=TESnum,xwin=False)
    print('\nplots for individual TES are found in directory: ./%s\n' % start_time.strftime('%Y/%Y%m%d'))



# and save the I-V fitting parameters for next time, so we don't have to re-run the processing
for asicobj in a.asic_list: asicobj.save_filter()

# now plot all the I-V curves
res = a.plot_iv_focalplane(xwin=False)
print('\nfocal plane plots are in the current directory: ./\n')

