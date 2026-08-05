'''
$Id: plot_options.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Wed 22 Jul 2026 13:50:46 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

options for various housekeeping plots
'''
import os,sys,subprocess,re
import datetime as dt
TZUTC = dt.timezone.utc

from satorchipy.datefunctions import str2dt, tstamp2dt, utcnow, utcfromtimestamp
from satorchipy.plotfunctions import nice_plot_colours as colours
from satorchipy.plotfunctions import nice_plot_markers as markers
from matplotlib import pyplot as plt

from ..utilities import hostname
from .utilities import qc_hk_dir

# parse command line arguments
plot_options = {}
plot_options['download'] = True
plot_options['qubic-central'] = 'qubic' # used for ssh connection
plot_options['estimate cold date'] = True
plot_options['show differences'] = True
plot_options['hk_dir'] = None
plot_options['tstart'] = None
plot_options['tend'] = None
plot_options['heaters'] = False
plot_options['pressure'] = False
plot_options['temperatures'] = None
plot_options['events'] = None

if os.path.isfile('plot_options.txt'):
    print('reading plot options from file: plot_options.txt')
    h = open('plot_options.txt','r')
    lines = h.read().split('\n')
    h.close()
    for line in lines:
        data_line = line.split('#')[0]
        if len(data_line)==0: continue
        optiondef = data_line.split('=')
        if len(optiondef)!=2: continue
        option = optiondef[0].strip()
        val_str = optiondef[1].strip()

        # special case for dates
        val_dt = str2dt(val_str)
        if val_dt is not None:
            val = val_dt.replace(tzinfo=TZUTC)
        else:
            try:
                val = eval(val_str)
            except:
                val = val_str

        # special case for list of values
        if val_str.find(',')>0:
            val = val_str.split(',')
            
        plot_options[option] = val
        
for arg in sys.argv:
    if arg=='--nodownload':
        plot_options['download'] = False
        continue

    if arg.find('--flagpos=')==0:
        plot_options['flagpos'] = float(arg.split('=')[-1])
        continue

    if arg.find('--tstart=')==0:
        plot_options['tstart'] = str2dt(arg.split('=')[-1])
        continue

    if arg.find('--tend=')==0:
        date_str = arg.split('=')[-1]
        plot_options['tend'] = str2dt(date_str)
        continue

    if arg.find('--noestimate')==0:
        plot_options['estimate cold date'] = False
        continue

    if arg.find('--nodiff')==0:
        plot_options['show differences'] = False
        continue

    if arg.find('--hk_dir=')==0:
        plot_options['hk_dir'] = arg.split('=')[-1]
        if not os.path.isdir(plot_options['hk_dir']):
            try:
                os.makedirs(plot_options['hk_dir'], exist_ok=True)
            except:
                print('PLOT_OPTIONS: could not make directory: %s' % plot_options['hk_dir'])
                plot_options['hk_dir'] = None
        continue

    if arg.find('--heater')==0:
        plot_options['heaters'] = True
        continue

    if arg.find('--pressure')==0:
        plot_options['pressure'] = True
        continue

    if arg.find('--log')==0:
        plot_options['log'] = True
        continue

    # generic
    match  = re.search('--(.*)=',arg)
    if match:
        arg_str = match.groups()[0]
        val_str = arg.split('=')[-1]
        
        # special case for list of values
        if val_str.find(',')>0:
            val_list = val_str.split(',')
            for idx,val in enumerate(val_list):
                val_list[idx] = val.strip()
            val_str = val_list
        try:
            val = eval(val_str)
        except:
            val = val_str
        plot_options[arg_str] = val
        
        continue

    match  = re.search('--(.*)',arg)
    if match:
        arg_str = match.groups()[0]
        val = True
        plot_options[arg_str] = val
        continue
    

# print plot_options
print('PLOT OPTIONS')
for key in plot_options.keys():
    val = plot_options[key]
    if isinstance(val,dt.datetime):
        print('%s = %s' % (key,val.strftime('%Y-%m-%d %H:%M:%S %Z')))
    else:
        print('%s = %s' % (key,plot_options[key]))
    
hk_dir = plot_options['hk_dir']
this_year = utcnow().strftime('%Y')
if hostname.find('qubic-central')<0:
    if hk_dir is None: hk_dir = '%s/hk' % this_year
else:
    hk_dir = qc_hk_dir

plot_options['hk_dir'] = hk_dir

# plot default parameters
plt.rcParams['figure.figsize'] = [24,12]
plt.rcParams['agg.path.chunksize'] = 10000
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['axes.labelsize'] = 20

boxprops = {}
boxprops['alpha'] = 0.99
boxprops['facecolor'] = 'white'
boxprops['boxstyle'] = 'round'
boxprops['edgecolor'] = 'black'
