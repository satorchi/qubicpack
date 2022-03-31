#!/usr/bin/env python3
'''
$Id: Vphi_analysis.py
$auth: Alejandro Almela
$auth: Guillaume Stankowiak
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Thu 24 Mar 2022 07:52:27 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

translation/update of SQUID analysis script by Guillaume Stankowiak and Alejandro Almela

'''
import os,sys
import numpy as np
import scipy as sp

from glob import glob
import datetime as dt

from qubicpack.utilities import Qubic_DataDir
from qubicpack.qubicfp import qubicfp
qubicfp.verbosity = 0

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
plt.ioff()
plt.rcParams['figure.figsize'] = [24,12]

# a help message
def help():
    msg = 'Analysis of QUBIC-SQUID data for optimisation' 
    msg += '\nusage: %s [options]' % sys.argv[0]
    msg += '\nOPTIONS:'
    msg += '\n --day=YYYY-MM-DD'
    msg += "\n       for example 2019-04-24 the SQUID data used in Guillaume's thesis"
    msg += '\n --hour=HH'
    msg += '\n        this can also be a pattern, for example 1[23] will find the data starting at 12:00 and 13:00'
    msg += '\n --datadir=/path/to/toplevel/data/directory'
    msg += '\n --squid_selection=<comma separated list of squid indexes for plotting, or "ALL">'
    msg += '\n                  WARNING:  If you select "ALL" you will have ~8000 plots'
    msg += '\n --help'
    msg += '\n'
    print(msg)
    return

do_plot = False
# squid_selection = np.arange(128)
squid_selection = [95,94,63,62] # for plotting

# choose the data
# day_str = '2019-04-24' # for comparison with previous data
# day_str = '2022-03-23'
# day_str = '2022-03-25'
# hour_pattern = '1[56]'

day_str = None
hour_pattern = ''


# parse command line arguments
data_dir = None
for arg in sys.argv:
    if arg.find('--day=')==0:
        day_str = arg.split('=')[-1]
        continue

    if arg.find('--hour=')==0:
        hour_pattern = arg.split('=')[-1]
        continue

    if arg.find('--datadir=')==0:
        data_dir = arg.split('=')[-1]
        continue

    if arg.find('--squid_selection=')==0:
        squid_selection_strlist = arg.split('=')[-1].split(',')
        if len(squid_selection)>1:
            squid_selection = []
            for val in squid_selection_strlist:
                squid_selection.append(int(val))
        else:
            keyword = squid_selection_strlist[0]
            if keyword.lower()=='all':
                squid_selection = np.arange(n_squids)
            else:
                squid_selection = [int(keyword)]
        
        do_plot = True
        continue

    if arg.find('help')>=0:
        help()
        quit()

if day_str is None:
    print('Enter a valid date to find the data.  Use the option --day=YYYY-MM-DD\n')
    help()
    quit()

# start processing
year_str = day_str.split('-')[0]

# use the command line to specifiy the data directory that is appropriate on your system
# Try to find the data
if data_dir is None:
    data_dir = os.sep.join([Qubic_DataDir(),year_str,day_str])
    if not os.path.isdir(data_dir):
        data_dir = os.sep.join([Qubic_DataDir(),day_str])
    if not os.path.isdir(data_dir):
        print('Could not find the data directory! %s' % data_dir)
        quit()
else:
    try_data_dir = os.path.join(data_dir,day_str)
    if not os.path.isdir(try_data_dir):
        print('Could not find data directory: %s' % try_data_dir)
        try_data_dir = os.path.join(year_str,data_dir,day_str)
        if not os.path.isdir(try_data_dir):
            print('Could not find data directory: %s' % try_data_dir)
            quit()
    data_dir = try_data_dir


# search for the datasets
pattern = '%s/*_%s*_opt_bias*' % (data_dir,hour_pattern)
dsets = glob(pattern)
if len(dsets)<1:
    print('No datasets found for pattern: %s' % pattern)
    quit()
    
dsets.sort()

# Find the dataset with bias index=7.  This one is used to calculate the offset value.
gotit = False
for idx,dset in enumerate(dsets):
    squid_index = dset.split('_')[-1].strip()
    if squid_index=='7':
        gotit = True
        idxoffset = idx
        break
# reorder the dataset list to do first the one which is used to calculate the Vsqoffset
if gotit:
    print('Found index 7')
    new_dsets = []
    new_dsets.append(dsets[idxoffset])
    for idx,dset in enumerate(dsets):
        if idx==idxoffset:continue
        new_dsets.append(dset)
    dsets = new_dsets

print('Vsqoffset will be calculated from dataset: %s' % dsets[0])


### some parameters ###
Min=1./10.4E-6 # Mutual inductance
gain = 70.*100
good_threshold = 5
n_indexes = 16 # number of possible bias values
n_squids = 128 # number of squid per ASIC
n_asics = 2

plotname_prefix = 'SQUID_%s' % day_str


# The value of the SQUID bias for each index
I = np.zeros(n_indexes)
I[0] = 0
I[1] = 5.1
I[2] = 7.65
I[3] = 10.20
I[4] = 12.75
I[5] = 15.30
I[6] = 17.85
I[7] = 20.41
I[8] = 22.96
I[9] = 25.51
I[10] = 28.06
I[11] = 30.61
I[12] = 33.16
I[13] = 35.71
I[14] = 38.26
I[15] = 40.81


# maximum length of a timeline
# this create the limit of the array, if not issues with some file of 4849 variable instead of 4950
lmt = 5000 

# Creation of empty array, fill with NaN
Vsquid = np.zeros((n_indexes,n_squids,lmt,n_asics))
Vsquid[:,:,:,:] = np.nan
VsquidSG = np.zeros((n_indexes,n_squids,lmt,n_asics))
VsquidSG[:,:,:,:] = np.nan

percent_good = np.zeros((n_indexes,n_asics))

histo = np.zeros((n_indexes,n_asics)) # create tab for peak to peak val
data = np.zeros((n_squids,n_indexes,n_asics)) # create a tab for each squid to keep all ptp value for each 
invdata = np.zeros((n_indexes,n_squids,n_asics))

# SQUID voltage offset
Vsqoffset = np.zeros((n_squids,n_asics))

# Now go through the data
for dset in dsets:
    squid_index_str = dset.split('_')[-1].strip()
    try:
        squid_index = int(squid_index_str)
    except:
        print('ERROR!  Is this a SQUID data acquisition? %s' % dset)
        continue

    # the first 8 index are too low to induce curves, so we start at 7 
    # if squid_index < 7:
    #     print('skipping SQUID index: %i' % squid_index)
    #     continue

    print('SQUID index: %i' % squid_index)
    
    # reading the data from qubicpack
    b = qubicfp()
    b.read_qubicstudio_dataset(dset)

    # bias is applied separately on each ASIC
    for asic in b.asic_list:

        ASICnum = asic.asic
        ASICidx = ASICnum - 1
        Rbias = asic.Rbias
        
        # Amplitude peak to peak of the sinus
        amp = asic.max_bias - asic.min_bias
        
        # offset of the sinus signal
        offset = 0.5*(asic.max_bias + asic.min_bias)
        Vbias = (offset+amp*asic.bias_phase()/2.)

        # this allow to avoid the noncoherent points in raw data for the flux
        sorted_index = np.argsort(Vbias)

        # convert to voltage (?)
        Vadjustment = 62.5/(70.*100)/(asic.nsamples-asic.n_masked())

        # timeline length limit
        lmt2 = asic.timeline(TES=1).shape[0]
        
        # go through the TES detectors
        for TESidx in range(n_squids):
            TESnum = TESidx + 1

            if squid_index == 7:
                Vsqoffset[TESidx,ASICidx] = np.mean( asic.timeline(TES=TESnum)*Vadjustment )

            timeline = asic.timeline(TES=TESnum) * Vadjustment - Vsqoffset[TESidx,ASICidx]

            
            # savitzky golay filter
            filt = sp.signal.savgol_filter(timeline, 51, 3) 
    
            Vsquid[squid_index,TESidx,0:lmt2,ASICidx] = timeline
            VsquidSG[squid_index,TESidx,0:lmt2,ASICidx] = filt

            histo[squid_index,ASICidx] = np.max(filt) - np.min(filt)
            data[TESidx,squid_index,ASICidx] = histo[squid_index,ASICidx]
            invdata[squid_index,TESidx,ASICidx] = histo[squid_index,ASICidx]


            ### plot parameters ###
            if do_plot and TESidx in squid_selection:
                fig = plt.figure()
                V0 = Vbias.min()
                # xpts = (0.5*Min/Rbias)*( Vbias[sorted_index] - V0 )
                # ypts = filt[sorted_index]
                xpts = (0.5*Min/Rbias)*( Vbias - V0 )
                ypts = filt
                plt.plot(xpts,ypts)
                plt.grid()
                plt.xlabel('Flux (in quantum of flux)')
                plt.ylabel('Voltage ($\mu$V)')
                plt.title('%s ASIC%i, SQUID number %i, bias index %i' % (day_str,ASICnum,TESnum,squid_index))

                figname = '%s_ASIC%02i_TES%03i_bias-index%02i_flux.png' % (plotname_prefix,ASICnum,TESnum,squid_index)
                plt.savefig(figname,format='png',dpi=100,bbox_inches='tight')
                plt.close(fig)
                

        # min, max, average
        Vmoy = np.nanmean(VsquidSG, axis=2)
        Vmin = np.nanmin(VsquidSG, axis=2)
        Vmax = np.nanmax(VsquidSG, axis=2)


        # Create variable , readable for the plot
        Vmoy2 = np.asarray([(-Vmoy[:,_,ASICidx]+Vmoy[0,_,ASICidx])*62.5/gain for _ in np.arange(n_squids)]) 
        Vmin2 = np.asarray([(-Vmin[:,_,ASICidx]+Vmoy[0,_,ASICidx])*62.5/gain for _ in np.arange(n_squids)]) 
        Vmax2 = np.asarray([(-Vmax[:,_,ASICidx]+Vmoy[0,_,ASICidx])*62.5/gain for _ in np.arange(n_squids)]) 


        ngood = (data[:,squid_index,ASICidx]>= good_threshold).sum()
        percent_good[squid_index,ASICidx] = 100*ngood/128

        
        ### plot all SQUIDS if requested ###
        if do_plot:
            for TESidx in squid_selection:
                TESnum = TESidx + 1

                fig = plt.figure()
                plt.plot(I,Vmin2[TESidx,:], label= "maximal Value" )
                plt.plot(I,Vmax2[TESidx,:], label= "minimal Value")
                plt.plot(I,Vmoy2[TESidx,:], label= "mean Value")

                plt.grid()
                plt.xlabel('Intensity ($\mu$A)')
                plt.ylabel('Voltage ($\mu$V)')
                plt.title('%s ASIC%i SQUID number %i, bias index %i' % (day_str,ASICnum,TESnum,squid_index)) 
                plt.legend(loc='upper right', bbox_to_anchor=(0.25,1 ))
                figname = '%s_ASIC%02i_TES%03i_bias-index%02i_IV.png' % (plotname_prefix,ASICnum,TESnum,squid_index)
                plt.savefig(figname,format='png',dpi=100,bbox_inches='tight')
                plt.close(fig)
    
## finished going through the datasets



# final plots
for ASICidx in range(n_asics):
    ASICnum = ASICidx + 1


    fig = plt.figure()
    plt.plot(np.arange(n_indexes),percent_good[:,ASICidx])
    plt.grid()
    plt.ylabel("Percentage working SQUID >%.1f$\mu$V ($\mu$V)" % good_threshold)
    plt.xlabel('Index')
    plt.title("%s ASIC%i Working SQUID by index" % (day_str,ASICnum))
    figname = '%s_ASIC%02i_good-percentage.png' % (plotname_prefix,ASICnum)
    plt.savefig(figname,format='png',dpi=100,bbox_inches='tight')
    

    fig = plt.figure()
    plt.plot(data[:,:,ASICidx])
    plt.grid()
    plt.ylabel("PtP value")
    plt.xlabel("Number of SQUID")
    plt.title('%s ASIC%i' % (day_str,ASICnum))
    figname = '%s_ASIC%02i_PtP.png' % (plotname_prefix,ASICnum)
    plt.savefig(figname,format='png',dpi=100,bbox_inches='tight')


    fig = plt.figure()
    plt.plot(invdata[:,:,ASICidx])
    plt.grid()
    plt.xlabel("Intensity (index of I)")
    plt.ylabel("SQUIDs")
    plt.title('%s ASIC%i' % (day_str,ASICnum))
    figname = '%s_ASIC%02i_inversePtP.png' % (plotname_prefix,ASICnum)
    plt.savefig(figname,format='png',dpi=100,bbox_inches='tight')

    # argmax take the position of the maxvalue for each squid

    fig = plt.figure()
    plt.hist(np.argmax(data[:,:,ASICidx], axis=1), range=[0,n_indexes], bins=16)
    plt.grid()
    plt.ylabel("Number of SQUIDs")
    plt.xlabel("Index of current")
    plt.title("%s Histogram of the optimum current for the SQUID response for ASIC %i" % (day_str,ASICnum))
    figname = "%s_ASIC%02i_Histogram.png" % (plotname_prefix,ASICnum)
    plt.savefig(figname,format='png',dpi=100,bbox_inches='tight')
    
    fig = plt.figure()
    plt.hist(data[:,9,ASICidx],range=[0,30], bins=30, alpha = 0.5, color='red' ,label="Isq = 25.5 $\mu$A")
    plt.hist(data[:,10,ASICidx],range=[0,30], bins=30, alpha = 0.5, color='blue',label="Isq = 28 $\mu$A ")
    plt.hist(data[:,11,ASICidx],range=[0,30], bins=30, alpha = 0.5, color='green', label="Isq = 30.6 $\mu$A")
    plt.legend()
    plt.grid()
    plt.xlabel("Voltage ($\mu$V)")
    plt.xlim(0,25)
    plt.ylabel('Number of SQUID')
    plt.title("%s ASIC %i histogram" % (day_str,ASICnum))
    figname = "%s_ASIC%02i_Histogram_multindex.png" % (plotname_prefix,ASICnum)
    plt.savefig(figname,format='png',dpi=100,bbox_inches='tight')

