"""
Script that reduces the size of a qubic science fits file
The input file should be a science file
The output will be a identical file but with a reduced number of time series points
The output file should point to a new directory with same format as original e.g. raws, hks, sums
Output file should have same name but in new directory since a.read_qubicstudio_dataset needs
filename to have asic information
Repeat this for each asic
startpoint and endpoint refer to time series points. e.g. if 10,000 time points for each TES
startpoint=0, and finishpoint=100 will produce a dataset with only the first 100 datapoints
"""
from __future__ import division, print_function

from astropy.io import fits
import matplotlib.pyplot as plt

def reduce_qubic_dataset(qfitsfile, outfile, startpoint, finishpoint):
    #load qubic science file
    print('loading and reducing {}'.format(qfitsfile))
    hdul = fits.open(qfitsfile, memmap=True)

    #select TOD data for reduction
    data = hdul[1].data
    #reduce data based on inputs
    print("starting data shape {}".format(data.shape))
    data = data[startpoint:finishpoint]
    print("final data shape {}".format(data.shape))
    hdul[1].data = data
    #plot reduced data
    plt.plot(data['GPSDate'], data['pixel93'], 'o-')
    plt.title('sample of reduced data for TES 93')
    #write data to new file
    hdul.writeto(outfile, overwrite=True)
    print('reduced file saved here: {}'.format(outfile))
    print('you should now copy the reduced fits files for asic 1&2 into a directory with correct structure with RAWS, HKS, SUMS. The 2 science files should have correct names too to load asic data properly')
    
    return