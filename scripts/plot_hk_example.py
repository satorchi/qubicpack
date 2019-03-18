#!/usr/bin/env python
'''
$Id: plot_hk_example.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Mon 31 Dec 2018 00:41:08 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

this is a simple example for reading the QUBIC housekeeping FITS file

'''
from __future__ import division, print_function
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import datetime as dt

# open the fits file
hdulist = fits.open('QUBIC_HK_20180612-103246.fits')

# look at the primary header
print(hdulist[0].header.values)

# the number of binary tables
nTables = len(hdulist) - 1

# the housekeeping labels
labels = []
for idx in range(nTables):
    table_number = idx + 1
    label = hdulist[table_number].header['TTYPE2']
    labels.append(label)


# find all the labels that have 300mK stage
labels_300mK = []
label_idx = []
for idx,label in enumerate(labels):
    if label.find('300mK')>=0:
        labels_300mK.append(label)
        label_idx.append(idx)



plt.ioff()

# make a plot with the 300mK data
for idx in label_idx:
    table_number = idx + 1
    label = hdulist[table_number].header['TTYPE2']
    tstamps = hdulist[table_number].data.field(0)
    temperatures = hdulist[table_number].data.field(1)

    # convert timestamps to human readable dates
    dates = []
    for tstamp in tstamps:
        dates.append(dt.datetime.fromtimestamp(tstamp))

    
    plt.plot(dates,temperatures,label=label)

plt.legend()
plt.show()
        


