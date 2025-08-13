#!/usr/bin/env python
'''
$Id: asd_analysis_example.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Mon 05 Mar 2018 15:11:53 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.


This is an example script for analysing the QUBIC noise measurement
You can run this in ipython, or Jupyter notebook.
'''
'''
Before beginning, you should have pystudio installed, and copy the example data files.

You will also need my datefunctions module

### INSTALLATION OF pystudio ###

get the package from github

$ git clone https://github.com/satorchi/pystudio.git
$ cd pystudio
$ sudo python ./setup.py install

$ cd ../
$ git clone https://github.com/satorchi/mypy.git
$ cd mypy
$ sudo python ./setup.py install

### Data files ###
You should copy the following data files into your working directory

QUBIC_timeline_20180301T142554UTC.fits
'''

# import necessary  stuff
import os,sys
if sys.argv[0].find('ipykernel')>=0:jupyter=True

if jupyter:        
    %matplotlib notebook

import matplotlib.pyplot as plt
from qubicpack import qubicpack as qp
from qubicpack.temperature_analysis import *
from qubicpack.plot_physical_layout import *
import datetime as dt


# create the qubicpack object and read the data
d0=qp()
if jupyter:d0.figsize=10,5
result=d0.read_fits('QUBIC_timeline_20180301T142554UTC.fits')

# The data contains multiple timelines at different temperatures
# you can print out a summary table
result=print_datlist(d0)

# we will work with the third timeline.
# (we set timeline_index=2)

# first of all, have a look at all timelines for all the TES
result=d0.plot_timeline_physical_layout(timeline_index=2)

# let's look at one TES and also look at the power spectrum
result=d0.plot_ASD(85,timeline_index=2)

# we can also bin the data.  Let's try with nbins=10
result=d0.plot_ASD(85,timeline_index=2,nbins=10)

# let's make all the plots, and get a list of the results
# we will make all the plots have the same ASD scale using the amin and amax keywords
# This will generate all the png files, but we won't show the plots on screen because there are too many.
asdresults=d0.plot_ASD_all(timeline_index=2,amin=1e-6,amax=1e-3,nbins=10)

# we can use the results list to generate the report document
texfile=d0.make_ASD_tex_report(reslist=asdresults,timeline_index=2)

# and we have to compile the latex file
result=os.system('pdflatex '+texfile)


