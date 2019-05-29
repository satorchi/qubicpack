#!/usr/bin/env python
'''
$Id: utilities.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Fri 24 May 2019 15:19:42 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

common utilities used in the qubicpack classes
(see also pix2tex.py)
'''
import datetime as dt

# on 6 Feb 2018, we reversed the wires for the ASICs
# so now QubicStudio and the dilution fridge use the same ASIC designation
asic_reversal_date=dt.datetime.strptime('2018-02-06 18:00','%Y-%m-%d %H:%M')

# number of pixels in the QUBIC detector matrix per ASIC
NPIXELS = 128

# maximum number of ASIC possible
NASIC = 2 # this will change to 16 or 32

# default figure size for plots
FIGSIZE = (12.80,6.40)


def TES_index(TES):
    '''
    return the index (counting from 0) given the TES number (counting from 1)
    '''
    if (not isinstance(TES,int))\
       or TES<1\
       or TES>NPIXELS:
        print('TES should have a value between 1 and %i' % NPIXELS)
        return None
    TES_idx=TES-1
    return TES_idx

def ASIC_index(asic):
    '''
    return the asic index (counting from zero) given the asic number
    the asic index is either 0 or 1 which is used for plotting the focal plane
    '''
    if (not isinstance(asic,int))\
       or asic<1\
       or asic>NASIC:
        print('asic should have a value between 1 and %i' % NASIC)
        return None
    asic_idx = asic-1
    remainder = asic_idx % 2
    if remainder == 0:
        return 0
    return 1

