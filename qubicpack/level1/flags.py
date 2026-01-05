'''
$Id: flags.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Mon 23 Jan 2023 11:48:12 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

definition and tools for flagging data
'''
import numpy as np
from scipy.interpolate import interp1d

# The maximum number of flags possible is 64
# Data is flagged by a 64-bit integer
# the more serious flags are at the higher bit positions
nflags = 64
flag_definition = []
for idx in range(nflags):
    flag_definition.append('available')
    
flag_definition[63] = 'saturation'
flag_definition[57] = 'cosmic ray'
flag_definition[51] = 'uncorrected flux jump'
flag_definition[45] = 'end of scan'
flag_definition[39] = 'bath temperature above 350mK'
flag_definition[38] = 'bath temperature above 340mK'
flag_definition[37] = 'bath temperature above 330mK'
flag_definition[36] = 'bath temperature rising'
flag_definition[33] = '1K temperature above 1.3K'
flag_definition[32] = '1K temperature above 1.2K'
flag_definition[31] = '1K temperature above 1.1K'
flag_definition[30] = '1K temperature rising'
flag_definition[27] = 'corrected flux jump'
flag_definition[17] = 'baseline adjusted'

# for ease-of-use, we make a dictionary for the bit numbers
flag_bit = {}
for idx,flagdef in enumerate(flag_definition):
    if flagdef!='available':
        flag_bit[flagdef] = idx

def set_flag(flagdef,flagval):
    '''
    set the flag bit for the given flag value
    '''
    bit_position = flag_bit[flagdef]
    flagnumber = np.uint(2**bit_position)
    flagval = np.uint(flagval)

    return flagval | flagnumber

def unset_flag(flagdef,flagval):
    '''
    unset the flag bit for a given flag value
    '''
    bit_position = flag_bit[flagdef]
    flagnumber = np.uint(2**bit_position)
    flagval = np.uint(flagval)
    return flagval & ~flagnumber

def clear_flag(flagdef,flag_array):
    '''
    clear all the flags in the flag_array for the given flag type and return the flag_array
    The flag_array is the an array the size of a TOD array
    '''
    for idx,flagval in enumerate(flag_array):
        flag_array[idx] = unset_flag(flagdef,flagval)

    return flag_array

def isset_flag(flagdef,flagval):
    '''
    return True or False if the given flag is set
    '''
    bit_position = flag_bit[flagdef]
    flagnumber = np.uint(2**bit_position)
    flagval = np.uint(flagval)

    return ( (flagnumber & flagval)==flagnumber )
    

def show_flags(flagval):
    '''
    show the flags which are set
    '''
    for flag in flag_bit.keys():
        if isset_flag(flag,flagval):print(flag)
    return

    
def interpolate_flags(timeaxis,t_tod,flag_array):
    '''
    interpolate flags onto TOD interpolated time axis
    best method determined by JCH in notebook FlaggingAndInterpolation
    https://github.com/qubicsoft/qubic/blob/master/qubic/doc/Time%20Domain%20Analysis/FlaggingAndInterpolation.ipynb
    '''
    
    # function to interpolate extrapolating from previous value
    interp_prev_func = interp1d(timeaxis,flag_array,kind='previous',fill_value='extrapolate')

    # function to interpolate extrapolating from next value
    interp_next_func = interp1d(timeaxis,flag_array,kind='next',fill_value='extrapolate')

    # take the worst flag of the two interpolations
    flag_interps = np.array( [interp_prev_func(t_tod), interp_next_func(t_tod)] )            
    flag_array_interp = np.max(flag_interps,axis=0)
    return flag_array_interp

