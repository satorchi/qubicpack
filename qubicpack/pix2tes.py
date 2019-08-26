'''
$Id: pix2tes.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Fri 21 Jul 2017 08:26:39 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

translate the TES number to the physical location on focal plane grid
'''
from __future__ import division, print_function
import numpy as np
import os
import pickle

from qubicpack.utilities import asic_reversal_date, NPIXELS, NASIC, TES_index, ASIC_index

TES2PIX = None

def assign_pix_grid():
    '''
    generate the layout of the TES array:  These are the pixel numbers
    zeros indicate empty grid positions
    '''
    ncols=[8,9,11,12,12,14,15,15,16,17,17,17,17,17,17,17,17]
    rows=[]
    prev_row_end=0
    for N in ncols:
        row=np.arange(N)+prev_row_end+1
        rows.append(row)
        prev_row_end=row[-1]

    full_grid=np.zeros((17,17), dtype=int)
    nrow=0
    for row in rows:
        Ntes=len(row)
        full_grid[nrow,0:Ntes]=row
        nrow+=1

    return full_grid

def assign_tes_grid(obsdate=None):
    '''
    make a 17x17 array with TES id's
    '''
    pix_grid = assign_pix_grid()
    TES2PIX = assign_pix2tes(obsdate)
    tes_grid = np.recarray(names='PIX,TES,ASIC',formats='int,int,int',shape=pix_grid.shape)    

    nrows = pix_grid.shape[0]
    ncols = pix_grid.shape[1]
    
    for row in range(nrows):
        for col in range(ncols):
            tes_grid[row,col].PIX = 0
            tes_grid[row,col].TES = 0
            tes_grid[row,col].ASIC = 0
            
            # the pixel identity associated with its physical location in the array
            physpix=pix_grid[row,col]
            tes_grid[row,col].PIX = physpix
            pix_index=physpix-1
            if physpix!=0:
                for asic_ctr in range(2):
                    asic = asic_ctr + 1
                    asic_idx = ASIC_index(asic)
                
                    if physpix in TES2PIX[asic_idx]:
                        TES=pix2tes(physpix)[0]
                        TES_idx = TES - 1
                        tes_grid[row,col].ASIC = asic
                        tes_grid[row,col].TES = TES
    return tes_grid                        

def assign_pix2tes(obsdate=None):
    '''
    translation from QubicStudio pixel number to physical location on the array
    and the opposite
    '''
    global TES2PIX

    # give the QubicStudio TES_index and get the physical pixel

    # asic "upper/left"
    asic_upper_left=[64,40,52,1004,28,39,51,63,120,104,91,76,50,62,90,119,89,103,27,38,61,118,75,135,49,88,151,152,74,17,26,37,102,117,134,1003,60,73,87,101,16,25,36,48,133,150,167,8,59,72,100,116,15,24,35,47,166,183,184,7,58,71,85,86,23,34,46,1002,182,199,6,14,57,132,149,165,13,22,33,45,84,99,115,5,44,56,69,70,4,12,21,32,181,198,215,216,43,55,164,1001,3,11,20,31,98,114,131,148,42,54,68,83,2,10,19,30,97,113,130,147,41,53,67,82,1,9,18,29]

    # asic "lower/right"
    asic_lower_right=[210,209,236,1008,214,213,212,211,223,222,221,220,227,226,225,224,231,230,229,228,240,239,238,237,244,243,242,241,248,247,246,245,161,160,159,1007,218,217,163,162,175,204,203,219,179,178,177,176,207,206,205,180,193,192,208,191,197,196,195,194,235,234,233,232,155,171,170,1006,124,140,156,172,128,127,126,125,187,186,185,129,141,157,189,188,145,144,143,142,202,201,200,146,158,174,190,173,77,66,65,1005,81,80,79,78,136,121,105,92,94,93,106,137,122,138,96,95,108,154,107,153,109,139,169,168,123,112,111,110]

    if obsdate is not None and obsdate < asic_reversal_date:
        TES2PIX = np.array([asic_upper_left,asic_lower_right])
    else:
        # we've switched the wiring for asic_1 and asic_2, and also the order of the columns
        asic_lower_right_row0=asic_lower_right[0:32]
        asic_lower_right_row1=asic_lower_right[32:64]
        asic_lower_right_row2=asic_lower_right[64:96]
        asic_lower_right_row3=asic_lower_right[96:128]
        new_asic_lower_right=np.concatenate([asic_lower_right_row3,asic_lower_right_row2,asic_lower_right_row1,asic_lower_right_row0])

        asic_upper_left_row0=asic_upper_left[0:32]
        asic_upper_left_row1=asic_upper_left[32:64]
        asic_upper_left_row2=asic_upper_left[64:96]
        asic_upper_left_row3=asic_upper_left[96:128]
        new_asic_upper_left=np.concatenate([asic_upper_left_row3,asic_upper_left_row2,asic_upper_left_row1,asic_upper_left_row0])

        TES2PIX = np.array([new_asic_lower_right,new_asic_upper_left])

    return TES2PIX

def tes2pix(TES,asic):
    global TES2PIX
    if TES2PIX is None: TES2PIX = assign_pix2tes()
    TES_idx = TES_index(TES)
    if TES_idx is None: return None
    if TES_idx == -1: return -1
    PIX = TES2PIX[ASIC_index(asic),TES_idx]
    return PIX

def pix2tes(PIX, asic=None):
    '''
    convert a pixel number to a TES number
    the ASIC number is not required, but is there for historical purposes
    so we don't have to change other scripts that use this function
    '''
    global TES2PIX
    if TES2PIX is None: TES2PIX = assign_pix2tes()
    pix_index=PIX-1
    asic_idx = None
    for idx in range(2):
        if PIX in TES2PIX[idx,:]:
            asic_idx = idx
            break
    if asic_idx is None:
        print('ERROR! invalid Pixel number request: %i' % PIX)
        return None

    idx_tuple = np.where(TES2PIX[asic_idx]==PIX)
    return (idx_tuple[0][0]+1, asic_idx+1)

