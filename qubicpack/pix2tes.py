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

def assign_pix_grid(self):
    '''
    generate the layout of the TES array
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

    self.pix_grid=full_grid
    return full_grid

def assign_pix2tes(self):
    '''
    translation from QubicStudio pixel number to physical location on the array
    and the opposite
    '''

    # give the QubicStudio TES_index and get the physical pixel

    # asic "upper/left"
    asic_upper_left=[64,40,52,1004,28,39,51,63,120,104,91,76,50,62,90,119,89,103,27,38,61,118,75,135,49,88,151,152,74,17,26,37,102,117,134,1003,60,73,87,101,16,25,36,48,133,150,167,8,59,72,100,116,15,24,35,47,166,183,184,7,58,71,85,86,23,34,46,1002,182,199,6,14,57,132,149,165,13,22,33,45,84,99,115,5,44,56,69,70,4,12,21,32,181,198,215,216,43,55,164,1001,3,11,20,31,98,114,131,148,42,54,68,83,2,10,19,30,97,113,130,147,41,53,67,82,1,9,18,29]

    # asic "lower/right"
    asic_lower_right=[210,209,236,1008,214,213,212,211,223,222,221,220,227,226,225,224,231,230,229,228,240,239,238,237,244,243,242,241,248,247,246,245,161,160,159,1007,218,217,163,162,175,204,203,219,179,178,177,176,207,206,205,180,193,192,208,191,197,196,195,194,235,234,233,232,155,171,170,1006,124,140,156,172,128,127,126,125,187,186,185,129,141,157,189,188,145,144,143,142,202,201,200,146,158,174,190,173,77,66,65,1005,81,80,79,78,136,121,105,92,94,93,106,137,122,138,96,95,108,154,107,153,109,139,169,168,123,112,111,110]

    if self.obsdate < self.asic_reversal_date:
        self.TES2PIX=np.array([asic_upper_left,asic_lower_right])
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

        self.TES2PIX=np.array([new_asic_lower_right,new_asic_upper_left])

    return

def tes2pix(self,TES):
    TES_index=self.TES_index(TES)
    PIX=self.TES2PIX[self.asic_index(),TES_index]
    return PIX

def pix2tes(self,PIX):
    pix_index=PIX-1
    if not PIX in self.TES2PIX[self.asic_index(),:]:
        self.printmsg('ERROR! invalid Pixel number request: %i' % PIX)
        return None

    TES_index=-1
    gotit=False
    while (not gotit) and TES_index<128:
        TES_index+=1
        _pix=self.TES2PIX[self.asic_index(),TES_index]
        if _pix==PIX:
            gotit=True

    if not gotit:
        self.printmsg('ERROR! how did we get here? PIX=%i' % PIX)
        return None
    
    TES=TES_index+1
    return TES    

def assign_lookup_table(self):
    '''
    make the lookup table with results for comparison

    this was only done for array P73, so don't ask for this table otherwise
    '''
    filename='TES_translation_table.pickle'
    cwd=os.getcwd()
    dirs=[cwd]
    if not isinstance(self.datadir,str):
        self.assign_datadir()
    if isinstance(self.datadir,str):
        dirs.append(self.datadir)
        
    gotit=False
    for d in dirs:
        filename_fullpath='%s/%s' % (d,filename)
        if os.path.exists(filename_fullpath):
            gotit=True
            break
    if not gotit:
        if self.detector_name=='P73':
            self.printmsg('WARNING! Cannot find translation table file: %s' % filename_fullpath,verbosity=2)
            self.printmsg('Open loop and Room Temperature tests will not be noted in plots etc.',verbosity=2)
        self.transdic=None
        return None

    h=open(filename_fullpath,'r')
    self.transdic=pickle.load(h)
    h.close()
    return self.transdic

def lookup_TEStable(self,key='PIX',value=100):
    if self.transdic is None:
        self.printmsg('No translation table.  Please load.')
        return None
    
    if not key in self.transdic[0].keys():
        self.printmsg('Please enter a valid key.  Choose from:',verbosity=0)
        for k in transdic[0].keys():
            print('    %s' % k)
        return None
    
    ret=[entry for entry in self.transdic if entry[key]==value]
    if len(ret)==1:
        ret=ret[0]
    return ret
