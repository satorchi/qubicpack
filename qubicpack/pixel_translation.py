'''
$Id: pixel_translation.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Mon 19 Aug 2019 01:15:24 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

translate various pixel identification numbers in the QUBIC focal plane
'''
from __future__ import division, print_function

import numpy as np
from matplotlib import pyplot as plt

from qubicpack.pix2tes import assign_tes_grid, tes2pix

quadrant_colour = ['blue','red','green','purple']
asic_colour = ['blue','darkblue','red','#cc0000','green','#00cc00','purple','#7210a7']
FPidentity = None


def plot_square(x,y,colour='black',label='null',labelcolour='white',ax=None,fontsize=10):
    '''
    plot a square with a label
    '''
    sidelen = 0.8
    if ax is None: ax = plt.gca()
    
    xpts = x + np.array([-0.5*sidelen,  0.5*sidelen,  0.5*sidelen,  -0.5*sidelen, -0.5*sidelen])
    ypts = y + np.array([-0.5*sidelen, -0.5*sidelen,  0.5*sidelen,   0.5*sidelen, -0.5*sidelen])

    ax.fill(xpts,ypts,color=colour)
    ax.text(x,y,label,va='center',ha='center',fontsize=fontsize,color=labelcolour)
    return


def plot_id_focalplane(figsize=(30,30)):
    '''
    plot all the different identity names of each pixel in the focal plane

    FPidentity is a recarray of shape 34*34
    '''
    global FPidentity
    if FPidentity is None: FPidentity = make_id_focalplane()

    scale_factor = figsize[0]
    title_fontsize = 0.67*scale_factor
    label_fontsize = 0.2*scale_factor

    fig = plt.figure(figsize=figsize)
    fig.canvas.set_window_title('plt: QUBIC Focal Plane ID Matrix')
    ax = fig.add_axes([0,0,1,1])
    ax.text(0.5,0.96,'QUBIC Focal Plane ID Matrix',
            ha='center',va='bottom',transform=ax.transAxes,fontsize=title_fontsize)
    ax.set_xlim(-1,35)
    ax.set_ylim(-1,35)
    ax.set_aspect('equal')
    
    for fp_idx in range(len(FPidentity)):
        txt = 'Q%i' % (FPidentity[fp_idx].quadrant)
        quadrant = FPidentity[fp_idx].quadrant
        asic = FPidentity[fp_idx].ASIC
        colour = asic_colour[asic-1]
        row = FPidentity[fp_idx].row
        col = FPidentity[fp_idx].col
        if FPidentity[fp_idx].TES==0:
            colour = 'black'
            txt += '\nFP%4i' % FPidentity[fp_idx].index
        else:
            txt += ' %s\nFP%4i\nPIX%03i\nASIC%i\nTES%03i'\
                % (FPidentity[fp_idx].matrix.decode('UTF-8'),
                   FPidentity[fp_idx].index,
                   FPidentity[fp_idx].PIX,
                   FPidentity[fp_idx].ASIC,
                   FPidentity[fp_idx].TES)
        plot_square(col,row,colour=colour,labelcolour='white',label=txt,fontsize=label_fontsize)
    return

def make_id_focalplane():
    '''
    make the matrix which has all the translation info for pixel identities
    '''
    global FPidentity
    tes_grid = assign_tes_grid()

    # initialize the matrix
    names = 'index,row,col,quadrant,matrix,TES,PIX,ASIC'
    fmts = 'int,int,int,int,a4,int,int,int'
    FPidentity = np.recarray(names=names,formats=fmts,shape=(34*34))

    fp_idx = 0
    for j in range(34):
        row = 33 - j
        for i in range(34):
            col = i
            if row < 17:
                if col < 17:
                    quadrant = 3
                    matrix = 'P87'
                    tes_y = 16 - col
                    tes_x = row
                else:
                    quadrant = 4
                    matrix = 'PXX'
                    tes_y = col - 17
                    tes_x = row
            else:
                if col < 17:
                    quadrant = 2
                    matrix = 'PXX'
                    tes_x = col
                    tes_y = row - 17
                else:
                    quadrant = 1
                    matrix = 'PXX'
                    tes_x = 33 - col
                    tes_y = row - 17
                    
                
            asic_no = tes_grid[tes_x,tes_y].ASIC
            TES_no = tes_grid[tes_x,tes_y].TES
            PIX = tes2pix(TES_no,asic_no)
            rotated_asic = 2*(quadrant-3) + asic_no
            if rotated_asic < 1:
                rotated_asic += 8
            if asic_no==0: rotated_asic = 0

            FPidentity[fp_idx].index = fp_idx
            FPidentity[fp_idx].quadrant = quadrant
            FPidentity[fp_idx].matrix = matrix
            FPidentity[fp_idx].TES = TES_no
            FPidentity[fp_idx].PIX = PIX
            FPidentity[fp_idx].ASIC =  rotated_asic
            FPidentity[fp_idx].row = row
            FPidentity[fp_idx].col = col
            fp_idx += 1
            
    return FPidentity


def tes2index(TES,ASIC):
    '''
    get the unique Focal Plane identifier for a given TES
    '''
    global FPidentity
    if FPidentity is None: FPidentity = make_id_focalplane()

    idx_range = np.where(FPidentity.TES==TES)
    TES_locations = FPidentity[idx_range]

    asic_idx = np.where(TES_locations.ASIC == ASIC)
    entry = TES_locations[asic_idx]
    return entry.index



            
