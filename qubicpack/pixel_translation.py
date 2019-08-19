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

from qubicpack.pix2tes import assign_tes_grid

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


def plot_id_focalplane(fpmatrix):
    '''
    plot all the different identity names of each pixel in the focal plane

    fpmatrix is a 34x34 matrix of dictionaries
    '''

    quadrant_colour = ['blue','red','green','purple']
    asic_colour = ['blue','darkblue','red','#cc0000','green','#00cc00','purple','#7210a7']

    fig = plt.figure(figsize=(20,20))
    fig.canvas.set_window_title('plt: QUBIC Focal Plane ID Matrix')
    ax = fig.add_axes([0,0,1,1])
    ax.text(0.5,0.96,'QUBIC Focal Plane ID Matrix', ha='center',va='bottom',transform=ax.transAxes,fontsize=20)
    ax.set_xlim(-1,35)
    ax.set_ylim(-1,35)
    ax.set_aspect('equal')
    
    
    for i in range(34):
        for j in range(34):
            txt = 'Q%i' % (fpmatrix[i][j]['quadrant'])
            quadrant = fpmatrix[i][j]['quadrant']
            asic = fpmatrix[i][j]['ASIC']
            colour = asic_colour[asic-1]
            if fpmatrix[i][j]['TES']==0:
                colour = 'black'
                txt += '\nFP%4i' % fpmatrix[i][j]['FP index']
            else:
                txt += ' %s\nFP%4i\nASIC%i\nTES%03i' % (fpmatrix[i][j]['matrix'],fpmatrix[i][j]['FP index'],fpmatrix[i][j]['ASIC'],fpmatrix[i][j]['TES'])
            plot_square(i,j,colour=colour,labelcolour='white',label=txt,fontsize=6)
    return

def make_id_focalplane():
    '''
    make the matrix which has all the translation info for pixel identities
    '''

    tes_grid = assign_tes_grid()

    # initialize the matrix
    fpmatrix = []
    for i in range(34):
        fpmatrix.append([])
        fpmatrix[i] = []
        for j in range(34):
            transdict = {}
            transdict['FP index'] = -1
            fpmatrix[i].append(transdict)

    # there should be a way to do this without a second double loop
    fp_idx = 0
    for j in range(34):
        row = 33 - j
        for i in range(34):
            col = i
            fpmatrix[col][row]['FP index'] = fp_idx
            fpmatrix[col][row]['PIX'] = -1
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
                    
                
            TES_parts = np.modf( tes_grid[tes_x,tes_y] )
            asic_no = int(round(10*TES_parts[0]))
            TES_no = int(round(TES_parts[1]))

            fpmatrix[col][row]['quadrant'] = quadrant
            fpmatrix[col][row]['matrix'] = matrix
            fpmatrix[col][row]['TES'] = TES_no
            if quadrant==3 or quadrant==4:
                rotated_asic = 2*(quadrant-3) + asic_no  # quadrant 3 & 4
            else:
                rotated_asic = 2*(quadrant+1) + asic_no  # quadrant 1 & 2                
            
            fpmatrix[col][row]['ASIC'] =  rotated_asic
                    
            fp_idx += 1
            
    return fpmatrix

            
