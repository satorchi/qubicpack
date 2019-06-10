'''
$Id: plot_fp.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Fri 24 May 2019 17:57:04 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

plot curves and background colour on the QUBIC focal plane
this replaces plot_physical_layout.py
'''
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt

from qubicpack.utilities import NPIXELS, NASIC, ASIC_index, FIGSIZE
from qubicpack.pix2tes import assign_pix_grid, assign_tes_grid

def mylut(v,vmin=3.0,vmax=9.0):
    vfractional=(v-vmin)/(vmax-vmin)
    colourmap = plt.cm.get_cmap('Spectral_r')
    rgb=colourmap(vfractional)
    return rgb

def plot_fp(args):
    '''
    plot curves and background colour on the QUBIC focal plane
    
    the argument is a dictionary with the curves and options
    valid keywords (n is the ASIC number):
      'ASIC<n>' : the array of NPIXELS curves
      'ASIC<n> x-axis' : x-axis (use this to plot bias curves)
      'ASIC<n> bg' : the value to determine the background colour
      'ASIC<n> good' : an array of bool for each TES (good or not)
      'obsdate' : observation date (datetime object)
      'title' : plot title
      'subtitle' : plot subtitle
      'lutmin' : min value for colour look up table
      'lutmax' : max value for colour look up table
      'pngname' : name of file for saving the plot

    '''
    
    # initialize stuff
    pix_grid = assign_pix_grid()
    nrows = pix_grid.shape[0]
    ncols = pix_grid.shape[1]

    if 'figsize' in args.keys():
        figsize = args['figsize']
    else:
        figsize = (FIGSIZE[0],FIGSIZE[0])
    fontsize=figsize[0]
    ttlfontsize=figsize[0]*1.2
        
    obsdate = None
    if 'obsdate' in args.keys():
        obsdate = args['obsdate']
    tes_grid = assign_tes_grid()

    if 'pngname' in args.keys():
        pngname = args['pngname']
    elif obsdate:
        pngname = 'QUBIC_focal_plane_%s.png' % obsdate.strftime('%Y%m%dT%H%M%S')
    else:
        pngname = 'QUBIC_focal_plane.png'


    if 'title' in args.keys():
        ttl = args['title']
    else:
        ttl = 'QUBIC TES array'

    subttl = None
    if 'subtitle' in args.keys():
        subttl = args['subtitle']

    lutmin = 0.0
    if 'lutmin' in args.keys():
        lutmin = args['lutmin']
    lutmax = 10.0
    if 'lutmax' in args.keys():
        lutmax = args['lutmax']

    label_boxprops = dict(boxstyle='round, pad=0.1', facecolor='white', alpha=1.0)
 
    face_colours = {}
    face_colours['ASIC0'] = 'black'
    face_colours['ASIC1'] = 'white'
    face_colours['ASIC2'] = 'white'

    curve_colours = {}
    curve_colours['ASIC0'] = 'black'
    curve_colours['ASIC1'] = 'black'
    curve_colours['ASIC2'] = 'blue'

    label_colours = {}
    label_colours['ASIC0'] = 'black'
    label_colours['ASIC1'] = 'black'
    label_colours['ASIC2'] = 'blue'
    
    plt.ion()
    fig,ax=plt.subplots(nrows,ncols,figsize=figsize)
    fig.text(0.5,0.985,ttl,ha='center',fontsize=ttlfontsize)
    fig.canvas.set_window_title('plt:  %s' % ttl)
    fig.suptitle(subttl,fontsize=ttlfontsize)
    

    for row in range(nrows):
        for col in range(ncols):
            ax[row,col].get_xaxis().set_visible(False)
            ax[row,col].get_yaxis().set_visible(False)

            # the pixel identity associated with its physical location in the array
            TES_parts = np.modf( tes_grid[row,col] )
            asic_str = '%.0f' % (10*TES_parts[0])
            TES_str = '%.0f' % TES_parts[1]
            TES = int(TES_str)
            TES_idx = TES - 1
            pix_label = TES_str
            asic_key = 'ASIC%s' % asic_str
            asicbg_key = '%s bg' % asic_key
            asicgood_key = '%s good' % asic_key
            xaxis_key = '%s x-axis' % asic_key
            
            face_colour = face_colours[asic_key]
            label_colour = label_colours[asic_key]
            curve_colour = curve_colours[asic_key]

            text_y=0.4
            text_x=0.9
            labelfontsize = ttlfontsize


            if asic_key in args.keys() and args[asic_key] is not None:
                curve_x = range(args[asic_key].shape[1])
                curve = args[asic_key][TES_idx]
                text_x = 0.9
                text_y = 0.1
                labelfontsize = fontsize
                if asicgood_key in args.keys() and not args[asicgood_key][TES_idx]:
                    face_colour='black'
                    label_colour='white'
                    curve_colour='white'
                elif asicbg_key in args.keys():
                    face_colour=mylut(args[asicbg_key][TES_idx],lutmin,lutmax)

                if xaxis_key in args.keys() and args[xaxis_key] is not None:
                    curve_x = args[xaxis_key]

                ax[row,col].plot(curve_x,curve,color=curve_colour)

            
                

            #print('(%i,%i) : facecolour=%s, labelcolour=%s' % (row,col,face_colour,label_colour))
            ax[row,col].set_facecolor(face_colour)
            label_boxprops['facecolor'] = face_colour
            ax[row,col].text(text_x,text_y,pix_label,
                             va='bottom',ha='right',color=label_colour,fontsize=labelfontsize,
                             bbox=label_boxprops,transform = ax[row,col].transAxes)
            
    plt.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    plt.show()

    return
