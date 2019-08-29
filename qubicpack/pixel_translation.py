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

import os
import datetime as dt
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits

import qubic
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

            
def extnames(hdulist):
    '''
    list the extension names of the HDU in a list of HDUs read from a FITS file
    '''
    for idx in range(len(hdulist)):
        if 'EXTNAME' in hdulist[idx].header.keys():
            print('%i:  %s' % (idx,hdulist[idx].header['EXTNAME']))

    return

def plot_fits_layout(filename):
    '''
    plot the QUBIC focal plane detector layout as found in the CalQubic_DetArray_...fits file
    '''
    basename = os.path.basename(filename)
    ttl = 'FITS Layout for %s' % basename
    
    hdulist = fits.open(filename)
    
    hdu_centre = hdulist[1]
    hdu_vertices = hdulist[2]
    hdu_removed = hdulist[3]
    hdu_index = hdulist[4]
    hdu_quadrant = hdulist[5]

    xmax = np.nanmax(hdu_vertices.data[:,:,:,0])
    xmin = np.nanmin(hdu_vertices.data[:,:,:,0])
    xspan = xmax - xmin
    plt_xmin = xmin - 0.1*xspan
    plt_xmax = xmax + 0.1*xspan
    
    ymax = np.nanmax(hdu_vertices.data[:,:,:,1])
    ymin = np.nanmin(hdu_vertices.data[:,:,:,1])
    yspan = ymax - ymin
    plt_ymin = ymin - 0.1*yspan
    plt_ymax = ymax + 0.1*yspan

    fig = plt.figure(figsize=(20,20))
    fig.canvas.set_window_title('plt: %s' % basename)
    ax = fig.add_axes([0,0,1,1])
    ax.set_aspect('equal')
    ax.set_xlim([plt_xmin,plt_xmax])
    ax.set_ylim([plt_ymin,plt_ymax])
    ax.text(0.5,1.0,ttl,va='top',ha='center',fontsize=20,transform=ax.transAxes)

    fp_idx = -1
    for i in range(34):
        for j in range(34):
            #fp_idx += 1
            fp_idx = 34*i + j
            quadrant_idx = hdu_quadrant.data[i,j]
            if quadrant_idx > 3: continue
            removed = hdu_removed.data[i,j]
            if removed==1: continue
            
            quadrant = quadrant_idx + 1

            fpindex = hdu_index.data[i,j]
            
            x = hdu_centre.data[i,j,0]
            y = hdu_centre.data[i,j,1]

            
            xpts = np.append(hdu_vertices.data[i,j,:,0],hdu_vertices.data[i,j,0,0])
            ypts = np.append(hdu_vertices.data[i,j,:,1],hdu_vertices.data[i,j,0,1])
            plt.fill(xpts,ypts,color=quadrant_colour[quadrant_idx])

            lbl = 'Q%i\n%04i\n%04i' % (quadrant,fpindex,fp_idx)
            ax.text(x,y,lbl,fontsize=6,color='white',ha='center',va='center')
            #ax.plot(x,y,ls='none',marker='s',markersize=10,color='black',fillstyle='none')


    return hdulist


def plot_instrument_layout(q):
    '''
    plot the layout of the focal plane according to the qubicsoft simulation software
    q is a QubicInstrument object 
    '''
    basename = os.path.basename(q.calibration.detarray)
    config = q.config
    ttl = 'Layout for %s Instrument' % config
    
    fig = plt.figure(figsize=(20,20))
    fig.canvas.set_window_title('plt: %s' % ttl)
    ax = fig.add_axes([0,0,1,1])
    ax.set_aspect('equal')
    q.plot()
    for idx,pos in enumerate(q.detector.center[:,0:2]):
        quadrant_idx = q.detector.quadrant[idx]
        quadrant = quadrant_idx + 1
        fpindex = q.detector.index[idx]
        lbl = 'Q%i\n%i' % (quadrant,fpindex)
        plt.text(pos[0],pos[1],lbl,va='center',ha='center',fontsize=12,color=quadrant_colour[quadrant_idx])

    ax.text(0.5,0.99,ttl,va='top',ha='center',fontsize=20,transform=ax.transAxes)

    return


        
def make_qubicsoft_detarray_fits(config='TD',detector=None,fp=None):
    '''
    make the file CalQubic_DetArray_TD.fits which is used in qubicsoft simulations
    fp is a qubicfp object which has a list of is_good based on I-V measurements
    '''
    global FPidentity
    if FPidentity is None: FPidentity = make_id_focalplane()

    
    
    # the primary header
    prihdr=fits.Header()
    prihdr['TELESCOP'] = ('QUBIC','Q & U Bolometric Interferometer for Cosmology')
    prihdr['AUTHOR'] = ('qubicpack','https://github.com/satorchi/qubicpack')
    prihdr['DATE'] = (dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC'),'date this file was created')
    prihdr['DESCRPN'] = ('QUBIC bolometer layout','description of this file')
    prihdr['format version'] = ('3.0','version of the organization of data')
    prihdr['COMMENT'] = 'Quadrant layout from Stefanos Marnieros and Claude Chapron'
    prihdr['COMMENT'] = ' In a quadrant, the spacing between detectors is 3000um'                 
    prihdr['COMMENT'] = ' The collecting surface of a bolometer is 2700umx2700um'                 
    prihdr['COMMENT'] = ' The extra gap between quadrants is 1800um'
    prihdu = fits.PrimaryHDU(header=prihdr)

    # there are 6 more HDUs
    # 1:  center
    # 2:  corner
    # 3:  removed
    # 4:  index
    # 5:  quadrant
    # 6:  efficiency

    gap = 0.0018 
    bol_length = 0.0027
    bol_sep = 0.003
    leftmost = -(16*bol_sep + (gap+bol_sep)/2)
    topmost = -leftmost
    
    centre = np.empty((34,34,2),dtype=np.float64)
    centre.fill(np.nan)

    corner = np.empty((34,34,4,2),dtype=np.float64)
    corner.fill(np.nan)
    
    removed = np.zeros((34,34),dtype=np.uint8)

    fpindex = -np.ones((34,34),dtype=np.int64)

    quadrant = np.empty((34,34),dtype=np.uint8)
    quadrant.fill(255)

    efficiency = np.empty((34,34),dtype=np.float64)
    efficiency.fill(0.8)    
    
    xgap = 0.0
    ygap = 0.0
    quad_choice = [0,1]
    fp_idx = 0
    for i in range(34):
        if i>16:
            ygap = gap
            quad_choice = [2,3]
        for j in range(34):
            quad = quad_choice[0]
            xgap = 0.0
            if j>16:
                xgap = gap
                quad = quad_choice[1]

            quadrant[i,j] = FPidentity[fp_idx].quadrant - 1 # quadrant 0-3 in the FITS file
            if FPidentity[fp_idx].PIX == -1:
                quadrant[i,j] = 255

            # for the Technical Demonstrator, we only have a quarter focal plane
            # and it's in quadrant 3 (quadrant_index = 2)
            if config=='TD' and quadrant[i,j]!=2:
                removed[i,j] = 1
            if fp is not None:
                asic = FPidentity[fp_idx].ASIC
                asic_idx = asic - 1
                TES = FPidentity[fp_idx].TES                
                
                if TES>0\
                   and asic<=len(fp.asic_list)\
                   and fp.asic_list[asic_idx] is not None\
                   and fp.asic_list[asic_idx].is_good_iv(TES):
                    removed[i,j] = 0
                    print('added: TES=%3i, asic=%i' % (TES,asic))
                else:
                    removed[i,j] = 1

                    
            x = leftmost + bol_sep*j + xgap
            y = topmost - bol_sep*i - ygap
            centre[i,j,:] = [x,y]

            vx = [x - 0.5*bol_length, x + 0.5*bol_length]
            vy = [y - 0.5*bol_length, y + 0.5*bol_length]
            corner[i,j,0,:] = [vx[0],vy[0]]
            corner[i,j,1,:] = [vx[0],vy[1]]
            corner[i,j,2,:] = [vx[1],vy[1]]
            corner[i,j,3,:] = [vx[1],vy[0]]

            #fp_idx = 34*i + j
            fpindex[i,j] = fp_idx
            fp_idx += 1


    cenhdr = fits.Header()
    cenhdr['EXTNAME'] = ('center','centre positions of the bolometers')
    cenhdr['BUNIT'] = ('m','dimensions in meters')                   
    cenhdu = fits.ImageHDU(header=cenhdr,data=centre)

    cornerhdr = fits.Header()
    cornerhdr['EXTNAME'] = ('corner','vertex positions of the bolometers')
    cornerhdr['BUNIT'] = ('m','dimensions in meters')                   
    cornerhdu = fits.ImageHDU(header=cornerhdr,data=corner)

    removedhdr = fits.Header()
    removedhdr['EXTNAME'] = ('removed','bolometers which are not used')
    removedhdr['BUNIT'] = ('','dimensionless')                   
    removedhdu = fits.ImageHDU(header=removedhdr,data=removed)

    indexhdr = fits.Header()
    indexhdr['EXTNAME'] = ('index','unique index for each location on the grid')
    indexhdr['BUNIT'] = ('','dimensionless')                   
    indexhdu = fits.ImageHDU(header=indexhdr,data=fpindex)

    quadranthdr = fits.Header()
    quadranthdr['EXTNAME'] = ('quadrant','0 to 3, counter-clockwise from top right')
    quadranthdr['BUNIT'] = ('','dimensionless')                   
    quadranthdu = fits.ImageHDU(header=quadranthdr,data=quadrant)

    efficiencyhdr = fits.Header()
    efficiencyhdr['EXTNAME'] = ('efficiency','bolometric efficiency')
    efficiencyhdr['BUNIT'] = ('','dimensionless')                   
    efficiencyhdu = fits.ImageHDU(header=efficiencyhdr,data=efficiency)
    
    
    # finally, make the whole FITS file
    hdulist = [prihdu,cenhdu,cornerhdu,removedhdu,indexhdu,quadranthdu,efficiencyhdu]
    thdulist = fits.HDUList(hdulist)

    if fp is not None:
        detector = fp.detector_name
    if detector is not None:
        filename = 'CalQubic_DetArray_%s_%s.fits' % (detector,config)
    else:
        filename = 'CalQubic_DetArray_%s.fits' % config
    thdulist.writeto(filename,overwrite=True)
    thdulist.close()

    
    return

