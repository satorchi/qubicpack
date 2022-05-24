#!/usr/bin/env python3
'''
$Id: mapmaking.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Tue 26 Apr 2022 13:01:09 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

Some basic map making:  

This is mainly for diagnostic purposes.  
Final mapmaking will be more sophisticated and will be implemented in qubicsoft rather than qubicpack

 - make a map from scan data with no modulation
 - plot maps
 - use simple sinusoidal projection for map images
'''
import numpy as np
import healpy as hp
from matplotlib import pyplot as plt
colour = None
def interpolate_azscan(az_axis,azpts,datapts,debug=False):
    '''
    interpolate an azimuth scan onto the given azimuth axis
    '''
    global colour
    if len(azpts)<3:return None


    # for debugging
    if debug==1:
        curve = plt.plot(azpts,datapts,ls='none',marker='^')
        colour = curve[0].get_color()
    elif debug==2:
        curve = plt.plot(azpts,datapts,ls='none',marker='v',color=colour)

    # np.interp only works for monotonically increasing x-axis, so do we have to flip?
    if azpts[-1]<azpts[0]:            
        datapts_interp = np.interp(az_axis,np.flip(azpts),np.flip(datapts))
    else:
        datapts_interp = np.interp(az_axis,azpts,datapts)

    # for debugging
    if debug:
        plt.plot(az_axis,datapts_interp,color=colour,ls='none',marker='.')
    if len(datapts_interp)!=len(az_axis):
        print("didn't get expected number of bins: ",len(datapts_interp))
        return None
    return datapts_interp

def project_azscan(phi_projected,azpts,elevation,datapts):
    '''
    simple projection, cos(elevation)
    '''
    phi = azpts*np.cos(np.radians(elevation))
    datapts_projected = np.interp(phi_projected,phi,datapts)
    unseen_indices = (phi_projected>phi.max()) | (phi_projected<phi.min())
    datapts_projected[unseen_indices] = np.nan
    return datapts_projected

def project_image(phi_projected,azpts,elpts,image):
    '''
    project the map onto the sky using simple cosine projection
    '''
    projected_image = np.empty(image.shape,dtype=float)
    projected_image[:] = np.nan
    
    for elidx,elevation in enumerate(elpts):
        azscan = image[elidx,:]
        projected_azscan = project_azscan(phi_projected,azpts,elevation,azscan)
        projected_image[elidx,:] = projected_azscan
        
    
    return projected_image

def project2healpix(img,az,el,nside=256):
    '''
    generate the Healpix map array for use with healpy tools

    nside = for healpix, to create the molleweide projection
    '''
    m = np.zeros(hp.nside2npix(nside))
    m[:] = hp.UNSEEN
    az_radians = np.radians(az)
    az_npts = az.size
    for elidx,elval in enumerate(el):
        thetaval = np.radians(90-elval)
        theta = np.zeros(az_npts)
        theta[:] = thetaval
        phi = az_radians*np.cos(theta)
        pixel_indexes = hp.ang2pix(nside, theta, phi)

        # we have to average over equal indices
        unique_indexes = np.unique(pixel_indexes)
        for pixidx in unique_indexes:
            img_indexes = pixel_indexes==pixidx
            pix_avg = img[elidx,img_indexes].mean()
            m[pixidx] = pix_avg

    return m
    

def make_map_no_modulation(self,
                           el_stepsize=0.3,
                           az_npts=1000,
                           asic=None,
                           TES=None,
                           backforth=True,
                           tau=None,
                           hk_timeoffset=None,
                           debug=False):
    '''
    make a map from scan data which does not use source modulation

    arguments:
         el_stepsize = how much we moved in elevation after each azimuth scan
                       (maybe we can determine this from the elevation data)
         az_npts = interpolate the azimuth scan onto this number of points
         asic = asic number
         TES = TES number
         backforth = True/False did we go back and forth in azimuth on the same elevation?
         tau = time constant (None if not applying the deconvolution)
         debug = True/False should we make a plot of the azimuth scans?

    '''
    args =self.args_ok(TES,asic,allow_multiple_TES=False)
    if args is None:return
    TES,asic = args
    
    retval = {} # dictionary with all kinds of stuff to return

    # azimuth and elevation
    el_raw = self.elevation()
    az_raw = self.azimuth()

    # get the data and time axes
    hk_timeaxis = self.timeaxis(datatype='hk')
    if hk_timeoffset is not None:
        hk_timeaxis += hk_timeoffset
    sci_timeaxis = self.timeaxis(datatype='sci',asic=asic)
    data_raw = self.timeline(TES=TES,asic=asic)

    # deconvolve the time constant, if requested
    retval['tau'] = tau
    if tau is not None:
        clock_diagnostic = self.timestamp_diagnostic(asic=asic)
        sampling_frequency = clock_diagnostic['mean_samples_per_pps']
        sampling_period = 1/sampling_frequency
        derivative = np.gradient(data_raw, sampling_period)
        ### Well known proxy for exponential deconvolution (JCH)
        signal = data_raw + derivative*tau
    else:
        signal = data_raw
        
    # interpolate science data to housekeeping time axis
    # data = np.interp(hk_timeaxis,sci_timeaxis,signal)

    # interpolate hk data to science data
    data = signal
    az = np.interp(sci_timeaxis,hk_timeaxis,az_raw)
    el = np.interp(sci_timeaxis,hk_timeaxis,el_raw)

    # find the indexes of the elevation which are the same
    el_index_list = []
    el0 = el.max()
    while el0>=el.min():
        el1 = el0 - el_stepsize
        el_index_list.append( (el<=el0) & (el>el1) )
        el0 -= el_stepsize
    el_npts = len(el_index_list)
    el_axis = np.zeros(el_npts)
    az_npts_scan1 = np.zeros(el_npts)
    az_npts_scan2 = np.zeros(el_npts)
    az_npts_fullscan = np.zeros(el_npts)
    az_crossover = np.zeros(el_npts)
    chi2 = np.zeros(el_npts) # difference between up-scan and down-scan per elevation
    

    # azimuth data will be interpolated to az_npts points
    az_min = az.min()
    az_max = az.max()
    az_range = az_max - az_min
    az_axis = az_min + np.arange(az_npts)*az_range/(az_npts-1)


    # make the two maps:  scan-up in azimuth, and scan-down in azimuth
    img_scan1 = np.zeros((el_npts,az_npts))
    img_scan2 = np.zeros((el_npts,az_npts))
    debug1 = False
    debug2 = False
    if debug:
        fig = plt.figure()
        fig.canvas.manager.set_window_title('plt: map debug')
        debug1 = 1
        debug2 = 2

    for idx,el_indexes in enumerate(el_index_list):
        el_axis[idx] = el[el_indexes].mean()
        if el_indexes.sum()<3:continue
        azpts_full = az[el_indexes]
        az_npts_fullscan[idx] = azpts_full.size
        
        az_idx_crossover = np.argmax(azpts_full)
        az_crossover[idx] = az_idx_crossover
        if not backforth: az_idx_crossover = 0

        azpts_scan1 = azpts_full[az_idx_crossover:]
        az_npts_scan1[idx] = azpts_scan1.size
        datapts_interp = interpolate_azscan(az_axis,azpts_scan1,-data[el_indexes][az_idx_crossover:],debug=debug1)
        if datapts_interp is not None:
            img_scan1[idx,:] = datapts_interp
            

        azpts_scan2 = azpts_full[:az_idx_crossover]
        az_npts_scan2[idx] = azpts_scan2.size
        datapts_interp = interpolate_azscan(az_axis,azpts_scan2,-data[el_indexes][:az_idx_crossover],debug=debug2)
        if datapts_interp is not None:
            img_scan2[idx,:] = datapts_interp


            

    # return values
    if azpts_scan1[0]<azpts_scan1[-1]:
        retval['image scan up'] = img_scan1
        retval['image scan down'] = img_scan2
        retval['az npts up'] = az_npts_scan1
        retval['az npts down'] = az_npts_scan2
    else:
        retval['image scan up'] = img_scan2
        retval['image scan down'] = img_scan1
        retval['az npts up'] = az_npts_scan2
        retval['az npts down'] = az_npts_scan1

    retval['az crossover'] = az_crossover
    retval['az npts per scan'] = az_npts_fullscan
    retval['az'] = az_axis
    retval['el'] = el_axis
    calinfo = self.calsource_info()
    if calinfo is None or calinfo['calsource']['status']=='OFF':
        calstr = 'calsource OFF'
    else:
        calstr = '%.3fGHz' % calinfo['calsource']['frequency']
    ttl = '%s %s, ASIC=%i, TES=%i' % (self.obsdate.strftime('%Y-%m-%d-%H:%M:%S'),calstr,asic,TES)
    if tau is not None:
        ttl += ' with time constant $\\tau=%.1f$ms' % (tau*1000)
    retval['title'] = ttl
    retval['obsdate'] = self.obsdate
    retval['asic'] = asic
    retval['TES'] = TES
    retval['calinfo'] = calinfo
    retval['el stepsize'] = el_stepsize
    retval['az npts'] = az_npts
    retval['el npts'] = el_npts
    retval['image combined'] = 0.5*(img_scan1+img_scan2)
    retval['image diff'] = img_scan1 - img_scan2


    # some stuff for the debug plot
    if debug:
        ax = fig.get_axes()[0]
        ax.text(0.5,1.01,ttl,ha='center',va='bottom',transform=ax.transAxes)
        ax.set_xlabel('azimuth / degrees')
            
    return retval

def plot_map(mapinfo,ax=None,plot_image=None,separate_figs=True,projected=True):
    '''
    plot the map using information in the mapinfo dictionary which was returned by make_map_no_modulation()

    The map types are:
    'scan up'      : scanning from negative azimuth to positive azimuth
    'scan down'    : scanning from positive azimuth to negative azimuth
    'scan combined': the linear combination of the up and down scans
    'scan diff'    : the difference map betwen the up and down scans
    'az scans'     : cross cuts at each elevation
 
    '''

    markers = {}
    markers['scan up'] = '^'
    markers['scan down'] = 'v'
    markers['combined scan'] = '+'
    markers['scan diff'] = '_'
    default_marker = '.'
    known_image_list = list(markers.keys()) + ['az scans']
    
    if plot_image is None or not (isinstance(plot_image,str) or isinstance(plot_image,list)):
        image_list = known_image_list
    else:
        if isinstance(plot_image,str):
            image_list = [plot_image]
        else:
            image_list = plot_image

    if ax is None:
        fig_h = plt.rcParams['figure.figsize'][1]
        figsize = (fig_h,fig_h)
        new_fig = True
    else:
        new_fig = False
        separate_figs = False

    az = mapinfo['az']
    el = mapinfo['el']
    azel_limits = (az[0],az[-1],el[0],el[-1])
    ttl = mapinfo['title']
    fileprefix = 'QUBIC_map_%s' % ttl.replace(', ','__').replace(' ','_')

    image = {}
    img_up = mapinfo['image scan up']
    img_dn = mapinfo['image scan down']
    if projected:
        img_up = project_image(mapinfo['az'],mapinfo['az'],mapinfo['el'],img_up)
        img_dn = project_image(mapinfo['az'],mapinfo['az'],mapinfo['el'],img_dn)
    img_combined = 0.5*(img_up + img_dn)
    image['scan up'] = img_up
    image['scan down'] = img_dn
    image['combined scan'] = img_combined
    image['scan diff'] = img_up - img_dn
    n_images = len(image_list)
    if new_fig and not separate_figs:
        fig = plt.figure()
        fig.canvas.manager.set_window_title('plt: map')
        fig.suptitle(ttl)

    errmsg = "Unknown image specification: %%s\n Known plots are: %s" % ', '.join(known_image_list)
    for idx,label in enumerate(image_list):      
        if label not in image.keys() and label!='az scans':
            print(errmsg % label)
            continue

        if label in image.keys():
            img = image[label]
            
        filename = '%s__%s.png' % (fileprefix,label.replace(' ','-'))
        if new_fig:
            if separate_figs:
                if label=='az scans':
                    fig = plt.figure()
                else:
                    fig = plt.figure(figsize=figsize)
                fig.canvas.manager.set_window_title('plt: map')
                ax = fig.add_axes((.1,.1,.88,.8))
            else:
                figwidth = 0.88/n_images
                ax = fig.add_axes((.1+idx*figwidth,.1,figwidth,.8))

        if label in image.keys():
            ax.imshow(img,extent=azel_limits)
        else:
            colour = []
            for elidx in range(el.size): colour.append(None)
            for img_type in image.keys():
                if img_type not in image_list:continue
                img = image[img_type]
                if img_type in markers.keys():
                    marker = markers[img_type]                    
                else:
                    marker = default_marker
                    
                for elidx in range(el.size):
                    if colour[elidx] is None:
                        curve = ax.plot(az,img[elidx,:],marker=marker,ls='none')
                        colour[elidx] = curve[0].get_color()
                    else:
                        curve = ax.plot(az,img[elidx,:],marker=marker,ls='none',color=colour[elidx])
                            
            
        ax.text(0.5,1.01,label,ha='center',va='bottom',transform=ax.transAxes)
        if projected:
            ax.set_xlabel('$\phi$ / degrees')
        else:
            ax.set_xlabel('azimuth / degrees')
        if idx==0 or separate_figs:
            if label=='az scans':
                ax.set_ylabel('ADU / arbitrary units')
            else:
                ax.set_ylabel('elevation / degrees')
        else:
            ax.get_yaxis().set_visible(False)
        if new_fig and separate_figs:
            fig.savefig(filename,format='png',dpi=100,bbox_inches='tight')
    if new_fig and not separate_figs:
        filename = '%s.png' % fileprefix
        fig.savefig(filename,format='png',dpi=100,bbox_inches='tight')
    return
