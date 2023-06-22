'''
$Id: demodulate.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Tue 30 Mar 2021 09:16:20 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

demodulate data
'''
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import fftconvolve

from qubicpack.utilities import figure_window_title, NPIXELS

def renorm(ar):
    return (ar - np.mean(ar)) / np.std(ar)

def get_index_interval(timeaxis,interval):
    '''
    translate time interval into array indexes
    '''
    idx_start = 0
    idx_end = -1
    tstart = interval[0]
    tend = interval[1]
    if timeaxis[0]<tstart:
        idxrange = np.where(timeaxis>=tstart)[0]
        idx_start = idxrange[0]
    if timeaxis[-1]>tend:
        idxrange = np.where(timeaxis>=tend)[0]
        idx_end = idxrange[0]
    
    return (idx_start,idx_end)

def sine_curve_model(x,a,period,timeshift,offset):
    '''
    model for fitting a sine curve
    '''
    y = a*np.sin( 2*np.pi*((x+timeshift)/period) ) + offset 
    
    return y

def fit_sine_curve(xpts,ypts,first_guess=None):
    '''
    fit given points to a sine curve
    '''
    retval = {}
    try:
        popt,pcov = curve_fit(sine_curve_model,xpts,ypts,p0=first_guess)
    except:
        retval['amplitude'] = None
        retval['period'] = None
        retval['timeshift'] = None
        retval['offset'] = None
        retval['popt'] = None
        retval['pcov'] = None
        return retval
    retval['amplitude'] = popt[0]
    retval['period'] = popt[1]
    retval['timeshift'] = popt[2]
    retval['offset'] = popt[3]
    retval['popt'] = popt
    retval['pcov'] = pcov
    return retval
                  
def baseline_by_period(xpts,ypts,period,verbosity=0):
    '''
    subtract a baseline which is calculated by binning on the expected period
    '''
    
    binsize = period
    binstart = xpts.min() + binsize
    nbins = (xpts.max() - xpts.min())/binsize
    nbins = int(nbins)+1
    binedges = binstart + np.arange(nbins)*binsize
    
    idx_bins = np.digitize(xpts,binedges)
    unique_idxbins = np.unique(idx_bins)
    nbins = len(unique_idxbins)

    baseline = np.empty(ypts.shape,dtype=float)

    ndims = len(ypts.shape)
    if ndims==1: # for a single TES
        for idx in range(nbins):
            idxrange = (idx_bins==idx)
            bin_mean = ypts[idxrange].mean()
            baseline[idxrange] = bin_mean
    else: # for a set of TES
        nTES = ypts.shape[0]
        for idx in range(nbins):
            idxrange = (idx_bins==idx)
            bin_mean = ypts[:,idxrange].mean(axis=1)
            for TESidx in range(nTES): # there must be a more pythonic way to do this
                baseline[TESidx,idxrange] = bin_mean[TESidx]
            
    return baseline

def fold_data(xpts,ypts,period,nbins=101,verbosity=0):
    '''
    fold data to the given period
    '''
    x_fold = xpts % period

    # now rebin
    binsize = (x_fold.max() - x_fold.min())/nbins
    binstart = x_fold.min() + binsize
    binedges = binstart + np.arange(nbins)*binsize

    # make sure the last bin edge is bigger than the maximum
    if binedges[-1] < x_fold.max():
        if verbosity>0:print('correcting for rounding error: %.9e < %.9e' % (binedges[-1],x_fold.max()))
        binedges[-1] = x_fold.max() + 1e-3 # make sure it's bigger
    
    idx_bins = np.digitize(x_fold,binedges)
    unique_idxbins = np.unique(idx_bins)

    new_nbins = len(unique_idxbins) # this should be the same!
    if new_nbins!=nbins and verbosity>0:
        print('fold_data: correcting for change of nbins! nbins=%i instead of given nbins=%i' % (new_nbins,nbins))
    x_bin = np.zeros(new_nbins,dtype=float)
    y_bin = np.zeros(new_nbins,dtype=float)
    err_bin = np.zeros(new_nbins,dtype=float)
    for idx,idx_bin in enumerate(unique_idxbins):
        idxrange = np.where(idx_bins==idx_bin)
        x_bin[idx] = x_fold[idxrange].mean()
        y_bin[idx] = ypts[idxrange].mean()
        err_bin[idx] = ypts[idxrange].std()

    return x_bin,y_bin,err_bin

def plot_folded(info,ax=None):
    '''
    make the plot of folded data
    '''
    savefig = False
    if ax is None:
        fig = plt.figure()
        figure_window_title(fig,'folded data')
        ax = fig.add_axes((0.07,0.05,0.88,0.88))
        savefig = True

    data = info['data']
    t_data = info['t_data']
    t0_data = info['t0 data']
    period = info['period']
    data_label = info['data label']
    folded_t_data = info['folded t_data']
    folded_data = info['folded data']
    folded_dataerr = info['folded data error']
    t_src = info['t_src']
    data_src = info['data src']
    t0_src = info['t0 source']
    use_calsource = info['use calsource']
    ax.plot((t_data-t0_data)%period, data,label='folded %s' % data_label,ls='none',marker='x')
    ax.plot(folded_t_data,folded_data,label='folded and averaged %s' % data_label,ls='solid',linewidth=3,color='red')
    ax.errorbar(folded_t_data,folded_data,yerr=folded_dataerr,ls='none',color='red',capthick=2,capsize=2)
    ax.text(0.5,1.0,'folding period = %.6f seconds' % period,ha='center',va='bottom',transform=ax.transAxes)
    ax.legend(loc='upper right',facecolor='wheat',framealpha=0.5)
    if use_calsource:
        folded_t_src = info['folded t_src']
        folded_src = info['folded src']        
    
        axcal = ax.twinx()
        axcal.plot((t_src-t0_data)%period, data_src, label='folded source',color='grey',ls='none',marker='+')
        axcal.plot(folded_t_src,folded_src,label='folded and averaged source',ls='solid',linewidth=3,color='black')
        axcal.legend(loc='lower left',facecolor='wheat',framealpha=0.5)

    if savefig:
        date_str = dt.datetime.utcfromtimestamp(t0_data).strftime('%Y%m%d-%H%M%S')
        if info['baselined']:
            baseline_str = '_baselined'
        else:
            baseline_str = ''
        fname = 'folded%s_%s_%s.png' % (baseline_str,info['data label'].replace(', ','_').replace(' ','-'),date_str)
        fig.savefig(fname,format='png',dpi=100,bbox_inches='tight')
        
    return    

def demodulate(self,
               asic=None,
               TES=None,
               offset=None,
               interval=None,
               calsource=True,
               period=None,
               align_clocks=False,
               timeaxistype='pps',
               remove_baseline=True,
               flip=True,
               doplot=True,
               xwin=True):
    '''
    Interpolating source data to the timestamps of the data

    asic: the asic number (1 to 16)

    TES: the TES number (1 to 128)
         keyword 'all' means average all TES
         an array of size 128 of bool with the TES to average together

    offset: the constant offset between calsource timestamps and data timestamps
    this is of the order of 0.1sec
    
    interval is the time interval in which to calculate the demodulation

    calsource: you can force to not use the calsource

    period: fold at the given period.  Default is to get it from the calsource information

    align_clocks:  you can force the timestamps of the calsource and the data to start at the same time

    flip: calsource is inverted compared to detector response.  Correct this with flip=True
          after 2023-02-01, use flip=True

    doplot: make plots

    xwin:  if making a plot, do not plot to screen if xwin==False
    '''
    args =self.args_ok(TES,asic,allow_multiple_TES=True)
    if args is None:return
    TES,asic = args
    
    if interval is None: # set interval to have no effect
        t0_interval = 0
        tend_interval = float((dt.datetime.utcnow() + dt.timedelta(days=7)).strftime('%s'))
        interval = (t0_interval,tend_interval)
    interval = np.array(interval,dtype=float)

    given_period = period

    errlist = []
    retval = {}
    retval['asic'] = asic
    retval['TES'] = TES
    retval['dataset'] = self.dataset_name
    retval['calsource info'] = self.calsource_info()

    if isinstance(TES,np.ndarray):
        TESstr = 'average of %i selected TES out of %i' % (TES.sum(),TES.size)
        TESfilenamestr = 'average-%i-TES' % TES.sum()
        if asic is None:
            ASICstr = 'all ASICs'
            ASICfilenamestr = 'all-ASICs'
            t_data_orig,adu = self.tod()
            data = adu[TES].mean(axis=0)
        else:            
            ASICstr = 'ASIC %i' % asic            
            ASICfilenamestr = 'ASIC%03i'
            adu = self.timeline_array(asic=asic)
            t_data_orig = self.timeaxis(datatype='sci',asic=asic,axistype=timeaxistype)
            data = adu[TES].mean(axis=0)
    else:
        TESstr = 'TES %i' % TES
        TESfilenamestr = 'TES%03i' % TES
        ASICstr = 'ASIC %i' % asic                    
        ASICfilenamestr = 'ASIC%03i' % asic                    
        data = self.timeline(asic=asic,TES=TES)
        if data is None: return retval
        t_data_orig = self.timeaxis(datatype='sci',asic=asic,axistype=timeaxistype)


    if flip: data = -data
    data_label = '%s, %s' % (ASICstr,TESstr)
    data_label_filenamestr = '%s_%s' % (ASICfilenamestr,TESfilenamestr)
    retval['data label'] = data_label
        
    t_data = t_data_orig.copy()
    t0_data = t_data[0]
    if t0_data > 1494453600: # if the timeaxis is made from the sampling time and starts at zero
        do_tzone_correction = True
    else:
        do_tzone_correction = False
    
    t_src_orig,v_src = self.calsource()
    if not calsource or v_src is None:
        use_calsource = False
        msg = 'No calsource.'
        self.printmsg(msg,verbosity=3)
        errlist.append(msg)
        t_src = t_data
        data_src = data
    else:
        # shift source to oscillate around zero
        use_calsource = True
        t_src = t_src_orig.copy()
        if flip:
            # source sampling inverted compared to data
            data_src = -v_src + v_src.mean()
        else:
            data_src =  v_src - v_src.mean()

    t0_src = t_src[0] # this will be modified after the offset is calculated below
    retval['use calsource'] = use_calsource
    
    # for some data, the source is in UT while the detector data is in localtime
    tz_offset = 0.0
    tz_count = 0
    if do_tzone_correction:
        delta_tz = np.abs(t0_data - t0_src)
        while delta_tz>3599:
            tz_offset += 3600
            delta_tz -= 3600
            tz_count += 1
        if (t0_data - t0_src) < 0:
            tz_offset = -tz_offset
    retval['time zone offset'] = tz_offset
    retval['time zone offset hours'] = tz_count
    self.printmsg('time zone offset %f' % tz_offset,verbosity=3)
    
    # put data in UT (assume calsource was in UT)
    t_data -= tz_offset
    t0_data = t_data[0]
    interval -= tz_offset

    # if forcing the alignment of the clocks, make t_src start at the same time as t_data
    if align_clocks:
        self.printmsg('forcing alignment of the clocks',verbosity=3)
        align_offset = t_src[0] - t_data[0]
        t_src -= align_offset
        t0_src = t_src[0]
        self.printmsg('t0_src=%f' % t0_src,verbosity=3)
    else:
        align_offset = 0
    retval['clock forced realignment offset'] = align_offset

    # truncate the data timeline to overlapping time or to the given interval
    t0_list = [t0_data,t0_src,interval[0]]
    self.printmsg('t0_list = %s' % t0_list,verbosity=3)
    tend_list = [t_data[-1],t_src[-1],interval[1]]
    self.printmsg('tend_list = %s' % tend_list,verbosity=3)
    tstart = max(t0_list)
    tend = min(tend_list)
    self.printmsg('tstart,tend = %.6f, %.6f' % (tstart,tend),verbosity=3)
    if tend<tstart:
        self.printmsg('Error!  The given interval does not have overlapping times with calsource and data.')
        self.printmsg('        Maybe try with option:  align_clocks=True')
        return retval
    idx_start,idx_end = get_index_interval(t_data,(tstart,tend))
    retval['data time start index'] = idx_start
    retval['data time end index'] = idx_end
    self.printmsg('data: idx_start, idx_end = %i, %i' % (idx_start,idx_end),verbosity=3)
    t_data = t_data[idx_start:idx_end]
    data = data[idx_start:idx_end]
    t0_data = t_data[0]
    retval['t0 data'] = t0_data
    t0_str = dt.datetime.utcfromtimestamp(t0_data).strftime('%Y-%m-%d %H:%M:%S.%f')
    t0_filenamestr = dt.datetime.utcfromtimestamp(t0_data).strftime('%Y%m%d-%H%M%S')

    # truncate the source timeline to overlapping time or to the given interval
    idxsrc_start,idxsrc_end = get_index_interval(t_src,(tstart,tend))
    self.printmsg('src: idx_start, idx_end = %i, %i' % (idxsrc_start,idxsrc_end),verbosity=3)
    t_src = t_src[idxsrc_start:idxsrc_end]
    data_src = data_src[idxsrc_start:idxsrc_end]
    
    # remove slow baseline if requested
    retval['baselined'] = remove_baseline
    if remove_baseline:
        baseline = baseline_by_period(t_data,data,period,verbosity=self.verbosity)
        data -= baseline
        baseline_str = '_baseline_removed'
    else:
        baseline_str = ''

    
    # fit data to sine curve
    amplitude = 0.5*(data.max() - data.min())
    calinfo = self.calsource_info()
    if given_period is None:
        if calinfo is not None and 'frequency' in calinfo['modulator'].keys() and calinfo['modulator']['frequency']!='none':
            period = 1.0/calinfo['modulator']['frequency']
        else:
            period = 1.0
    
    timeshift = 0.0
    amplitude_offset = data.mean()
    first_guess = (amplitude,period,timeshift,amplitude_offset)
    data_fit = fit_sine_curve(t_data-t0_data,data,first_guess=first_guess)
    retval['data fit'] = data_fit
    model_data = sine_curve_model(t_data-t0_data,
                                  data_fit['amplitude'],
                                  data_fit['period'],
                                  data_fit['timeshift'],
                                  data_fit['offset'])
    # find the first time of max (when sin(2pi((x+timeshift)/period) = 1)
    t_data_max = data_fit['period']/4 - data_fit['timeshift']
    # sometimes, curve_fit returns a negative amplitude, which is equivalent to a half-period phase shift
    if data_fit['amplitude']<0:
        t_data_max += 0.5*data_fit['period']
    retval['t_data max'] = t_data_max
    

    # fit source to sine curve
    if use_calsource:
        # period is either given as an argument, or taken from calsource_info (see above)
        amplitude = 0.5*(data_src.max() - data_src.min())
        timeshift = 0.0
        amplitude_offset = 0.0
        first_guess = (amplitude,period,timeshift,amplitude_offset)
        source_fit = fit_sine_curve(t_src-t0_data,data_src,first_guess=first_guess)
        # find the first time of max (when sin(2pi((x+timeshift)/period) = 1)
        t_src_max = source_fit['period']/4 - source_fit['timeshift']
        # sometimes, curve_fit returns a negative amplitude, which is equivalent to a half-period phase shift
        if source_fit['amplitude']<0:
            t_src_max += 0.5*source_fit['period']
    else:
        source_fit = data_fit
        t_src_max = t_data_max
        data_src = data
        t_src = t_data

    retval['source fit'] = source_fit
    retval['t_src max'] = t_src_max
    if given_period is None:
        period = source_fit['period']
    else:
        period = given_period
    retval['period'] = period # this is the period we use for demodulation

    # find the constant timestamp offset between source and data
    # for more info: http://qubic.in2p3.fr/wiki/pmwiki.php/TD/Demodulation
    if offset is None:
        offset = t_data_max - t_src_max
    retval['source-data time offset'] = offset
    self.printmsg('source-data time offset %f' % offset,verbosity=3)
    
    # adjust the source timeline by the constant offset
    # go back to the original source data before applying the offset
    # and then find the appropriate interval
    if use_calsource:
        t_src_orig,v_src = self.calsource()
        t_src = t_src_orig.copy()
        data_src = -v_src + v_src.mean()

        t_src += (offset-align_offset)

        t0_src = t_src[0]
        # truncate the source timeline to overlapping time or to the given interval
        idxsrc_start,idxsrc_end = get_index_interval(t_src,(tstart,tend))
        retval['source time start index'] = idxsrc_start
        retval['source time end index'] = idxsrc_end
        self.printmsg('src: idx_start, idx_end = %i, %i' % (idxsrc_start,idxsrc_end),verbosity=3)
        t_src = t_src[idxsrc_start:idxsrc_end]
        data_src = data_src[idxsrc_start:idxsrc_end]
        model_src = sine_curve_model(t_src-t0_data,
                                     source_fit['amplitude'],
                                     source_fit['period'],
                                     source_fit['timeshift']-offset,
                                     source_fit['offset'])

    # some values to return
    npts = len(data)
    retval['t_data'] = t_data
    retval['data'] = data
    retval['n data points'] = npts
    retval['interval'] = (t_data[0],t_data[-1])

    retval['t0 source'] =  t0_src
    retval['t_src'] = t_src
    retval['data src'] = data_src

    ### Interpolating source data to the timestamps of the data
    ### and making the product of the detector and the source
    ### the source is renormalized, but the data is only shifted to oscillate around 0
    if use_calsource:
        self.printmsg('demodulate: number of points for t_data, t_src, data_src: %i, %i, %i' % (len(t_data),len(t_src),len(data_src)),verbosity=3)
        data_src_interp = np.interp(t_data-t0_data, t_src-t0_data, data_src)
    else:
        data_src_interp = data
    product = renorm(data_src_interp) * (data - data.mean())
    retval['calsource interpolated to data time axis'] = data_src_interp

    # demodulate, smoothed over a period
    freq_sampling = 1. / ( (t_data[-1] - t_data[0])/npts )
    size_period = int(freq_sampling * period) + 1
    filter_period = np.ones(size_period) / size_period
    try:
        demodulated = fftconvolve(product, filter_period, mode='same')
    except:
        msg = 'ERROR calculating demodulation'
        errlist.append(msg)
        self.printmsg(msg,verbosity=3)
        demodulated = np.zeros(size_period)
        
    retval['freq_sampling'] = freq_sampling
    retval['size_period'] = size_period
    retval['filter_period'] = filter_period
    retval['demodulated'] =  demodulated
    retval['n demodulated points'] = len(demodulated)

    # rebin the result by period bins
    period_index = ((t_data - t_data[0]) / period).astype(int)
    allperiods = np.unique(period_index)
    binned_t = np.zeros(len(allperiods))
    binned_demodulated = np.zeros((len(allperiods)))
    binned_stdev = np.zeros((len(allperiods)))
    binned_npts = np.zeros((len(allperiods)))
    for idx,p_index in enumerate(allperiods):
        idxrange = period_index == p_index
        binned_t[idx] = np.mean(t_data[idxrange] - t_data[0])
        binned_demodulated[idx] = np.mean(demodulated[idxrange])
        binned_stdev[idx] = np.std(demodulated[idxrange]) # / np.sqrt(idxrange.sum()) really?
        binned_npts[idx] = idxrange.sum()
    retval['binned t'] = binned_t
    retval['binned demodulated'] = binned_demodulated
    retval['binned stdev'] = binned_stdev
    retval['demodulated signal'] = binned_demodulated[1:-1].mean()
    retval['demodulated signal error']  = binned_stdev[1:-1].mean()
    retval['npts per bin'] = binned_npts

    # fold the data and the calsource
    folded_t_data,folded_data,folded_dataerr = fold_data(t_data-t0_data,data,period,verbosity=self.verbosity)
    retval['folded t_data'] = folded_t_data
    retval['folded data'] = folded_data
    retval['folded data error'] = folded_dataerr
    if use_calsource:
        folded_t_src, folded_src, folded_srcerr  = fold_data(t_src-t0_data,data_src,period,verbosity=self.verbosity)
        retval['folded t_src'] = folded_t_src
        retval['folded src'] = folded_src
        retval['folded src error'] = folded_srcerr

    # if not plotting, exit now
    if not doplot:
        retval['axes'] = None
        retval['fig'] = None
        return retval

    if xwin: plt.ion()
    else: plt.ioff()
    fig = plt.figure()
    retval['fig'] = fig
    figure_window_title(fig,'data/calsource interpolation')
    axes = []

    # first column of plots
    axes.append(fig.add_axes((0.03,0.68,0.62,0.3)))
    axes[-1].plot((t_data-t0_data), renorm(data),label=data_label,ls='none',marker='x')
    axes[-1].plot((t_data-t0_data), renorm(model_data),label='data sine model',color='orange')
    axes[-1].plot([t_data_max,t_data_max],[-2,2],color='orange',linewidth=3,label='first data peak')
    if use_calsource:
        axes[-1].plot((t_src-t0_data), renorm(data_src), label='source',ls='none',marker='+')
        axes[-1].plot((t_src-t0_data), renorm(model_src),label='source sine model',color='red')
        axes[-1].plot([t_src_max,t_src_max],[-2,2],color='red',linewidth=3,label='first source peak')
    axes[-1].legend()
    axes[-1].text(0.5,1.01,self.infotext(),ha='center',va='bottom',transform=axes[-1].transAxes)
    axes[-1].set_xlim(0,t_data[-1]-t0_data)
    axes[-1].tick_params(labelbottom=False)
    
    axes.append(fig.add_axes((0.03,0.35,0.62,0.3)))
    axes[-1].plot((t_data-t0_data), renorm(data), label=data_label,marker='v',ls='none')
    if use_calsource:
        axes[-1].plot((t_data-t0_data), renorm(data_src_interp), label='SRC (interp on data)',marker='^',ls='none')
    axes[-1].set_xlim(0,t_data[-1]-t0_data)
    axes[-1].legend()
    axes[-1].tick_params(labelbottom=False)
        
    axes.append(fig.add_axes((0.03,0.03,0.62,0.3)))
    axes[-1].plot((t_data-t0_data), product, label='data $\\times$ source',marker='o',markersize=1,ls='none')
    axes[-1].legend()        
    axes[-1].set_xlabel('date / seconds since %s' % t0_str)
    axes[-1].set_xlim(0,t_data[-1]-t0_data)

    # second column of plots

    # folded
    axes.append(fig.add_axes((0.70,0.68,0.28,0.3)))
    plot_folded(retval,ax=axes[-1])

    # demodulated
    axes.append(fig.add_axes((0.70,0.35,0.28,0.3)))
    axes[-1].plot((t_data-t0_data),demodulated, label='demodulated',ls='solid')
    axes[-1].errorbar(binned_t,binned_demodulated, yerr=binned_stdev,label='binned',
                      marker='x',markersize=3,ls='none',color='red',capsize=3,capthick=3)
    axes[-1].legend(loc='lower left')
    axes[-1].tick_params(labelbottom=False)
    # try:
    #     expo = int(np.log10(retval['demodulated signal']))
    #     val = retval['demodulated signal']/10**expo
    #     sig_str = '$%.6f\\times10^{%i}$' % (val,expo)
    #     expo = int(np.log10(retval['demodulated signal error']))
    #     val = retval['demodulated signal error']/10**expo
    #     err_str = '$%.2f\\times10^{%i}$' % (val,expo)
    # except:
    sig_str = '%.6f' % retval['demodulated signal']
    err_str = '%.2f' % retval['demodulated signal error']
    txt = 'demodulated signal = %s $\\pm$ %s' % (sig_str,err_str)
    if len(errlist)>0:
        txt = '\n'.join([txt]+errlist)               
    boxprops = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    axes[-1].text(0.5,0.05,txt,
                  ha='center',va='bottom',transform=axes[-1].transAxes,bbox=boxprops)
    xlim = axes[-1].axis()[:2]
    
    # number of points per bin
    axes.append(fig.add_axes((0.70,0.03,0.28,0.3)))
    axes[-1].bar(binned_t, binned_npts, align='center')
    axes[-1].set_xlim(xlim)
    axes[-1].set_xlabel('period / seconds')
    axes[-1].set_ylabel('number of points per bin')

    
    pngname = 'demodulation_diagnostic%s_%s_%s.png' % (baseline_str,data_label_filenamestr,t0_filenamestr)    
    fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    retval['pngname'] = pngname
        
    retval['axes'] = axes
    if not xwin:
        retval['fig'] = None
        retval['axes'] = None
        plt.close(fig)
    return retval
