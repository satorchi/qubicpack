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

def renorm(ar):
    return (ar - np.mean(ar)) / np.std(ar)

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
    popt,pcov = curve_fit(sine_curve_model,xpts,ypts,p0=first_guess)
    retval['amplitude'] = popt[0]
    retval['period'] = popt[1]
    retval['timeshift'] = popt[2]
    retval['offset'] = popt[3]
    retval['popt'] = popt
    retval['pcov'] = pcov
    return retval
                  

def fold_data(xpts,ypts,period,nbins=100):
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
        print('correcting for rounding error: %.9e < %.9e' % (binedges[-1],x_fold.max()))
        binedges[-1] = x_fold.max() + 1e-3 # make sure it's bigger
    
    idx_bins = np.digitize(x_fold,binedges)
    unique_idxbins = np.unique(idx_bins)

    new_nbins = len(unique_idxbins) # this should be the same!
    if new_nbins!=nbins:
        print('correcting for increased nbins! nbins=%i' % new_nbins)
    x_bin = np.zeros(new_nbins,dtype=np.float)
    y_bin = np.zeros(new_nbins,dtype=np.float)
    for idx in unique_idxbins:
        idxrange = np.where(idx_bins==idx)
        x_bin[idx] = x_fold[idxrange].mean()
        y_bin[idx] = ypts[idxrange].mean()

    return x_bin,y_bin

def demodulate(self,asic=None,TES=None,offset=None,interval=None,calsource=True,align_clocks=False,doplot=True,xwin=True):
    '''
    Interpolating source data to the timestamps of the data

    offset is the constant offset between calsource timestamps and data timestamps
    this is of the order of 0.1sec
    
    interval is the time interval in which to calculate the demodulation

    calsource: you can force to not use the calsource

    align_clocks:  you can force the timestamps of the calsource and the data to start at the same time

    doplot: make plots

    xwin:  if making a plot, do not plot to screen if xwin==False
    '''
    if asic is None:
        print("\nPlease give an asic number\n")
        return None
    if TES is None:
        print("\nPlease give a TES number\n")
        return None
    if interval is None: # set interval to have no effect
        t0_interval = 0
        tend_interval = float((dt.datetime.utcnow() + dt.timedelta(days=7)).strftime('%s'))
        interval = (t0_interval,tend_interval)
    interval = np.array(interval)

    errlist = []
    retval = {}
    retval['asic'] = asic
    retval['TES'] = TES
    retval['dataset'] = self.dataset_name
    retval['calsource info'] = self.calsource_info()

    
    data = self.timeline(asic=asic,TES=TES)
    t_data_orig = self.timeaxis(datatype='sci',asic=asic)
    t_data = t_data_orig.copy()
    t0_data = t_data[0]
    
    t_src_orig,v_src = self.calsource()
    if v_src is None or not calsource:
        use_calsource = False
        msg = 'No calsource.'
        self.printmsg(msg,verbosity=3)
        errlist.append(msg)
        t_src = t_data
        data_src = np.ones(t_data.shape,dtype=np.float)
    else:
        # shift source to oscillate around zero
        use_calsource = True
        t_src = t_src_orig.copy()
        data_src = -v_src + v_src.mean()# source sampling inverted compared to data
        
    t0_src = t_src[0] # this will be modified after the offset is calculated below
    retval['t0 source'] =  t0_src

    # for some data, the source is in UT while the detector data is in localtime
    delta_tz = np.abs(t0_data - t0_src)
    tz_offset = 0.0
    tz_count = 0
    while delta_tz>3599:
        tz_offset += 3600
        delta_tz -= 3600
        tz_count += 1
    if (t0_data - t0_src) < 0:
        tz_offset = -tz_offset
    retval['time zone offset'] = tz_offset
    retval['time zone offset hours'] = tz_count
    
    # put data in UT (assume calsource was in UT)
    t_data -= tz_offset
    t0_data = t_data[0]
    interval -= tz_offset

    # if forcing the alignment of the clocks, make t_src start at the same time as t_data
    if align_clocks:
        align_offset = t_src[0] - t_data[0]
        t_src -= align_offset
        t0_src = t_src[0]
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
    idx_start = 0
    idx_end = -1
    if t0_data<tstart:
        idxrange = np.where(t_data>=tstart)[0]
        idx_start = idxrange[0]
    if t_data[-1]>tend:
        idxrange = np.where(t_data>=tend)[0]
        idx_end = idxrange[0]
    retval['data time start index'] = idx_start
    retval['data time end index'] = idx_end
    self.printmsg('idx_start, idx_end = %i, %i' % (idx_start,idx_end),verbosity=3)
    t_data = t_data[idx_start:idx_end]
    data = data[idx_start:idx_end]
    t0_data = t_data[0]
    retval['t0 data'] = t0_data
    t0_str = dt.datetime.utcfromtimestamp(t0_data).strftime('%Y-%m-%d %H:%M:%S.%f')
    t0_filenamestr = dt.datetime.utcfromtimestamp(t0_data).strftime('%Y%m%d-%H%M%S')

    # truncate the source timeline to overlapping time or to the given interval
    idxsrc_start = 0
    idxsrc_end = -1
    if t0_src<tstart:
        idxrange = np.where(t_src>=tstart)[0]
        idxsrc_start = idxrange[0]
    if t_src[-1]>tend:
        idxrange = np.where(t_src>=tend)[0]
        idxsrc_end = idxrange[0]        
    retval['source time start index'] = idxsrc_start
    retval['source time end index'] = idxsrc_end
    t_src = t_src[idxsrc_start:idxsrc_end]
    data_src = data_src[idxsrc_start:idxsrc_end]
    
    
    # fit data to sine curve
    amplitude = 0.5*(data.max() - data.min())
    if self.calsource_info() is not None:
        period = self.calsource_info()['modulator']['frequency']
    else:
        period = 1.0
    # timeshift = 0.5*period # there's a half period shift:  see http://qubic.in2p3.fr/wiki/pmwiki.php/TD/Demodulation
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
        amplitude = 0.5*(data_src.max() - data_src.min())
        if self.calsource_info() is not None:
            period = self.calsource_info()['modulator']['frequency']
        else:
            period = 1.0
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

    retval['source fit'] = source_fit
    retval['t_src max'] = t_src_max
    period = source_fit['period']
    retval['period'] = period # this is the period we use for demodulation

    # find the constant timestamp offset between source and data
    # for more info: http://qubic.in2p3.fr/wiki/pmwiki.php/TD/Demodulation
    if offset is None:
        offset = t_data_max - t_src_max
    retval['source-data time offset'] = offset
    
    
    # adjust the source timeline by the constant offset
    # go back to the original source data before applying the offset
    # and then find the appropriate interval
    if use_calsource:
        t_src_orig,v_src = self.calsource()
        t_src = t_src_orig.copy()
        data_src = -v_src + v_src.mean()
    t_src += offset

    #t_src_max += offset
    t0_src = t_src[0]
    # truncate the source timeline to overlapping time or to the given interval
    if t0_src<tstart:
        idxrange = np.where(t_src>=tstart)[0]
        idxsrc_start = idxrange[0]
    if t_src[-1]>tend:
        idxrange = np.where(t_src>=tend)[0]
        idxsrc_end = idxrange[0]        
    retval['source time start index'] = idxsrc_start
    retval['source time end index'] = idxsrc_end
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

    ### Interpolating source data to the timestamps of the data
    ### and making the product of the detector and the source
    ### the source is renormalized, but the data is only shifted to oscillate around 0
    data_src_interp = np.interp(t_data-t0_data, t_src-t0_data, data_src)
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
    folded_t_data,folded_data = fold_data(t_data-t0_data,data,period)
    folded_t_src, folded_src = fold_data(t_src-t0_data,data_src,period)
    retval['folded t_data'] = folded_t_data
    retval['folded data'] = folded_data
    retval['folded t_src'] = folded_t_src
    retval['folded src'] = folded_src

    # if not plotting, exit now
    if not doplot:
        retval['axes'] = None
        retval['fig'] = None
        return retval

    if xwin: plt.ion()
    else: plt.ioff()
    fig = plt.figure()
    retval['fig'] = fig
    fig.canvas.set_window_title('plt: data/calsource interpolation')
    axes = []

    # first column of plots
    
    axes.append(fig.add_axes((0.03,0.68,0.62,0.3)))
    axes[-1].plot((t_data-t0_data), renorm(data),label='TES {} ASIC {}'.format(TES,asic),ls='none',marker='x')
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
    axes[-1].plot((t_data-t0_data), renorm(data), label='TES {} ASIC {}'.format(TES,asic),marker='v',ls='none')
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
    axes[-1].plot((t_data-t0_data)%period, renorm(data),label='folded TES {} ASIC {}'.format(TES,asic),ls='none',marker='x')
    axes[-1].plot(folded_t_data,renorm(folded_data),label='folded and averaged TES {} ASIC {}'.format(TES,asic),ls='solid',linewidth=3,color='yellow')
    if use_calsource:
        axes[-1].plot((t_src-t0_data)%period, renorm(data_src), label='folded source',ls='none',marker='+')
        axes[-1].plot(folded_t_src,renorm(folded_src),label='folded and averaged source',ls='solid',linewidth=3,color='red')
    axes[-1].text(0.5,1.0,'folding period = %.6f seconds' % period,ha='center',va='bottom',transform=axes[-1].transAxes)
    axes[-1].legend()

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

    pngname = 'demodulation_diagnostic_ASIC%i_TES%03i_%s.png' % (asic,TES,t0_filenamestr)
    fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    retval['pngname'] = pngname
        
    retval['axes'] = axes
    if not xwin:
        retval['fig'] = None
        retval['axes'] = None
        plt.close(fig)
    return retval
