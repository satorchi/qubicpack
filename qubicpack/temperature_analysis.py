#!/usr/bin/env python
'''
$Id: temperature_analysis.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Fri 18 Aug 2017 18:57:44 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

analyze I-V curves at different temperatures to get the NEP
'''
import os
from pprint import pprint
import matplotlib.pyplot as plt
import numpy as np
from iminuit import Minuit
import pickle
from glob import glob
import datetime as dt

from qubicpack.qubicfp import qubicfp
from qubicpack.pix2tes import tes2pix
from qubicpack.utilities import asic_reversal_date, NPIXELS, TES_index, ASIC_index, figure_window_title, fmt4latex
NASIC = qubicfp.NASIC

# some constants and values required
kBoltzmann=1.3806485279e-23
temperature_precision = 0.005 # close enough for temperature

# text box in plots
boxprops = {}
boxprops['alpha'] = 0.4
boxprops['color'] = 'grey'
boxprops['boxstyle'] = 'round'
    
def print_datlist(datlist,obsdate=None,temperature=None):
    '''
    check if we're using a list of qubicfp objects or qubicasic/qubicpack objects
    '''
    if not isinstance(datlist,list):datlist=[datlist]
    print(' idx        array ASIC date                temp')

    if datlist[0].__object_type__=='qubicpack' or datlist[0].__object_type__=='qubicasic':
        return print_asic_datlist(datlist,obsdate,temperature)

    if not datlist[0].__object_type__=='qubicfp':
        print('ERROR! This is not a list of qubicfp or qubicpack objects')
        return None

    for idx,go in enumerate(datlist):
        print('[%2i]' % idx)
        print_asic_datlist(go.asic_list,obsdate,temperature)

    return    
    

def print_asic_datlist(datlist,obsdate=None,temperature=None):
    '''
    print some of the main parameters of all the data in a list of qubicfp objects
    select members of the list according to obsdate and/or temperature
    '''
    datstr='    [%2i][%2i] %s ASIC%i %s %.3fmK'
    for idx,go in enumerate(datlist):
        if go is None: continue
        
        printit=True
        if go.exist_timeline_data():
            ntimelines=go.ntimelines()
            for t_idx in range(ntimelines):
                printit=True
                if 'TES_TEMP' in go.tdata[t_idx].keys():
                    T = go.tdata[t_idx]['TES_TEMP']
                else:
                    T = None
                if 'BEG-OBS' in go.tdata[t_idx].keys():
                    d = go.tdata[t_idx]['BEG-OBS']
                else:
                    d = go.tdata[t_idx]['DATE-OBS']
                if obsdate is not None and d!=obsdate:
                    printit = False
                if T is None:
                    printit = False
                elif temperature is not None and not (T<temperature+0.001 and T>temperature-0.001):
                    printit = False
                if printit: print(datstr % (idx,t_idx,go.detector_name,go.asic,d,1000*T))
        else:
            T = go.temperature
            d = go.obsdate
            t_idx = 0
            if obsdate is not None and d!=obsdate:
                printit = False
            if temperature is not None and not (T<temperature+0.001 and T>temperature-0.001):
                printit = False
            if printit: print(datstr % (idx,t_idx,go.detector_name,go.asic,d,1000*T))
            
    return

def is_300mK_measurement(go):
    '''
    check if this is a measurement taken at 300mK
    the argument "go" is an object of type qubicfp, qubicasic, or qubicpack
    '''
    if go.temperature>=0.3-temperature_precision\
       and go.temperature<=0.3+temperature_precision:
        return True
    return False

def read_data_from_20170804():
    '''
    read data from the test campaign of 2017-08-04
    I-V curves of P73 at different bath temperatures
    '''
    files=[]
    datlist=[]

    files.append('QUBIC_TES_20170804T134238UTC.fits')
    files.append('QUBIC_TES_20170804T144916UTC.fits')
    files.append('QUBIC_TES_20170804T150120UTC.fits')
    files.append('QUBIC_TES_20170804T151319UTC.fits')
    files.append('QUBIC_TES_20170804T152431UTC.fits')
    files.append('QUBIC_TES_20170804T153704UTC.fits')
    files.append('QUBIC_TES_20170804T155105UTC.fits')
    files.append('QUBIC_TES_20170804T160111UTC.fits')

    idx=0
    for F in files:
        datlist.append(qp())
        datlist[idx].read_fits(F)
        idx+=1
    return datlist

def read_data_from_20170905():
    '''
    read the data from the measurement campaign of 5/6 Sept 2017
    '''
    fitslist=['QUBIC_TES_20170905T082847UTC.fits',
              'QUBIC_TES_20170905T091628UTC.fits',
              'QUBIC_TES_20170905T095550UTC.fits',
              'QUBIC_TES_20170905T145150UTC.fits',
              'QUBIC_TES_20170905T154040UTC.fits',
              'QUBIC_TES_20170905T165531UTC.fits',
              'QUBIC_TES_20170905T172101UTC.fits',
              'QUBIC_TES_20170906T075624UTC.fits',
              'QUBIC_TES_20170906T085757UTC.fits',
              'QUBIC_TES_20170906T091909UTC.fits',
              'QUBIC_TES_20170906T094433UTC.fits',
              'QUBIC_TES_20170906T101411UTC.fits',
              'QUBIC_TES_20170906T113411UTC.fits',
              'QUBIC_TES_20170906T121954UTC.fits',
              'QUBIC_TES_20170906T123554UTC.fits']
    datlist=[]
    for F in fitslist:
        datlist.append(qp())
        datlist[-1].read_fits(F)
    
    return datlist

def verify_temperature_arguments(fplist,TES,asic):
    # make sure TES is a valid selection
    try:
        TES = int(TES)
        asic = int(asic)
    except:
        print('ERROR! Please enter a valid TES and asic number.')
        return False

    if TES<1 or TES>NPIXELS:
        print('ERROR! Please enter a valid TES number between 1 and %i.' % NPIXELS)
        return False
        
    # make sure we have a list of qubicpack.qubicfp objects
    if not isinstance(fplist, list):
        print('ERROR!  Please provide a list of qubicpack.qubicfp objects')
        return False

    if asic<1 or asic>NASIC:
        print('ERROR! Please enter a valid asic number: %i' % asic)
        return False
    asic_idx = asic - 1

    detector_name = None
    for fp in fplist:
        if not isinstance(fp,qubicfp):
            print('ERROR! The list should contain qubicpack.qubicfp objects')
            return False
        
        if fp.asic_list[asic_idx] is None:
            continue

        if detector_name is None:
            detector_name = fp.asic_list[asic_idx].detector_name
            
        if fp.asic_list[asic_idx].detector_name != detector_name:
            print('ERROR! These data are not for the same detector array.')
            return False
        

    if detector_name is None: return False

    return True

def get_temperature_info(fplist,TES,asic):
    '''
    return some information to label the plots
    '''
    if not verify_temperature_arguments(fplist,TES,asic): return None
    
    asic_idx = asic - 1
    detector_name = None
    temps_list = []
    temps_index_list = []
    
    turnover_list = []
    turnover_index_list = []
    turnover_temps_list = []
    
    obsdates = []

    for idx,fp in enumerate(fplist):
        go = fp.asic_list[asic_idx]
        if go is None: continue
        if go.temperature is None: continue
        if go.temperature=='': continue
        
        if detector_name is None:
            detector_name = go.detector_name

        temps_index_list.append(idx)
        temps_list.append(go.temperature)                        
        obsdates.append(go.obsdate)
            
        if go.turnover(TES) is not None:
            turnover_index_list.append(idx)
            turnover_list.append(go.turnover(TES))
            turnover_temps_list.append(go.temperature)

    temps_list = np.array(temps_list)
    temps_index_list = np.array(temps_index_list)
    turnover_list = np.array(turnover_list)
    turnover_temps_list = np.array(turnover_temps_list)
    
    sorted_turnover_temps_index = sorted(range(len(turnover_temps_list)), key=lambda i: turnover_temps_list[i])
    sorted_turnover = turnover_list[sorted_turnover_temps_index]
    sorted_turnover_temps = turnover_temps_list[sorted_turnover_temps_index]

    sorted_temps_index = sorted(range(len(temps_list)), key=lambda i: temps_list[i])
    sorted_temps = temps_list[sorted_temps_index]
    sorted_index = temps_index_list[sorted_temps_index]

    startdate = min(obsdates)
    enddate = max(obsdates)
    ymd_start = (startdate.year,startdate.month,startdate.day)
    ymd_end = (enddate.year,enddate.month,enddate.day)
    if ymd_start==ymd_end:
        datadate_str = '%s to %s' % (startdate.strftime('%Y-%m-%d %H:%M'),enddate.strftime('%H:%M'))
        fname_datestr = '%s-%s' % (startdate.strftime('%Y%m%dT%H%M%S'),enddate.strftime('%H%M%S'))
    else:
        datadate_str = '%s to %s' % (startdate.strftime('%Y-%m-%d %H:%M'),enddate.strftime('%Y-%m-%d %H:%M'))
        fname_datestr = '%s-%s' % (startdate.strftime('%Y%m%dT%H%M%S'),enddate.strftime('%Y%m%dT%H%M%S'))

    ret = {}
    ret['ASIC'] = asic
    ret['TES'] = TES
    ret['detector_name'] = detector_name
    ret['temps_index_list'] = np.array(temps_index_list)
    ret['sorted_index'] = np.array(sorted_index)
    ret['sorted_temps'] = np.array(sorted_temps)
    
    ret['turnover_list'] = np.array(turnover_list)
    ret['turnover_index_list'] = np.array(turnover_index_list)
    ret['sorted_turnover'] = np.array(sorted_turnover)
    ret['sorted_turnover_temps'] = np.array(sorted_turnover_temps)
    
    ret['obsdates'] = obsdates
    ret['datadate_str'] = datadate_str
    ret['fname_datestr'] = fname_datestr
    
    return ret

def plot_TES_turnover_temperature(fplist,TES,asic,xwin=True):
    '''
    plot the turnover point as a function of temperature for a given TES
    '''
    info = get_temperature_info(fplist,TES,asic)
    if info is None: return
    
    datadate_str = info['datadate_str']
    fname_datestr = info['fname_datestr']
    detector_name = info['detector_name']
    
    sorted_temps = info['sorted_turnover_temps']
    sorted_turnover = info['sorted_turnover']
    
    pngname='QUBIC_Array-%s_TES%03i_ASIC%i_Turnover_Temperature_%s.png' % (detector_name,TES,asic,fname_datestr)
    xlabel='T$_{bath}$ / mK'
    ylabel='V$_{turnover}$ / V'

    ttl='QUBIC Array %s, ASIC %i, TES #%i: Turnover at Different Temperatures' % (detector_name,asic,TES)
    subttl = 'Measurements of %s' % datadate_str
    
    if xwin:plt.ion()
    else:
        plt.close('all')
        plt.ioff()
        
    fig = plt.figure()
    figure_window_title(fig,ttl)

    plt.suptitle(ttl+'\n'+subttl)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.plot(sorted_temps,sorted_turnover,linestyle='none',marker='D')
    plt.plot(sorted_temps,sorted_turnover,color='green')

    if len(sorted_turnover)==0:
        print('No turnover found at any temperature for TES %i on ASIC %i' % (TES,asic))
        xmax = info['sorted_temps'].max()
        xmin = info['sorted_temps'].min()
    else:
        xmax = sorted_temps.max()
        xmin = sorted_temps.min()
    
    span=xmax-xmin
    plot_xlim=(xmin-0.05*span,xmax+0.1*span)
    
    ax=plt.gca()
    ax.set_xlim(plot_xlim)
    
    plt.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    if xwin: plt.show()
    else: plt.close('all')
    
    return pngname


def plot_TES_temperature_curves(fplist,TES,asic,plot='I',xwin=True):
    '''
    plot the I-V, P-V, R-P curves for each temperature
    '''
    info = get_temperature_info(fplist,TES,asic)
    if info is None: return

    detector_name = info['detector_name']
    fname_datestr = info['fname_datestr']
    datadate_str = info['datadate_str']
    sorted_index = info['sorted_index']

    plot_type='I'
    if plot.upper()[0]=='R':
        plot_type='R'
        pngname='QUBIC_Array-%s_TES%03i_ASIC%i_R-V_Temperatures_%s.png' % (detector_name,TES,asic,fname_datestr)
        #xlabel='P$_{TES}$ / $p$W'
        xlabel='V$_{bias}$ / V'
        ylabel='$\\frac{R_\mathrm{TES}}{R_\mathrm{normal}}$ / %'
    elif plot.upper()[0]=='P':
        plot_type='P'
        pngname='QUBIC_Array-%s_TES%03i_ASIC%i_P-V_Temperatures_%s.png' % (detector_name,TES,asic,fname_datestr)
        xlabel='V$_{bias}$ / V'
        ylabel='P$_{TES}$ / $p$W'
    else:
        plot_type='I'
        pngname='QUBIC_Array-%s_TES%03i_ASIC%i_I-V_Temperatures_%s.png' % (detector_name,TES,asic,fname_datestr)
        xlabel='V$_{bias}$ / V'
        ylabel='I$_{TES}$ / $\mu$A'
        
    ttl='QUBIC Array %s ASIC %i TES%03i at Different Temperatures\nmeasurements from %s' % (detector_name,asic,TES,datadate_str)
    if xwin:plt.ion()
    else:
        plt.close('all')
        plt.ioff()
        
    fig = plt.figure()
    figure_window_title(fig,ttl)

    plt.title(ttl)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    min_bias=1000.
    max_bias=-1000.
    min_Rn_ratio=1000.
    max_Rn_ratio=-1000.
    min_P=1000.
    max_P=-1000.
    for idx in sorted_index:
        go = fplist[idx].asic(asic)
        if go is None: continue
        
        lbl = '%.0f mK' % (1000*go.temperature)

        startend = go.selected_iv_curve(TES)
        if startend is None:
            print('ERROR! index=%i, Tbath=%s: Could not get selected curve for TES#%i' % (idx,lbl,TES))
            continue
        istart,iend = startend

        Iadjusted = go.adjusted_iv(TES)
        I = Iadjusted[istart:iend]
        bias = go.bias_factor*go.vbias[istart:iend]


        # power calculation
        Ites=go.Ites(TES)[istart:iend]
        Vtes=go.Vtes(TES)[istart:iend]
        Ptes=go.Ptes(TES)[istart:iend]

        if plot_type=='R':
            # plot normal resistance vs power
            if not go.R1(TES) is None:
                Rn_ratio=go.Rn_ratio(TES)[istart:iend]
                if min(Rn_ratio)<min_Rn_ratio:min_Rn_ratio=min(Rn_ratio)
                if max(Rn_ratio)>max_Rn_ratio:max_Rn_ratio=max(Rn_ratio)
                Pbias=go.Pbias(TES)
                lbl+=', P$_\mathrm{bias}=$%.2f pW' % Pbias                
                #plt.plot(Ptes,Rn_ratio,label=lbl)
                plt.plot(bias,Rn_ratio,label=lbl)
        elif plot_type=='P':
            # plot power vs bias
            plt.plot(bias,Ptes,label=lbl)
        else:
            # plot current vs bias
            plt.plot(bias,I,label=lbl)

        if min(bias)<min_bias:min_bias=min(bias)
        if max(bias)>max_bias:max_bias=max(bias)        
        if min(Ptes)<min_P:min_P=min(Ptes)
        if max(Ptes)>max_P:max_P=max(Ptes)        

    if plot_type=='R':
        xmin=min_P
        xmax=max_P
        xmin=min_bias
        xmax=max_bias
        # draw the line at 90% Rn
        plt.plot([xmin,xmax],[90,90],color='black',linestyle='dashed',linewidth=3)
        
    else:
        xmin=min_bias
        xmax=max_bias
        
    span=xmax-xmin
    plot_xlim=(xmin-0.05*span,xmax+0.35*span)
    
    ax=plt.gca()
    ax.set_xlim(plot_xlim)
        
    plt.legend()
    plt.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    if xwin: plt.show()
    else: plt.close('all')
    
    return pngname


def P_bath_function(Tbath,K,T0,n):
    '''
    TES power should follow this as a function of bath temperature
    see Perbost 2016, PhD thesis, eq. 2.11 and 6.16
    '''
    P = K*(T0**n - Tbath**n)
    return P

def least_squares(K,T0,n):
    global P_data, T_data
    Pmean = P_data.mean()
    ls = np.sum( ((P_data - P_bath_function(T_data,K,T0,n))/Pmean)**2 )
    return ls
        
def redo_fit_Pbath(T_pts, P_pts, p0=None,ftol=1e-8,verbosity=0):
    '''
    run through fit multiple times giving a better first guess each time
    '''

    ret = {}
    K_prev = 1.0
    Kdiff_prev = 1e9
    Klist = []
    
    ctr = 1
    while True:
    
        res = fit_Pbath(T_pts, P_pts, p0=p0,ftol=ftol,verbosity=verbosity)
        if res is None:
            ret['fit'] = None
            ret['tolerance'] = ftol
            ret['counter'] = ctr
            return ret

        K = res['K']
        T0 = res['T0']
        n = res['n']
        Klist.append(K)
        Kdiff = np.abs(1 - K/K_prev)
        Kdiffdiff = Kdiff - Kdiff_prev

        ret = res
        ret['counter'] = ctr
        ret['Kdiff'] = Kdiff
        ret['Kdiffdiff'] = Kdiffdiff
        ret['Klist'] = Klist

        if Kdiff<1e-10:
            ret['comment'] = 'K converged'
            return ret

        if K>9e-10 and ctr>2:
            ret['comment'] = 'K big'
            return ret

        if K<0 and ctr>2:
            ret['comment'] = 'K negative'
            return ret
        
        if Kdiffdiff>0 and ctr>21:
            ret['comment'] = 'K diverging'
            return ret

        if ctr==2000:
            ret['comment'] = 'max iterations'
            return ret

        Kdiff_prev = Kdiff
        K_prev = K
        p0 = np.array([K,T0,n])
        ctr += 1

        if verbosity>0:
            print('retrying fit with better first guess')

    return None

def fit_Pbath(T_pts, P_pts, p0=None,ftol=1e-8,verbosity=0):
    '''
    find best fit to the P_bath function
    '''
    global T_data, P_data
    start_ftol = ftol
    ret = {}
    ret['success'] = False
    ret['fit'] = None

    # first guess based on results on array P87 measured 2020-Jan-06
    K = 3.2e-07
    T0 = 0.412
    n = 3.0
    # if p0 is None: p0 = np.array([1E-10, 0.5,4.5]) #MP
    if p0 is None:
        p0 = np.array([K,T0,n])
    else:
        (K,T0,n) = p0

    ret['K'] = K
    ret['T0'] = T0
    ret['n'] = n
        

    # make sure T_pts and P_pts are 1d arrays
    npts = T_pts.size
    if npts<3:
        if verbosity>0:print('insufficient number of points for fit_Pbath: %i' % npts)
        return None
    T_data = np.array(T_pts).reshape(npts)
    P_data = np.array(P_pts).reshape(npts)

    try:
        m = Minuit(least_squares,
                   K=K,
                   T0=T0,
                   n=n,
                   errordef=1,
                   error_K=0.2*K,
                   error_T0=0.2*T0,
                   error_n=0.2*n,
                   limit_T0=(0,2),
                   limit_n=(0,10))
        
        m.migrad()
        m.get_fmin()
        ans = m.get_param_states()
        for idx,parm in enumerate(ans):
            ret[parm['name']] = parm['value']
            if verbosity>0:
                print('%i: %s=%.6e fixed=%s, constant=%s' % (idx,parm['name'],parm['value'],parm['is_fixed'],parm['is_const']))
        Chisq = least_squares(ret['K'],ret['T0'],ret['n'])
        ret['Chi square'] = Chisq
        ret['fit'] = ((K,T0,n),None) # retaining compatibility with curve_fit output
        ret['success'] = True
        ret['minuit'] = True
        return ret
    
    except:
        ret['minuit'] = False
        m = None
        if verbosity>0:
            print('Minuit unsuccessful')
            
    
    
    # try to fit the curve with the default tolerance (1e-8)
    # if unsuccessful, relax the tolerance by a factor 10 and try again

    # bounds are only applicable to dogbox and trf, but those don't work properly
    for ctr in range(9):
        try:
            fit = curve_fit(P_bath_function,T_data,P_data,p0=p0,ftol=ftol,method='lm')
            ret['fit'] = fit
            ret['K'],ret['T0'],ret['n'] = fit[0][0:3]
            ret['tolerance'] = ftol
            ret['curve fit attempts'] = ctr+1
            ret['success'] = True
            return ret
        except:
            ret['tolerance'] = ftol
            ret['curve fit attempts'] = ctr+1
            ftol *= 10
            if verbosity>0: print('retrying fit with relaxed tolerance: %.1e' % ftol)
       
    return ret

def calculate_TES_NEP(fplist,TES,asic,p0=None,T0_limit=0.7,n_limit=8,mean_istart=0,mean_iend=10,verbosity=0):
    '''
    make the list of temperatures and associated P0
    and calculate the NEP
    '''
    if not verify_temperature_arguments(fplist,TES,asic):return None
    asic_idx = asic - 1
    detector_name = None
    temps_list=[]
    for fp in fplist:
        go = fp.asic_list[asic_idx]
        if go is None: continue
        if detector_name is None:
            detector_name = go.detector_name
        temps_list.append(go.temperature)
    temps_list=np.array(temps_list)
    sorted_index=sorted(range(len(temps_list)), key=lambda i: temps_list[i])
    sorted_temps=temps_list[sorted_index]
    
    ret={}
    ret['DET_NAME'] = detector_name
    ret['TES'] = TES
    ret['ASIC'] = asic
    ret['T0_limit'] = T0_limit
    ret['Tmin'] = temps_list.min()
    ret['Tmax'] = temps_list.max()
    ret['mean_istart'] = []
    ret['mean_iend'] = []
    
    # make the arrays of Power and T_bath
    P = []
    T = []
    all_T = []
    all_P = []
    obsdates = []
    for idx in sorted_index:
        fp = fplist[idx]
        go = fp.asic_list[asic_idx]
        if go is None: continue
        if verbosity>0:
            print('calculating NEP for dataset %2i: ASIC%i, TES%3i, Tbath=%.1fmK, %s'\
                  % (idx,go.asic,TES,1000*go.temperature,go.obsdate))
        if not go.is_good_iv(TES):
            if verbosity>0: print('  rejected.  bad I-V')
            continue
        if go.turnover(TES) is None:
            if verbosity>0: print('   rejected.  no turnover')
            continue
        
        filterinfo = go.filterinfo(TES)
        
        istart,iend = go.selected_iv_curve(TES)
        npts = iend-istart

        Iadjusted = go.adjusted_iv(TES)
        I = Iadjusted[istart:iend]
        if 'Iturnover' in filterinfo['fit'].keys() and filterinfo['fit']['Iturnover'] is not None:
            Iturnover = 1e-6*filterinfo['fit']['Iturnover']
        else:
            Iturnover = 1e-6*I.min()
            
        Tbath = go.temperature
        Ptes = go.Ptes(TES)[istart:iend] #MP
        Vbias = go.vbias[istart:iend]
        Vtes_turnover = go.Rshunt*(go.turnover(TES)/go.Rbias-Iturnover)
        all_P.append(Ptes.mean()*1e-12)
        all_T.append(Tbath)
        obsdates.append(go.obsdate)

        # take the mean of the Pbias in the superconducting region
        Vsuper = filterinfo['fit']['Vsuper']
        super_idx = np.where(Vbias<Vsuper)[0]
        if len(super_idx)==0:
            ret['mean_istart'] = None
            ret['mean_iend'] = None
            continue
        
        Pbeg = np.mean(Ptes[super_idx])
        ret['mean_istart'] = istart+super_idx[0]
        ret['mean_iend'] = istart+super_idx[-1]

        '''
        if Ptes[-1]<Ptes[0]:
            curve_mean_istart = npts - mean_istart
            curve_mean_iend = npts - mean_iend
        else:
            curve_mean_istart = mean_istart
            curve_mean_iend = mean_iend
            
            
        ret['mean_istart'].append(curve_mean_istart)
        ret['mean_iend'].append(curve_mean_iend)
        
        Pbeg = np.mean(Ptes[curve_mean_istart:curve_mean_iend])
        #if ((Pbeg > 5) and (Pbeg < 40)):
        '''
        
        if Pbeg>0.0:
            P.append(Pbeg*1e-12)
            T.append(Tbath)

    P = np.array(P)
    T = np.array(T)
    all_P = np.array(all_P)
    all_T = np.array(all_T)
    ret['P'] = P
    ret['T'] = T
    ret['all_temperatures'] = all_T
    ret['all_P'] = all_P
    ret['comment'] = ''
    comments = []
    obsdates.sort()
    ret['obsdates'] = obsdates 
        
    ret['is_good'] = None
    npts = len(P)
    temperature_fit = {}
    if npts<3:
        temperature_fit['fit'] = None
    else:
        temperature_fit = redo_fit_Pbath(T,P,p0=p0,verbosity=verbosity)
        ret['fit points'] = '$<P>$ determined from superconducting region'
        ret['fit npts'] = npts

    if temperature_fit['fit'] is None:
        # try again with the full range but only for P>0
        idx_range = np.where(all_P>0)[0]
        ret['fit npts'] = len(idx_range)
        if len(idx_range)!=len(all_P):
            ret['fit points'] = '$<P>$ determined from full range, temperatures filtered for P values>0'
        else:
            ret['fit points'] = '$<P>$ determined from full range, all available temperatures'

        if len(idx_range)<3:
            comments.append('insufficient data: Npts=%i' % len(idx_range))
        else:
            temperature_fit = redo_fit_Pbath(all_T[idx_range],all_P[idx_range],p0=p0,verbosity=verbosity)

        
    
    if temperature_fit['fit'] is None:
        ret['K'] = None
        ret['T0'] = None
        ret['n'] = None
        ret['NEP'] = None
        ret['G'] = None
        ret['gamma'] = None
        ret['is_good'] = False
        if verbosity>0:
            print('insufficient data for curve fit:  ASIC%i, TES=%i' % (ret['ASIC'],ret['TES']))
        return ret
    
        
    K = temperature_fit['K']
    T0 = temperature_fit['T0']
    n = temperature_fit['n']
    ret['K'] = K
    ret['T0'] = T0
    ret['n'] = n
    for key in ['fit tolerance','counter','Kdiff','Kdiffdiff','Klist','Chi square']:
        if key in temperature_fit.keys():
            ret[key] = temperature_fit[key]
    ret['fit comment'] = temperature_fit['comment']

    if n<0:
        ret['is_good'] = False
        comments.append('$n<0$')

    if n>n_limit:
        ret['is_good'] = False
        comments.append('$n>%.1f$' % n_limit)
    
    if T0>T0_limit:
        ret['is_good'] = False
        comments.append('T$_0 > %.1f$ mK' % (1000*T0_limit))

    ret['comment'] = '\n'.join(comments)

    G = n*K*(T0**(n-1))
    ret['G'] = G
    Tratio = 0.35/T0
    # gamma is defined in Perbost PhD eq. 2.72 (page 82)
    gamma = (n/(2*n+1)) * (1-Tratio**(2*n+1))/(1-Tratio**n)
    ret['gamma'] = gamma
    discr = gamma*kBoltzmann*G
    if discr<0.0:
        if verbosity>0: print('ERROR! Imaginary NEP!  TES=%i' % TES)
        NEP = -2*T0*np.sqrt(-discr)
        ret['is_good'] = False
    else:
        NEP = 2*T0*np.sqrt(discr)
    ret['NEP'] = NEP    

    if np.isnan(NEP):
        ret['is_good'] = False
        
    if ret['is_good'] is None:
        ret['is_good'] = True

    return ret

def make_TES_NEP_resultslist(fplist,p0=None,verbosity=0):
    '''
    make a list of NEP calculation results, one for each TES
    '''
    results=[]
    for asic_idx in range(NASIC):
        asic = asic_idx + 1
        for idx in range(NPIXELS):
            TES = 1 + idx
            if verbosity>0: print('calculating NEP for ASIC%i, TES%03i' % (asic,TES))
            res = calculate_TES_NEP(fplist,TES,asic,p0=p0,verbosity=verbosity)
            results.append(res)
            
    return results

def plot_TES_NEP(fplist=None,TES=None,asic=None,result=None,xwin=True,p0=None,mean_istart=0,mean_iend=10):
    '''
    plot the P vs. Temperature for a TES
    '''

    if result is None:
        result = calculate_TES_NEP(fplist,TES,asic,p0=p0,mean_istart=mean_istart,mean_iend=mean_iend)
        
    if result is None:return None

    TES = result['TES']
    asic = result['ASIC']
    detector_name = result['DET_NAME']

    all_T = result['all_temperatures']
    all_P = result['all_P']
    P = result['P']
    T = result['T']

    NEP = result['NEP']
    K = result['K']
    T0 = result['T0']
    n = result['n']
    G = result['G']
    npts = result['fit npts']
        
    Tmin = result['Tmin']
    Tmax = result['Tmax']
    T_span = Tmax-Tmin
    # add 30% to the plot edges
    plot_T_min = Tmin - T_span
    plot_T_max = Tmax + T_span
    T_span = plot_T_max - plot_T_min
    T_stepsize = T_span/1000

    if NEP is None:
        txt = 'NEP estimate is not possible.'
        if len(all_P)<3:
            txt += '\nInsufficient number of points: %i' % len(all_P)
        
        if len(all_P)==1:
            plot_P_min=all_P[0]-0.5*all_P[0]
            plot_P_max=all_P[0]+0.5*all_P[0]
        elif len(all_P)==0:
            plot_P_min=0.0
            plot_P_max=5e-11
        else:
            P_span=np.nanmax(all_P)-np.nanmin(all_P)
            plot_P_min=np.nanmin(all_P)-0.5*P_span
            plot_P_max=np.nanmax(all_P)+0.5*P_span
        if np.isnan(plot_P_min):
            plot_P_min = 0.0
        if np.isnan(plot_P_max):
            plot_P_max = 5e-11
    else:
        fit_T = np.arange(plot_T_min,plot_T_max,T_stepsize)
        fit_P = P_bath_function(fit_T,K,T0,n)

        P_span = np.nanmax(fit_P) - np.nanmin(fit_P)
        plot_P_min = np.nanmin(fit_P)-0.5*P_span
        plot_P_max = np.nanmax(fit_P)+0.5*P_span
        
        txt = '$\kappa=$%s W/K$^n$' % fmt4latex(K,2)
        txt += '\nT$_c$=%.1f mK' % (1000*T0)
        txt += '\nn=%.3f' % n
        txt += '\nG=%s W/K' % fmt4latex(G,2)
        txt += '\nNEP=%s at T$_{bath}$=350mK' % fmt4latex(NEP,2)
        txt += '\nN$_\mathrm{fit\,\,points}$ = %i' % npts
        txt += '\n%s' % result['fit points']
        txt += '\n$\chi^2 = $%s' % fmt4latex(result['Chi square'],2)
        if result['is_good']:
            txt += '\ndetector is GOOD'
        else:
            txt += '\ndetector is BAD\n'+result['comment']

        P_span=plot_P_max-plot_P_min
        if len(P)==0:
            dat_min = np.nanmin(all_P)
            dat_max = np.nanmax(all_P)
        else:
            dat_min = np.nanmin( [np.nanmin(P),np.nanmin(all_P)] )
            dat_max = np.nanmax( [np.nanmax(P),np.nanmax(all_P)] )
        dat_span = dat_max - dat_min
        span = max([dat_span,P_span])
        if plot_P_min>dat_min:
            plot_P_min = dat_min-0.5*span
        
        if plot_P_max<dat_max:
            plot_P_max = dat_max+0.5*span

    obsdate_list = result['obsdates']
    if obsdate_list:
        ymd_start = (obsdate_list[0].year,obsdate_list[0].month,obsdate_list[0].day)
        ymd_end = (obsdate_list[-1].year,obsdate_list[-1].month,obsdate_list[-1].day)
        if ymd_start==ymd_end:
            datadate_str = '%s to %s' % (obsdate_list[0].strftime('%Y-%m-%d %H:%M'),obsdate_list[-1].strftime('%H:%M'))
            fname_datestr = '%s-%s' % (obsdate_list[0].strftime('%Y%m%dT%H%M%S'),obsdate_list[-1].strftime('%H%M%S'))
        else:
            datadate_str = '%s to %s' % (obsdate_list[0].strftime('%Y-%m-%d %H:%M'),obsdate_list[-1].strftime('%Y-%m-%d %H:%M'))
            fname_datestr = '%s-%s' % (obsdate_list[0].strftime('%Y%m%dT%H%M%S'),obsdate_list[-1].strftime('%Y%m%dT%H%M%S'))
    else:
        datadate_str = 'insufficient data'
        fname_datestr = '_'

    pngname='QUBIC_Array-%s_TES%03i_ASIC%i_NEP_%s.png' % (detector_name,TES,asic,fname_datestr)
    result['pngname'] = pngname
    ttl='QUBIC Array %s, ASIC %i, TES #%i: NEP (%s)' % (detector_name,asic,TES,datadate_str)

    if xwin:plt.ion()
    else:
        plt.close('all')
        plt.ioff()
    fig = plt.figure()
    figure_window_title(fig,ttl)

    ax = plt.gca()
    ax.set_xlim(plot_T_min,plot_T_max)
    ax.set_ylim(1e12*plot_P_min,1e12*plot_P_max)
    plt.title(ttl)
    ax.set_xlabel('T$_\mathrm{bath}$ / K')
    ax.set_ylabel('Power / pWatt')
    if result['fit points'].find('full range')>0:
        ax.plot(all_T,1e12*all_P,ls='none',marker='D')
    else:
        ax.plot(T,1e12*P,ls='none',marker='D')
    if not NEP is None: plt.plot(fit_T,1e12*fit_P,color='red')
    ax.text(0.99,0.96,txt,va='top',ha='right',fontsize=plt.rcParams['legend.fontsize'],transform=ax.transAxes,bbox=boxprops)
    fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    if xwin:plt.show()
    else:plt.close('all')
    return result

def plot_NEP_histogram(NEPresults,xwin=True,nbins=10):
    '''
    plot the histogram of the NEP calculations

    NEPresults is a list of dictionaries returned by make_TES_NEP_resultslist()
    '''
    retval = {}

    T0_limit = NEPresults[0]['T0_limit']

    NEP_estimate = []
    TESlist = []
    Klist = []
    T0list = []
    nlist = []
    Glist = []
    good_idx = []
    bad_idx = []
    Chilist = []

    # check if these results are all for the same ASIC, and detector array
    # and make lists of results
    asic = -1
    for idx,res in enumerate(NEPresults):
        if res is None: continue
        
        if asic==-1:
            asic = res['ASIC']
            detector_name = res['DET_NAME']
            obsdate_list = res['obsdates']
            
        if res['ASIC']!=asic:
            asic = None
        if res['DET_NAME']!=detector_name:
            detector_name = None

        if res['is_good']:
            NEP_estimate.append(res['NEP'])
            TESlist.append(res['TES'])
            Klist.append(res['K'])
            T0list.append(res['T0'])
            nlist.append(res['n'])
            Glist.append(res['G'])
            Chilist.append(res['Chi square'])
            good_idx.append(idx)
            for obsdate in res['obsdates']:
                if obsdate not in obsdate_list:
                    obsdate_list.append(obsdate)
        else:
            bad_idx.append(idx)
                
    if asic==-1:
        print('ERROR!  No valid results in the dataset.')
        return None

    if len(good_idx)==0:
        print('ERROR!  No NEP estimates possible in this dataset')
        return None
        
    retval['ASIC'] = asic
    retval['DET_NAME'] = detector_name
    obsdate_list.sort()
    retval['obsdates'] = obsdate_list

    good_idx = np.array(good_idx)        
    bad_idx = np.array(bad_idx)  
    NEP_estimate=np.array(NEP_estimate)
    nNEP=len(NEP_estimate)
    NEPmean=NEP_estimate.mean()

    retval['good_idx'] = good_idx
    retval['ngood'] = len(good_idx)
    retval['bad_idx'] = bad_idx
    retval['NEP mean'] = NEPmean
    retval['NEP'] = NEP_estimate
    retval['TES'] = np.array(TESlist)
    retval['K'] = np.array(Klist)
    retval['T0'] = np.array(T0list)
    retval['n'] = np.array(nlist)
    retval['G'] = np.array(Glist)
    retval['Chi square'] = np.array(Chilist)

    ymd_start = (obsdate_list[0].year,obsdate_list[0].month,obsdate_list[0].day)
    ymd_end = (obsdate_list[-1].year,obsdate_list[-1].month,obsdate_list[-1].day)
    if ymd_start==ymd_end:
        datadate_str = '%s to %s' % (obsdate_list[0].strftime('%Y-%m-%d %H:%M'),obsdate_list[-1].strftime('%H:%M'))
        fname_datestr = '%s-%s' % (obsdate_list[0].strftime('%Y%m%dT%H%M%S'),obsdate_list[-1].strftime('%H%M%S'))
    else:
        datadate_str = '%s to %s' % (obsdate_list[0].strftime('%Y-%m-%d %H:%M'),obsdate_list[-1].strftime('%Y-%m-%d %H:%M'))
        fname_datestr = '%s-%s' % (obsdate_list[0].strftime('%Y%m%dT%H%M%S'),obsdate_list[-1].strftime('%Y%m%dT%H%M%S'))
    
    txt = 'NEP$_\mathrm{mean}=%.4f \\times 10^{-17}\mathrm{W}/\sqrt{\mathrm{Hz}}$' % (1e17*NEPmean)
    txt += '\nG$_\mathrm{mean}=%.4f \\times 10^{-10}\mathrm{W}/{\mathrm{K}}$' % (1e10*retval['G'].mean())
    txt += '\nn$_\mathrm{mean}=%.4f$' % (retval['n'].mean())
    txt += '\nT$_\mathrm{0\,mean}=%.4f$ mK' % (1e3*retval['T0'].mean())
    txt += '\n%i good TES out of %i. yield=%.1f%%' % (nNEP,len(NEPresults),100*nNEP/len(NEPresults))
    #txt += '\nT$_0$ limit = %.1f mK' % (1000*T0_limit)

    
    if asic is not None:
        pngname = 'QUBIC_TES_ASIC%i_KEYVAL_histogram_%s.png' % (asic,fname_datestr)
        ttl = 'QUBIC TES ASIC %i KEYVAL' % asic
    else:
        pngname = 'QUBIC_TES_KEYVAL_histogram_%s.png' % fname_datestr
        ttl = 'QUBIC TES KEYVAL'
    
    ttl+='\nhistogram of %i TES\ndata taken %s' % (nNEP,datadate_str)
    if xwin:plt.ion()
    else:
        plt.close('all')
        plt.ioff()


    xlabel = {}
    xlabel['NEP'] = 'NEP  /  ${W}/\sqrt{Hz}$'
    xlabel['K'] = 'K / $W/K^{n}$'
    xlabel['G'] = 'G / $W/K$'
    xlabel['T0'] = 'T$_0$ / mK'
    xlabel['n'] = '$n$'
    xlabel['Chi square'] = '$\chi^2$'
    for keyval in xlabel.keys():        
        retval['pngname %s' % keyval] = pngname.replace('KEYVAL',keyval).replace(' ','_')
        figttl = ttl.replace('KEYVAL',keyval)
        
        fig = plt.figure()
        figure_window_title(fig,figttl)

        ax = plt.gca()
        plt.title(figttl)
        ax.set_xlabel(xlabel[keyval])
        ax.set_ylabel('Number per bin')

        hist,binedge=np.histogram(retval[keyval],bins=nbins)
        # np.histogram returns the bin edges.  change this to the bin centres
        bincentre = (binedge[:-1] + binedge[1:]) / 2
        width = 0.7 * (binedge[1] - binedge[0])
        ax.bar(bincentre, hist, align='center',width=width)

        bin_span=binedge[-1] - binedge[0]
        plot_bin_min=binedge[0]-0.1*bin_span
        plot_bin_max=binedge[-1]+0.1*bin_span
        ax.set_xlim(plot_bin_min,plot_bin_max)

        n_max=max(hist)
        ax.set_ylim(0,n_max+1)

        ax.text(0.985,0.97,txt,ha='right',va='top',fontsize=plt.rcParams['legend.fontsize'],transform=ax.transAxes,bbox=boxprops)
        fig.savefig(retval['pngname %s' % keyval],format='png',dpi=100,bbox_inches='tight')
    
        if xwin:fig.show()
        else:plt.close('all')

    
    return retval
            
def make_TES_NEP_tex_report(fplist,NEPresults=None,refresh=True):
    '''
    make a LaTeX source file for a test report of NEP estimates
    the input arguments are a list of qubicpack objects with the I-V data at different temperatures
    and a list of results from plot_TES_NEP()
    '''
    if NEPresults is None:
        print('Please enter a list of NEP results.  You can use NEPresults=make_TES_NEP_resultslist(fplist)')
        return None
    
    
    # find the data at 300mK
    go300 = {}
    datelist = ''
    T300_diff = {1:1e9,2:1e9}
    asic_list = []
    obsdates = []
    for fp in fplist:
        for go in fp.asic_list:
            if go is None:continue
            obsdates.append(go.obsdate)
            asic_list.append(go.asic)
            delta = np.abs(go.temperature - 0.3)
            if delta < T300_diff[go.asic]:
                go300[go.asic] = go
                T300_diff[go.asic] = delta
            datelist += '\n\\item %.3fmK on %s'\
                % (1000*go.temperature,go.obsdate.strftime('%Y-%m-%d %H:%M:%S'))

    obsdates.sort()
    ymd_start = (obsdates[0].year,obsdates[0].month,obsdates[0].day)
    ymd_end = (obsdates[-1].year,obsdates[-1].month,obsdates[-1].day)
    if ymd_start==ymd_end:
        datadate_str = '%s to %s' % (obsdates[0].strftime('%Y-%m-%d %H:%M'),obsdates[-1].strftime('%H:%M'))
        fname_datestr = '%s-%s' % (obsdates[0].strftime('%Y%m%dT%H%M%S'),obsdates[-1].strftime('%H%M%S'))
    else:
        datadate_str = '%s to %s' % (obsdates[0].strftime('%Y-%m-%d %H:%M'),obsdates[-1].strftime('%Y-%m-%d %H:%M'))
        fname_datestr = '%s-%s' % (obsdates[0].strftime('%Y%m%dT%H%M%S'),obsdates[-1].strftime('%Y%m%dT%H%M%S'))



    for asic in go300.keys():
        go = go300[asic]
        print('300mK results for ASIC %i from %s, Tbath=%.1fmK'\
              % (go.asic,go.obsdate.strftime('%Y-%m-%d %H:%M:%S'),go.temperature*1000))
        observer = go.observer.replace('<','$<$').replace('>','$>$')
        detector_name = go.detector_name

        if go.transdic is None:
            show_extra_columns = False
        else:
            show_extra_columns = True

    asic_list = sorted(set(asic_list))
    if len(asic_list)==1:
        asic = asic_list[0]
        texfilename = 'QUBIC_Array-%s_ASIC%i_NEP_%s.tex' % (detector_name,asic,fname_datestr)
    else:
        asic = None
        texfilename = 'QUBIC_Array-%s_NEP_%s.tex' % (detector_name,fname_datestr)

    NEP_estimate=[]
    TESlist=[]
    for res in NEPresults:
        NEP=res['NEP']
        if res['is_good']:
            NEP_estimate.append(NEP)
            TESlist.append(res['TES'])
    nNEP=len(TESlist)
    NEP_estimate=np.array(NEP_estimate)
    NEPmean=NEP_estimate.mean()
    NPIXELS = len(NEPresults)
    
    h=open(texfilename,'w')
    h.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    h.write('%%%%% WARNING!  Automatically generated file.  Do not edit! %%%%%\n')
    h.write('%%%%% This file could be overwritten                        %%%%%\n')
    h.write(dt.datetime.utcnow().strftime('%%%%%%%%%% File generated %Y-%m-%d %H:%M:%S UTC                %%%%%%%%%%\n'))
    h.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    h.write('\\documentclass[a4paper,12pt]{article}\n')
    h.write('\\usepackage{graphicx}\n')
    h.write('\\usepackage{hyperref}\n')
    h.write('\\usepackage{longtable}\n')
    h.write('\\usepackage{setspace}\n')
    h.write('\\oddsidemargin=0.0cm\n')
    h.write('\\evensidemargin=0.0cm\n')
    h.write('\\newcommand{\\comment}[1]{\n\\begin{minipage}[t]{20ex}\n\\setstretch{0.5}\\flushleft\\noindent\n#1\n\\vspace*{1ex}\n\\end{minipage}}\n')
    h.write('\\newlength{\\openlooplen}\n')
    h.write('\\settowidth{\\openlooplen}{ooop}\n')
    h.write('\\newlength{\\cflen}\n')
    h.write('\\settowidth{\\cflen}{carbon}\n')
    h.write('\\newcommand{\\openloopheading}{\n\\begin{minipage}[t]{\\openlooplen}\nopen\\\\\nloop\n\\end{minipage}\n}\n')
    h.write('\\newcommand{\\cfheading}{\n\\begin{minipage}[t]{\\cflen}\ncarbon\\\\\nfibre\n\\end{minipage}\n}\n')

    cfibre_ack = 'The carbon fibre measurements are from Sophie Henrot Versill\\a\'e, see \\url{http://qubic.in2p3.fr/wiki/pmwiki.php/TD/P73TestWithACarbonFiberSource}.'
    openloop_ack = 'Results of the open loop test and the room temperature measurements are from Damien Pr\\^ele.'
    
    
    h.write('\\begin{document}\n')
    h.write('\\begin{center}\n')
    h.write('QUBIC TES Report\\\\\n')
    h.write('NEP estimates\\\\\n')
    h.write('data from %s\\\\\n' % datadate_str)
    h.write('compiled by %s\\\\\nusing QubicPack: \\url{https://github.com/satorchi/qubicpack}\\\\\n' % observer)
    h.write(dt.datetime.utcnow().strftime('this report compiled %Y-%m-%d %H:%M UTC\\\\\n'))
    h.write('\\end{center}\n')

    h.write('\\vspace*{3ex}\n')
    h.write('\\noindent Summary:\n')
    h.write('NEP are estimated according to the method described by \\mbox{C.~Perbost}, Ph.D, Section~6.3\n\n')
    h.write('\\noindent\\begin{itemize}\n')
    h.write('\\item Array %s\n' % detector_name)
    if asic is not None:h.write('\\item ASIC %i\n' % asic)
    h.write('\\item NEP estimated for %i TES out of %i\n' % (nNEP,len(NEPresults)))
    h.write('\\item average NEP=%.2f $\\times10^{-17}$ W / $\\sqrt{\\mathrm{Hz}}$\n' % (NEPmean*1e17))
    h.write('\\item data from:\n\\begin{itemize}\n')
    h.write(datelist)
    h.write('\n\\end{itemize}\n')
    h.write('\\end{itemize}\n')
    
    h.write('\n\\vspace*{3ex}\n\\noindent This document includes the following:\n')
    h.write('\\begin{itemize}\n')
    h.write('\\item Histograms of NEP, T$_0$, $n$, and $G$ values for the array\n')
    h.write('\\item Summary Table of NEP estimate for each TES\n')
    h.write('\\item Plot of I-V curves at the different bath temperatures\n')
    h.write('\\item Plot of P-V curves at the different bath temperatures\n')
    h.write('\\item Plot of R-V curves at the different bath temperatures\n')
    h.write('\\item Plot of turnover voltage at the different bath temperatures\n')
    h.write('\\item Plot of P-Temperature curves with NEP estimate\n')
    h.write('\\end{itemize}\n\\clearpage\n')

    reshisto = plot_NEP_histogram(NEPresults=NEPresults,xwin=False)
    histoparms = ['NEP','T0','n','G']
    for parm in histoparms:
        png = reshisto['pngname %s' % parm]
        if os.path.exists(png):
            h.write('\n\n\\noindent\\includegraphics[width=0.9\\linewidth,clip]{%s}' % png)
    h.write('\n\\clearpage\n')
       
    nrows = nNEP
    if show_extra_columns:
        colfmt='|r|r|r|r|r|r|l|l|l|r|'
        headline = '\\multicolumn{1}{|c|}{ASIC} & '\
            '\\multicolumn{1}{|c|}{TES} & '\
            '\\multicolumn{1}{|c|}{pix} & '\
            '\\multicolumn{1}{c|}{V$_\\mathrm{turnover}$} & '\
            '\\multicolumn{1}{c|}{R$_1$} & '\
            '\\multicolumn{1}{c|}{R$_\\mathrm{300K}$} & '\
            '\\multicolumn{1}{c|}{\\openloopheading} &'\
            '\\multicolumn{1}{c|}{\\cfheading} &'\
            '\\multicolumn{1}{c|}{comment} & '\
            '\\multicolumn{1}{c|}{NEP}'
    else:
        colfmt='|r|r|r|r|r|l|r|'
        headline = '\\multicolumn{1}{|c|}{ASIC} & '\
            '\\multicolumn{1}{|c|}{TES} & '\
            '\\multicolumn{1}{|c|}{pix} & '\
            '\\multicolumn{1}{c|}{V$_\\mathrm{turnover}$} & '\
            '\\multicolumn{1}{c|}{R$_1$} & '\
            '\\multicolumn{1}{c|}{comment} & '\
            '\\multicolumn{1}{c|}{NEP}'
    h.write('\\noindent\\begin{longtable}{%s}\n' % colfmt)
    h.write('\\caption{Summary Table for TES\\\\\n')
    for asic in go300.keys():
        h.write('V$_\\mathrm{turnover}$ and R$_1$ taken from measurement at T$_\mathrm{bath}$=%.1fmK for ASIC %i\\\\'\
                % (1000*go300[asic].temperature,asic))
    if show_extra_columns:
        h.write(cfibre_ack)
        h.write('\\\\\n')
        h.write(openloop_ack)
    h.write('}\\\\\n\\hline\n')
    h.write(headline+'\\\\ \n')
    h.write('\\hline\\endhead\n')
    h.write('\\hline\\endfoot\n')
    for result in NEPresults:
        NEP = result['NEP']
        TES = result['TES']
        asic = result['ASIC']
        asic_str = '%i &' % asic
        rowstr = asic_str+go300[asic].iv_tex_table_entry(TES)
        if not result['is_good']:
            if result['comment']:
                rowstr = rowstr.replace('good',result['comment'])
            else:
                rowstr = rowstr.replace('good','bad')
        if NEP is None:
            rowstr += ' & - \\\\\n'
        else:
            rowstr += ' & %.2f \\\\\n' % (1e17*NEP)
            
        h.write(rowstr)
    h.write('\\hline\n')
    h.write('\\end{longtable}\n')
    h.write('\n\\clearpage')
        

    for result in NEPresults:
        TES = result['TES']
        asic = result['ASIC']
        h.write('\n\\clearpage')
        pngIV = 'QUBIC_Array-%s_TES%03i_ASIC%i_I-V_Temperatures.png' % (detector_name,TES,asic)
        if refresh or not os.path.exists(pngIV):
            pngIV = plot_TES_temperature_curves(fplist,TES,asic,plot='I',xwin=False)

        pngPV = 'QUBIC_Array-%s_TES%03i_ASIC%i_P-V_Temperatures.png' % (detector_name,TES,asic)
        if refresh or not os.path.exists(pngPV):
            pngPV = plot_TES_temperature_curves(fplist,TES,asic,plot='P',xwin=False)

        pngRP = 'QUBIC_Array-%s_TES%03i_ASIC%i_R-V_Temperatures.png' % (detector_name,TES,asic)
        if refresh or not os.path.exists(pngRP):
            pngRP = plot_TES_temperature_curves(fplist,TES,asic,plot='R',xwin=False)

        pngTurnover='QUBIC_Array-%s_TES%03i_ASIC%i_Turnover_Temperature.png' % (detector_name,TES,asic)
        if refresh or not os.path.exists(pngTurnover):
            pngTurnover = plot_TES_turnover_temperature(fplist,TES,asic,xwin=False)
            
        pngNEP = 'QUBIC_Array-%s_TES%03i_ASIC%i_NEP.png' % (detector_name,TES,asic)
        if refresh or not os.path.exists(pngNEP):
            res = plot_TES_NEP(result=result,xwin=False)
            pngNEP = res['pngname']
        
        pngFiles=[pngIV,pngPV,pngRP,pngTurnover,pngNEP]
        for png in pngFiles:
            if png is not None and os.path.exists(png):
                h.write('\n\\noindent\\includegraphics[width=0.7\\linewidth,clip]{%s}\\\\' % png)
    
    h.write('\n\n\\end{document}\n')
    h.close()
    return texfilename

def rt_analysis(fplist,TES,asic,xwin=True):
    '''
    do the analysis to find the critical temperature (resistance vs. temperature)
    datlist is a list of qubicpack.qubicfp objects containing timeline data
    '''
    if not verify_temperature_arguments(fplist,1,asic):return None
    asic_idx = asic - 1

    # get all the results, including for multiple timelines in a single qp object
    detector_name=None
    reslist=[]
    for fp in fplist:
        go = fp.asic_list[asic_idx]
        if not go.exist_timeline_data():continue
        if detector_name is None:detector_name=go.detector_name
            
        ntimelines=go.ntimelines()
        for idx in range(ntimelines):
            res=go.plot_timeline(TES,timeline_index=idx,fit=True,xwin=False)
            reslist.append(res)

    plot_rt_analysis(reslist,xwin)
    return reslist

def plot_rt_analysis(reslist,xwin=True):
    '''
    plot the results of the R-T analysis
    we need to run fit_timeline() for this!
    '''
    TES=reslist[0]['TES']
    detector_name=reslist[0]['DET_NAME']
    asic=reslist[0]['ASIC']

    PIX=tes2pix(TES,asic)

    Tbath=[]
    R=[]
    dates=[]
    for res in reslist:
        if res['R amplitude'] is not None:
            Tbath.append(1000*res['Tbath'])
            R.append(res['R amplitude'])
            if isinstance(res['date'],list):
                date = res['date'][0]
            else:
                date = res['date']
            dates.append(date)
                     
    ntemps=len(Tbath)
    sorted_index=sorted(range(ntemps), key=lambda i: Tbath[i])

    Tsorted=np.array(Tbath)[sorted_index]
    Rsorted=np.array(R)[sorted_index]

    # show the dates of the measurements (just the day, not the time)
    date_txt='Measurements taken:'
    dates.sort()
    d_prev=''
    for d in dates:
        d_str=d.strftime('%Y-%m-%d')
        if d_str!=d_prev: date_txt+='\n   %s' % d_str
        d_prev=d_str
    
    ttl='Array %s: Critical Temperature for TES %i (PIX %i) on ASIC %i' % (detector_name,TES,PIX,asic)
    pngname='QUBIC_Array-%s_TES%03i_ASIC%i_Tcritical.png' % (detector_name,TES,asic)

    if xwin: plt.ion()
    else: plt.ioff()
    fig = plt.figure()
    figure_window_title(fig,ttl)
    ax = plt.gca()
    plt.title(ttl)
    ax.plot(Tsorted,1e6*Rsorted,ls='none',marker='D',color='blue')
    ax.set_xlabel('T$_\mathrm{bath}$ / mK')
    ax.set_ylabel('R$_\mathrm{TES}$ / $m\Omega$')

    ax.text(0.98,0.02,date_txt,va='bottom',ha='right',fontsize=plt.rcParams['legend.fontsize'],transform=ax.transAxes,bbox=boxprops)

    plt.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    if xwin:plt.show()
    else: plt.close('all')

    retval = {}
    retval['Tbath'] = Tsorted
    retval['R'] = Rsorted
    return retval

