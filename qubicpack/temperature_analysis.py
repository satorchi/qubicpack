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
from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
import pickle
from math import sqrt
from glob import glob
import pickle
import os
import datetime as dt
from scipy.optimize import curve_fit

from qubicpack.qubicfp import qubicfp
from qubicpack.pix2tes import tes2pix
from qubicpack.utilities import FIGSIZE, asic_reversal_date, NPIXELS, NASIC, TES_index, ASIC_index

# some constants and values required
kBoltzmann=1.3806485279e-23
temperature_precision = 0.005 # close enough for temperature

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
                T=go.tdata[t_idx]['TES_TEMP']
                if 'BEG-OBS' in go.tdata[t_idx].keys():
                    d=go.tdata[t_idx]['BEG-OBS']
                else:
                    d=go.tdata[t_idx]['DATE-OBS']
                if not obsdate is None and d!=obsdate:
                    printit=False
                if not temperature is None and not (T<temperature+0.001 and T>temperature-0.001):
                    printit=False
                if printit: print(datstr % (idx,t_idx,go.detector_name,go.asic,d,1000*T))
        else:
            T=go.temperature
            d=go.obsdate
            t_idx=0
            if obsdate is not None and d!=obsdate:
                printit=False
            if temperature is not None and not (T<temperature+0.001 and T>temperature-0.001):
                printit=False
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

    if asic<1 or asic>NASIC:
        print('ERROR! Please enter a valid asic number: %i' % asic)
    asic_idx = asic - 1

    if TES<1 or TES>NPIXELS:
        print('ERROR! Please enter a valid TES number: %i' % TES)
        
    # make sure we have a list of qubicpack.qubicfp objects
    if not isinstance(fplist, list):
        print('ERROR!  Please provide a list of qubicpack.qubicfp objects')
        return False

    if not isinstance(fplist[0], qubicfp):
        print('ERROR!  Please provide a list of qubicpack.qubicfp objects')
        return False

    detector_name=fplist[0].asic_list[asic_idx].detector_name
    for fp in fplist:
        if not isinstance(fp,qubicfp):
            print('ERROR! The list should contain qubicpack.qubicfp objects')
            return False
        if TES<1 or TES>NPIXELS:
            print('ERROR! Please enter a valid TES number between 1 and %i.' % NPIXELS)
            return False
        if fp.asic_list[asic_idx] is None:
            print('ERROR! There is no data for the requested ASIC.')
            return False
        if fp.asic_list[asic_idx].detector_name != detector_name:
            print('ERROR! These data are not for the same detector array.')
            return False
        

    return True

def plot_TES_turnover_temperature(fplist,TES,asic,xwin=True):
    '''
    plot the turnover point as a function of temperature for a given TES
    '''
    if not verify_temperature_arguments(fplist,TES,asic):return None
    asic_idx = asic - 1
    detector_name=fplist[0].asic_list[asic_idx].detector_name
    
    temps_list=[]
    turnover_list=[]
    for fp in fplist:
        go = fp.asic_list[asic_idx]
        if not go.turnover(TES) is None:
            temps_list.append(go.temperature)
            turnover_list.append(go.turnover(TES))


    if len(turnover_list)==0:
        print('no turnover for TES %i on ASIC %i' % (TES,asic))
        return None
    
    temps_list=np.array(temps_list)
    turnover_list=np.array(turnover_list)

    sorted_index=sorted(range(len(temps_list)), key=lambda i: temps_list[i])
    sorted_temps=temps_list[sorted_index]
    sorted_turnover=turnover_list[sorted_index]
    
    pngname='QUBIC_Array-%s_TES%03i_ASIC%i_Turnover_Temperature.png' % (detector_name,TES,asic)
    xlabel='T$_{bath}$ / mK'
    ylabel='V$_{turnover}$ / V'

    ttl='QUBIC Array %s, ASIC %i, TES #%i: Turnover at Different Temperatures' % (detector_name,asic,TES)
    subttl=fplist[0].asic_list[asic_idx].obsdate.strftime('Measurements of %Y-%m-%d')
    
    if xwin:plt.ion()
    else:
        plt.close('all')
        plt.ioff()
        
    fig=plt.figure(figsize=FIGSIZE)
    fig.canvas.set_window_title('plt: '+ttl)

    plt.suptitle(ttl+'\n'+subttl)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.plot(temps_list,turnover_list,linestyle='none',marker='D')
    plt.plot(sorted_temps,sorted_turnover,color='green')

    xmax=temps_list.max()
    xmin=temps_list.min()
    
    span=xmax-xmin
    plot_xlim=(xmin-0.05*span,xmax+0.1*span)
    
    ax=plt.gca()
    ax.set_xlim(plot_xlim)
    
    plt.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    if xwin: plt.show()
    else: plt.close('all')
    
    return


def plot_TES_temperature_curves(fplist,TES,asic,plot='I',xwin=True):
    '''
    plot the I-V, P-V, R-P curves for each temperature
    '''
    if not verify_temperature_arguments(fplist,TES,asic):return None
    asic_idx = asic - 1
    detector_name = fplist[0].asic_list[asic_idx].detector_name

    temps_list=[]
    for fp in fplist:
        go = fp.asic_list[asic_idx]
        temps_list.append(go.temperature)
        
    temps_list=np.array(temps_list)
    sorted_index=sorted(range(len(temps_list)), key=lambda i: temps_list[i])
    sorted_temps=temps_list[sorted_index]

    plot_type='I'
    if plot.upper()[0]=='R':
        plot_type='R'
        pngname='QUBIC_Array-%s_TES%03i_ASIC%i_R-V_Temperatures.png' % (detector_name,TES,asic)
        xlabel='P$_{TES}$ / $p$W'
        ylabel='$\\frac{R_\mathrm{TES}}{R_\mathrm{normal}}$ / %'
    elif plot.upper()[0]=='P':
        plot_type='P'
        pngname='QUBIC_Array-%s_TES%03i_ASIC%i_P-V_Temperatures.png' % (detector_name,TES,asic)
        xlabel='V$_{bias}$ / V'
        ylabel='P$_{TES}$ / $p$W'
    else:
        plot_type='I'
        pngname='QUBIC_Array-%s_TES%03i_ASIC%i_I-V_Temperatures.png' % (detector_name,TES,asic)
        xlabel='V$_{bias}$ / V'
        ylabel='I$_{TES}$ / $\mu$A'
        
    ttl='QUBIC Array %s ASIC %i TES%03i at Different Temperatures' % (detector_name,asic,TES)
    if xwin:plt.ion()
    else:
        plt.close('all')
        plt.ioff()
        
    fig=plt.figure(figsize=FIGSIZE)
    fig.canvas.set_window_title('plt: '+ttl)

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
        go = fplist[idx].asic_list[asic_idx]
        lbl = '%.0f mK' % (1000*go.temperature)

        startend = go.selected_iv_curve(TES)
        if startend is None:
            print('ERROR! index=%i, Tbath=%s: Could not get selected curve for TES#%i' % (idx,lbl,TES))
            continue
        istart,iend = startend

        Iadjusted=go.adjusted_iv(TES)
        I=Iadjusted[istart:iend]
        bias=go.vbias[istart:iend]


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
                plt.plot(Ptes,Rn_ratio,label=lbl)
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
        # draw the line at 90% Rn
        plt.plot([xmin,xmax],[90,90],color='black',linestyle='dashed',linewidth=3)
        
    else:
        xmin=min_bias
        xmax=max_bias
        
    span=xmax-xmin
    plot_xlim=(xmin-0.05*span,xmax+0.35*span)
    
    ax=plt.gca()
    ax.set_xlim(plot_xlim)
        
    plt.legend(fontsize=10)
    plt.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    if xwin: plt.show()
    else: plt.close('all')
    
    return


def P_bath_function(Tbath,K,T0,n):
    '''
    TES power should follow this as a function of bath temperature
    see Perbost 2016, PhD thesis, eq. 2.11 and 6.16
    '''
    P=K*(T0**n - Tbath**n)
    return P

def fit_Pbath(T_pts, P_pts, p0=None,ftol=1e-8):
    '''
    find best fit to the P_bath function
    '''
    if p0 is None: p0=np.array([1E-10, 0.5,4.5]) #MP

    # make sure T_pts and P_pts are 1d arrays
    npts = T_pts.size
    T=np.array(T_pts).reshape(npts)
    P=np.array(P_pts).reshape(npts)

    # try to fit the curve with the default tolerance (1e-8)
    # if unsuccessful, relax the tolerance by a factor 100 and try again
    for ctr in range(5):
        try:
            fit = curve_fit(P_bath_function,T,P,p0=p0,ftol=ftol)
            ret = {}
            ret['fit'] = fit
            ret['tolerance'] = ftol
            return ret
        except:
            ret=None
        ftol *= 100
        #print('retrying fit with relaxed tolerance: %.1e' % ftol)
        
    return ret

def calculate_TES_NEP(fplist,TES,asic,p0=None):
    '''
    make the list of temperatures and associated P0
    and calculate the NEP
    '''
    if not verify_temperature_arguments(fplist,TES,asic):return None
    asic_idx = asic - 1
    detector_name=fplist[0].asic_list[asic_idx].detector_name

    temps_list=[]
    for fp in fplist:
        go = fp.asic_list[asic_idx]
        temps_list.append(go.temperature)    
    temps_list=np.array(temps_list)
    sorted_index=sorted(range(len(temps_list)), key=lambda i: temps_list[i])
    sorted_temps=temps_list[sorted_index]
    
    ret={}
    ret['DET_NAME'] = detector_name
    ret['TES'] = TES
    ret['ASIC'] = asic
    ret['Tmin'] = temps_list.min()
    ret['Tmax'] = temps_list.max()

    # make the arrays of Power and T_bath
    P = []
    T = []
    all_T = []
    all_P = []
    for idx in sorted_index:
        fp = fplist[idx]
        go = fp.asic_list[asic_idx]
        if go.turnover(TES) is None or go.turnover(TES)<0.0:continue
        
        filterinfo = go.filterinfo(TES)
        
        istart,iend = go.selected_iv_curve(TES)

        Iadjusted = go.adjusted_iv(TES)
        I = Iadjusted[istart:iend]
        if 'Iturnover' in filterinfo['fit'].keys():
            Iturnover = 1e-6*filterinfo['fit']['Iturnover']
        else:
            Iturnover = 1e-6*I.min()

        Tbath = go.temperature
        Ptes = go.Ptes(TES)[istart:iend] #MP
        Vtes_turnover = go.Rshunt*(go.turnover(TES)/go.Rbias-Iturnover)
        all_P.append(Ptes.mean()*1e-12)
        all_T.append(Tbath)
        Pbeg = np.mean(Ptes[0:10])
        if ((Pbeg > 5) and (Pbeg < 40)):
            P.append(Pbeg*1e-12)
            T.append(Tbath)

    P = np.array(P)
    T = np.array(T)
    all_P = np.array(all_P)
    all_T = np.array(all_T)
    ret['P'] = P
    ret['T'] = T
    ret['all temperatures'] = all_T
    ret['all P'] = all_P
    
    temperature_fit = fit_Pbath(T,P,p0=p0)
    ret['fit points'] = 'filter by first points'
    if temperature_fit is None:
        # try again with the full range
        temperature_fit = fit_Pbath(all_T,all_P,p0=p0)
        ret['fit points'] = 'all'

    
    if temperature_fit is None:
        ret['K'] = None
        ret['T0'] = None
        ret['n'] = None
        ret['NEP'] = None
        ret['G'] = None
        ret['gamma'] = None
        print ('insufficient data for curve fit:  TES=%i' % TES)
        return ret
    
        
    K = temperature_fit['fit'][0][0]
    T0 = temperature_fit['fit'][0][1]
    n = temperature_fit['fit'][0][2]
    ret['K'] = K
    ret['T0'] = T0
    ret['n'] = n
    ret['fit tolerance'] = temperature_fit['tolerance']

    G = n*K*(T0**(n-1))
    ret['G'] = G
    Tratio = 0.35/T0
    # gamma is defined in Perbost PhD eq. 2.72 (page 82)
    gamma = (n/(2*n+1)) * (1-Tratio**(2*n+1))/(1-Tratio**n)
    ret['gamma'] = gamma
    discr = gamma*kBoltzmann*G
    if discr<0.0:
        print('ERROR! Imaginary NEP!  TES=%i' % TES)
        NEP = -2*T0*sqrt(-discr)
    else:
        NEP = 2*T0*sqrt(discr)
    ret['NEP'] = NEP    


    return ret

def make_TES_NEP_resultslist(fplist,asic):
    '''
    make a list of NEP calculation results, one for each TES
    '''
    if not verify_temperature_arguments(fplist,1,asic):return None
    results=[]
    for idx in range(NPIXELS):
        TES = 1 + idx
        res = calculate_TES_NEP(fplist,TES,asic)
        results.append(res)
            
    return results

def plot_TES_NEP(fplist,TES,asic,xwin=True,p0=None):
    '''
    plot the P vs. Temperature for a TES
    '''

    result = calculate_TES_NEP(fplist,TES,asic,p0=p0)
    if result is None:return None

    TES = result['TES']
    asic = result['ASIC']
    detector_name = result['DET_NAME']

    all_T = result['all temperatures']
    all_P = result['all P']
    P = result['P']
    T = result['T']

    NEP = result['NEP']
    K = result['K']
    T0 = result['T0']
    n = result['n']
    G = result['G']
    
    Tmin = result['Tmin']
    Tmax = result['Tmax']
    T_span = Tmax-Tmin
    # add 30% to the plot edges
    plot_T_min = Tmin - 0.3*T_span
    plot_T_max = Tmax + 0.3*T_span
    T_stepsize = 1.1*T_span/100

    if NEP is None:
        txt='NEP estimate is not possible'
        if len(P)==1:
            plot_P_min=P[0]-0.2*P[0]
            plot_P_max=P[0]+0.2*P[0]
        elif len(P)==0:
            plot_P_min=0.0
            plot_P_max=5e-11
        else:
            P_span=max(P)-min(P)
            plot_P_min=min(P)-0.1*P_span
            plot_P_max=max(P)+0.1*P_span
    else:
        fit_T=np.arange(plot_T_min,plot_T_max,T_stepsize)
        fit_P=P_bath_function(fit_T,K,T0,n)

        P_span=max(fit_P)-min(fit_P)
        plot_P_min=min(fit_P)-0.1*P_span
        plot_P_max=max(fit_P)+0.1*P_span
        
        txt='K=%.4e' % K
        txt+='\nT$_0$=%.1f mK' % (1000*T0)
        if T0<0.3:
            txt+=' ERROR!  T$_0$<300 mK'
        txt+='\nn=%.3f' % n
        txt+='\nG=%.4e' % G
        txt+='\nNEP=%.4e at T$_{bath}$=350mK' % NEP
            
    

    pngname='QUBIC_Array-%s_TES%03i_ASIC%i_NEP.png' % (detector_name,TES,asic)
    ttl='QUBIC Array %s, ASIC %i, TES #%i: NEP' % (detector_name,asic,TES)

    boxprops = {}
    boxprops['alpha'] = 0.8
    boxprops['color'] = 'wheat'
    boxprops['boxstyle'] = 'round'
    
    if xwin:plt.ion()
    else:
        plt.close('all')
        plt.ioff()
    fig=plt.figure(figsize=FIGSIZE)
    fig.canvas.set_window_title('plt: '+ttl)

    ax=plt.gca()
    ax.set_xlim(plot_T_min,plot_T_max)
    ax.set_ylim(plot_P_min,plot_P_max)
    plt.title(ttl)
    ax.set_xlabel('T$_\mathrm{bath}$ / K')
    ax.set_ylabel('Power / Watt')
    P_span=plot_P_max-plot_P_min
    text_y=plot_P_min+0.5*P_span
    if result['fit points'] == 'all':
        ax.plot(all_T,all_P,ls='none',marker='D')
    else:
        ax.plot(T,P,ls='none',marker='D')
    if not NEP is None: plt.plot(fit_T,fit_P,color='red')
    #ax.text(Tmin,text_y,txt,fontsize=14)
    ax.text(0.95,0.95,txt,va='top',ha='right',fontsize=14,transform=ax.transAxes,bbox=boxprops)
    plt.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    if xwin:plt.show()
    else:plt.close('all')
    return result

def plot_NEP_histogram(fplist,asic,NEPresults=None,xwin=True):
    '''
    plot the histogram of the NEP calculations
    '''
    if not verify_temperature_arguments(fplist,1,asic):return None
    asic_idx = asic - 1
    
    # find the data at 300mK
    go300=None
    for fp in fplist:
        go = fp.asic_list[asic_idx]
        if go300 is None and is_300mK_measurement(go):
            go300=go
    if go300 is None:
        go300=fplist[-1].asic_list[asic_idx]
        
    datadate=fplist[0].asic_list[asic_idx].obsdate
    detector_name=go300.detector_name

    # generate the results if not already done
    if NEPresults is None:NEPresults=make_TES_NEP_resultslist(qplist,asic)

    NEP_estimate=[]
    TESlist=[]
    for res in NEPresults:
        NEP=res['NEP']
        if not NEP is None:
            NEP_estimate.append(NEP)
            TESlist.append(res['TES'])
    nNEP=len(NEP_estimate)
    NEP_estimate=np.array(NEP_estimate)
    NEPmean=NEP_estimate.mean()
    txt='NEP$_\mathrm{mean}=%.4f \\times 10^{-17}W/\sqrt{\mathrm{Hz}}$' % (1e17*NEPmean)
    
    pngname='QUBIC_TES_ASIC%i_NEP_histogram.png' % asic
    ttl='QUBIC TES ASIC %i NEP' % asic
    ttl+='\nhistogram of %i TES\ndata taken %s' % (nNEP,datadate.strftime('%Y-%m-%d'))
    if xwin:plt.ion()
    else:
        plt.close('all')
        plt.ioff()
    fig=plt.figure(figsize=FIGSIZE)
    fig.canvas.set_window_title('plt: '+ttl)

    ax=plt.gca()
    plt.title(ttl)
    ax.set_xlabel('NEP  /  ${W}/\sqrt{Hz}$')
    ax.set_ylabel('Number per bin')

    hist,binedge=np.histogram(NEP_estimate,bins=10)
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

    ax.text(0.01,0.99,txt,ha='left',va='top',fontsize=16,transform=ax.transAxes)
    fig.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    if xwin:fig.show()
    else:plt.close('all')
    return
            
def make_TES_NEP_tex_report(fplist,asic,NEPresults=None,refresh=True):
    '''
    make a LaTeX source file for a test report of NEP estimates
    the input arguments are a list of qubicpack objects with the I-V data at different temperatures
    and a list of results from plot_TES_NEP()
    '''
    if not verify_temperature_arguments(fplist,1,asic):return None
    asic_idx = asic - 1
    
    # find the data at 300mK
    go300=None
    datelist=''
    for fp in fplist:
        go = fp.asic_list[asic_idx]
        datelist+='\n\\item %.3fmK on %s'\
            % (1000*go.temperature,go.obsdate.strftime('%Y-%m-%d %H:%M:%S'))
        if go300 is None and is_300mK_measurement(go):
            go300=go
    if go300 is None:
        go300=fplist[-1].asic_list[asic_idx]

    print('300mK results from %s, Tbath=%.1fmK'\
          % (go300.obsdate.strftime('%Y-%m-%d %H:%M:%S'),go300.temperature*1000))
    observer=go300.observer.replace('<','$<$').replace('>','$>$')
    detector_name=go300.detector_name
    show_extra_columns = True
    if go300.transdic is None:
        show_extra_columns = False
    
    # generate the plots if not already done
    if NEPresults is None:NEPresults=make_TES_NEP_resultslist(fplist,asic)

    NEP_estimate=[]
    TESlist=[]
    for res in NEPresults:
        NEP=res['NEP']
        if not NEP is None:
            NEP_estimate.append(NEP)
            TESlist.append(res['TES'])
    nNEP=len(TESlist)
    NEP_estimate=np.array(NEP_estimate)
    NEPmean=NEP_estimate.mean()
    
    texfilename=str('QUBIC_Array-%s_ASIC%i_NEP.tex' % (detector_name,asic))
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
    h.write(go300.obsdate.strftime('data from %Y-%m-%d\\\\\n'))
    h.write('compiled by %s\\\\\nusing QubicPack: \\url{https://github.com/satorchi/qubicpack}\\\\\n' % observer)
    h.write(dt.datetime.utcnow().strftime('this report compiled %Y-%m-%d %H:%M UTC\\\\\n'))
    h.write('\\end{center}\n')

    h.write('\\vspace*{3ex}\n')
    h.write('\\noindent Summary:\n')
    h.write('NEP are estimated according to the method described by C.~Perbost, Ph.D, Section~6.3\n\n')
    h.write('\\noindent\\begin{itemize}\n')
    h.write('\\item Array %s\n' % detector_name)
    h.write('\\item ASIC %i\n' % asic)
    h.write('\\item NEP estimated for %i TES out of %i\n' % (nNEP,len(NEPresults)))
    h.write('\\item average NEP=%.2f $\\times10^{-17}$ W / $\\sqrt{\\mathrm{Hz}}$\n' % (NEPmean*1e17))
    h.write('\\item data from:\n\\begin{itemize}\n')
    h.write(datelist)
    h.write('\n\\end{itemize}\n')
    h.write('\\end{itemize}\n')
    
    h.write('\n\\vspace*{3ex}\n\\noindent This document includes the following:\n')
    h.write('\\begin{itemize}\n')
    h.write('\\item Summary Table of NEP estimate for each TES, where the data permits\n')
    h.write('\\item Histogram of NEP values for the array\n')
    h.write('\\item Plot of I-V curves at the different bath temperatures\n')
    h.write('\\item Plot of P-V curves at the different bath temperatures\n')
    h.write('\\item Plot of P-Temperature curves with NEP estimate, where possible\n')
    h.write('\\end{itemize}\n\\clearpage\n')

    png='QUBIC_TES_ASIC%i_NEP_histogram.png' % asic
    if refresh or not os.path.exists(png):
        plot_NEP_histogram(fplist,asic,NEPresults,xwin=False)                           

    if os.path.exists(png):
       h.write('\n\\noindent\\includegraphics[width=0.9\\linewidth,clip]{%s}' % png)
       h.write('\n\\clearpage\n')
       
    nrows=nNEP
    if show_extra_columns:
        colfmt='|r|r|r|r|r|l|l|l|r|'
        headline='\\multicolumn{1}{|c|}{TES} & '\
            '\\multicolumn{1}{|c|}{pix} & '\
            '\\multicolumn{1}{c|}{V$_\\mathrm{turnover}$} & '\
            '\\multicolumn{1}{c|}{R$_1$} & '\
            '\\multicolumn{1}{c|}{R$_\\mathrm{300K}$} & '\
            '\\multicolumn{1}{c|}{\\openloopheading} &'\
            '\\multicolumn{1}{c|}{\\cfheading} &'\
            '\\multicolumn{1}{c|}{comment} & '\
            '\\multicolumn{1}{c|}{NEP}'
    else:
        colfmt='|r|r|r|r|l|r|'
        headline='\\multicolumn{1}{|c|}{TES} & '\
            '\\multicolumn{1}{|c|}{pix} & '\
            '\\multicolumn{1}{c|}{V$_\\mathrm{turnover}$} & '\
            '\\multicolumn{1}{c|}{R$_1$} & '\
            '\\multicolumn{1}{c|}{comment} & '\
            '\\multicolumn{1}{c|}{NEP}'
    h.write('\\noindent\\begin{longtable}{%s}\n' % colfmt)
    h.write('\\caption{Summary Table for TES\\\\\n')
    if show_extra_columns:
        h.write(cfibre_ack)
        h.write('\\\\\n')
        h.write(openloop_ack)
    h.write('}\\\\\n\\hline\n')
    h.write(headline+'\\\\ \n')
    h.write('\\hline\\endhead\n')
    h.write('\\hline\\endfoot\n')
    for result in NEPresults:
        NEP=result['NEP']
        TES=result['TES']
        if not NEP is None:
            rowstr=go300.iv_tex_table_entry(TES)
            rowstr+=' & %.2f \\\\\n' % (1e17*NEP)
            h.write(rowstr)
    h.write('\\hline\n')
    h.write('\\end{longtable}\n')
    h.write('\n\\clearpage')
        

    for TES in np.arange(1,129):
        h.write('\n\\clearpage')
        pngIV ='QUBIC_Array-%s_TES%03i_ASIC%i_I-V_Temperatures.png' % (detector_name,TES,asic)
        if refresh or not os.path.exists(pngIV):
            res=plot_TES_temperature_curves(fplist,TES,asic,plot='I',xwin=False)

        pngPV ='QUBIC_Array-%s_TES%03i_ASIC%i_P-V_Temperatures.png' % (detector_name,TES,asic)
        if refresh or not os.path.exists(pngPV):
            res=plot_TES_temperature_curves(fplist,TES,asic,plot='P',xwin=False)

        pngRP ='QUBIC_Array-%s_TES%03i_ASIC%i_R-V_Temperatures.png' % (detector_name,TES,asic)
        if refresh or not os.path.exists(pngRP):
            res=plot_TES_temperature_curves(fplist,TES,asic,plot='R',xwin=False)

        pngTurnover='QUBIC_Array-%s_TES%03i_ASIC%i_Turnover_Temperature.png' % (detector_name,TES,asic)
        if refresh or not os.path.exists(pngTurnover):
            res=plot_TES_turnover_temperature(fplist,TES,asic,xwin=False)
            
        pngNEP='QUBIC_Array-%s_TES%03i_ASIC%i_NEP.png' % (detector_name,TES,asic)
        if refresh or not os.path.exists(pngNEP):
            res=plot_TES_NEP(fplist,TES,asic,xwin=False)

        
        pngFiles=[pngIV,pngPV,pngRP,pngTurnover,pngNEP]
        for png in pngFiles:
            if os.path.exists(png):
                h.write('\n\\noindent\\includegraphics[width=0.7\\linewidth,clip]{%s}\\\\' % png)
    
    h.write('\n\n\\end{document}\n')
    h.close()
    return texfilename

def rt_analysis(datlist,TES,asic,xwin=True):
    '''
    do the analysis to find the critical temperature (resistance vs. temperature)
    datlist is a list of qubicpack.qubicfp objects containing timeline data
    '''
    if not verify_temperature_arguments(fplist,1,asic):return None
    asic_idx = asic - 1

    # get all the results, including for multiple timelines in a single qp object
    detector_name=None
    reslist=[]
    for fs in datlist:
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
        Tbath.append(1000*res['Tbath'])
        R.append(res['R amplitude'])
        dates.append(res['date'])
                     
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
    fig=plt.figure(figsize=FIGSIZE)
    fig.canvas.set_window_title('plt: '+ttl)
    ax=plt.gca()
    plt.title(ttl)
    ax.plot(Tsorted,1e6*Rsorted,marker='D',color='blue')
    ax.set_xlabel('T$_\mathrm{bath}$ / mK')
    ax.set_ylabel('R$_\mathrm{TES}$ / $\mu\Omega$')

    text_x=0.98
    text_y=0.02
    boxprops = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(text_x,text_y,date_txt,va='bottom',ha='right',fontsize=10,transform = ax.transAxes,bbox=boxprops)

    plt.savefig(pngname,format='png',dpi=100,bbox_inches='tight')
    if xwin:plt.show()
    else: plt.close('all')
    return
