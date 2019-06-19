#!/usr/bin/env python
"""
$Id: iv.py

$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Wed 05 Jul 2017 14:39:42 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

methods to plot and analyse I-V curves

"""
from __future__ import division, print_function
import numpy as np
import sys,os,time
import datetime as dt
import matplotlib.pyplot as plt
from glob import glob
import math
import pickle
from scipy.optimize import curve_fit
#from PIL import Image

from qubicpack.utilities import TES_index
from qubicpack.pix2tes import assign_pix2tes,pix2tes,tes2pix
from qubicpack.plot_fp import plot_fp

def exist_iv_data(self):
    '''
    check if we have I-V data
    '''
    if not isinstance(self.adu,np.ndarray):
        self.printmsg('No I-V data!',verbosity=2)
        return False

    if not isinstance(self.vbias,np.ndarray):
        self.printmsg('No bias data!',verbosity=2)
        return False

    if self.vbias.shape[0]!=self.adu.shape[1]:
        self.printmsg("I and V don't match! npts V = %i, npts I = %i" % (self.vbias.shape[0],self.adu.shape[1]))
        return False
    return True

def wait_a_bit(self,pausetime=None):
    if pausetime is None:
        pausetime=self.pausetime
    else:
        self.assign_pausetime(pausetime)
        
    print("waiting %.3f seconds" % pausetime)
    time.sleep(pausetime)
    return

def lut(self,v,vmin=3.0,vmax=9.0):
    '''
    a colour look up table for showing Vturnover in the I-V plots
    '''
    vfractional=(v-vmin)/(vmax-vmin)
    colourmap = plt.cm.get_cmap('Spectral_r')
    rgb=colourmap(vfractional)
    return rgb

def n_masked(self):
    '''
    find the number of masked samples per integration
    '''
    if self.rawmask is None:
        rawmask_added_date=dt.datetime.strptime('2018-02-14 11:30','%Y-%m-%d %H:%M')
        if self.obsdate < rawmask_added_date: return 8
        return 0

    n_masked=0
    for val in self.rawmask:
        for idx in range(8):
            masked = val & 1
            n_masked+=masked
            val = val >> 1
    return n_masked
        

def ADU2I(self,ADU, offset=None, R1adjust=1.0):
    ''' 
    This is the magic formula to convert the measured output of the TES to current
    the voltage (ADU) returned by the TES is converted to a current in uA        
    '''
    q_ADC = 20./(2**16-1)
    G_FLL = (10.4 / 0.2) * self.Rfeedback
    n_masked=self.n_masked()

    # since QubicStudio 3.5 (2019-01-29), and perhaps before, the masking is taken into account by QubicStudio
    # I = 1e6 * (ADU / 2**7) * (q_ADC/G_FLL) * (self.nsamples - n_masked) * R1adjust
    I = 1e6 * (ADU / 2**7) * (q_ADC/G_FLL) * R1adjust
    if not offset is None: return I+offset
    return I

def setup_plot_Vavg(self,axes=None):
    ttl=str('Average Current per TES with different Vbias')
    plt.ion()
    fig=plt.figure(figsize=self.figsize)
    fig.canvas.set_window_title('plt: '+ttl) 
    fig.suptitle(ttl,fontsize=16)
    if isinstance(axes,list) or isinstance(axes,np.ndarray): plt.axis(axes)
    plt.xlabel('TES number')
    plt.ylabel('I  /  $\mu$A')
    return fig

def plot_Vavg(self,Vavg,Vbias,offset=None,axes=None):
    Iavg=self.ADU2I(Vavg,offset)
    
    lbl=str('V$_{bias}$ = %.2fV' % Vbias)
    plt.cla()
    if isinstance(axes,list) or isinstance(axes,np.ndarray): plt.axis(axes)
    plt.xlabel('TES number')
    plt.ylabel('I  /  $\mu$A')
    # plot markers with no lines
    plt.plot(Iavg,marker='D',drawstyle='steps-mid',linestyle='none',color='green',label=lbl)
    # plot bars up to the markers
    tes_axis=np.arange(self.NPIXELS)-0.25
    plt.bar(tes_axis,height=Iavg,color='pink',width=0.5)
    plt.legend()
    plt.pause(0.01)
    return

def plot_iv_all(self,selection=None,xwin=True):
    if not isinstance(self.vbias,np.ndarray):
        self.vbias=make_Vbias()
    ttl=str('QUBIC I-V curves (%s)' % (self.obsdate.strftime('%Y-%b-%d %H:%M UTC')))
    if isinstance(selection,list):
        nselection=0
        for val in selection:
            if val: nselection+=1
    else:
        nselection=self.NPIXELS
                
    subttl=str('plotting curves for %i TES out of %i' % (nselection,self.NPIXELS))
    if xwin: plt.ion()
    else: plt.ioff()
    fig=plt.figure(figsize=self.figsize)
    if xwin: fig.canvas.set_window_title('plt: '+ttl) 
    fig.suptitle(ttl+'\n'+subttl,fontsize=16)
    plt.xlabel('Bias Voltage  /  V')
    plt.ylabel('Current  /  $\mu$A')

    nbias=self.adu.shape[1]
    
    offset=[]
    colour_idx=0
    ncolours=len(self.colours)
    for TES_idx in range(self.NPIXELS):
        TES=TES_idx+1

        if (not isinstance(selection,list)) or (selection[TES_idx]):
            istart,iend=self.selected_iv_curve(TES)
            Iadjusted=self.adjusted_iv(TES)[istart:iend]
            bias=self.vbias[istart:iend]
            if colour_idx >= ncolours:colour_idx=0
            plt.plot(bias,Iadjusted,color=self.colours[colour_idx])
            colour_idx+=1

    pngname=str('TES_IV_array-%s_ASIC%i_all_%s.png' % (self.detector_name,self.asic,self.obsdate.strftime('%Y%m%dT%H%M%SUTC')))
    pngname_fullpath=self.output_filename(pngname)
    if isinstance(pngname_fullpath,str): plt.savefig(pngname_fullpath,format='png',dpi=100,bbox_inches='tight')
    if xwin: plt.show()
    else: plt.close('all')
    return fig

def setup_plot_iv_multi(self,nrows=16,ncols=8,xwin=True):
    if not isinstance(self.vbias,np.ndarray): self.vbias=make_Vbias()
    ttl=str('QUBIC I-V curves (%s)' % (self.obsdate.strftime('%Y-%b-%d %H:%M UTC')))

    nbad=0
    for val in self.is_good_iv():
        if not val:nbad+=1
    ttl+=str('\n%i flagged as bad pixels' % nbad)
    
    if xwin: plt.ion()
    else: plt.ioff()
    fig,axes=plt.subplots(nrows,ncols,sharex=True,sharey=False,figsize=self.figsize)
    if xwin: fig.canvas.set_window_title('plt: '+ttl)
    fig.suptitle(ttl,fontsize=16)
    plt.xlabel('Bias Voltage  /  V')
    plt.ylabel('Current  /  $\mu$A')
    return fig,axes

def plot_iv_multi(self, xwin=True):
    '''
    plot all TES I-V curves on a grid
    the optional list "is_good" is boolean for each TES
    if present, a label will be shown on each curve indicating 
    whether that TES is considered good or not
    '''
    ngood=self.ngood()
    nrows=16
    ncols=8
    fig,axes=self.setup_plot_iv_multi(nrows,ncols,xwin)
    
    TES_idx=0
    for row in range(nrows):
        for col in range(ncols):
            TES=TES_idx+1
            
            axes[row,col].get_xaxis().set_visible(False)
            axes[row,col].get_yaxis().set_visible(False)

            Iadjusted=self.adjusted_iv(TES)
            self.draw_iv(Iadjusted,colour='blue',axis=axes[row,col])
            text_y=min(Iadjusted)
            axes[row,col].text(self.bias_factor*self.max_bias,text_y,str('%i' % (TES_idx+1)),va='bottom',ha='right',color='black')

            if (not self.is_good_iv() is None)\
               and (not self.is_good_iv()[TES_idx]):
                axes[row,col].set_facecolor('red')

            TES_idx+=1

    pngname=str('TES_IV_array-%s_ASIC%i_thumbnail_%s.png' % (self.detector_name,self.asic,self.obsdate.strftime('%Y%m%dT%H%M%SUTC')))
    pngname_fullpath=self.output_filename(pngname)
    if isinstance(pngname_fullpath,str): plt.savefig(pngname_fullpath,format='png',dpi=100,bbox_inches='tight')
    if xwin: plt.show()
    else: plt.close('all')
    
    return fig

def plot_iv_physical_layout(self,xwin=True):
    '''
    plot the I-V curves in thumbnails mapped to the physical location of each detector
    '''
    if not self.exist_iv_data():return None

    args= {}
    key = 'ASIC%i' % self.asic
    
    args['title'] = 'QUBIC Focal Plane I-V curves: %s' % self.obsdate.strftime('%Y-%b-%d %H:%M UTC')
    subttl_list = []
    subttl_list.append(self.infotext())
    args['obsdate'] = self.obsdate

    bias,adu = self.best_iv_curve()
    args[key] = adu

    keyx = '%s x-axis' % key
    args[keyx] = bias

    keygood = '%s good' % key
    args[keygood] = self.is_good_iv()

    keybg = '%s bg' % key
    args[keybg] = self.turnover()

    ngood = self.ngood()
    if ngood is not None:
        subttl_list.append('%i flagged as bad pixels : yield = %.1f%%' %
                           (self.NPIXELS-ngood,100.0*ngood/self.NPIXELS))
    args['subtitle'] = '\n'.join(subttl_list)
    
    pngname='TES_IV_array-%s_ASIC%i_%s.png' % (self.detector_name,self.asic,self.obsdate.strftime('%Y%m%dT%H%M%SUTC'))
    pngname_fullpath=self.output_filename(pngname)
    args['pngname'] = pngname_fullpath

    plot_fp(args)
    return args



def make_line(self,pt1,pt2,xmin,xmax):
    '''
    make a line with extents xmin,xmax which includes the two points pt1, and pt2
    y=mx+b
    '''
    m = (pt2[1]-pt1[1]) / (pt2[0]-pt1[0])
    b = pt1[1] - (m*pt1[0])

    ymin=m*xmin+b
    ymax=m*xmax+b
    print('straight line should go through the origin: b=%.5e =  0?' % b)
    return [ymin,ymax]

def draw_tangent(self,TES):
    '''
    make a tangent line of the I-V curve fit at the maximum bias
    '''
    R1=self.R1(TES)
    if R1 is None: return None,None

    offset=self.offset(TES)
    
    # tangent is determined by the fit
    slope=1/R1

    Iinfinity=self.Vinfinity
    # The line is described by: y=slope*x + I0 + offset
    I0=Iinfinity - slope*self.Vinfinity - offset

    V1=self.bias_factor*self.max_bias
    I1=slope*V1 + I0 + offset

    #V2=self.bias_factor*self.min_bias
    V2=0.0
    I2=slope*V2 + I0 + offset

    xpts=[V2,V1]
    ypts=[I2,I1]
    plt.plot(xpts,ypts,linestyle='dashed',color='green',label='model normal region')
    
    return I0

def fitted_iv_curve(self,TES):
    '''
    make a curve from the fit parameters
    '''
    filterinfo=self.filterinfo(TES)
    if filterinfo is None:return None

    offset=self.offset(TES)

    fit=filterinfo['fit']
    self.TES=TES # this is required for the "mixed" and "combined" models

    istart,iend=self.selected_iv_curve(TES)
    bias=self.bias_factor*self.vbias[istart:iend]

    # polynomial fit
    if 'fitfunction' not in fit.keys() or fit['fitfunction']=='POLYNOMIAL':
        func=np.poly1d(fit['fitinfo'][0]) + offset
        f=func(bias)
        return bias,f

    # combined polynomial fit
    Vsuper,Vnormal,a0,a1,b0,b1,b2,b3,c0,c1=fit['fitinfo'][0]
    f=self.model_iv_combined(bias,Vsuper,Vnormal,a0,a1,b0,b1,b2,b3,c0,c1) + offset
    return bias,f

def filter_jumps(self,I,jumplimit=2.0):
    '''
    filter out big jumps and return the remaining curve for fitting
    '''
    self.debugmsg('filtering jumps')
    npts_curve=len(I)

    # don't do anything if jumplimit not given
    if not (isinstance(jumplimit,float) or isinstance(jumplimit,int)):
        self.debugmsg('no jump filtering.  returning full range.')
        return (0,npts_curve-1)

    # find the step size between subsequent bias points
    # we make this relative to the lower point in the step
    # so that a big step is a big number
    steps=[]
    for idx in np.arange(1,npts_curve):
        stepsize=I[idx]-I[idx-1]
        if stepsize > self.zero:
            rel_stepsize = abs(stepsize/I[idx-1])
        elif stepsize < -self.zero:
            rel_stepsize = abs(stepsize/I[idx])
        else:
            rel_stepsize=1/self.zero
        steps.append(rel_stepsize)
        self.debugmsg('%i) stepsize: %.4f' % (idx,rel_stepsize))
    steps=np.array(steps)

    # find big steps 
    xpts=[]
    # don't forget to include the first and last point!
    xpts.append(0)
    meanval=steps.mean()
    for idx in range(len(steps)):
        sigma=steps[idx]/meanval
        self.debugmsg('%i) sigma stepsize: %.4f' % (idx+1,sigma))
        if sigma>jumplimit: xpts.append(idx)
    xpts.append(npts_curve-1)
    msg=''
    for idx in xpts:
        msg=str('%s %i,' % (msg,idx))
    msg='found big jumps at: '+msg
    self.debugmsg(msg)
    # we will return the largest span of non-jump data found in the curve
    maxspan=0
    maxspan_idx1=0
    maxspan_idx2=npts_curve
    for idx in range(len(xpts)-1):
        idx_next=idx+1
        span=abs(xpts[idx_next]-xpts[idx])
        if span>maxspan:
            maxspan=span
            maxspan_idx1=xpts[idx]
            maxspan_idx2=xpts[idx_next]

    # we only return the indexes of the good span
    self.debugmsg('best span of points in the curve is %i:%i' % (maxspan_idx1,maxspan_idx2))
    return maxspan_idx1,maxspan_idx2

def polynomial_fit_parameters(self,fit):
    '''
    determine the TES characteristics from the polynomial fit
    this is called from fit_iv()
    '''
    TES=fit['TES']
    TES_idx=TES_index(TES)
    R1adjust=self.R1adjust(TES)
    I=self.ADU2I(self.adu[TES_idx,:],R1adjust=R1adjust)
    npts=len(I)

    # Vinfinity is the virtual point where we calculate the offset
    
    # the coefficients of the polynomial fit
    if fit['fitfunction']=='POLYNOMIAL':
        # one polynomial was used for the entire I-V curve
        a3=fit['fitinfo'][0][0]
        a2=fit['fitinfo'][0][1]
        a1=fit['fitinfo'][0][2]
        a0=fit['fitinfo'][0][3]
    else:
        # the polynomial fit is only for the "mixed" region
        a0=fit['fitinfo'][0][4]
        a1=fit['fitinfo'][0][5]
        a2=fit['fitinfo'][0][6]
        a3=fit['fitinfo'][0][7]


    # find turning where polynomial tangent is zero (i.e. first derivative is zero)
    t1=-a2/(3*a3)
    discriminant=a2**2 - 3*a1*a3
    if discriminant<0.0:
        fit['turning']=[None,None]
        fit['concavity']=[None,None]
        fit['turnover']=None
        fit['Iturnover']=None
        fit['inflection']=None

    else:
        t2=np.sqrt(discriminant)/(3*a3)
        x0_0=t1+t2
        x0_1=t1-t2
        fit['turning']=[x0_0,x0_1]

        # check the concavity of the turning: up or down (+ve is up)
        fit['concavity']=[2*a2 + 6*a3*x0_0, 2*a2 + 6*a3*x0_1]

        
    # check if we have a valid turnover point within the range
    found_turnover=False
    idx=0
    for V0 in fit['turning']:
        concavity=fit['concavity'][idx]
        if not V0 is None:
            if (not concavity is None) and concavity>0:
                found_turnover=True
                fit['turnover']=V0
        idx+=1
        
    if not found_turnover:
        fit['turnover']=None
        fit['Iturnover']=None

    n_turnings_within_range=0
    for V0 in fit['turning']:
        if (not V0 is None) and V0>self.bias_factor*self.min_bias and V0<self.bias_factor*self.max_bias:
            n_turnings_within_range+=1     
    fit['turnings within range']=n_turnings_within_range

    # find the inflection point between the turning points
    inflection_V=-a2/(3*a3)
    fit['inflection']=inflection_V

    # if we're using the combined fit, then we can exit now
    if fit['fitfunction']=='COMBINED':return fit

    fit['Vsuper'] = None
    fit['Vnormal'] = None
    # the inflection point greater than turnover is the Vnormal
    if fit['turnover'] is not None and inflection_V>fit['turnover']:
        fit['Vnormal'] = inflection_V
    
    # if the inflection is between the turnover and the max bias,
    # then we fit a straight line to the final points
    # instead of using the fit all the way through
    if fit['Vnormal'] is not None and (inflection_V<self.bias_factor*self.max_bias):
        # find the corresponding points to fit
        istart=fit['curve index']*fit['npts_curve']
        iend=istart+fit['npts_curve']
        xpts=self.bias_factor*self.vbias[istart:iend]
        gotit=False
        dv=xpts[1]-xpts[0]
        idx=0
        while not gotit and idx<fit['npts_curve']:
            if (dv>0.0 and xpts[idx]>=inflection_V)\
               or (dv<0.0 and xpts[idx]<=inflection_V): 
                gotit=True
            idx+=1

        idx-=1
        if (dv>0.0):
            ibeg=istart+idx
            istop=iend
        else:
            ibeg=istart
            istop=istart+idx+1

        xpts=self.bias_factor*self.vbias[ibeg:istop]
        ypts=I[ibeg:istop]
        fit['linefit xpts']=xpts
        fit['linefit ypts']=ypts
        self.debugmsg('polynomial_fit_parameters(%i): ibeg=%i, istop=%i' % (TES,ibeg,istop))
        linefit=np.polyfit(xpts,ypts,1,full=True)
        slope=linefit[0][0]
        b=linefit[0][1]
        if abs(slope)>self.zero:
            R1=1/slope
        else:
            R1=1/self.zero
        # offset forces the line to have I(max_bias)=max_bias (i.e. R=1 Ohm)
        Imax=slope*self.Vinfinity + b
        offset=self.Vinfinity-Imax
        fit['R1']=R1
        self.debugmsg('polynomial_fit_parameters (2): setting R1 to %.3f' % R1)
        fit['offset']=offset

        # get the current, according to the model, at the region transitions: Vsuper, Vnormal
        Vsuper=fit['Vsuper']
        if Vsuper is not None:
            fit['Isuper']=a0 + a1*Vsuper + a2*Vsuper**2 + a3*Vsuper**3 + offset
        Vnormal=fit['Vnormal']
        if Vnormal is not None:
            fit['Inormal']=a0 + a1*Vnormal + a2*Vnormal**2 + a3*Vnormal**3 + offset

        if found_turnover:
            V0=fit['turnover']
            fit['Iturnover']=a0 + a1*V0 + a2*V0**2 + a3*V0**3 + offset
        return fit
        

        
    # if the above didn't work, we use the original curve fit 
    # we shift the fit curve up/down to have I(max_bias)=max_bias
    # which puts the max bias position at a point on the R=1 Ohm line
    Imax=a0 + a1*self.Vinfinity + a2*(self.Vinfinity**2) + a3*(self.Vinfinity**3)
    offset=self.Vinfinity-Imax
    fit['offset']=offset
    if found_turnover:
        V0=fit['turnover']
        fit['Iturnover']=a0 + a1*V0 + a2*V0**2 + a3*V0**3 + offset
    
    # find the tangent line of the fit to the I-V curve at the maximum bias
    # this should be equivalent to a circuit with resistance 1 Ohm
    # tangent is the first derivative of the polyfit
    slope=a1 + 2*a2*self.Vinfinity + 3*a3*(self.Vinfinity**2)
        
    # The line should have a slope equal to a circuit with resistance 1 Ohm
    if abs(slope)>self.zero:
        R1 = 1/slope
    else:
        R1=1/self.zero
    fit['R1']=R1
    self.debugmsg('polynomial_fit_parameters (3): setting R1 to %.3f' % R1)
    return fit

def combined_fit_parameters(self,fit):
    '''
    determine the TES characteristics from the fit of multiple polynomials
    this is called from fit_iv()
    '''
    self.debugmsg('calling combined_fit_parameters')
    # first of all, get the parameters from the mixed region
    fit=self.polynomial_fit_parameters(fit)
    Vturnover=fit['turnover']

    # Vinfinity is the virtual point where we calculate the offset
    
    # superconducting region and normal region
    Vsuper=fit['fitinfo'][0][0]
    Vnormal=fit['fitinfo'][0][1]
    fit['Vsuper']=Vsuper
    fit['Vnormal']=Vnormal
    
    # the normal resistance comes from the fit to the normal region
    c0=fit['fitinfo'][0][8]
    c1=fit['fitinfo'][0][9]
    R1=1/c1
    fit['R1']=R1
    self.debugmsg('combined_fit_parameters (1): setting R1 to %.3f' % R1)
    
    # find offset that puts current equal to max bias voltage (R=1 at Vbias=Vinfinity)
    Imax=self.model_iv_normal(self.Vinfinity,c0,c1)
    offset=self.Vinfinity-Imax
    fit['offset']=offset

    # For comparison, we calculate the intersection of the two models: super/normal
    # We expect intuitively that this is the Vturnover, but it doesn't work out that way
    a0=fit['fitinfo'][0][2]
    a1=fit['fitinfo'][0][3]
    b=a0-c0
    discr=b**2 + 4*c1*a1
    if discr<0.0:
        # if no intersection between the models,
        Vsupernormal=None
    else:
        v01= 0.5*(b+np.sqrt(discr))/c1
        v02= 0.5*(b-np.sqrt(discr))/c1
        fit['model intersection']=(v01,v02)
        if v01>=Vsuper and v01<=Vnormal:
            Vsupernormal=v01
        else:
            Vsupernormal=v02
    fit['supernormal']=Vsupernormal
            
    b0=fit['fitinfo'][0][4]
    b1=fit['fitinfo'][0][5]
    b2=fit['fitinfo'][0][6]
    b3=fit['fitinfo'][0][7]
    if isinstance(Vturnover,float):
        fit['Iturnover']=self.model_iv_mixed(Vturnover,b0,b1,b2,b3) + offset
    else:
        fit['Iturnover']=None
        self.debugmsg('combined_fit_parameters: setting Iturnover to None')

    # get the current, according to the model, at the region transitions: Vsuper, Vnormal
    fit['Isuper'] =self.model_iv_mixed(Vsuper ,b0,b1,b2,b3) + offset
    fit['Inormal']=self.model_iv_mixed(Vnormal,b0,b1,b2,b3) + offset
        
    return fit

def do_linefit(self,bias,curve):
    '''
    fit a straight line
    '''
    self.debugmsg('I-V linefit.  A straight line through the whole I-V curve.')
    npts=len(bias)
    polyfit=np.polyfit(bias,curve,1,full=True)
    residual=polyfit[1][0]/npts
    ret={}
    ret['fitinfo']=polyfit
    ret['residual']=residual
    return ret

def do_polyfit(self,bias,curve):
    '''
    fit I-V curve to polynomial degree 3
    normalize the residual to the number of points in the fit
    '''
    self.debugmsg('I-V polyfit.  Single polynomial for the whole I-V curve.')
    npts=len(bias)
    polyfit=np.polyfit(bias,curve,3,full=True)
    residual=polyfit[1][0]/npts
    ret={}
    ret['fitinfo']=polyfit
    ret['residual']=residual
    return ret

def do_combinedfit(self,TES,bias,curve,Vsuper=None,Vnormal=None):
    '''
    fit I-V curve to a combined polynomial 
    (see model_iv_combined() below)
    '''
    self.debugmsg('calling do_combinedfit')
    # dictionary for return
    fit={}

    '''
    self.TES=TES # required for model_iv_combined()
    Vturnover=self.turnover(TES)
    '''
    
    Vmax=max(bias)
    Vmin=min(bias)
    Vturnover_idx=np.argmin(curve)
    Vturnover=bias[Vturnover_idx]
    if Vnormal is None:
        Vnormal=Vturnover + 0.5*(Vmax-Vturnover)
    if Vsuper is None:
        Vsuper=Vturnover - 0.5*(Vturnover-Vmin)
    self.debugmsg('I-V do_combined_fit, using Vsuper=%.2f and Vnormal=%.2f' % (Vsuper,Vnormal))

    a0=1.0
    a1=1.0
    b0=1.0
    b1=1.0
    b2=1.0
    b3=1.0
    c0=1.0
    c1=1.0
    p0=[Vsuper,Vnormal,a0,a1,b0,b1,b2,b3,c0,c1]

    popt,pcov=curve_fit(self.model_iv_combined,bias,curve,p0=p0)
    Vsuper,Vnormal,a0,a1,b0,b1,b2,b3,c0,c1=popt
    fitinfo=[popt,pcov]
    fit['fitinfo']=fitinfo
    fit['Vsuper']=Vsuper
    fit['Vnormal']=Vnormal
    fit['R1']=1.0/c1
    self.debugmsg('do_combinedfit: setting R1 to %.3f' % (1.0/c1))

    # calculate a performance measure
    npts=len(bias)
    Vfit=self.model_iv_combined(bias,Vsuper,Vnormal,a0,a1,b0,b1,b2,b3,c0,c1)
    sigma2=((curve-Vfit)/Vfit)**2
    residual = np.sqrt( np.sum(sigma2) )/npts
    fit['residual']=residual

    return fit


def model_iv_super(self,V,coeff0,coeff1):
    '''
    the model of the superconducting portion of the I-V curve
    '''
    self.printmsg('DEBUG: call to model_iv_super()',verbosity=5)
    I=coeff0 + coeff1/V
    return I

def model_iv_normal(self,V,coeff0,coeff1):
    '''
    the model of the normal portion of the I-V curve
    '''
    self.printmsg('DEBUG: call to model_iv_normal()',verbosity=5)
    I=coeff0 + coeff1*V
    return I

def model_iv_mixed(self,V,coeff0,coeff1,coeff2,coeff3):
    '''
    mixed region modeled with polynomial

    We impose the condition that the turnover occurs at the intersection of the super and normal models
    Turnover occurs when the tangent to the curve is zero (first derivative=0)
    0 = coeff1*Vturnover + 2*coeff2*Vturnover + 3*coeff3*Vturnover**2
    with the above we impose a relation between coeff3 and coeff2 and coeff1, knowing Vturnover
    '''
    self.printmsg('DEBUG: call to model_iv_mixed()',verbosity=5)

    # Vturnover sneaks in via a class variable.  self.TES should be assigned in fit_iv()
    #Vturnover=self.turnover(self.TES)
    
    #coeff3 = - (coeff1 + 2*coeff2*Vturnover)/(3*Vturnover**2)
    I=coeff0 + coeff1*V + coeff2*V**2 + coeff3*V**3
    return I

def model_iv_combined(self,V,Vsuper,Vnormal,a0,a1,b0,b1,b2,b3,c0,c1):
    '''function to fit I-V curve in three parts.  Damien Prele suggests some physics...
      1. 1/V : TES in super conducting state
      2. polynomial 2nd order to join the super to the normal state
      3. straight line in normal region

    there are 10 parameters to fit!

    BUG: The Vsuper and Vnormal remain fixed at the first guess. (or not?)

    '''
    self.printmsg('DEBUG: call to model_iv_combined()',verbosity=5)
    
    if isinstance(V,np.ndarray):
        I=np.empty(len(V))
        for idx,v in enumerate(V):
            if v<Vsuper:
                I[idx]=self.model_iv_super(v,a0,a1)
            elif v<Vnormal:
                I[idx]=self.model_iv_mixed(v,b0,b1,b2,b3)
            else:
                I[idx]=self.model_iv_normal(v,c0,c1)
        return I

    
    if V<Vsuper:
        return self.model_iv_super(V,a0,a1)
    if V<Vnormal:
        return self.model_iv_mixed(v,b0,b1,b2,b3)
    
    return self.model_iv_normal(V,c0,c1)

def fit_iv(self,TES,
           jumplimit=None,
           curve_index=None,
           fitfunction='COMBINED',
           Vsuper=None,
           Vnormal=None,
           istart=None,
           iend=None,
           R1adjust=1.0):
    '''
    fit the I-V curve to a polynomial

    if we're cycling, we always go down-up-down or up-down-up for a single cycle
    so a single IV curve is 1/(2*ncycles)

    we work directly with the uncalibrated data.
    The fit will be used to make the final adjustments

    optional arguments: 
       jumplimit:    this is the smallest step considered to be a jump in the data
       curve_index:  force the fit to use a particular curve in the cycle, and not simply the "best" one
       fitfunction:  use a 3rd degree polynomial, or a combination of polynomials
       istart, iend: window for fitting (start and end indices)
       Vsuper:       Bias voltage below which the TES is superconducting
       Vnormal:      Bias voltage above which the TES is normal
    '''
    if not self.exist_iv_data():return None

    self.debugmsg('calling fit_iv with Vsuper,Vnormal = %s, %s' % (str(Vsuper),str(Vnormal)))
    self.TES=TES # we need this global variable for the "mixed model".  see above in model_iv_mixed().
    TES_idx=TES_index(TES)
    self.debugmsg('fit_iv() : R1adjust=%f' % R1adjust)
    override_istart=istart
    override_iend=iend
    
    # return is a dictionary with various info
    fit={}
    fit['TES']=TES
    fit['fitfunction']=fitfunction.upper()

    I=self.ADU2I(self.adu[TES_idx,:],R1adjust=R1adjust)
    npts=len(I)

    if self.cycle_vbias:
        ncurves=self.nbiascycles*2
    else:
        ncurves=self.nbiascycles
    npts_curve=int(npts/ncurves)
    fit['ncurves']=ncurves
    fit['npts_curve']=npts_curve
    self.debugmsg('number of curves: %i' % ncurves)
    self.debugmsg('npts per curve: %i' % npts_curve)
    
    # fit for each measured curve and find the best one
    best_residual=1./self.zero
    best_curve_index=0
    allfits=[]
    fitranges=[]
    residuals=[]
    istart=0
    for idx in range(ncurves):
        iend=istart+npts_curve
        ypts=I[istart:iend]
        self.debugmsg('cycle %i: fitting curve istart=%i, iend=%i' % ((idx+1),istart,iend))
        xpts=self.bias_factor*self.vbias[istart:iend] # should maybe use Vtes here? but there's a chance of recursion with Voffset.

        # filter out the big jumps
        # the return is the range of indexes of the acceptable points
        good_start,good_end=self.filter_jumps(ypts,jumplimit)
        # check if we are overriding the jump algorithm with start and end indices
        if isinstance(override_istart,int):
            good_start=override_istart
            self.debugmsg('override of fit window.  good_start=%i' % good_start)
        if isinstance(override_iend,int):
            good_end=override_iend
            self.debugmsg('override of fit window.  good_end=%i' % good_end)
        self.debugmsg('fit_iv:  good_start,good_end=%i,%i' % (good_start,good_end))
        npts_span=good_end-good_start
        if npts_span<11:
            self.debugmsg('couldn\'t find a large span without jumps! Fitting the whole curve...')
            good_start=0
            good_end=len(xpts)
            npts_span=npts_curve
        curve=ypts[good_start:good_end]
        bias=xpts[good_start:good_end]

        if fitfunction=='POLYNOMIAL':
            ivfit=self.do_polyfit(bias,curve)
        else:
            ivfit=self.do_combinedfit(TES,bias,curve,Vsuper,Vnormal)

        #for key in ivfit.keys(): fit[key]=ivfit[key] # it's a hack.    
        residual=ivfit['residual']
        residuals.append(residual)
        allfits.append(ivfit['fitinfo'])
        fitranges.append((good_start,good_end))
        if abs(residual)<best_residual:
            best_residual=abs(residual)
            best_curve_index=idx
        istart+=npts_curve

    # from now on we use the best curve fit
    # unless there is request to override with the curve_index option
    fit['best curve index']=best_curve_index
    if not curve_index is None:
        if not isinstance(curve_index,int) or curve_index>=ncurves or curve_index<0:
            print('Invalid option for curve index:  Please give an integer between 0 and %i' % (ncurves-1))
            print('Using default:  best curve index=%i' % best_curve_index)
            curve_index=best_curve_index
    else:
        curve_index=best_curve_index
    fit['curve index']=curve_index
    fitinfo=allfits[curve_index]
    fit['fitinfo']=fitinfo
    fit['fit range']=fitranges[curve_index]
    fit['residual']=residuals[curve_index]

    if fitfunction=='POLYNOMIAL':
        fit=self.polynomial_fit_parameters(fit)
    else:
        fit=self.combined_fit_parameters(fit)

    keys=''
    for key in fit.keys():keys+=key+', '
    self.debugmsg('returning from fit_iv() with keys: %s' % keys)
    return fit


def draw_iv(self,I,bias=None,colour='blue',axis=None,label=None):
    '''
    draw an individual I-V curve
    '''
    if axis is None: axis = plt.gca()
    npts=len(I)
    if npts<len(self.vbias) and npts>0:
        # this is a partial curve, it should come with a bias axis
        if bias is None: # assume it's the first part of Vbias during an ongoing acquisition
            bias = self.bias_factor*self.vbias[0:npts]      
            plt.cla()
            axis.set_xlim([self.bias_factor*self.min_bias,self.bias_factor*self.max_bias])
            axis.plot(bias,I,color=colour)

            # we mark the last point to show the progressive plot
            axis.plot(bias[-1],I[-1],color='red',marker='D',linestyle='none')

            # and the temperature during the acquisition
            if self.temperature is None:
                tempstr = 'unknown'
            else:
                tempstr = '%0.2f mK' % (1000*self.temperature)            
            tempstr = 'T$_\mathrm{bath}$=%s' % tempstr
            axis.text(0.05,0.95,tempstr,va='top',ha='left',fontsize=12,transform=axis.transAxes)
            plt.pause(0.01)
            return

        # if bias is given, just plot it
        axis.plot(bias,I,color=colour)
        return

    if self.cycle_vbias:
        # plot down and up voltage with different linestyles
        mid=int(len(self.bias_factor*self.vbias)/2)
        axis.plot(self.bias_factor*self.vbias[0:mid],I[0:mid],linestyle='solid', color=colour)
        axis.plot(self.bias_factor*self.vbias[mid:-1], I[mid:-1], linestyle='dashed',color=colour)
        return
    
    axis.plot(self.bias_factor*self.vbias,I,color=colour)
    return

def setup_plot_iv(self,TES,xwin=True):
    ttl=str('QUBIC I-V curve for TES#%3i (%s)' % (TES,self.obsdate.strftime('%Y-%b-%d %H:%M UTC')))
    if self.temperature is None:
        tempstr='unknown'
    else:
        tempstr=str('%.0f mK' % (1000*self.temperature))
    subttl=str('Array %s, ASIC #%i, Pixel #%i, Temperature %s' % (self.detector_name,self.asic,tes2pix(TES,self.asic),tempstr))
    if xwin: plt.ion()
    else: plt.ioff()
    fig=plt.figure(figsize=self.figsize)
    fig.canvas.set_window_title('plt: '+ttl) 
    fig.suptitle(ttl+'\n'+subttl,fontsize=16)
    ax=plt.gca()
    ax.set_xlabel('Bias Voltage  /  V')
    ax.set_ylabel('Current  /  $\mu$A')
    bias_range=self.bias_factor*(self.max_bias-self.min_bias)
    bias_plot_window=[self.bias_factor*self.min_bias - 0.1*bias_range,self.bias_factor*self.max_bias + 0.1*bias_range]
    ax.set_xlim(bias_plot_window)
    # Y plot limits calculated in plot_iv()
    return fig,ax

def adjusted_iv(self,TES):
    '''
    return the adjusted I-V curve
    '''
    offset=self.offset(TES)
    R1adjust=self.R1adjust(TES)
    Iadjusted=self.ADU2I(self.adu[TES_index(TES),:],offset=offset,R1adjust=R1adjust)
    return Iadjusted

def oplot_iv(self,TES,label=None,best=True,axis=None):
    Iadjusted = self.adjusted_iv(TES)
    bias = None
    if best:
        istart,iend=self.selected_iv_curve(TES)
        Iadjusted = Iadjusted[istart:iend]
        bias = self.bias_factor*self.vbias[istart:iend]
    return self.draw_iv(Iadjusted,bias=bias,label=label,axis=axis)

def plot_iv(self,TES=None,multi=False,xwin=True,best=True):
    '''
    plot the I-V curve.
      keywords:
        TES: the TES number to plot.  If "None" then plot all in the plot of the focal plane
        multi: plot all the TES I-V curves in a multi-window plot
        xwin: plot to the screen (or not)
        best: plot only the "best" curve if there are multiple I-V curves
    '''
    filterinfo=self.filterinfo(TES)
    if filterinfo is None:return None
    
    if multi:return self.plot_iv_multi()
    if TES is None:return self.plot_iv_physical_layout()
    if not isinstance(TES,int): return self.plot_iv_physical_layout()

    self.debugmsg('plotting I-V curve for TES %i' % TES)
    self.TES=TES
    TES_idx=TES_index(TES)
    fit=filterinfo['fit']
    
    fig,ax=self.setup_plot_iv(TES,xwin)

    # identify which fit function was used
    if 'fitfunction' not in fit.keys():
        fit['fitfunction']='POLYNOMIAL'
    txt='fit function = %s' % fit['fitfunction']
    
    # normalize the Current so that R=1 Ohm at the highest Voffset
    offset=self.offset(TES)
    txt+=str('\noffset=%.4e' % offset)
    Iadjusted=self.adjusted_iv(TES)
    Itop=max(Iadjusted)
    Ibot=min(Iadjusted)
    Irange=Itop-Ibot
    I_plot_window=[Ibot-0.5*Irange,Itop+0.1*Irange]
    ax.set_ylim(I_plot_window)
    Vturnover=self.turnover(TES)
    self.oplot_iv(TES,axis=ax,best=best)
        
    # draw a line tangent to the fit at the highest Vbias
    # I0 here is the current extrapolated to Vbias=0
    I0=self.draw_tangent(TES)
    
    R1=self.R1(TES)
    if not R1 is None: txt+=str('\ndynamic normal resistance:  R$_1$=%.4f $\Omega$' % R1)

    # draw a fit to the I-V curve
    txt+=str('\nfit residual: %.4e' % filterinfo['residual'])
    bias,f=self.fitted_iv_curve(TES)
    ax.plot(bias,f,linestyle='dashed',color='red',label='model')

    # draw the curve fitting the super conducting region
    if fit['fitfunction']=='COMBINED':
        a0=fit['fitinfo'][0][2]
        a1=fit['fitinfo'][0][3]
        Isuper=self.model_iv_super(bias,a0,a1)+self.offset(TES)
        ax.plot(bias,Isuper,linestyle='dashed',color='magenta',label='model superconducting region')

        # redefine Ibot as the intersection of the super/normal models.
        Vsupernormal=fit['supernormal']
        if isinstance(Vsupernormal,float):
            Ibot=self.model_iv_super(Vsupernormal,a0,a1)+self.offset(TES)

        Vsuper=fit['Vsuper']
        Vnormal=fit['Vnormal']
        ax.plot([Vsuper,Vsuper],[Ibot,Itop],linestyle='dashed',color='cyan')
        plt.text(Vsuper,Itop,'Superconducting  \nregion  ',ha='right',va='top',fontsize=12)
        ax.plot([Vnormal,Vnormal],[Ibot,Itop],linestyle='dashed',color='cyan')
        plt.text(Vnormal,Itop,'  Normal\n  region',ha='left',va='top',fontsize=12)
        ax.plot([Vsupernormal,Vsupernormal],[Ibot,Itop],linestyle='dashed',color='cyan')

    # draw vertical lines to show the range used for the fit
    if 'fit range' in fit.keys():
        fit_istart,fit_iend = fit['fit range']
        istart,iend = self.selected_iv_curve(TES)
        bias = self.bias_factor*self.vbias[istart:iend]
        fit_vstart = bias[fit_istart]
        fit_vend = bias[fit_iend]
        ax.plot([fit_vstart,fit_vstart],[Ibot,Itop],color='red',linestyle='dashed')
        ax.plot([fit_vend,fit_vend],[Ibot,Itop],color='red',linestyle='dashed')
    

    # note the turnover point
    if Vturnover is None:
        txt+='\nNo turnover!'
    else:
        ax.plot([Vturnover,Vturnover],[Ibot,Itop],linestyle='dashed',color='green')
        txt+=str('\nturnover Vbias=%.2fV' % Vturnover)
        # Iturnover is the current at Vturnover
        if 'Iturnover' in fit.keys():
            Iturnover=fit['Iturnover']
        else:
            Iturnover=Ibot
        txt+=str('\nI$_\mathrm{turnover}$=%.2f $\mu$A' % Iturnover)
        

    # add room temp results, if loaded
    if not self.transdic is None:
        PIX=tes2pix(TES,self.asic)
        # self.debugmsg('table lookup for PIX=%i' % PIX)
        entry=self.lookup_TEStable(key='PIX',value=PIX)
        R300=entry['R300']
        if isinstance(R300,float):
            R300str='%.2f $\Omega$' % R300
        else:
            R300str=R300
        txt+='\nRoom Temperature Resistance: %s' % R300str
        openloop=entry['OpenLoop']
        txt+='\nOpen Loop Test:  %s' % openloop
    
    is_good=self.is_good_iv(TES)
    comment=filterinfo['comment']
    if not is_good:txt+=str('\nFlagged as BAD:  %s' % comment)

    # write out the comments
    text_x=0.98
    text_y=0.02
    boxprops = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(text_x,text_y,txt,va='bottom',ha='right',fontsize=10,transform = ax.transAxes,bbox=boxprops)
    ax.legend(loc='lower left',bbox_to_anchor=(0.02, 0.02),fontsize=10)
    pngname=str('TES%03i_IV_array-%s_ASIC%i_%s.png' % (TES,self.detector_name,self.asic,self.obsdate.strftime('%Y%m%dT%H%M%SUTC')))
    pngname_fullpath=self.output_filename(pngname)
    if isinstance(pngname_fullpath,str): plt.savefig(pngname_fullpath,format='png',dpi=100,bbox_inches='tight')
    if xwin: plt.show()
    else: plt.close('all')
    return fig


def plot_pv(self,TES,xwin=True):
    '''
    plot the power vs voltage curve for the TES
    '''
    ttl=str('QUBIC P-V curve for TES#%3i (%s)' % (TES,self.obsdate.strftime('%Y-%b-%d %H:%M UTC')))
    if self.temperature is None:
        tempstr='unknown'
    else:
        tempstr=str('%.0f mK' % (1000*self.temperature))
    subttl=str('Array %s, ASIC #%i, Pixel #%i, Temperature %s' % (self.detector_name,self.asic,tes2pix(TES,self.asic),tempstr))
    if xwin: plt.ion()
    else: plt.ioff()
    fig,ax=plt.subplots(1,1,figsize=self.figsize)
    fig.canvas.set_window_title('plt: '+ttl) 
    fig.suptitle(ttl+'\n'+subttl,fontsize=16)
    ax.set_xlabel('Bias Voltage  /  V')
    ax.set_ylabel('P$_\mathrm{TES}$  /  $p$A')
    ax.set_xlim([self.bias_factor*self.min_bias,self.bias_factor*self.max_bias])

    istart,iend=self.selected_iv_curve(TES)
    Ptes=self.Ptes(TES)[istart:iend]
    bias=self.bias_factor*self.vbias[istart:iend]
    ax.plot(bias,Ptes)
    
    pngname=str('TES%03i_PV_array-%s_ASIC%i_%s.png' % (TES,self.detector_name,self.asic,self.obsdate.strftime('%Y%m%dT%H%M%SUTC')))
    pngname_fullpath=self.output_filename(pngname)
    if isinstance(pngname_fullpath,str): plt.savefig(pngname_fullpath,format='png',dpi=100,bbox_inches='tight')
    if xwin: plt.show()
    else: plt.close('all')
    return fig,ax
    
def plot_rp(self,TES,xwin=True):
    if self.R1(TES) is None:
        print('No normal resistance estimate.')
        return None

    istart,iend=self.selected_iv_curve(TES)
    
    Rn_ratio=self.Rn_ratio(TES)[istart:iend]
    Ptes=self.Ptes(TES)[istart:iend]
    Pbias=self.Pbias(TES)
    lbl='P$_\mathrm{bias}=$%.2f pW' % Pbias

    Rmin=min(Rn_ratio)
    Rmax=max(Rn_ratio)
    Rspan=Rmax-Rmin
    plot_Rmin=Rmin-0.2*Rspan
    plot_Rmax=100.
    
    Pmin=min(Ptes)
    Pmax=max(Ptes)
    Pspan=Pmax-Pmin
    plot_Pmin=Pmin-0.05*Pspan
    plot_Pmax=Pmax+0.2*Pspan
    
    ttl=str('QUBIC R-P curve for TES#%3i (%s)' % (TES,self.obsdate.strftime('%Y-%b-%d %H:%M UTC')))
    if self.temperature is None:
        tempstr='unknown'
    else:
        tempstr=str('%.0f mK' % (1000*self.temperature))
    subttl=str('Array %s, ASIC #%i, Pixel #%i, Temperature %s' % (self.detector_name,self.asic,tes2pix(TES,self.asic),tempstr))
    if xwin: plt.ion()
    else: plt.ioff()
    fig,ax=plt.subplots(1,1,figsize=self.figsize)
    fig.canvas.set_window_title('plt: '+ttl) 
    fig.suptitle(ttl+'\n'+subttl,fontsize=16)
    ax.set_xlabel('P$_\mathrm{TES}$  /  pW')
    ax.set_ylabel('$\\frac{R_\mathrm{TES}}{R_\mathrm{normal}}$ / %')

    ax.plot(Ptes,Rn_ratio)
    ax.plot([Pbias,Pbias],[0,90],linestyle='dashed',color='green')
    ax.plot([plot_Pmin,Pbias],[90,90],linestyle='dashed',color='green')
    ax.set_xlim([plot_Pmin,plot_Pmax])
    ax.set_ylim([plot_Rmin,plot_Rmax])

    text_x=0.98
    text_y=0.02
    boxprops = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(text_x,text_y,lbl,va='bottom',ha='right',fontsize=10,transform = ax.transAxes,bbox=boxprops)
    
    pngname=str('TES%03i_RP_array-%s_ASIC%i_%s.png' % (TES,self.detector_name,self.asic,self.obsdate.strftime('%Y%m%dT%H%M%SUTC')))
    pngname_fullpath=self.output_filename(pngname)
    if isinstance(pngname_fullpath,str): plt.savefig(pngname_fullpath,format='png',dpi=100,bbox_inches='tight')
    if xwin: plt.show()
    else: plt.close('all')
    return fig,ax


def Vbias2Vtes(self,V,Ites):
    '''
    return the Vtes for a given bias voltage.
    the current must also be given
    '''
    Vtes=self.Rshunt*(V/self.Rbias-Ites)
    return Vtes

def responsivity_func(self,Vtes,I,f_prime):
    '''
    the responsivity is dI/dP

    f_prime is the derivative of the model I=f(Vbias), dI/dVbias = f'(Vbias) at the point I,Vbias
    for more details see plot_responsivity() below
    '''
    responsivity = 1.0/( I*self.Rshunt * (1.0/(self.Rbias*f_prime) - 1) + Vtes )
    return responsivity

def conductance_func(self,f_prime):
    '''
    the conducance G0

    f_prime is the derivative of the model I=f(Vbias), dI/dVbias = f'(Vbias) at the point I,Vbias
    for more details see plot_responsivity() below
    '''
    G0=1.0/( self.Rshunt*(1.0/(self.Rbias*f_prime) - 1) )
    return G0

def plot_responsivity(self,TES,xwin=True,npts_region=500,window_size=51,filter_sigma=10,
                      plot_model=True,rmax=None,rmin=None):
    '''
    plot the responsivity of the TES

    The formula is:
      Si = -1/2Vbias * (Z0-R0)/(Z0+RL)                       [Eq.1]

      with
        Z0 = dVtes/dI = the slope of the curve at the point (Vtes,I)
        R0 = Vtes/I (at the point (Vtes,I)
        RL = Rshunt + Rparasitic = Rshunt + (approx) 10e-3 Ohm

        Vtes = Rshunt * (Vbias/Rbias-I)                      [Eq.2]
        dVtes/dP = Rshunt * ( (dVbias/dP)/Rbias - dI/dP )    [Eq.3]
        dVtes/dI = Rshunt * ( (dVbias/dI)/Rbias - 1 )        [Eq.4]

    # NOTE: at turnover, dV/dI becomes infinite (slope zero in the I-V curve)
    We can rewrite the formula in terms of Conductance (G0) rather than Resistance (Z0)

    Si = -1/2Vbias * (1-R0*G0)/(1-RL*G0)                     [Eq.5]

    with
       G0 = dI/dVtes = the slope of the curve at the point (Vtes,I)



    Is the above formula correct?
    Shih-Fu Lee et al (1998) give the definition:
    
    Si = dI/dP
    I  = P/Vtes
    dI/dP = 1/Vtes * ( 1 - (P/Vtes) dVtes/dP )               [Eq.6]

    and then substituting [Eq.3]

    dI/dP = 1/(Vtes-1) * (1 - I * Rshunt * (dVbias/dP)/Rbias) [Eq.7]


    in each region, we have 
      I=f(Vbias)
    so
      dI/dP = dVbias/dP f'(Vbias)

    and substituting into [Eq.7] and rearranging gives:

      dI/dP = [ IRshunt (1/Rbias f'(Vbias) - 1) + Vtes ]^-1    [Eq.8]
      

    # NOTE: the modelling gives the current in uA (crap programming, sorry)
    in the super region,
        I = f(Vbias)  = C0 + C1/Vbias
            f'(Vbias) = -C1/Vbias**2

        C0 has units of A
        C1 has units of V * uA
    
    in the overlap region
        I = f(Vbias)  = C0 + C1*Vbias + C2*Vbias**2 + C3*Vbias**3
            f'(Vbias) = C1 + 2*C2*Vbias + 3*C3*Vbias**2

        C0 has units of A
        C1 has units of A/V
        C2 has units of A/V^2
        C3 has units of A/V^3

     in the normal region
        I = f(Vbias)  = C0 + C1Vbias
            f'(Vbias) = C1

        C0 has units of A
        C1 has units of A/V

    '''
    filterinfo=self.filterinfo(TES)
    if filterinfo is None:return None

    if not filterinfo['fit']['fitfunction']=='COMBINED':
        print("I need to have the COMBINED model.  Please rerun the filter and select COMBINED for the fitfunction.")
        return None

    RL = self.Rshunt + 10e-3


    # some parameters
    fitparms=filterinfo['fit']['fitinfo'][0]
    turnover_bias=self.turnover(TES)
    Iturnover_TES=filterinfo['fit']['Iturnover']
    if not Iturnover_TES is None:
        Iturnover_TES*=1e-6
        Vturnover_TES=self.Vbias2Vtes(turnover_bias,Iturnover_TES)
    else:
        Vturnover_TES=None
        
    Ioffset_TES=1e-6*self.offset(TES)
    Vsuper =filterinfo['fit']['Vsuper']
    Vnormal=filterinfo['fit']['Vnormal']

    # take npts_region points in each region: super,combined,normal
    dVscale=1./npts_region
    responsivity=[]
    voltage=[]
    Z0_model=[]
    G0_model=[]
    Imodel=[]

    # super region
    range_min=self.bias_factor*self.min_bias
    range_max=self.bias_factor*Vsuper
    dV=dVscale*(range_max-range_min)
    vbias=np.arange(range_min,range_max,dV)
    V=range_min
    C0=fitparms[2] * 1e-6 # units of Amps
    C1=fitparms[3] * 1e-6 # units of Amps * Volt
    Isuper_TES=self.model_iv_super(Vsuper,C0,C1) + Ioffset_TES
    Vsuper_TES=self.Vbias2Vtes(Vsuper,Isuper_TES)
    for V in vbias:
        Ites=self.model_iv_super(V,C0,C1) + Ioffset_TES
        Imodel.append(Ites)
        Vtes=self.Vbias2Vtes(V,Ites)
        f_prime=-C1/V**2
        R0 = Vtes/Ites
        G0 = self.conductance_func(f_prime)
        Si = self.responsivity_func(Vtes,Ites,f_prime)
        responsivity.append(Si)
        voltage.append(Vtes)
        G0_model.append(G0)

    # mixed region
    range_min=self.bias_factor*Vsuper
    range_max=self.bias_factor*Vnormal
    dV=dVscale*(range_max-range_min)
    V=range_min
    vbias=np.arange(range_min,range_max,dV)
    C0=fitparms[4] * 1e-6 # units of Amps
    C1=fitparms[5] * 1e-6 # units of Amps/Volt
    C2=fitparms[6] * 1e-6 # units of Amps/Volt^2
    C3=fitparms[7] * 1e-6 # units of Amps/Volt^3
    Inormal_TES=self.model_iv_mixed(V,C0,C1,C2,C3) + Ioffset_TES
    Vnormal_TES=self.Vbias2Vtes(Vnormal,Inormal_TES)
    for V in vbias:
        Ites=self.model_iv_mixed(V,C0,C1,C2,C3) + Ioffset_TES
        Imodel.append(Ites)
        Vtes=self.Vbias2Vtes(V,Ites)
        f_prime = C1 + 2*C2*V + 3*C3*V**2
        R0 = Vtes/Ites
        G0 = self.conductance_func(f_prime)
        Si = self.responsivity_func(Vtes,Ites,f_prime)
        responsivity.append(Si)
        voltage.append(Vtes)
        G0_model.append(G0)
        
    # normal region
    range_min=self.bias_factor*Vnormal
    range_max=self.bias_factor*self.max_bias
    dV=dVscale*(range_max-range_min)
    V=range_min
    vbias=np.arange(range_min,range_max,dV)
    C0=fitparms[8] * 1e-6 # units of Amps
    C1=fitparms[9] * 1e-6 # units of Amps/Volt
    for V in vbias:
        Ites=self.model_iv_normal(V,C0,C1) + Ioffset_TES
        Imodel.append(Ites)
        Vtes=self.Vbias2Vtes(V,Ites)
        f_prime = C1
        Si = self.responsivity_func(Vtes,Ites,f_prime)
        R0 = Vtes/Ites
        G0 = self.conductance_func(f_prime)
        Si = self.responsivity_func(Vtes,Ites,f_prime)
        responsivity.append(Si)
        voltage.append(Vtes)
        G0_model.append(G0)

    # convert to numpy arrays
    responsivity=np.array(responsivity)
    

    # Now make the calculation based on the smoothed measurement, instead of using the model
    istart,iend=self.selected_iv_curve(TES)
    I=self.Ites(TES)[istart:iend]
    V=self.Vtes(TES)[istart:iend]
    Vbias=self.vbias[istart:iend]
    P=I*V
    
    # make sure window_size is an odd number
    if window_size % 2 == 0:window_size+=1
    window=np.hanning(window_size)
    Ismooth=np.convolve(window/window.sum(),I,mode='valid')
    Vsmooth=V[window_size//2:-window_size//2]
    if Vsmooth.size != Ismooth.size:
        npts_diff=Vsmooth.size-Ismooth.size
        self.debugmsg('applying correction for unequal windows: npts_diff=%i' % npts_diff)
        Vsmooth=V[npts_diff+window_size//2:-window_size//2]
    npts_smooth=Vsmooth.size
    G0 = np.gradient(Ismooth,Vsmooth)
    R0 = Vsmooth/Ismooth
    Psmooth = Ismooth*Vsmooth
    meas_responsivity = np.gradient(Ismooth,Psmooth)
    
    self.debugmsg('min(G0_model),max(G0_model)=%.4e,%.4e' % (min(G0_model),max(G0_model)))
    bias,Imodel_intrinsic=self.fitted_iv_curve(TES)
    Imodel_intrinsic*=1.0e-6
    Pmodel=Imodel_intrinsic*V
    
    ttl=str('QUBIC Responsivity curve for TES#%3i (%s)' % (TES,self.obsdate.strftime('%Y-%b-%d %H:%M UTC')))
    if self.temperature is None:
        tempstr='unknown'
    else:
        tempstr=str('%.0f mK' % (1000*self.temperature))
    subttl=str('Array %s, ASIC #%i, Pixel #%i, Temperature %s' % (self.detector_name,self.asic,tes2pix(TES,self.asic),tempstr))
    if xwin: plt.ion()
    else: plt.ioff()
    fig,ax=plt.subplots(1,1,figsize=self.figsize)
    fig.canvas.set_window_title('plt: '+ttl) 
    fig.suptitle(ttl+'\n'+subttl,fontsize=16)
    ax.set_xlabel('TES Voltage  /  $\mu$V')
    ax.set_ylabel('S$_\mathrm{i}$  /  A$\cdot\mathrm{W}^{-1}$')

    # rescale V axis to uV
    voltage=1e6*np.array(voltage)
    Vsuper_TES*=1.0e6
    Vnormal_TES*=1.0e6
    plot_xrange=max(voltage)-min(voltage)
    ax.set_xlim([min(voltage)-0.1*plot_xrange,max(voltage)+0.1*plot_xrange])

    # filter out the noisy first part
    resp_middle=meas_responsivity[npts_smooth//5:-npts_smooth//5]
    std_dev=resp_middle.std()
    start_idx=0
    for idx,val in enumerate(meas_responsivity[:-npts_smooth//5]):
        if abs(val/std_dev) > filter_sigma:start_idx=idx
    ax.plot(1e6*Vsmooth[start_idx:],meas_responsivity[start_idx:],label='responsivity from data')
    
    # y axis limits
    top=max(meas_responsivity[start_idx:])
    bot=min(meas_responsivity[start_idx:])
    if not rmin is None:bot=rmin
    if not rmax is None:top=rmax
    plot_yrange=top-bot
    ax.set_ylim([bot-0.1*plot_yrange,top+0.1*plot_yrange])

    # draw vertical lines showing the three regions
    ax.plot([Vsuper_TES,Vsuper_TES],[bot,top],linestyle='dashed',color='green')
    plt.text(Vsuper_TES,top,'Superconducting  \nregion  ',ha='right',va='top',fontsize=12,color='green')
    ax.plot([Vnormal_TES,Vnormal_TES],[bot,top],linestyle='dashed',color='green')
    plt.text(Vnormal_TES,top,'  Normal\n  region',ha='left',va='top',fontsize=12,color='green')


    if plot_model:
        ax.plot(voltage,responsivity,label='responsivity from model')

    self.debugmsg('len(Vsmooth)=%i' % len(Vsmooth))
    self.debugmsg('len(Ismooth)=%i' % len(Ismooth))
    self.debugmsg('window_size=%i' % window_size)

    ax.plot(Vbias,0.5/V,label='$\dfrac{1}{2\mathrm{V}_\mathrm{TES}}$',color='red')

    # draw a vertical line at the turnover
    if not Vturnover_TES is None:
        Vturnover=1e6*Vturnover_TES
        turnover_label='V$_\mathrm{turnover}$=%.3f $\mu$V' % Vturnover
        ax.plot([Vturnover,Vturnover],[bot,top],linestyle='dashed',color='blue',label=turnover_label)

    ax.legend()

    
    pngname=str('TES%03i_responsivity_Array-%s_ASIC%i_%s.png' % (TES,self.detector_name,self.asic,self.obsdate.strftime('%Y%m%dT%H%M%SUTC')))
    pngname_fullpath=self.output_filename(pngname)
    if isinstance(pngname_fullpath,str): plt.savefig(pngname_fullpath,format='png',dpi=100,bbox_inches='tight')
    if xwin: plt.show()
    else: plt.close('all')

    retval={}
    retval['responsivity']=responsivity
    retval['meas_responsivity']=meas_responsivity
    
    return retval



def plot_ip(self,TES,xwin=True):
    '''
    plot Current vs Power for the TES
    '''

    if not self.exist_iv_data():return None
    filterinfo=self.filterinfo(TES)
    if filterinfo is None:return None    
    
    istart,iend=self.selected_iv_curve(TES)

    # data
    Ites=1e6*self.Ites(TES)[istart:iend] # uA
    Ptes=self.Ptes(TES)[istart:iend] # pW

    # model
    bias,Imodel=self.fitted_iv_curve(TES)
    Imodel=1.0e-6*Imodel[istart:iend]
    Vmodel=self.Vbias2Vtes(bias[istart:iend],Imodel)
    Pmodel=Imodel*Vmodel*1e12 # pW
    Imodel*=1e6 # uA

    # turnover
    Pturnover=None
    Vturnover=self.turnover(TES)
    Iturnover=self.Iturnover(TES)
    if not Vturnover is None and not Iturnover is None:
        Vturnover_TES=self.Vbias2Vtes(Vturnover,Iturnover*1e-6)
        Pturnover=Iturnover*Vturnover_TES*1e6 # pW

    # bias voltage at boundary between super region and overlap region
    Vsuper = filterinfo['fit']['Vsuper']
    Isuper = filterinfo['fit']['Isuper']
    Psuper=None
    if not Vsuper is None and not Isuper is None:
        Vsuper_TES=self.Vbias2Vtes(Vsuper,Isuper*1e-6)
        Psuper=Isuper*Vsuper_TES*1e6 # pW

    # bias voltage at boundary between overlap region and normal region
    Vnormal=filterinfo['fit']['Vnormal']
    Inormal = filterinfo['fit']['Inormal']
    Pnormal=None
    if not Vnormal is None and not Inormal is None:
        Vnormal_TES=self.Vbias2Vtes(Vnormal,Inormal*1e-6)
        Pnormal=Inormal*Vnormal_TES*1e6 # pW

    
    Imin=min(Ites)
    Imax=max(Ites)
    Ispan=Imax-Imin
    plot_Imin=Imin-0.2*Ispan
    plot_Imax=Imax+0.2*Ispan
    
    Pmin=min(Ptes)
    Pmax=max(Ptes)
    Pspan=Pmax-Pmin
    plot_Pmin=Pmin-0.05*Pspan
    plot_Pmax=Pmax+0.2*Pspan
    
    ttl=str('QUBIC I-P curve for TES#%3i (%s)' % (TES,self.obsdate.strftime('%Y-%b-%d %H:%M UTC')))
    if self.temperature is None:
        tempstr='unknown'
    else:
        tempstr=str('%.0f mK' % (1000*self.temperature))
    subttl=str('Array %s, ASIC #%i, Pixel #%i, Temperature %s' % (self.detector_name,self.asic,tes2pix(TES,self.asic),tempstr))
    if xwin: plt.ion()
    else: plt.ioff()
    fig,ax=plt.subplots(1,1,figsize=self.figsize)
    fig.canvas.set_window_title('plt: '+ttl) 
    fig.suptitle(ttl+'\n'+subttl,fontsize=16)
    ax.set_xlabel('P$_\mathrm{TES}$  /  pW')
    ax.set_ylabel('I$_\mathrm{TES}$  /  $\mu$A')

    ax.plot(Ptes,Ites,label='data')
    ax.plot(Pmodel,Imodel,label='model')

    # show turnover
    if not Pturnover is None:
        turnover_label='I$_\mathrm{turnover}$=%.3f$\mu$A, P$_\mathrm{turnover}$=%.3fpW' % (Iturnover,Pturnover)
        ax.plot([Pturnover,Pturnover],[plot_Imin,Iturnover],linestyle='dashed',color='red',label=turnover_label)
        ax.plot([plot_Pmin,Pturnover],[Iturnover,Iturnover],linestyle='dashed',color='red')

    # show boundaries
    if not Psuper is None:
        ax.plot([Psuper,Psuper],[plot_Imin,plot_Imax],linestyle='dashed',color='green')
        plt.text(Psuper,plot_Imax,'Superconducting  \nregion  ',ha='right',va='top',fontsize=10,color='green')
    if not Pnormal is None:
        ax.plot([Pnormal,Pnormal],[plot_Imin,plot_Imax],linestyle='dashed',color='green')
        plt.text(Pnormal,plot_Imax,'  Normal\n  region',ha='left',va='top',fontsize=10,color='green')
        
    ax.set_xlim([plot_Pmin,plot_Pmax])
    ax.set_ylim([plot_Imin,plot_Imax])
    ax.legend()
    
    pngname=str('TES%03i_IP_array-%s_ASIC%i_%s.png' % (TES,self.detector_name,self.asic,self.obsdate.strftime('%Y%m%dT%H%M%SUTC')))
    pngname_fullpath=self.output_filename(pngname)
    if isinstance(pngname_fullpath,str): plt.savefig(pngname_fullpath,format='png',dpi=100,bbox_inches='tight')
    if xwin: plt.show()
    else: plt.close('all')
    return fig,ax

                
def make_Vbias(self,cycle=True,ncycles=2,vmin=0.5,vmax=3.0,dv=0.002,lowhigh=True):
    '''
    the bias voltage values used during the I-V curve measurement
    '''

    if vmax<0.0:
        vmax=np.abs(vmax)
        print('No negative values for bias! Setting max bias to %.2f V' % vmax)

    if vmax>self.max_permitted_bias:
        print('It is dangerous to set the bias voltage greater than %.2f V.' % self.max_permitted_bias)
        print('Setting maximum bias to %.2f V' % self.max_permitted_bias)
        vmax=self.max_permitted_bias

    max_offset=self.DAC2V * 2**15
    if vmax>max_offset:
        print('WARNING! Cannot set bias offset greater than %.3f V.' % max_offset)
        print('Setting maximum bias to %.2f V' % max_offset)
        vmax=max_offset

    if vmin<0.0:
        print('No negative values! Setting minimum bias to 0 V')
        vmin=0.0

    if ncycles<1:
        print('You need at least one cycle! Setting ncycles=1')
        ncycles=1
    
    going_up=np.arange(vmin,vmax+dv,dv)
    going_dn=np.flip(going_up,0)

    if cycle:
        if lowhigh:
            onecycle=np.concatenate((going_up,going_dn),0)
        else:
            onecycle=np.concatenate((going_dn,going_up),0)
    else:
        if lowhigh:
            onecycle=going_up
        else:
            onecycle=going_dn

    self.cycle_vbias=cycle
    self.nbiascycles=ncycles

    vbias=onecycle
    for n in range(ncycles-1):
        vbias=np.concatenate((vbias,onecycle),0)
    
    self.vbias=vbias
    self.min_bias=min(self.vbias)
    self.max_bias=max(self.vbias)
    self.max_bias_position=np.argmax(self.vbias)
    return vbias


def filter_iv(self,TES,
              R1adjust=1.0,
              R1_limit=10.0,
              residual_limit=3.0,
              abs_amplitude_limit=0.01,
              rel_amplitude_limit=0.1,
              bias_margin=0.2,
              jumplimit=None,
              curve_index=None,
              fitfunction='COMBINED',
              Vsuper=None,
              Vnormal=None,
              istart=None,
              iend=None):
    '''
    determine if this is a good TES from the I-V curve
    '''
    self.TES=TES # in case we use the model which requires known Vturnover (see fit_iv)
    TES_idx=TES_index(TES)

    # reset the filterinfo to avoid confusion regarding offset used in Ites and Vtes
    self.assign_filterinfo(TES,None)
    
    # dictionary to return stuff
    ret={}
    ret['R1adjust']=R1adjust
    self.debugmsg('filter_iv() : R1adjust=%f' % R1adjust)
    ret['TES']=TES
    ret['is_good']=True
    ret['comment']='no comment'

    # fit to the chosen model. The fit will be for the best measured curve if it's cycled bias
    fit = self.fit_iv(TES,jumplimit,curve_index,fitfunction,Vsuper,Vnormal,istart,iend,R1adjust)
    ret['fit'] = fit
    residual = fit['residual']
    ret['residual'] = residual
    offset = fit['offset']
    ret['offset'] = offset
    ADU = self.adu[TES_idx,:]
    Iadjusted = self.ADU2I(ADU,offset=offset,R1adjust=R1adjust)
    ret['turnover'] = fit['turnover']
    R1 = fit['R1']
    ret['R1'] = R1

    ##########################
    # filter out bad detectors
    ##########################

    # negative normal resistance is bad
    if R1<0.0:
        ret['is_good']=False
        ret['comment']='negative normal resistance'
        return self.assign_filterinfo(TES,ret)

    # big normal resistance is bad
    if R1>R1_limit:
        ret['is_good']=False
        ret['comment']='high normal resistance'
        return self.assign_filterinfo(TES,ret)
    
    # is it a good fit?
    if residual>residual_limit:
        ret['is_good']=False
        ret['comment']='bad poly fit'
        return self.assign_filterinfo(TES,ret)
    
    # small amplitude is rejected
    # we use the best curve as determined by the fit unless curve_index was specified
    curve_index=fit['curve index']
    npts=fit['npts_curve']
    istart=npts*curve_index
    iend=npts*(curve_index+1)
    meanval=Iadjusted[istart:iend].mean()
    maxval=max(Iadjusted[istart:iend])
    minval=min(Iadjusted[istart:iend])
    spread=abs(maxval-minval)
    self.debugmsg('maxval=%f, minval=%f, abs amplitude=%f' % (maxval,minval,spread))
    ret['abs_amplitude']=spread
    if spread<abs_amplitude_limit:
        ret['is_good']=False
        ret['comment']='current too low'
        return self.assign_filterinfo(TES,ret)
        
    # peak to peak amplitude
    rel_amplitude=abs(spread/meanval)
    ret['rel_amplitude']=rel_amplitude
    if rel_amplitude<rel_amplitude_limit:
        ret['is_good']=False
        ret['comment']='current peak-to-peak too small'
        return self.assign_filterinfo(TES,ret)
    
    # do we find a valid turnover for the Vbias?
    ret['turnover']=fit['turnover']
    if fit['turning'] is None or fit['turnover'] is None:
        ret['is_good']=False
        ret['comment']='no turnover'
        return self.assign_filterinfo(TES,ret)

    # is the operational point (the turnover) within the acceptable range?
    if ret['turnover']<self.bias_factor*self.min_bias+bias_margin or ret['turnover']>self.bias_factor*self.max_bias-bias_margin:
        ret['is_good']=False
        ret['comment']='operation point outside acceptable range'
        return self.assign_filterinfo(TES,ret)

    # do we have both turning points within the bias range?
    # maybe I should delete this filter
    #if fit['turnings within range']>1:
    #    ret['is_good']=False
    #    ret['comment']='bad I-V profile'
    #    return self.assign_filterinfo(TES,ret)
    
    
    # we only get this far if it's a good I-V
    return self.assign_filterinfo(TES,ret)

def filter_iv_all(self,
                  R1adjust=1.0,
                  R1_limit=10,
                  residual_limit=3.0,
                  abs_amplitude_limit=0.01,
                  rel_amplitude_limit=0.1,
                  bias_margin=0.2,
                  jumplimit=None,
                  fitfunction='COMBINED',
                  Vsuper=None,
                  Vnormal=None,
                  istart=None,
                  iend=None):
    '''
    find which TES are good
    '''
    if not self.exist_iv_data():
        print('No data!  Please read a file, or run a measurement.')
        return None

    # return a list with the filter info for each TES
    filtersummary=[]

    # assign the R1 adjustment for each TES
    R1adjust_vector=np.ones(self.NPIXELS)
    if R1adjust is None:R1adjust=1.0
    if isinstance(R1adjust,float) or isinstance(R1adjust,int):
        for idx in range(self.NPIXELS):R1adjust_vector[idx]=R1adjust
    elif len(R1adjust)==self.NPIXELS:R1adjust_vector=R1adjust
        
    # go through each filter.  Jump out and examine the next I-V curve as soon as a bad I-V is found
    for TES_idx in range(self.NPIXELS):
        TES=TES_idx+1
        self.debugmsg('running filter on TES %03i' % TES)
        filterinfo=self.filter_iv(TES,
                                  R1adjust_vector[TES_idx],
                                  R1_limit,
                                  residual_limit,
                                  abs_amplitude_limit,
                                  rel_amplitude_limit,
                                  bias_margin,
                                  jumplimit,
                                  curve_index=None,
                                  fitfunction=fitfunction,
                                  Vsuper=Vsuper,
                                  Vnormal=Vnormal,
                                  istart=istart,
                                  iend=iend)
        filtersummary.append(filterinfo)
        
    self.filtersummary=filtersummary
    return filtersummary

def save_filter(self):
    '''
    save the filter to a picke file
    '''
    datefmt='%Y%m%dT%H%M%SUTC'
    datestr=self.obsdate.strftime(datefmt)
    picklename=str('QUBIC_TES_%s.filter.pickle' % datestr)
    h=open(picklename,'w')
    pickle.dump(self.filtersummary,h)
    h.close()
    return

def read_filter(self):
    '''
    read the filter from a pickle file
    '''
    datefmt='%Y%m%dT%H%M%SUTC'
    datestr=self.obsdate.strftime(datefmt)
    picklename=str('QUBIC_TES_%s.filter.pickle' % datestr)
    if not os.path.exists(picklename):
        print('No previously saved filter information: %s' % picklename)
        return None

    print('Reading previously saved filter information: %s' % picklename)
    h=open(picklename,'r')
    filtersummary=pickle.load(h)
    h.close()
    self.filtersummary=filtersummary
    return True

def read_ADU_file(self,filename):
    '''
    legacy:  a few files were produced in this format before the FITS file was defined
    '''
    if not os.path.exists(filename):
        print("file not found: ",filename)
        return None

    # try to get date from filename
    self.assign_obsdate(self.read_date_from_filename(filename))

    handle=open(filename,'r')
    raw=handle.read()
    handle.close()
    lines=raw.split('\n')
    nlines=0
    X=[]
    for line in lines:
        if line=='':
            continue
        nlines+=1
        line_list=line.split('\t')
        row=np.empty(len(line_list))
        i=0
        for strval in line_list:
            val=eval(strval)
            row[i]=val
            i+=1
        X.append(row)
        adu=np.array(X)
    self.assign_ADU(adu)
    return adu

def iv_tex_table_entry(self,TES):
    TES_idx=TES_index(TES)
    PIX=tes2pix(TES,self.asic)
    if self.turnover(TES) is None:
        turnover='-'
    else:
        turnover=str('%.2f V' % self.turnover(TES))

    R1=self.R1(TES)
    if R1 is None or R1>10000:
        R1str='-'
    else:
        if abs(R1)<100:
            R1str=str('%.2f $\Omega$' % R1)
        else:
            R1str=str('%.2e $\Omega$' % R1)

    comment=self.filtersummary[TES_idx]['comment']
    if comment=='no comment': comment='good'

    if self.transdic is None:
        R300str='--'
        openloop='--'
        cf='--'
    else:
        # self.debugmsg('table lookup for PIX=%i' % PIX)
        entry=self.lookup_TEStable(key='PIX',value=PIX)
        R300=entry['R300']
        if isinstance(R300,float):
            if abs(R300)<10000:
                R300str='%.1f $\Omega$' % R300
            else:
                R300str='%.2e $\Omega$' % R300
        else:
            R300str=R300

        openloop=entry['OpenLoop']
        cf=entry['CarbonFibre']

    comment_entry=str('\\comment{%s}' % comment)
    rowstr='%3i & %3i & %s & %s & %s & %s & %s & %s' % (TES, PIX, turnover, R1str, R300str, openloop, cf, comment_entry)
    return rowstr

def make_iv_tex_report(self,tableonly=False):
    '''
    make a report in LaTeX.  
    This relies on the data in self.filtersummary.  See self.filter_iv_all() above
    '''
    if not self.exist_iv_data():return None
    
    thumbnailplot=str('TES_IV_array-%s_ASIC%i_%s.png'     % (self.detector_name,self.asic,self.obsdate.strftime('%Y%m%dT%H%M%SUTC')))
    allplot      =str('TES_IV_array-%s_ASIC%i_all_%s.png' % (self.detector_name,self.asic,self.obsdate.strftime('%Y%m%dT%H%M%SUTC')))
    pattern      =str('TES???_IV_array-%s_ASIC%i_%s.png'  % (self.detector_name,self.asic,self.obsdate.strftime('%Y%m%dT%H%M%SUTC')))

    # do the globbing in working directory
    cwd=os.getcwd()
    subdir=self.data_subdir()
    if isinstance(subdir,str):
        workdir='%s/%s' % (self.datadir,subdir)
    else:
        workdir=self.datadir
    os.chdir(workdir) # move to data directory
    iv_plots=glob(pattern)
    os.chdir(cwd) # and return to previous directory

    iv_plots.sort()

    if len(iv_plots)<self.NPIXELS:
        print('WARNING: Did not find all the I-V plots!')

    observer=self.observer.replace('<','$<$').replace('>','$>$')
    
    texfilename=str('TES_IV_array-%s_ASIC%i_%s.tex' % (self.detector_name,self.asic,self.obsdate.strftime('%Y%m%dT%H%M%SUTC')))
    texfilename_fullpath=self.output_filename(texfilename)
    if not isinstance(texfilename_fullpath,str):
        print('ERROR! Not possible to write tex file.')
        return None

    cfibre_ack = 'The carbon fibre measurements are from Sophie Henrot Versill\\a\'e, see \\url{http://qubic.in2p3.fr/wiki/pmwiki.php/TD/P73TestWithACarbonFiberSource}.\\\\\n'
    openloop_ack = 'Results of the open loop test and the room temperature measurements are from Damien Pr\\^ele\\\\\n'
    
    h=open(texfilename_fullpath,'w')
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
    h.write('\\newcommand{\\comment}[1]{\n\\begin{minipage}[t]{20ex}\n\\setstretch{0.5}\\flushleft\\noindent\n#1\n\\vspace*{1ex}\n\\end{minipage}}\n')
    h.write('\\newlength{\\openlooplen}\n')
    h.write('\\settowidth{\\openlooplen}{ooop}\n')
    h.write('\\newlength{\\cflen}\n')
    h.write('\\settowidth{\\cflen}{carbon}\n')
    h.write('\\newcommand{\\openloopheading}{\n\\begin{minipage}[t]{\\openlooplen}\nopen\\\\\nloop\n\\end{minipage}\n}\n')
    h.write('\\newcommand{\\cfheading}{\n\\begin{minipage}[t]{\\cflen}\ncarbon\\\\\nfibre\n\\end{minipage}\n}\n')
    
    h.write('\\begin{document}\n')
    h.write('\\begin{center}\n')
    h.write('QUBIC TES Report\\\\\n')
    h.write(self.obsdate.strftime('data from %Y-%m-%d %H:%M UTC\\\\\n'))
    h.write('compiled by %s\\\\\nusing QubicPack: \\url{https://github.com/satorchi/qubicpack}\n' % observer)
    h.write(dt.datetime.utcnow().strftime('this report compiled %Y-%m-%d %H:%M UTC\\\\\n'))
    h.write('\\end{center}\n')

    h.write('\\vspace*{3ex}\n')
    h.write('\\noindent Summary:\n')
    h.write('\\noindent\\begin{itemize}\n')
    h.write('\\item Array %s\n' % self.detector_name)
    h.write('\\item ASIC %i\n' % self.asic)
    if self.temperature is None:
        tempstr='unknown'
    else:
        tempstr=str('%.0f mK' % (1000*self.temperature))
    h.write('\\item TES physical temperature: %s\n' % tempstr)

    h.write('\\item %i pixels are flagged as bad.\n\\item %.1f\\%s of the array is good\n'\
            % ( self.NPIXELS-self.ngood(), 100.0*self.ngood()/self.NPIXELS, '%' ))
    h.write('\\end{itemize}\n')
    
    h.write('\n\\vspace*{3ex}\n\\noindent This document includes the following:\n')
    h.write('\\begin{itemize}\n')
    h.write('\\item Summary Table including turnover points and other parameters for each TES\n')
    h.write('\\item Plot of all the I-V curves, each in its corresponding location in the focal plane\n')
    h.write('\\item Plot of all the good I-V curves on a single plot\n')
    h.write('\\item Plot of each TES I-V curve (%i plots)\n' % self.NPIXELS)
    h.write('\\end{itemize}\n\\clearpage\n')

    ncols=1
    nrows=int(self.NPIXELS/ncols)
    colfmt='|r|r|r|r|r|l|l|l|'
    headline1='\\multicolumn{1}{|c|}{TES} & '\
               '\\multicolumn{1}{|c|}{pix} & '\
               '\\multicolumn{1}{c|}{V$_{\\rm turnover}$} & '\
               '\\multicolumn{1}{c|}{R$_1$} & '\
               '\\multicolumn{1}{c|}{R$_{\\rm 300K}$} & '\
               '\\multicolumn{1}{c|}{\\openloopheading} &'\
               '\\multicolumn{1}{c|}{\\cfheading} &'\
               '\\multicolumn{1}{c|}{comment}'
    headline=''
    headline+=headline1
    if ncols>1:
        for j in range(ncols-1):
            colfmt+='|||r|r|'
            headline+=' & '+headline1 
    h.write('\\noindent\\begin{longtable}{%s}\n' % colfmt)
    h.write('\\caption{Summary Table for TES\\\\\n')
    if self.transdic is not None:
        h.write(cfibre_ack)
        h.write(openloop_ack)
    h.write('}\\\\\n\\hline\n')
    h.write(headline+'\\\\ \n')
    h.write('\\hline\\endhead\n')
    h.write('\\hline\\endfoot\n')
    for i in range(nrows):
        for j in range(ncols):
            TES_idx=i+j*nrows
            TES=TES_idx+1
            rowstr=self.iv_tex_table_entry(TES)
            h.write(rowstr)
            if j<ncols-1: h.write(' &')
            else: h.write('\\\\\n')
    h.write('\\hline\n')
    h.write('\\end{longtable}\n\\clearpage\n')


    # make a table of disagreement
    if self.transdic is not None:
        h.write('\\noindent\\begin{longtable}{%s}\n' % colfmt)
        h.write('\\caption{Table of Disagreement\\\\\n')
        h.write(cfibre_ack)
        h.write(openloop_ack)
        h.write('}\\hline\n')
        h.write(headline+'\\\\ \n')
        h.write('\\hline\\endhead\n')
        h.write('\\hline\\endfoot\n')
        for TES_idx in range(self.NPIXELS):
            TES=TES_idx+1
            PIX=tes2pix(TES,self.asic)
            entry=self.lookup_TEStable(key='PIX',value=PIX)

            if entry['CarbonFibre']=='good':
                is_good_CF=True
            else:
                is_good_CF=False

            if entry['OpenLoop']=='good':
                is_good_OL=True
            else:
                is_good_OL=False

            R300=entry['R300']
            is_good_R300=False
            if isinstance(entry['R300'],float):
                if abs(R300)<10000: is_good_R300=True
                
            '''
            if ((self.is_good_iv(TES) and (not is_good_CF
                                           or not is_good_OL
                                           or not is_good_R300))
                or (not self.is_good_iv(TES) and (is_good_CF
                                                  or is_good_OL
                                                  or is_good_R300))):
            '''
            if ((self.is_good_iv(TES) and not is_good_CF)\
                or (not self.is_good_iv(TES) and is_good_CF)):
                
                rowstr=self.iv_tex_table_entry(TES)
                h.write(rowstr)
                h.write('\\\\\n')
                
        h.write('\\hline\n')
        h.write('\\end{longtable}\n\\clearpage\n')
    

    if tableonly:
        h.write('\n\n\\end{document}\n')
        h.close()
        return texfilename_fullpath

    h.write('\n\\noindent')

    #includefmt = '\n\\includegraphics[width=0.8\\linewidth,natwidth=%i,natheight=%i,clip]{%s}\\\\'
    includefmt = '\n\\includegraphics[width=0.8\\linewidth,clip]{%s}\\\\'
    #img = Image.open(self.output_filename(thumbnailplot))
    #natwidth,natheight = img.size
    #img.close()
    #h.write(includefmt % (natwidth,natheight,thumbnailplot))
    h.write(includefmt % thumbnailplot)

    #img = Image.open(self.output_filename(allplot))
    #natwidth,natheight = img.size
    #img.close()
    #h.write(includefmt % (natwidth,natheight,allplot))
    h.write(includefmt % allplot)
    
    h.write('\n\\clearpage\n\\noindent')
    for png in iv_plots:
        #img = Image.open(self.output_filename(png))
        #natwidth,natheight = img.size
        #img.close()
        #h.write(includefmt % (natwidth,natheight,png))
        h.write(includefmt % png)
    
    h.write('\n\n\\end{document}\n')
    h.close()
    return texfilename_fullpath


def make_iv_report(self):
    '''
    do all the business to generate the I-V report document
    '''
    if not self.exist_iv_data():return None

    # plot all the I-V in the focal-plane map
    self.figsize=(14,14)
    self.plot_iv_physical_layout(xwin=False)

    # plot all the good I-V curves on a single plot
    self.plot_iv_all(selection=self.is_good_iv(),xwin=False)

    # plot each I-V curve
    self.figsize=(16,12)
    for TES_idx in range(self.NPIXELS):
        self.plot_iv(TES_idx+1,xwin=False)

    # generate the LaTeX file
    texname=self.make_iv_tex_report()
    if not isinstance(texname,str):return None

    # process the LaTeX file a couple of times
    cwd=os.getcwd()
    subdir=self.data_subdir()
    if isinstance(subdir,str):
        workdir='%s/%s' % (self.datadir,subdir)
    else:
        workdir=self.datadir
    os.chdir(workdir) # move to data directory
    cmd='pdflatex %s' % texname
    os.system(cmd)
    os.system(cmd)
    os.chdir(cwd) # and return to previous directory
    pdfname=texname.replace('.tex','.pdf')
    return

def iv2txt(self,TES):
    '''
    extract the I-V data from a given TES to a text file with two columns
    '''
    if not self.exist_iv_data():return None
    
    fname='QUBIC_TES%03i_array-%s_ASIC%i_%.0fmK_IV_%s.txt' % (TES,self.detector_name,self.asic,1000*self.temperature,self.obsdate.strftime('%Y%m%dT%H%M%S'))
    h=open(fname,'w')
    Ites=self.Ites(TES)
    if not isinstance(Ites,np.ndarray):return None
    
    Vtes=self.Vtes(TES)
    for idx in range(len(Ites)):
        h.write('%.6e %.6e\n' % (Vtes[idx],Ites[idx]))
    h.close()
    return fname


###################################################
### helper functions to return info from the filter
###################################################

def filterinfo(self,TES=None):
    '''
    return the filterinfo for a given TES
    '''
    if not self.exist_iv_data():return None

    # if no TES is specified, return the whole list
    if TES is None:
        # if filter has not been run, run it with defaults
        for TES_idx in range(self.NPIXELS):
            f=self.filtersummary[TES_idx]
            if f is None: f=self.filter_iv(TES_idx+1)
        return self.filtersummary

    # if not a valid TES, return None
    if not isinstance(TES,int) or TES<1 or TES>self.NPIXELS:
        print('please enter a valid TES number between 1 and %i.' % self.NPIXELS)
        return None

    # if filter has not been run, run it with defaults
    f=self.filtersummary[TES_index(TES)]
    if f is None: f=self.filter_iv(TES)
    return self.filtersummary[TES_index(TES)]

def assign_filterinfo(self,TES,filterinfo):
    '''
    assign the dictionary of filter info to the filtersummary list
    '''
    self.filtersummary[TES_index(TES)]=filterinfo
    return filterinfo
    
def is_good_iv(self,TES=None):
    '''
    return the judgement about a TES
    if TES is None, return a list of all TES determinations
    '''

    filterinfo=self.filterinfo(TES)
    if filterinfo is None:return False

    if TES is None:
        filtersummary=filterinfo
        is_good=[]
        for finfo in filtersummary:
            is_good.append(finfo['is_good'])
        return is_good
    return filterinfo['is_good']

def good_index(self):
    '''
    return a list of indexes corresponding to the good TES
    '''
    good_index=[]
    for TES_idx in range(self.NPIXELS):
        TES=TES_idx+1
        if self.is_good_iv(TES):good_index.append(TES_idx)
    return good_index
    
def ngood(self):
    '''
    return the number of good TES
    '''    
    ngood=0
    for TES_idx in range(self.NPIXELS):
        TES=TES_idx+1
        if self.is_good_iv(TES):ngood+=1
    return ngood

def turnover(self,TES=None):
    '''
    return the turnover (operation) voltage for the TES
    if TES is None, return a list for all the TES
    '''
    filterinfo=self.filterinfo(TES)
    if filterinfo is None:return None
        
    if TES is None:
        filtersummary=filterinfo
        turnover=[]
        for finfo in filtersummary:
            turnover.append(finfo['fit']['turnover'])
        return turnover
    return filterinfo['fit']['turnover']

def Iturnover(self,TES=None):
    '''
    return the current at turnover (operation) voltage for the TES
    '''
    filterinfo=self.filterinfo(TES)
    if filterinfo is None:return None
        
    if TES is None:return None
    if 'Iturnover' in filterinfo['fit'].keys():                                  
        return filterinfo['fit']['Iturnover']
    return None

def offset(self,TES=None):
    '''
    return the offset current for the TES
    if TES is None, return a list for all the TES
    '''
    filterinfo=self.filterinfo(TES)
    if filterinfo is None:return 0.0
        
    if TES is None:
        filtersummary=filterinfo
        offset=[]
        for finfo in filtersummary:
            offset.append(finfo['fit']['offset'])
        return offset
    return filterinfo['fit']['offset']

def R1adjust(self,TES=None):
    '''
    return the R1 adjustment if any
    '''
    if not self.exist_iv_data():return None

    # if no TES is specified, return the whole list
    if TES is None:
        R1adjust_vector=np.ones(self.NPIXELS)
        # if no R1adjust defined, return 1.0 for each
        for TES_idx in range(self.NPIXELS):
            f=self.filtersummary[TES_idx]
            if (not f is None) and ('R1adjust' in f.keys()):
                R1adjust_vector[TES_idx]=f['R1adjust']
        return R1adjust_vector

    # if not a valid TES, return None
    if not isinstance(TES,int) or TES<1 or TES>self.NPIXELS:
        print('please enter a valid TES number between 1 and %i.' % self.NPIXELS)
        return None

    # if filter has not been run, return 1.0
    f=self.filtersummary[TES_index(TES)]
    if f is None: return 1.0

    if 'R1adjust' in f.keys():
        return self.filtersummary[TES_index(TES)]['R1adjust']
    return 1.0


def R1(self,TES=None):
    '''
    return the dynamic normal resistance for the TES
    if TES is None, return a list for all the TES
    '''
    filterinfo=self.filterinfo(TES)
    if filterinfo is None:return None
    if TES is None:
        filtersummary=filterinfo
        R1=[]
        for finfo in filtersummary:
            R1.append(finfo['fit']['R1'])
        return R1
    return filterinfo['fit']['R1']

def Rn(self,TES=None):
    '''
    this is an alias for R1
    '''
    return self.R1(TES)

def selected_iv_curve(self,TES):
    '''
    return the index end points which selects the I-V cycle of the measurement used in the fit
    '''
    filterinfo=self.filterinfo(TES)
    if filterinfo is None:return None
    
    if 'curve index' in filterinfo['fit'].keys():
        curve_index=filterinfo['fit']['curve index']
    else:
        curve_index=filterinfo['fit']['best curve index']

    npts_curve=filterinfo['fit']['npts_curve']
    istart=curve_index*npts_curve
    iend=istart+npts_curve
    return (istart,iend)

def best_iv_curve(self,TES=None):
    '''
    return the best I-V curve for a TES
    if TES is None, return all best I-V curves
    '''
    filterinfo=self.filterinfo(TES)
    if filterinfo is None:return None
    
    if TES is not None:
        TES_idx = TES_index(TES)
        istart,iend=self.selected_iv_curve(TES)
        bias = self.bias_factor*self.vbias[istart:iend]
        adu = self.adu[TES_idx,istart:iend]
        return bias, adu

    adu_best = []
    bias_best = []
    for TES_idx in range(self.NPIXELS):
        TES = TES_idx + 1
        fit = filterinfo[TES_idx]['fit']
        istart,iend=self.selected_iv_curve(TES)
        bias = self.bias_factor*self.vbias[istart:iend]
        bias_best.append(bias)
        adu_best.append(self.adu[TES_idx,istart:iend])
        
    return bias_best,adu_best
        
def Ites(self,TES):
    '''
    return the TES current in Amps (not in microAmps)
    '''
    filterinfo=self.filterinfo(TES)
    if filterinfo is None:return None
    
    Ites=self.adjusted_iv(TES)*1e-6 # Amps
    return Ites

def Vtes(self,TES):
    '''
    return the Vtes
    '''
    Ites=self.Ites(TES)
    if not isinstance(Ites,np.ndarray):return None
    
    Vtes=self.Vbias2Vtes(self.vbias*self.bias_factor,Ites)
    return Vtes

def Ptes(self,TES):
    '''
    return the power on the TES as a function of bias
    '''
    filterinfo=self.filterinfo(TES)
    if filterinfo is None:return None
    
    Ptes=self.Ites(TES)*self.Vtes(TES)*1e12 # pW
    return Ptes


def Rn_ratio(self,TES):
    '''
    return the ratio of TES resistance to Normal resistance
    '''
    if self.R1(TES) is None:return None

    Rn=self.Vtes(TES)/self.Ites(TES)
    Rn_ratio=100*Rn/self.R1(TES) # percent

    return Rn_ratio

def Pbias(self,TES):
    '''
    find the Pbias at 90% Rn
    '''    
    filterinfo=self.filterinfo(TES)
    if filterinfo is None:return None

    Rn_ratio=self.Rn_ratio(TES)
    if not isinstance(Rn_ratio,np.ndarray):return None

    istart,iend=self.selected_iv_curve(TES)

    Rn_ratio=Rn_ratio[istart:iend]
    Ptes=self.Ptes(TES)
    Ptes=Ptes[istart:iend]
    
    # check that Rn_ratio is increasing
    increasing=np.diff(Rn_ratio).mean()
    if increasing<0:
        Pbias=np.interp(90., np.flip(Rn_ratio,0), np.flip(Ptes,0))
    else:
        Pbias=np.interp(90., Rn_ratio, Ptes)

    return Pbias

