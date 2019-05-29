'''
$Id: fpmethods.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Tue 28 May 2019 12:43:24 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

methods for the QUBIC focal plane class qubicfp
'''
from __future__ import division, print_function

import numpy as np
import datetime as dt
import sys,os,time

from qubicpack.qubicasic import qubicasic
from qubicpack.utilities import NASIC
from qubicpack.plot_fp import plot_fp

def assign_defaults(self):
    '''default values for object variables
    '''
    self.obsdate = None
    self.logfile = None
    self.figsize = (12.80,6.40)
    self.colours = ['blue','green','red','cyan','magenta','yellow','black']
    self.observer = 'APC LaboMM'
    self.tdata = [{}]
    self.temperature = None
    self.hk  =  {}
    return

def infotext(self):
    '''
    some basic info to put on plots
    '''
    txt = []
    for asic_obj in self.asic_list:
        if asic_obj is not None:
            txt.append(asic_obj.infotext())
    infotxt = '\n'.join(txt)
    return infotxt

    
def read_qubicstudio_science_fits(self,hdu):
    '''
    read the science data for an ASIC
    The HDU passed here as the argument should already have been identified as the Science HDU
    '''
    self.printmsg('DEBUG: read_qubicstudio_science_fits object type is %s' % self.__object_type__,verbosity=3)
        
    asic_no = hdu.header['ASIC_NUM']
    asic_ctr = asic_no - 1
    qubicasic.verbosity = self.verbosity
    self.asic_list[asic_ctr] = qubicasic()
    self.asic_list[asic_ctr].read_qubicstudio_science_fits(hdu)
    obsdate = self.asic_list[asic_ctr].obsdate
    if self.obsdate is None:
        self.printmsg('DEBUG: setting obsdate which was None: %s' % obsdate.strftime('%Y-%m-%d %H:%M:%S.%f'),verbosity=3)
        self.obsdate = obsdate
    if self.obsdate<>obsdate:
        self.printmsg('PROBLEM! Observation date does not correspond between ASIC data:',verbosity=3)
        for idx,asic_obj in enumerate(self.asic_list):
            if asic_obj is not None:
                self.printmsg('ASIC%i: %s' % (idx+1,asic_obj.obsdate.strftime('%Y-%m-%d %H:%M:%S.%f')),verbosity=3)

    self.printmsg('Observation date: %s' % obsdate.strftime('%Y-%m-%d %H:%M:%S.%f'))
    return

def read_qubicstudio_asic_fits(self,hdulist):
    '''
    read the data giving the ASIC configuration
    The HDU passed here as the argument should already have been identified as the ASIC HDU

    we should read the science data first, and then read the corresponding ASIC table
    so we read this file for each of the defined asic objects
    '''
    for asic_obj in self.asic_list:
        if asic_obj is not None:
            asic_obj.read_qubicstudio_asic_fits(hdulist)
    return


def args_ok(self,TES=None,asic=None):
    '''
    check if arguments are okay for the wrapper
    '''
    if asic is None:
        self.printmsg('Please give an asic number')
        return False

    if TES is None:
        self.printmsg('Please give a TES number')
        return False

    return True
    

#### timeline methods
def sample_period(self,asic=None):
    '''
    wrapper to get the sample period for an asic
    '''
    if asic is None:
        self.printmsg('Please enter an asic number')
        return None

    asic_idx = asic-1
    return self.asic_list[asic_idx].sample_period()

def timeline_array(self,asic=None):
    '''
    wrapper to get the timeline array for an asic
    '''
    if asic is None:
        self.printmsg('Please enter an asic number')
        return None

    asic_idx = asic-1
    return self.asic_list[asic_idx].timeline_array()
    

def plot_timeline(self,TES=None,asic=None):
    '''
    wrapper to plot timeline of the asic object
    '''

    if not self.args_ok(TES,asic):return
    
    asic_idx = asic - 1
    if self.asic_list[asic_idx] is None:
        self.printmsg('No data for ASIC %i' % asic)
        return

    return self.asic_list[asic_idx].plot_timeline(TES=TES)
    
def plot_timeline_focalplane(self):
    '''
    plot all the timelines in the focal plane
    '''

    args= {}
    args['title'] = 'QUBIC Focal Plane: %s' % self.dataset_name
    subttl_list = []
    obsdates = []
    for idx,asic_obj in enumerate(self.asic_list):
        if asic_obj is not None:
            key = 'ASIC%i' % (idx+1)
            
            subttl_list.append(asic_obj.infotext())

            args[key] = asic_obj.timeline_array()
            obsdates.append(asic_obj.obsdate)
            
    args['subtitle'] = '\n'.join(subttl_list)
    args['obsdate'] = min(obsdates)

    return plot_fp(args)


#### I-V methods 
def plot_iv(self,TES=None,asic=None):
    '''
    wrapper to plot I-V of the asic object
    '''

    if not self.args_ok(TES,asic):return
    
    asic_idx = asic - 1
    if self.asic_list[asic_idx] is None:
        self.printmsg('No data for ASIC %i' % asic)
        return

    return self.asic_list[asic_idx].plot_iv(TES=TES)

def plot_iv_focalplane(self):
    '''
    plot all the I-V curves in the focal plane
    '''

    args= {}
    args['title'] = 'QUBIC Focal Plane I-V curves: %s' % self.dataset_name
    subttl_list = []
    obsdates = []
    for idx,asic_obj in enumerate(self.asic_list):
        if asic_obj is not None:
            key = 'ASIC%i' % (idx+1)
            
            subttl_list.append(asic_obj.infotext())

            args[key] = asic_obj.adu
            obsdates.append(asic_obj.obsdate)

            keyx = '%s x-axis' % key
            args[keyx] = asic_obj.vbias

            keygood = '%s good' % key
            args[keygood] = asic_obj.is_good_iv()

            keybg = '%s bg' % key
            args[keybg] = asic_obj.turnover()
            
    args['subtitle'] = '\n'.join(subttl_list)
    args['obsdate'] = min(obsdates)

    return plot_fp(args)

    
def filter_iv_all(self,
                  R1adjust=1.0,
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

    for asic_obj in self.asic_list:
        if asic_obj is not None:
            asic_obj.filter_iv_all(R1adjust,residual_limit,abs_amplitude_limit,rel_amplitude_limit,
                                   bias_margin,jumplimit,fitfunction,Vsuper,Vnormal,istart,iend)

    return
    
