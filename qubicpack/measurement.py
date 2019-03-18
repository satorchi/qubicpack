#!/usr/bin/env python
'''
$Id: measurement.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Tue 20 Feb 2018 10:15:44 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

some functions used to setup a TES measurement
'''
from __future__ import division, print_function
from qubicpack import qubicpack as qp
from qubicpack.dummy_client import dummy_client
import datetime as dt
import subprocess,os,sys,time
import numpy as np
from copy import copy
from satorchipy.datefunctions import *
reload(sys)
sys.setdefaultencoding('utf8')

def get_from_keyboard(msg,default=None):
    ''''
    get interactive input from the keyboard
    '''
    prompt='%s (default: %s) ' % (msg,str(default))
    ans=raw_input(prompt)
    if ans=='':return default
    if type(default)==str:
        return ans
    
    try:
        x=eval(ans)
    except:
        print('invalid response.')
        return None
    return x
    
def read_bath_temperature(qpobject):
    '''
    read the current bath temperature and write it out to the log file
    '''
    qpobject.writelog('reading temperature')
    Tbath=qpobject.oxford_read_bath_temperature()
    if Tbath is None:
        qpobject.writelog('ERROR! Could not read bath temperature.')
        Tbath=qpobject.temperature
    qpobject.writelog('Tbath=%.2f mK' % (1000*Tbath))
    return Tbath

def measurement_defaultparameters(meastype='IV'):
    '''
    the default parameters for a measurement
    '''
    params={}
    return params

def measurement_parseargs(argv):
    '''
    parse the argument list and return a dictionary with parameters for the measurement
    '''
    params={}
    
    for arg in argv:
        if arg.upper()=='--TESTMODE':
            params['TESTMODE']=True
            continue

        if arg.upper().find('--ASIC=')==0:
            params['asic']=int(eval(arg.split('=')[1]))
            continue

        if arg.upper().find('--ARRAY=')==0:
            params['detname']=arg.split('=')[1]
            continue

        if arg.upper().find('--MEASTYPE=')==0:
            params['meastype']=arg.split('=')[1]
            continue

        if arg.upper().find('--MIN_BIAS=')==0:
            params['min_bias']=eval(arg.split('=')[1])
            continue
    
        if arg.upper().find('--MAX_BIAS=')==0:
            params['max_bias']=eval(arg.split('=')[1])
            continue

        if arg.upper().find('--START_TEMP=')==0:
            params['start_temp']=eval(arg.split('=')[1])
            continue
    
        if arg.upper().find('--END_TEMP=')==0:
            params['end_temp']=eval(arg.split('=')[1])
            continue

        if arg.upper().find('--STEP_TEMP=')==0:
            params['step_temp']=eval(arg.split('=')[1])
            continue

        if arg.upper().find('--MONITOR_TES=')==0:
            params['monitor_TES']=int(eval(arg.split('=')[1]))
            continue

        if arg.upper().find('--TIMEOUT=')==0:
            temp_timeout_secs=eval(arg.split('=')[1])
            params['temp_timeout']=dt.timedelta(seconds=temp_timeout_secs)
            continue

        if arg.upper().find('--MINWAIT=')==0:
            temp_minwait_secs=eval(arg.split('=')[1])
            params['temp_minwait']=dt.timedelta(seconds=temp_minwait_secs)
            continue

        if arg.upper().find('--WAIT=')==0:
            temp_wait_secs=eval(arg.split('=')[1])
            params['temp_wait']=dt.timedelta(seconds=temp_wait_secs)
            continue
    
        if arg.upper().find('--PRECISION=')==0:
            params['temp_precision']=eval(arg.split('=')[1])
            continue

        if arg.upper().find('--CYCLE')==0:
            params['cycle_temp']=True
            continue

    return params


def measurement_askargs(params=None):
    '''
    get parameters from keyboard if not already specified
    '''
    if params is None:params={}
    
    if not 'detname' in params.keys() or params['detname'] is None:
        detname=get_from_keyboard('Which array is it? ','P90')

    if not 'asic' in params.keys() or params['asic'] is None:
        params['asic']=get_from_keyboard('Which ASIC?  ',2)

    # setup bias voltage range
    # don't ask if it's an R-T measurement
    if (not 'meastype' in params.keys() or not params['meastype']=='RT')\
       and (not 'min_bias' in params.keys() or params['min_bias'] is None):
        params['min_bias']=get_from_keyboard('minimum bias voltage ',3.5)

    
if not meastype=='RT' and max_bias is None:
    max_bias=get_from_keyboard('maximum bias voltage ',max_possible_bias)
    if max_bias is None:quit()
go.max_bias=max_bias

# setup temperature range

# cycle temperature back to start?
if cycle_temp is None:
    ans=get_from_keyboard('cycle temperature back to start (y/n)?','y')
    if ans.upper()=='Y':
        cycle_temp=True
    else:
        cycle_temp=False

if start_temp is None:
    start_temp=get_from_keyboard('start bath temperature ',0.6)
    if start_temp is None:quit()
if end_temp is None:    
    end_temp=get_from_keyboard('end bath temperature ',0.3)
    if end_temp is None:quit()
if step_temp is None:
    step_temp_default=(end_temp-start_temp)/8.
    step_temp=get_from_keyboard('temperature steps',step_temp_default)
    if step_temp is None:quit()

# make sure steps are negative if we're going down in temperature
if start_temp>end_temp:
    if step_temp>0:step_temp=-step_temp
else:
    if step_temp<0:step_temp=-step_temp

Tbath_target=list(np.arange(start_temp,end_temp,step_temp))
if cycle_temp:
    T_return=copy(Tbath_target)
    T_return.reverse()
    Tbath_target=Tbath_target+T_return
    
if monitor_TES is None:    
    monitor_TES=get_from_keyboard('which TES would you like to monitor during the measurement? ',monitor_TES_default)
    if monitor_TES is None:quit()

# if running in test mode, use a random generated result
if TESTMODE:
    if meastype=='IV': go.adu=np.random.rand(go.NPIXELS,len(go.vbias))
    go.temperature=0.3
    go.nsamples=100
    go.chunk_size=300
    go.OxfordInstruments_ip='127.0.0.1'


# parameters specific to measurement type
params_IV={}
params_IV['meastype']='IV'
params_IV['timeline_period']=240.0
params_IV['frequency']=99
params_IV['PID_I']=20
params_IV['min_bias']=min_bias
params_IV['max_bias']=max_bias

params_RT={}
params_RT['meastype']='RT'
params_RT['timeline_period']=60.0
params_RT['frequency']=10.0
params_RT['PID_I']=50
params_RT['min_bias']=-0.5
params_RT['max_bias']= 0.5


# check if we need a second qubicpack object
if meastype is None:
    meastype=get_from_keyboard('Which type of measurement (IV or RT)? ','IV')
meastype=meastype.upper()    
if meastype=='BOTH':
    meastype_str='I-V and R-T'
    
