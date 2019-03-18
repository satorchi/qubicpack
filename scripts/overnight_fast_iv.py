#!/usr/bin/env python
'''
$Id: overnight_fast_iv.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Tue 31 Oct 2017 16:42:05 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

wrapper script to run several I-V measurements, changing the bath temperature for each one
we use the "fast I-V" method with the Bias modulated by a sine curve
'''
from __future__ import division, print_function
import os
if 'SSH_CLIENT' in os.environ.keys():
    import matplotlib
    matplotlib.use('Agg')
from qubicpack import qubicpack as qp
from qubicpack.dummy_client import dummy_client
import matplotlib.pyplot as plt
import datetime as dt
import subprocess,os,sys,time
import numpy as np
from copy import copy
reload(sys)
sys.setdefaultencoding('utf8')
from satorchipy.datefunctions import *

def read_bath_temperature(qpobject):
    qpobject.writelog('reading temperature')
    Tbath=qpobject.oxford_read_bath_temperature()
    if Tbath is None:
        qpobject.writelog('ERROR! Could not read bath temperature.')
        Tbath=qpobject.temperature
    qpobject.writelog('Tbath=%.2f mK' % (1000*Tbath))
    return Tbath

# create the  qubicpack object
go=qp()
# set debuglevel to 1 if you want lots of messages on the screen
go.debuglevel=0

# set TESTMODE to False for a real measurement (default)
TESTMODE=False

# default is to ask for parameters
# command line arguments will suppress questions for the corresponding parameters
asic=None
detname=None
meastype=None

# Mon 22 Jan 2018 08:46:16 CET: we have removed the 5x bias factor
go.max_permitted_bias=10.0
min_bias=None
max_possible_bias=go.DAC2V * 2**15
max_bias=None
shape=None

# temperature start, end, and stepsize
start_temp=None
end_temp=None
step_temp=None
cycle_temp=None

# temperature waiting for temperature to settle:
#  timeout = max time to wait
#  minwait = minimum time to wait
#  wait    = time to wait between temperature measurements
temp_timeout=None
temp_minwait=None
temp_wait=None
temp_precision=None
Tcrit=0.375

monitor_TES=None
monitor_TES_default=34


PID_I=None

argv=[]
if len(sys.argv)>1:argv=sys.argv[1:]
for arg in argv:
    if arg.upper()=='--TESTMODE':
        TESTMODE=True
        continue

    if arg.upper().find('--ASIC=')==0:
        asic=int(eval(arg.split('=')[1]))
        continue

    if arg.upper().find('--ARRAY=')==0:
        detname=arg.split('=')[1]
        continue

    if arg.upper().find('--MEASTYPE=')==0:
        meastype=arg.split('=')[1]
        continue

    if arg.upper().find('--MIN_BIAS=')==0:
        min_bias=eval(arg.split('=')[1])
        continue
    
    if arg.upper().find('--MAX_BIAS=')==0:
        max_bias=eval(arg.split('=')[1])
        continue

    if arg.upper().find('--START_TEMP=')==0:
        start_temp=eval(arg.split('=')[1])
        continue
    
    if arg.upper().find('--END_TEMP=')==0:
        end_temp=eval(arg.split('=')[1])
        continue

    if arg.upper().find('--STEP_TEMP=')==0:
        step_temp=eval(arg.split('=')[1])
        continue

    if arg.upper().find('--MONITOR_TES=')==0:
        monitor_TES=int(eval(arg.split('=')[1]))
        continue

    if arg.upper().find('--TIMEOUT=')==0:
        temp_timeout_secs=eval(arg.split('=')[1])
        temp_timeout=dt.timedelta(seconds=temp_timeout_secs)
        continue

    if arg.upper().find('--MINWAIT=')==0:
        temp_minwait_secs=eval(arg.split('=')[1])
        temp_minwait=dt.timedelta(seconds=temp_minwait_secs)
        continue

    if arg.upper().find('--WAIT=')==0:
        temp_wait_secs=eval(arg.split('=')[1])
        temp_wait=dt.timedelta(seconds=temp_wait_secs)
        continue
    
    if arg.upper().find('--PRECISION=')==0:
        temp_precision=eval(arg.split('=')[1])
        continue

    if arg.upper().find('--CYCLE')==0:
        cycle_temp=True
        continue

    if arg.upper().find('--NO-CYCLE')==0:
        cycle_temp=False
        continue

    if arg.upper().find('--shape=')==0:
        shape=eval(arg.split('=')[1])
        continue
    

'''
get parameters from keyboard if not already specified
'''

# precision required for bath temperature
if temp_precision is None:temp_precision=0.002 # in Kelvin

# timeout for waiting for temperature to settle
if temp_minwait is None:temp_minwait=dt.timedelta(seconds=30)
if temp_timeout is None:temp_timeout=dt.timedelta(minutes=10)
if temp_wait    is None:temp_wait=dt.timedelta(seconds=15)
wait_msg='waiting %.0f seconds for temperature to settle' % tot_seconds(temp_wait)
    
if detname is None:
    detname=go.get_from_keyboard('Which array is it? ','P90')
go.assign_detector_name(detname)

# can I get ASIC from QubicStudio?
if asic is None:
    asic=go.get_from_keyboard('Which ASIC?  ',2)
    if asic is None:quit()
ret=go.assign_asic(asic)

# setup bias voltage range
# don't ask if it's an R-T measurement
if not meastype=='RT' and min_bias is None:
    min_bias=go.get_from_keyboard('minimum bias voltage ',3.5)
    if min_bias is None:quit()
go.min_bias=min_bias

if not meastype=='RT' and max_bias is None:
    max_bias=go.get_from_keyboard('maximum bias voltage ',max_possible_bias)
    if max_bias is None:quit()
go.max_bias=max_bias

if shape is None:shape=0 # default to sinusoid without asking
if not shape in [0,1,2]:
    shape=go.get_from_keyboard('enter bias modulation shape:  0) sinusoid 1) triangular 2) continuous',0)
    go.bias_mode=shape
    

# setup temperature range

# cycle temperature back to start?
if cycle_temp is None:
    ans=go.get_from_keyboard('cycle temperature back to start (y/n)?','y')
    if ans.upper()=='Y':
        cycle_temp=True
    else:
        cycle_temp=False

if start_temp is None:
    start_temp=go.get_from_keyboard('start bath temperature ',0.6)
    if start_temp is None:quit()
if end_temp is None:    
    end_temp=go.get_from_keyboard('end bath temperature ',0.3)
    if end_temp is None:quit()
if step_temp is None:
    step_temp_default=(end_temp-start_temp)/8.
    step_temp=go.get_from_keyboard('temperature steps',step_temp_default)
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
    monitor_TES=go.get_from_keyboard('which TES would you like to monitor during the measurement? ',monitor_TES_default)
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
params_IV['shape']=0

params_RT={}
params_RT['meastype']='RT'
params_RT['timeline_period']=60.0
params_RT['frequency']=10.0
params_RT['PID_I']=50
params_RT['min_bias']=-0.5
params_RT['max_bias']= 0.5
params_RT['shape']=0

# check if we need a second qubicpack object
if meastype is None:
    meastype=go.get_from_keyboard('Which type of measurement (IV or RT)? ','IV')
meastype=meastype.upper()    
if meastype=='BOTH':
    params=params_IV # first measurement is I-V
    meastype_str='I-V and R-T'
    go2=qp()
    go2.debuglevel=1
    go2.assign_asic(go.asic)
    go2.assign_detector_name(go.detector_name)
    go2.min_bias=params_RT['min_bias']
    go2.max_bias=params_RT['max_bias']
    if TESTMODE:
        go2.client=dummy_client()
        go2.temperature=0.3
        go2.nsamples=100
        go2.chunk_size=300
        go2.OxfordInstruments_ip='127.0.0.1'
    measurement_period=params_IV['timeline_period']+params_RT['timeline_period']
    
elif meastype=='RT':
    params=params_RT
    meastype_str='R-T'
    measurement_period=params['timeline_period']
else:
    params=params_IV
    meastype_str='I-V'
    measurement_period=params['timeline_period']

# check if we can connect
if TESTMODE:
    go.client=dummy_client()
else:
    ret=go.verify_QS_connection()
    if not ret:quit()

# make a log file
go.assign_logfile('temperature_%s' % meastype)
if meastype=='BOTH':go2.logfile=go.logfile
go.writelog('starting %s measurements at different temperatures using the timeline (fast) method' % meastype_str)
go.writelog('ASIC=%i' % go.asic)
go.writelog('start temperature=%.3f K' % start_temp)
go.writelog('end temperature=%.3f K' % end_temp)
go.writelog('temperature step=%.3f K' % step_temp)
nsteps=len(Tbath_target)
go.writelog('number of temperatures=%i' % nsteps)
    
# estimated time: temperature settle time plus measurement time for I-V
duration_estimate=(nsteps+1)*(temp_minwait+dt.timedelta(seconds=measurement_period))
endtime_estimate=dt.datetime.utcnow()+duration_estimate
go.writelog(endtime_estimate.strftime('estimated end at %Y-%m-%d %H:%M:%S'))

# run the measurement
for T in Tbath_target:
    # set the desired bath temperature
    cmdret=go.oxford_set_point(T)
    # make sure the set point was accepted
    Tsetpt=go.oxford_read_set_point()
    if Tsetpt is None:
        go.writelog('ERROR! Could not read set point temperature.')
        Tsetpt=T
    go.writelog('Temperature set point = %.2f mK' % (1000*Tsetpt))
    Tbath=read_bath_temperature(go)
    Tbath_previous=Tbath
    delta=np.abs(Tbath - T)
    delta_step=np.abs(Tbath - Tbath_previous)
    start_waiting=dt.datetime.utcnow()
    end_waiting=start_waiting+temp_timeout
    min_endtime=start_waiting+temp_minwait
    
    while (delta>temp_precision
          or delta_step>temp_precision\
          or dt.datetime.utcnow()<min_endtime)\
          and dt.datetime.utcnow()<end_waiting:

        go.writelog(wait_msg)
        time.sleep(tot_seconds(temp_wait))
        Tbath=read_bath_temperature(go)
        delta_step=np.abs(Tbath - Tbath_previous)
        Tbath_previous=Tbath
        delta=np.abs(Tbath - T)

        # check heater percentage
        heatpercent=go.oxford_read_heater_level()
        if heatpercent>99:
            go.writelog('We need to increase the maximum current to the heater')
            cmdret=go.oxford_increase_heater_range()
            heater=go.oxford_read_heater_range()
            go.writelog('heater range: %f mA' % heater)
        
    go.writelog('starting %s measurement' % meastype)
    if delta>temp_precision:
        go.writelog('WARNING! Did not reach target temperature!')
        go.writelog('Tbath=%0.2f mK, Tsetpoint=%0.2f mK' % (1000*Tbath,1000*T))

    # reset FLL and re-compute the offsets before measurement
    if not TESTMODE:
        go.assign_integration_time(1.0) # int time 1sec for offset calculation
        go.compute_offsets()
        go.feedback_offsets()
        
        go.configure_PID(I=params['PID_I']) # feedback_offsets() configured the PID.  Now we set it to what we want.
        go.assign_integration_time(params['timeline_period']) 
        go.writelog('minimum bias=%.2f V' % params['min_bias'])
        go.writelog('maximum bias=%.2f V' % params['max_bias'])
        do_measurement=go.get_iv_timeline(vmin=params['min_bias'],
                                          vmax=params['max_bias'],
                                          frequency=params['frequency'],
                                          shape=params['shape'])
        if do_measurement is None:
            go.writelog('ERROR! Did not successfully acquire a timeline!')
        else:
            go.write_fits()
        go.writelog('end %s measurement' % meastype)

        # if we're doing both measurements, the second one is always the RT measurement
        if meastype=='BOTH' and Tbath>Tcrit:
            go2.configure_PID(I=params_RT['PID_I']) 
            go2.assign_integration_time(params_RT['timeline_period'])
            go.writelog('minimum bias=%.2f V' % params_RT['min_bias'])
            go.writelog('maximum bias=%.2f V' % params_RT['max_bias'])
            do_measurement=go2.get_iv_timeline(vmin=params_RT['min_bias'],vmax=params_RT['max_bias'],frequency=params_RT['frequency'])
            if do_measurement is None:
                go.writelog('ERROR! Did not successfully acquire a timeline for %s measurement!' % params_RT['meastype'])
            else:
                go2.write_fits()
            go.writelog('end %s measurement' % params_RT['meastype'])
            
        '''
        # do this in post processing instead
        if meastype=='IV' or meastype=='BOTH':        
            # generate the test document
            go.timeline2adu(monitor_TES)
            go.writelog('generating test document')
            pdfname=go.make_iv_report()
            go.writelog('test document generated')
        '''
    # reset data
    if not TESTMODE:go.adu=None



# finally, switch off the temperature control loop
go.oxford_pidoff()

