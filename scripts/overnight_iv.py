#!/usr/bin/env python
'''
$Id: overnight_iv.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Mon 25 Sep 2017 16:59:04 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

wrapper script to run several I-V measurements, changing the bath temperature for each one
'''
from __future__ import division, print_function
import os
if 'SSH_CLIENT' in os.environ.keys():
    import matplotlib
    matplotlib.use('Agg')
from qubicpack import qubicpack as qp
import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
import subprocess,os,sys,time
reload(sys)
sys.setdefaultencoding('utf8')
from satorchipy.datefunctions import *

# set TESTMODE to False for a real measurement (default)
TESTMODE=False
if len(sys.argv)>1:
    if sys.argv[1].upper()=='--TESTMODE':
        TESTMODE=True

# precision required for bath temperature
temp_precision=0.005 # in Kelvin

# timeout for waiting for temperature to settle
if TESTMODE:
    temp_minwait=dt.timedelta(seconds=30)
    temp_timeout=dt.timedelta(seconds=60)
    temp_wait=dt.timedelta(seconds=1)
    wait_msg='waiting %.0f seconds for temperature to settle' % tot_seconds(temp_wait)
else:
    temp_minwait=dt.timedelta(minutes=30)
    temp_minwait=dt.timedelta(minutes=5)
    temp_timeout=dt.timedelta(minutes=60)
    temp_wait=dt.timedelta(minutes=1)
    wait_msg='waiting %.1f minutes for temperature to settle' % (tot_seconds(temp_wait)/60.)


def read_bath_temperature(qpobject,logfile):
    Tbath=qpobject.oxford_read_bath_temperature()
    if Tbath is None:
        qpobject.writelog(logfile,'ERROR! Could not read bath temperature.')
        Tbath=qpobject.temperature
    return Tbath

go=qp()
figsize=go.figsize

# set debuglevel to 1 if you want lots of messages on the screen
# go.debuglevel=1


'''
get parameters
'''
detname=go.get_from_keyboard('Which array is it? ','P87')
go.assign_detector_name(detname)

# can I get ASIC from QubicStudio?
asic=go.get_from_keyboard('Which ASIC?  ',2)
if asic is None:quit()
ret=go.assign_asic(asic)

# verify that we can get stuff from QubicStudio
if not TESTMODE:
    ret=go.verify_QS_connection()
    if not ret:quit()


# setup bias voltage range
# Mon 22 Jan 2018 08:46:16 CET: we have removed the 5x bias factor
go.max_permitted_bias=10.0

min_bias=go.get_from_keyboard('minimum bias voltage ',3.5)
if min_bias is None:quit()
max_bias=go.get_from_keyboard('maximum bias voltage ',9.0)
if max_bias is None:quit()
default_dv=(max_bias-min_bias)/300.0
dv=go.get_from_keyboard('bias step size ',default_dv)
cycle=go.get_from_keyboard('cycle bias up/down? ','y')
if cycle is None:quit()
if cycle.upper()=='N':
    cyclebias=False
else:
    cyclebias=True
ncycles=go.get_from_keyboard('number of bias cycles ',1)
if ncycles is None:quit()
monitor_TES=go.get_from_keyboard('which TES would you like to monitor during the measurement? ',82)
if monitor_TES is None:quit()

go.make_Vbias(vmin=min_bias,vmax=max_bias,cycle=cyclebias,ncycles=ncycles,dv=dv)

tinteg=go.get_from_keyboard('integration time ',1.0)
go.assign_integration_time(tinteg)

# setup temperature range
start_temp=go.get_from_keyboard('start bath temperature ',0.6)
if start_temp is None:quit()
end_temp=go.get_from_keyboard('end bath temperature ',0.3)
if end_temp is None:quit()
step_temp=go.get_from_keyboard('temperature steps',0.025)
if step_temp is None:quit()

# make sure steps are negative if we're going down in temperature
if start_temp>end_temp:
    if step_temp>0:step_temp=-step_temp
else:
    if step_temp<0:step_temp=-step_temp

Tbath_target=np.arange(start_temp,end_temp,step_temp)

# if running in test mode, use a random generated result
if TESTMODE:
    go.adu=np.random.rand(go.NPIXELS,len(go.vbias))
    go.temperature=0.3
    go.nsamples=100
    go.OxfordInstruments_ip='127.0.0.1'

# make a log file
logfile=dt.datetime.utcnow().strftime('temperature_IV_logfile_%Y%m%dT%H%M%SUTC.txt')
logfile_fullpath=go.output_filename(logfile)

go.writelog(logfile_fullpath,'starting I-V measurements at different temperatures using the stepped bias (slow) method')
go.writelog(logfile_fullpath,'ASIC=%i' % go.asic)
go.writelog(logfile_fullpath,'minimum bias=%.2f V' % min_bias)
go.writelog(logfile_fullpath,'maximum bias=%.2f V' % max_bias)
go.writelog(logfile_fullpath,'start temperature=%.3f K' % start_temp)
go.writelog(logfile_fullpath,'end temperature=%.3f K' % end_temp)
go.writelog(logfile_fullpath,'temperature step=%.3f K' % step_temp)
nsteps=len(Tbath_target)
go.writelog(logfile_fullpath,'number of temperatures=%i' % nsteps)

# estimated time: half hour for temperature to settle, 25 minutes for I-V measurement
duration_estimate=nsteps*(temp_minwait+dt.timedelta(minutes=25))
endtime_estimate=dt.datetime.utcnow()+duration_estimate
go.writelog(logfile_fullpath,endtime_estimate.strftime('estimated end at %Y-%m-%d %H:%M:%S'))


# run the measurement
for T in Tbath_target:
    # set the desired bath temperature
    cmdret=go.oxford_set_point(T)
    # make sure the set point was accepted
    Tsetpt=go.oxford_read_set_point()
    if Tsetpt is None:
        go.writelog(logfile_fullpath,'ERROR! Could not read set point temperature.')
        Tsetpt=T
    go.writelog(logfile_fullpath,'Temperature set point = %.2f mK' % (1000*Tsetpt))
    Tbath=read_bath_temperature(go,logfile_fullpath)
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

        go.writelog(logfile_fullpath,wait_msg)
        time.sleep(tot_seconds(temp_wait))
        go.writelog(logfile_fullpath,'reading temperature')
        Tbath=read_bath_temperature(go,logfile_fullpath)
        delta_step=np.abs(Tbath - Tbath_previous)
        Tbath_previous=Tbath
        delta=np.abs(Tbath - T)
        go.writelog(logfile_fullpath,'Tbath=%0.2f mK' %  (1000*go.temperature))

        # check heater percentage
        heatpercent=go.oxford_read_heater_level()
        if heatpercent>99:
            go.writelog(logfile_fullpath,'We need to increase the maximum current to the heater')
            cmdret=go.oxford_increase_heater_range()
            heater=go.oxford_read_heater_range()
            go.writelog(logfile_fullpath,'heater range: %f mA' % heater)
        
    go.writelog(logfile_fullpath,'starting I-V measurement')
    if delta>temp_precision:
        go.writelog(logfile_fullpath,'WARNING! Did not reach target temperature!')
        go.writelog(logfile_fullpath,'Tbath=%0.2f mK, Tsetpoint=%0.2f mK' % (1000*Tbath,1000*T))

    if not TESTMODE:
        # reset FLL before measurement
        go.configure_PID()

        # recalculate the offsets
        go.compute_offsets()

        # and the feedback offsets
        go.feedback_offsets()
        

    # get the I-V curve
    go.get_iv_data(TES=monitor_TES,replay=TESTMODE)
    go.writelog(logfile_fullpath,'end I-V measurement')
    plt.close('all')

    # generate the test document
    go.writelog(logfile_fullpath,'generating test document')
    if not TESTMODE: pdfname=go.make_iv_report()
    go.writelog(logfile_fullpath,'test document generated')

    # reset the plotting figure size
    go.figsize=figsize

    # reset data
    if not TESTMODE:go.adu=None



# finally, switch off the temperature control loop
go.oxford_pidoff()

