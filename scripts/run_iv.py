#!/usr/bin/env python
"""
$Id: run_iv.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Fri 07 Jul 2017 13:45:13 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

wrapper script to run the I-V curve data gathering
"""
from __future__ import division, print_function
from qubicpack import qubicpack as qp
import matplotlib.pyplot as plt
import subprocess,os,sys
reload(sys)
sys.setdefaultencoding('utf8')

go=qp()

# set debuglevel to 1 if you want lots of messages on the screen
go.debuglevel=1

detname=go.get_from_keyboard('Which array is it? ','P87')
go.assign_detector_name(detname)

asic=go.get_from_keyboard('Which ASIC?  ',2)
if asic is None:quit()
ret=go.assign_asic(asic)

# verify that we can get stuff from QubicStudio
ret=go.verify_QS_connection()
if not ret:quit()

# Fri 19 Jan 2018 13:45:44 CET: we have removed the 5x bias factor
#min_bias=go.get_from_keyboard('minimum bias voltage ',0.5)
go.max_permitted_bias=10.0
min_bias=go.get_from_keyboard('minimum bias voltage ',3.0)
if min_bias is None:quit()
#max_bias=go.get_from_keyboard('maximum bias voltage ',3.0)
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
ncycles=go.get_from_keyboard('number of bias cycles ',3)
if ncycles is None:quit()
monitor_TES=go.get_from_keyboard('which TES would you like to monitor during the measurement? ',70)
if monitor_TES is None:quit()

go.make_Vbias(vmin=min_bias,vmax=max_bias,cycle=cyclebias,ncycles=ncycles,dv=dv)

# run the measurement
go.get_iv_data(TES=monitor_TES)

# generate the test document
pdfname=go.make_iv_report()

# find pdf viewer
viewers=['xpdf','evince','okular','acroread']
use_viewer=None
for viewer in viewers:
    cmd='which %s' % viewer
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    if not out=='':
        use_viewer=out.strip()
        break
        
cmd='%s %s' % (use_viewer,pdfname)
if not pdfname is None and os.path.exists(pdfname) and not use_viewer is None:
    os.system(cmd)

    

