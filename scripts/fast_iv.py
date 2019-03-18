#!/usr/bin/env python
'''
$Id: fast_iv.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Mon 30 Oct 2017 08:31:17 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

wrapper script to run an I-V measurement using the "fast" method
'''
from __future__ import division, print_function
from qubicpack import qubicpack as qp
import matplotlib.pyplot as plt
import datetime as dt
import subprocess,os,sys,time
import numpy as np
reload(sys)
sys.setdefaultencoding('utf8')
from satorchipy.datefunctions import *

go=qp()
go.debuglevel=1

detname=go.get_from_keyboard('Which array is it? ','P90')
go.assign_detector_name(detname)

# can I get ASIC from QubicStudio?
asic=go.get_from_keyboard('Which ASIC?  ',2)
if asic is None:quit()
ret=go.assign_asic(asic)

# setup bias voltage range
# Mon 22 Jan 2018 08:46:16 CET: we have removed the 5x bias factor
go.max_permitted_bias=10.0
min_bias=go.get_from_keyboard('minimum bias voltage ',3.0)
if min_bias is None:quit()
max_possible_bias=go.DAC2V * 2**15
max_bias=go.get_from_keyboard('maximum bias voltage ',max_possible_bias)
if max_bias is None:quit()

monitor_TES=go.get_from_keyboard('which TES would you like to monitor during the measurement? ',34)
if monitor_TES is None:quit()


go.configure_PID()
go.assign_integration_time(1.0) # int time 1sec for offset calculation
if not go.compute_offsets():quit()
if not go.feedback_offsets():quit()
go.assign_integration_time(240.0) # int time 4 minutes for I-V from timeline
go.get_iv_timeline(vmin=min_bias,vmax=max_bias)

go.timeline2adu(monitor_TES)

go.write_fits()
#go.plot_iv(monitor_TES)
