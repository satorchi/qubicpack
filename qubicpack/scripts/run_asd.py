#!/usr/bin/env python
"""
$Id: run_asd.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Mon 17 Jul 2017 20:20:33 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

wrapper script to run the timeline / ASD function
"""
from qubicpack import qubicpack as qp
import matplotlib.pyplot as plt
go=qp()

'''
if ntimelines is a positive integer, the timeline data will be saved
otherwise, the plots are run in an infinite loop in "monitor mode" and no data is saved
to exit from monitor mode, press Ctrl-C at the command line where the script was run (not on the plot)
'''
go.get_ASD(tinteg=1.,TES=70,ntimelines=-1)
