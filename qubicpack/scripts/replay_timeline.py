#!/usr/bin/env python
'''
$Id: replay_timeline.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Wed 19 Jul 2017 16:21:01 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

wrapper script to replay the QUBIC TES timeline data
read the data from a FITS  file
'''
from qubicpack import qubicpack as qp
import matplotlib.pyplot as plt
go=qp()
go.read_fits('QUBIC_timeline_20170713T151340UTC.fits')
go.plot_ASD(replay=True,TES=71)

raw_input('Hit return to exit. ')
