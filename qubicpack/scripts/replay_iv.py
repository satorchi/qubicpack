#!/usr/bin/env python
'''
$Id: replay_iv.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Wed 19 Jul 2017 16:26:01 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

example wrapper script to replay the QUBIC I-V data
data is read from a FITS file
'''
from qubicpack import qubicpack as qp
import matplotlib.pyplot as plt
go=qp()

go.read_fits('QUBIC_TES_20170712T154242UTC.fits')
go.get_iv_data(replay=True,TES=70)
raw_input('Hit return to exit. ')
