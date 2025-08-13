#!/usr/bin/env python
'''
$Id: run_plot_both_asics.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Thu 31 Aug 2017 08:00:02 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

plot all I-V curves on a map of the focal plane from both ASICs
'''
from qubicpack import qubicpack as qp
from qubicpack.plot_physical_layout import *

a1=qp()
a1.read_fits('QUBIC_TES_20170711T151016UTC.fits')
a2=qp()
a2.read_fits('QUBIC_TES_20170712T154242UTC.fits')

ret=plot_physical_layout(a1=a1,a2=a2)
