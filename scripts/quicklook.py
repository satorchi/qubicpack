#!/usr/bin/env python3
'''
$Id: quicklook.py<scripts>
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Mon 23 Dec 2019 08:59:23 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

find the latest data and make a quicklook plot
'''
from glob import glob
from qubicpack.qubicfp import qubicfp
qubicfp.verbosity=1
datadir = '/archive'

datasets = glob('%s/20??-??-??/*' % datadir)
datasets.sort()

latest_dset = datasets[-1]

a = qubicfp()
a.read_qubicstudio_dataset(latest_dset)
a.quicklook()

prmpt = 'Press return to exit. '
ans = input(prmpt)
