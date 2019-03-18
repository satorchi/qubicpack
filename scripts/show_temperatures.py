#!/usr/bin/env python
from qubicpack import qubicpack as qp
import sys,os,time
import datetime as dt

go=qp()
now=dt.datetime.utcnow()
onehour=dt.timedelta(hours=1)
endtime=now+8*onehour
while now<endtime:
    print now.strftime('%Y-%m-%d %H:%M:%S UTC\n')
    print go.oxford_read_all_temperatures()
    time.sleep(10)
    now=dt.datetime.utcnow()

    
