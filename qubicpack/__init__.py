"""
$Id: __init__.py<qubicpack>
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>

$created: Mon 10 Jul 2017 11:55:24 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

a class with tools which are generally useful for scripts using pystudio
"""
import numpy as np
import datetime as dt
import sys,os,time

class qubicpack:
    verbosity=0 # class variable.  You can change this before instantiating an object.
    __object_type__ = 'qubicpack'

    
    def __init__(self):
        print('WARNING!  The qubicpack object is DEPRECATED.  Please use qubicfp() instead.')
        return

