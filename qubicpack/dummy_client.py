'''
$Id: dummy_client.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Sun 18 Feb 2018 08:11:29 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

a dummy pystudio client to use for testing
'''
from __future__ import division, print_function
import numpy as np

class dummy_client:

    def __init__(self):
        self.connected=True
        return

    def connect_QubicStudio(self):
        return True
    
    def fetch(self,param_str):
        print('dummy fetch: %s' % param_str)
        return np.ones((1))

    def request(self,param_str):
        print('dummy request: %s' % param_str)
        return np.ones((1))

    def sendSetTESDAC(self, asic_index, shape, frequency, DACamplitude, DACoffset):
        return True

    def sendActivatePID(self,asic_index,stat):
        return True

    def sendConfigurePID(self,asic_index, P,I,D):
        return True
    
    def sendSetFeedbackRelay(self, asic_index,stat):
        return True

    def sendSetOffsetTable(self, asic_index, offsets):
        return True

    
