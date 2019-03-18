'''
$Id: squids.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>

     modified from the squids2.py by Michel Piat


$created: Mon 11 Sep 2017 19:24:29 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

test the squid
'''
from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np

def squid_test(self,vmin=0.0,vmax=15.0,dv=1.0,tinteg=None):
    '''
    test the squids
    '''

    client = self.connect_QubicStudio()
    if client is None: return None
    
    # ask for return values to be voltages
    client.sendSetScientificDataTfUsed(1)

    self.assign_integration_time(tinteg)

    voltages=np.arange(vmin,vmax+dv,dv)
    
    mean_SQUIDs = np.zeros((self.NPIXELS, voltages.size))
    min_SQUIDs  = np.zeros((self.NPIXELS, voltages.size))
    max_SQUIDs  = np.zeros((self.NPIXELS, voltages.size))
    for idx, bias in enumerate(voltages):
        client.sendSetAsicVicm(self.QS_asic_index, 3)
        client.waitMs(500)
        client.sendSetAsicSpol(self.QS_asic_index, bias)
        client.waitMs(500)
        timeline = self.integrate_scientific_data()
        mean_SQUIDs[:, idx] = timeline.mean(axis=-1)
        min_SQUIDs[:, idx] = timeline.min(axis=-1)
        max_SQUIDs[:, idx] = timeline.max(axis=-1)


    delta_SQUIDs = max_SQUIDs - min_SQUIDs

    # Recherche max
    indexes_max = np.argmax(delta_SQUIDs, axis=-1)
    for idx, maxidx in enumerate(indexes_max):
        print('===>Pixel {:3}: delta max={}mV at bias={}'.format(idx, delta_SQUIDs[idx, maxidx]*1000, voltages[maxidx]))


    '''
    for idx in range(len(delta_SQUIDs)):
        plt.plot(delta_SQUIDs[idx])
    '''


    return delta_SQUIDs
