'''
$Id: compute_offsets.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Fri 15 Dec 2017 17:58:35 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

this is a translation of the Dscript: Comp_offset_ASIC1b.dscript
once tested, it will be moved into qubicpack acquisition.py

'''
from qubicpack import qubicpack as qp
import numpy as np


go=qp()
go.assign_asic(1)
go.assign_integration_time(0.5)

client=go.connect_QubicStudio()
if client is None:quit()


count=10
consigne=0

# first switch off the loop
client.sendActivatePID(go.QS_asic_index,0)

# make sure relay=10kOhm  val=1 -> 10kOhm, val=0 -> 100kOhm
client.sendSetFeedbackRelay(go.QS_asic_index,1)

# set sampling frequency 400Hz
freq=400.
# set sampling amplitude 0.1V
amplitude=1.0
# set sampling offset 6V
bias=6.0
# set shape to sinus
shape=0
go.set_VoffsetTES(bias, amplitude, freq, shape)

# to begin, assign zero offset
offsets = np.zeros(go.NPIXELS)
client.sendSetOffsetTable(go.QS_asic_index, offsets)
go.wait_a_bit()

# to begin, get the current offsets
#parameter='QUBIC_OffsetDACValues_%i' % go.QS_asic_index
#offsets=client.fetch(parameter)

# set the running average to zero
data_avg=np.zeros(go.NPIXELS)

k=1.0 # the first step is big
for counter in range(count):

    print('count: %i/%i' % (counter+1,count))
    timeline = go.integrate_scientific_data()
    this_data_avg=timeline.mean(axis=-1)
    prev_offsets=offsets
    offsets=-k*(this_data_avg-consigne)+prev_offsets
    client.sendSetOffsetTable(go.QS_asic_index, offsets)
    go.wait_a_bit()

    data_avg+=this_data_avg
    k=0.5 # and subsequent steps are smaller
data_avg=data_avg/count

# finally, switch back on the FLL
#client.sendActivatePID(go.QS_asic_index,1)
#go.wait_a_bit()


# Damien algo
'''
timeline = go.integrate_scientific_data()
data_avg_0=timeline.mean(axis=-1)
offsets=-(data_avg_0-consigne)

client.sendSetOffsetTable(go.QS_asic_index, offsets)
go.wait_a_bit()

timeline = go.integrate_scientific_data()
data_avg_1=timeline.mean(axis=-1)

offsets=offsets - (data_avg_0-data_avg_1)
client.sendSetOffsetTable(go.QS_asic_index, offsets)
go.wait_a_bit()
'''

#data_avg_0=data_avg_1



