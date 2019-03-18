'''
$Id: feedback_offsets.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Tue 23 Jan 2018 15:56:06 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

compute the feedback offsets and apply them
once tested, this will be moved into qubicpack acquisition.py
'''
from qubicpack import qubicpack as qp
import numpy as np
import time

go=qp()
go.assign_asic(1)
go.assign_integration_time(1.0)

client=go.connect_QubicStudio()
if client is None:quit()

count=10
consigne=0

## switch off the feedback loop
client.sendActivatePID(go.QS_asic_index,0)

# make sure relay=10kOhm  val=1 -> 10kOhm, val=0 -> 100kOhm
client.sendSetFeedbackRelay(go.QS_asic_index,1)


# set sampling frequency 10Hz
freq=10.
# set sampling amplitude 0.0V
amplitude=0.0
# set sampling offset 6V
bias=6.0
# set shape to sinus
shape=0
go.set_VoffsetTES(bias, amplitude, freq, shape)

# to begin, assign zero offset
offsets = np.zeros(go.NPIXELS)
client.sendSetFeedbackTable(go.QS_asic_index, offsets)
go.wait_a_bit(1.0)

## switch on the feedback loop
go.configure_PID(P=0,I=10,D=0)
go.wait_a_bit(5.0)
go.assign_pausetime(0.5)

# set the running average to zero
data_avg=np.zeros(go.NPIXELS)

# correction direction changes with ASIC
if go.QS_asic_index==0:
    correction_direction = 1
else:
    correction_direction = -1
    

k=1.0 # the first step is big
for counter in range(count):

    print('count %i/%i: integrating...' % (counter+1,count))
    timeline = go.integrate_scientific_data()
    print('count %i/%i: finished integrating' % (counter+1,count))
    this_data_avg=timeline.mean(axis=-1)
    prev_offsets=offsets
    offsets = correction_direction*k*(this_data_avg-consigne)+prev_offsets
    print('count %i/%i: applying feedback offsets...' % (counter+1,count))
    client.sendSetFeedbackTable(go.QS_asic_index, offsets)
    print('count %i/%i: feedback offsets applied.' % (counter+1,count))
    go.wait_a_bit()

    data_avg+=this_data_avg
    print('count %i/%i: data for TES 37: %.5e' % (counter+1,count,data_avg[36]))
    k=0.2 # and subsequent steps are smaller
data_avg=data_avg/count
