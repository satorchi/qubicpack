#!/usr/bin/env python
'''
$Id: oxford.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Mon 25 Sep 2017 16:26:59 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

control and read output from the Oxford Instruments Triton dilution fridge

see document: "Triton manual issue 3.7.pdf" (page 100)
I put it on the QUBIC wiki.  Here's the link:
http://qubic.in2p3.fr/wiki/uploads/DetectorWorkingGroup/MeasurementsDetectors/Triton%20manual%20issue%203.7.pdf
'''
from __future__ import division, print_function
import numpy as np
import socket,time

def oxford_assign_temperature_labels(self):
    labels=[]
    labels.append('4K Head')
    labels.append('4K Plate')
    labels.append('Still RuO2')
    labels.append('MC Plate cernox')
    labels.append('MC Plate RuO2')
    labels.append('100mK Plate RuO2')
    labels.append('NOT USED')
    labels.append('NOT USED')
    labels.append('70K Head')
    labels.append('70K Plate')
    self.oxford_temperature_labels=labels
    return

def oxford_assign_heater_ranges(self):
    ranges=[]
    ranges.append(0.1)
    ranges.append(0.316)
    ranges.append(1.0)
    ranges.append(3.16)
    self.oxford_heater_ranges=ranges
    return

def oxford_send_cmd(self, cmd=None):
    '''
    send a command to the Oxford Instruments control computer for the dilution fridge
    '''
    if cmd is None:
        cmd='READ:SYS:TIME\n'

    if not isinstance(cmd,str):
        print('please enter a valid command')
        return None

    s=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        s.connect((self.OxfordInstruments_ip, 33576))
    except:
        print('ERROR! Oxford Instruments is not available.')
        return None
    s.send(cmd)

    d=''
    b=''
    count=0
    while not b=='\n':
        b=s.recv(1)
        d+=b
        count+=1

    #print('return length: %i' % count)
    s.close()
    return d.strip('\n').split(':')


def oxford_init(self):
    '''
    initialize the Oxford Instruments computer so it accepts commands
    '''
    return self.oxford_send_cmd('SET:SYS:USER:NORM\n')

def oxford_pidoff(self):
    '''
    switch off the temperature control loop
    '''
    pidoff='SET:DEV:T5:TEMP:LOOP:MODE:OFF\n'
    d=self.oxford_init()
    return self.oxford_send_cmd(pidoff)

def oxford_set_point(self, T=None, heater=None, ramp=0.1):
    '''
    configure the loop set point for bath temperature, and activate the loop
    '''
    if (not isinstance(T,float)) and (not isinstance(T,int)):
        print('ERROR! invalid temperature')
        return None

    # first initialize Oxford Inst.
    d=self.oxford_init()

    # disable all the other sensors
    for idx,label in enumerate(self.oxford_temperature_labels):
        if not (label=='NOT USED' or label=='MC Plate RuO2'):
            cmd='SET:DEV:T%i:TEMP:MEAS:ENAB:OFF\n' % (idx+1)
            d=self.oxford_send_cmd(cmd)

    # enable the bath temperature sensor
    cmd='SET:DEV:T5:TEMP:MEAS:ENAB:ON\n'
    d=self.oxford_send_cmd(cmd)

    # activate the loop:  This must be done first!
    # and then configure the temperature set-point
    cmd ='SET:DEV:T5:TEMP:LOOP:MODE:ON\n'        # ON/OFF
    cmd+='SET:DEV:T5:TEMP:LOOP:TSET:%0.3f\n' % T # K
    d=self.oxford_send_cmd(cmd)

    # wait a second and then activate the heater
    time.sleep(1)
    d=self.oxford_set_heater_range(heater)

    # set the ramp rate for temperature control
    cmdramp='SET:DEV:T5:TEMP:LOOP:RAMP:RATE:%f\n' % ramp # K/min
    d=self.oxford_send_cmd(cmdramp)

    # enable the ramp
    time.sleep(1)
    rampenable='SET:DEV:T5:TEMP:LOOP:RAMP:ENAB:ON\n'
    d=self.oxford_send_cmd(cmdramp)
    
    return d

def oxford_read_set_point(self):
    '''
    read the loop set point for bath temperature
    '''
    cmd='READ:DEV:T5:TEMP:LOOP:TSET\n'
    d=self.oxford_send_cmd(cmd)
    try:
        T=eval(d[-1].replace('K',''))
    except:
        print('ERROR! could not read set point temperature: %s' % d)
        return None
    return T

def oxford_read_temperature(self,chan=5):
    '''
    read the temperature from one of the thermometers
    '''
    if not isinstance(chan,int):
        print('ERROR! invalid thermometer channel.  Enter a number from 1 to 10.')
        return None
    
    cmd='READ:DEV:T%i:TEMP:SIG:TEMP\n' % chan
    d=self.oxford_send_cmd(cmd)
    try:
        T=eval(d[-1].replace('K',''))
    except:
        print('ERROR! could not read temperature: %s' % d)
        return None
    return T

def oxford_check_calibration(self,chan=5):
    '''
    check if a temperature sensor is calibrated
    '''
    if not isinstance(chan,int):
        print('ERROR! invalid thermometer channel.  Enter a number from 1 to 10.')
        return None

    cmd='READ:DEV:T%i:TEMP:CAL:CHK\n' % chan
    d=self.oxford_send_cmd(cmd)
    if d is None:return False
    ok=d[-1]
    if ok=='OK':return True
    return False
    

def oxford_read_bath_temperature(self):
    '''
    read the bath temperature of the dilution fridge
    '''
    cmd='READ:DEV:T5:TEMP:SIG:TEMP\n'
        
    d=self.oxford_send_cmd(cmd)
    try:
        T=eval(d[-1].replace('K',''))
    except:
        print('ERROR! could not read bath temperature: %s' % d)
        return None

    return self.assign_temperature(T)

def oxford_read_all_temperatures(self):
    '''
    read all the temperatures from the dilution fridge
    '''
    temperature_table=''
    for idx,label in enumerate(self.oxford_temperature_labels):
        if not label=='NOT USED':
            chan=idx+1
            val=self.oxford_read_temperature(chan)
            if val is None:
                tempstr='INACCESSIBLE'
            else:
                tempstr='%8.3f K' % val
            calmsg=''
            '''
            calok=self.oxford_check_calibration(chan)
            if calok:
                calmsg=''
            else:
                calmsg=' == UNCALIBRATED =='
            '''
            temperature_table+='T%02i) %s -- %s%s\n' % (chan,tempstr,label,calmsg)

    return temperature_table

def oxford_read_heater_level(self):
    '''
    read the percent level of the heater.  
    This is a percentage of the max current specified in the set heater command
    '''

    # first we read the heater output power in microWatts
    cmd='READ:DEV:H1:HTR:SIG:POWR\n'
    d=self.oxford_send_cmd(cmd)    
    try:
        P=1e-6*eval(d[-1].replace('uW',''))
    except:
        print('ERROR! could not read heater power: %s' % d)
        return None
    self.debugmsg('heater power: %f W' % P)
    
    # next we read the heater resistance in Ohms    
    cmd='READ:DEV:H1:HTR:RES\n'
    d=self.oxford_send_cmd(cmd)
    try:
        R=eval(d[-1])
    except:
        print('ERROR! could not read heater resistance: %s' % d)
        return None
    self.debugmsg('heater resistance: %f Ohm' % R)

    # and we calculate the current with Ohm's law: P=I^2 R
    I=np.sqrt(P/R)
    self.debugmsg('heater current: %f A' % I)

    # and we read the maximum range that we provided in the set_point command
    Imax=self.oxford_read_heater_range()
    if Imax is None:return None
    # convert Imax to Amps
    Imax=0.001*Imax
    
    htrpercent=100.0*I/Imax
    return htrpercent

def oxford_set_heater_range(self,heater=None):
    '''
    set the heater maximum current level
    '''    
    # determine heater level to apply
    if (not isinstance(heater,float)) and (not isinstance(heater,int)):
        heater=self.oxford_determine_best_heater_level()
        if heater is None:return None

    self.debugmsg('setting heater range: %f mA' % heater)
    cmdheat='SET:DEV:T5:TEMP:LOOP:RANGE:%f\n' % heater # mA
    d=self.oxford_send_cmd(cmdheat)
    return d

def oxford_read_heater_range(self):
    '''
    read the maximum current range for the heater
    NOTE: the return is in mA.
    '''
    cmd='READ:DEV:T5:TEMP:LOOP:RANGE\n'
    d=self.oxford_send_cmd(cmd)
    try:
        val=eval(d[-1].replace('mA','').replace('uA',''))
        if d[-1].find('mA')>0:
            f=1e-3
        elif d[-1].find('uA')>0:
            f=1e-6
        else:
            f=1
        Imax=val*f
    except:
        print('ERROR! could not read the heater power range maximum')
        return None

    # return value in mA rather than Amps
    return Imax*1000
    
def oxford_determine_best_heater_level(self):
    '''
    determine the correct heater maximum level based on the set point
    NOTE: this is returned in mA
    '''
    Tsetpt=self.oxford_read_set_point()
    if Tsetpt is None:return None

    if Tsetpt>=0.3:
        heater=1.0
    elif Tsetpt>=0.2:
        heater=0.316
    else:
        heater=0.1

    return heater

def oxford_increase_heater_range(self,testidx=-1):
    '''
    increase by one level the maximum heater range
    '''
    # first read the current setting
    heater=self.oxford_read_heater_range()
    if heater is None:
        if testidx<0:return None
        heater=self.oxford_heater_ranges[testidx]

    if heater>=self.oxford_heater_ranges[-1]:
        print('WARNING! Already at maximum heater range!')
        return heater
    
    self.debugmsg('current heater range: %f mA' % heater)
    # check which level number
    for heater_range in self.oxford_heater_ranges:
        self.debugmsg('next heater range: %f mA' % heater_range)
        if heater_range>heater:break

    d=self.oxford_set_heater_range(heater_range)
    return d
