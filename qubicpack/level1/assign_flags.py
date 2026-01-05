'''
$Id: assign_flags.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Fri 02 Jan 2026 18:42:44 CET
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

assign the flags for the TOD
'''
import numpy as np
from .flags import interpolate_flags, set_flag
from .saturation import assign_saturation_flags


def assign_temperature_flags(self,timeaxis,flag_array):
    '''
    Assign the flags for temperature exceeding certain limits
    We assign flags for the 300mK stage and the 1K stage above the following limits (see flags.py):
    330mK, 340mK, 350mK and 1.1K, 1.2K, 1.3K
    
    see details on the wiki: https://qubic.in2p3.fr/wiki/TD/ProcessedDataFormat

    We interpolate the Housekeeping timeaxis to the science data timeaxis, and determine the flags
    '''

    # flags for bath temperature
    timeaxis_Tbath,Tbath = self.Tbath
    if Tbath is not None:
        Tbath_interp = np.interp(timeaxis,timeaxis_Tbath,Tbath)

        for T_limit in [0.350,0.340,0.320]:
            mask = Tbath_interp>T_limit
            flagname = 'bath temperature above %.0fmK' % (1000*T_limit)
            flag_array[mask] = set_flag(flagname,flag_array[mask])


    # flags for 1K Stage
    T_1K = self.get_hk('1K Stage')
    timeaxis_1K = self.timeaxis(datatype='1K Stage')
    if T_1K is not None:
        T_interp = np.interp(timeaxis,timeaxis_1K,T_1K)
        
        for T_limit in [1.3,1.2,1.1]:
            mask = T_interp>T_limit
            flagname = '1K temperature above %.1fK' % (T_limit)
            flag_array[mask] = set_flag(flagname,flag_array[mask])
    
    
    return flag_array
    

def assign_flags(self,indextype='TES',t_tod=None,flag_array=None):
    '''
    assign the flags
    '''

    if indextype.upper().find('QS')==0 or indextype.upper().find('QUBICSOFT')>=0:
        is_QSindex = True
    else:
        is_QSindex = False

    if t_tod is None:
        do_interpolate = False
    else:
        do_interpolate = True
        npts_interp = len(t_tod)
    
        
    # some flags must be assigned for each ASIC separately in the original time axis
    flag_list = []
    asic_ctr = 0
    for asicobj in self.asic_list:
        if asicobj is None: continue
        asic_ctr += 1
        timeline_array = asicobj.timeline_array()
        timeaxis = asicobj.timeaxis()
        
        # assign flags for saturation
        asic_flag_array = assign_saturation_flags(timeline_array)

        # assign flags for flux jumps
        # asic_flag_array = assign_fluxjump_flags(timeline_array,asic_flag_array)

        # assign flags for cosmic ray events
        # asic_flag_array = assign_cosmicray_flags(timeline_array,asic_flag_array)

        if do_interpolate:
            flag_list.append( interpolate_flags(timeaxis,t_tod,asic_flag_array) )
        else:
            flag_list.append(asic_flag_array)
            
    # assign flags for temperature limits


    return
