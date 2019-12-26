'''
$Id: qubicfp.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Mon 20 May 2019 19:15:41 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

class to read/plot data from the entire QUBIC focal plane
'''
from __future__ import division, print_function

import numpy as np
import datetime as dt
import sys,os,time

from qubicpack.qubicasic import qubicasic
from qubicpack.utilities import NASIC

class qubicfp:
    verbosity = 0 # class variable.  You can change this before instantiating an object
    __object_type__ = 'qubicfp'

    from .tools import\
        debugmsg,\
        printmsg,\
        read_fits_field,\
        read_fits,\
        find_calsource,\
        find_hornswitch,\
        read_qubicstudio_dataset,\
        read_calsource_fits,\
        read_qubicstudio_fits,\
        read_qubicstudio_hkextern_fits,\
        read_qubicstudio_hkfits,\
        read_bins,\
        output_filename,\
        data_subdir,\
        get_from_keyboard,\
        writelog,\
        pps2date,\
        RawMask,\
        get_hk,\
        gps,\
        pps,\
        azimuth,\
        elevation,\
        hwp_position,\
        calsource,\
        qubicstudio_filetype_truename,\
        qubicstudio_hk_truename,\
        azel_etc

    from .assign_variables import\
        assign_constants,\
        assign_fitsblurbs,\
        assign_observer,\
        assign_datadir,\
        assign_obsdate,\
        assign_detector_name,\
        guess_detector_name,\
        assign_logfile

    from .timeline import\
        timeaxis

    from .timestamping_diagnostic import\
        timestamp_diagnostic,\
        plot_timestamp_diagnostic,\
        plot_pps,\
        plot_pps_nsamples,\
        plot_timestamp_diagnostic_fig1,\
        plot_timestamp_diagnostic_fig2,\
        lost_packets

    from .fpmethods import\
        infotext,\
        calsource_info,\
        calsource_infotext,\
        assign_defaults,\
        assign_temperature,\
        read_qubicstudio_science_fits,\
        read_qubicstudio_asic_fits,\
        read_qubicpack_fits,\
        args_ok,\
        bias_phase,\
        timeline_vbias,\
        Rfeedback,\
        feedback_on,\
        sample_period,\
        timeline,\
        timeline_array,\
        plot_timeline,\
        plot_timeline_focalplane,\
        plot_iv,\
        plot_responsivity,\
        plot_iv_focalplane,\
        filter_iv_all,\
        ngood,\
        make_iv_tex_report,\
        make_iv_report


    from .iv import\
        ADU2I

    from .quicklook import\
        plot_calsource,\
        plot_temperatures,\
        plot_300mKtemperatures,\
        plot_1Ktemperatures,\
        plot_switchstatus,\
        plot_azel,\
        plot_hwp,\
        quicklook
    
    def __init__(self):

        self.asic_list = []
        for idx in range(NASIC):
            self.asic_list.append(None)

        self.assign_defaults()
        return None


