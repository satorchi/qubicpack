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
import numpy as np
import datetime as dt
import sys,os,time

from qubicpack.qubicasic import qubicasic

class qubicfp:
    verbosity = 1 # class variable.  You can change this before instantiating an object
    __object_type__ = 'qubicfp'

    # maximum number of ASIC possible
    NASIC = 16 # this will change to 16 or 32

    from .tools import\
        debugmsg,\
        printmsg,\
        read_fits_field,\
        read_fits,\
        find_calsource,\
        find_hornswitch,\
        assign_bath_temperature,\
        read_qubicstudio_dataset,\
        read_calsource_fits,\
        read_qubicstudio_fits,\
        read_qubicstudio_hkextern_fits,\
        read_qubicstudio_hkfits,\
        output_filename,\
        data_subdir,\
        get_from_keyboard,\
        writelog,\
        RawMask,\
        get_hk,\
        gps,\
        pps,\
        assign_pointing_data,\
        azimuth_redmount,\
        elevation_redmount,\
        azimuth_cortex,\
        elevation_cortex,\
        azimuth,\
        elevation,\
        rotation,\
        hwp_position,\
        calsource,\
        qubicstudio_filetype_truename,\
        qubicstudio_hk_truename
    
    from .assign_variables import\
        assign_constants,\
        assign_fitsblurbs,\
        assign_observer,\
        assign_datadir,\
        assign_obsdate,\
        assign_detector_name,\
        guess_detector_name,\
        assign_logfile,\
        assign_temperature_labels

    from .timeline import\
        timeaxis,\
        pointing_timeaxis

    from .timestamping_diagnostic import\
        timestamp_diagnostic,\
        plot_timestamp_diagnostic,\
        plot_pps,\
        plot_pps_nsamples,\
        plot_timestamp_diagnostic_fig1,\
        plot_timestamp_diagnostic_fig2,\
        lost_packets

    from .fpmethods import\
        assign_verbosity,\
        infotext,\
        calsource_oldinfo,\
        calsource_info,\
        calsource_infotext,\
        assign_defaults,\
        assign_temperature,\
        read_qubicstudio_science_fits,\
        read_qubicstudio_raw_fits,\
        read_qubicstudio_asic_fits,\
        read_qubicpack_fits,\
        exist_data,\
        nasics,\
        args_ok,\
        asic,\
        bias_phase,\
        timeline_vbias,\
        Rfeedback,\
        FLL_State,\
        relay_heater,\
        sample_period,\
        timeline,\
        timeline_array,\
        tod,\
        plot_timeline,\
        plot_timeline_focalplane,\
        plot_iv,\
        plot_responsivity,\
        plot_iv_focalplane,\
        filter_iv_all,\
        ngood,\
        is_good,\
        make_iv_tex_report,\
        make_iv_report


    from .iv import\
        ADU2I

    from .quicklook import\
        assign_imagename,\
        plot_calsource,\
        plot_temperatures,\
        plot_300mKtemperatures,\
        plot_1Ktemperatures,\
        plot_switchstatus,\
        plot_azel,\
        plot_hwp,\
        quicklook

    from .demodulate import\
        demodulate

    from .mapmaking import\
        make_map_no_modulation,\
        make_all_maps,\
        plot_maps_focalplane
    
    from .hwp_analysis import\
        assign_hwp_chunkinfo,\
        chunkify_by_position

    from .ASD import\
        plot_powerspectrum,\
        plot_powerspectrum_focalplane
    

    from .qubiciv import\
        write_iv_fits

    from .level1.filetools import\
        write_level1
    
    def __init__(self):

        self.asic_list = []
        for idx in range(self.NASIC):
            self.asic_list.append(None)

        self.assign_defaults()
        return None


