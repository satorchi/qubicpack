'''
$Id: qubicasic.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Wed 15 May 2019 10:50:03 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

class for the data/methods associated with the QUBIC asic (1/8th of the focal plane)
'''
import numpy as np
import datetime as dt
import sys,os,time

    
class qubicasic:
    verbosity=0 # class variable.  You can change this before instantiating an object.
    __object_type__ = 'qubicasic'

    from .assign_variables import\
        assign_constants,\
        assign_defaults,\
        assign_fitsblurbs,\
        assign_observer,\
        assign_asic,\
        asic_index,\
        assign_integration_time,\
        assign_ADU,\
        assign_pausetime,\
        assign_temperature,\
        assign_datadir,\
        assign_obsdate,\
        assign_bias_factor,\
        assign_detector_name,\
        guess_detector_name,\
        assign_logfile,\
        assign_lookup_table,\
        lookup_TEStable


    from .tools import\
        debugmsg,\
        printmsg,\
        read_date_from_filename,\
        keyvals,\
        write_fits,\
        read_fits_field,\
        read_fits,\
        read_qubicpack_fits,\
        read_qubicstudio_science_fits,\
        read_qubicstudio_asic_fits,\
        read_qubicstudio_hkfits,\
        read_bins,\
        output_filename,\
        data_subdir,\
        get_from_keyboard,\
        pps2date,\
        RawMask,\
        FLL_State,\
        get_hk,\
        gps,\
        pps,\
        bias_phase,\
        infotext,\
        qubicstudio_filetype_truename,\
        qubicstudio_hk_truename

    from .iv import\
        exist_iv_data,\
        wait_a_bit,\
        lut,\
        n_masked,\
        ADU2I,\
        setup_plot_Vavg,\
        plot_Vavg,\
        plot_iv_all,\
        setup_plot_iv_multi,\
        plot_iv_multi,\
        plot_iv_physical_layout,\
        make_line,\
        filter_jumps,\
        polynomial_fit_parameters,\
        combined_fit_parameters,\
        do_linefit,\
        do_polyfit,\
        do_combinedfit,\
        model_iv_super,\
        model_iv_mixed,\
        model_iv_normal,\
        model_iv_combined,\
        fit_iv,\
        draw_tangent,\
        fitted_iv_curve,\
        draw_iv,\
        setup_plot_iv,\
        adjusted_iv,\
        Vbias2Vtes,\
        oplot_iv,\
        plot_iv,\
        responsivity_func,\
        conductance_func,\
        etf_func,\
        calculate_responsivity,\
        plot_fom,\
        plot_ip,\
        plot_pv,\
        plot_rp,\
        make_Vbias,\
        filter_iv,\
        filter_iv_all,\
        read_filter,\
        save_filter,\
        filterinfo,\
        assign_filterinfo,\
        is_good_iv,\
        good_index,\
        ngood,\
        turnover,\
        Iturnover,\
        offset,\
        R1adjust,\
        R1,\
        Rn,\
        selected_iv_curve,\
        best_iv_curve,\
        Vtes,\
        Ites,\
        Ptes,\
        Rn_ratio,\
        Pbias,\
        read_ADU_file,\
        iv_tex_table_entry,\
        make_iv_tex_report,\
        make_iv_report,\
        iv2txt

    from .timeline import\
        exist_timeline_data,\
        ntimelines,\
        timeline,\
        timeline_array,\
        amplitude2DAC,\
        bias_offset2DAC,\
        sample_period,\
        timeline_npts,\
        timeline_timeaxis,\
        timeline_computertime,\
        timeaxis,\
        determine_bias_modulation,\
        timeline2adu,\
        plot_timeline,\
        plot_timeline_physical_layout,\
        model_timeline,\
        fit_timeline

    from .ASD import\
        plot_ASD,\
        plot_ASD_all,\
        plot_ASD_physical_layout,\
        make_ASD_tex_report
    
    from .timestamping_diagnostic import\
        plot_timestamp_diagnostic,\
        lost_packets
    
    def __init__(self):
        self.assign_defaults()
        return

    #### some dummy defs to be compatible with the old qubicpack
    def oxford_assign_temperature_labels(self):
        '''
        we don't have the Oxford in this class, but I want to use the generic "assign_defaults"
        '''
        self.printmsg('DEBUG:  Call to dummy oxford_assign_temperature_labels.',verbosity=3)
        return None

    def oxford_assign_heater_ranges(self):
        '''
        we don't have the Oxford in this class, but I want to use the generic "assign_defaults"
        '''
        self.printmsg('DEBUG:  Call to dummy oxford_assign_heater_ranges.',verbosity=3)
        return None


    
