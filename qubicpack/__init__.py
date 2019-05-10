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
from __future__ import division, print_function
import numpy as np
import datetime as dt
import sys,os,time

try:
    import pystudio
    with_pystudio=True
except:
    with_pystudio=False

    
class qubicpack:
    from .assign_variables import\
        assign_defaults,\
        assign_observer,\
        assign_asic,\
        asic_index,\
        TES_index,\
        assign_integration_time,\
        assign_ADU,\
        assign_pausetime,\
        assign_temperature,\
        assign_datadir,\
        assign_obsdate,\
        assign_bias_factor,\
        assign_detector_name,\
        guess_detector_name,\
        assign_logfile

    from .pix2tes import\
        assign_pix_grid,\
        assign_pix2tes,\
        tes2pix,\
        pix2tes,\
        assign_lookup_table,\
        lookup_TEStable
    
    from .tools import\
        read_date_from_filename,\
        keyvals,\
        write_fits,\
        read_fits_field,\
        read_fits,\
        read_qubicpack_fits,\
        read_qubicstudio_dataset,\
        read_calsource_fits,\
        read_qubicstudio_fits,\
        read_qubicstudio_science_fits,\
        read_qubicstudio_asic_fits,\
        read_qubicstudio_hkextern_fits,\
        read_qubicstudio_hkfits,\
        read_bins,\
        output_filename,\
        data_subdir,\
        get_from_keyboard,\
        writelog,\
        pps2date,\
        RawMask,\
        gps,\
        pps,\
        azimuth,\
        elevation,\
        bias_phase,\
        calsource,\
        infotext,\
        qubicstudio_filetype_truename,\
        azel_etc

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
        plot_responsivity,\
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

    from .oxford import\
        oxford_assign_temperature_labels,\
        oxford_assign_heater_ranges,\
        oxford_send_cmd,\
        oxford_init,\
        oxford_pidoff,\
        oxford_set_point,\
        oxford_read_set_point,\
        oxford_read_temperature,\
        oxford_read_bath_temperature,\
        oxford_read_all_temperatures,\
        oxford_check_calibration,\
        oxford_read_heater_level,\
        oxford_read_heater_range,\
        oxford_set_heater_range,\
        oxford_determine_best_heater_level,\
        oxford_increase_heater_range
    
    from .timestamping_diagnostic import\
        plot_timestamp_diagnostic,\
        lost_packets
    
    if with_pystudio:
        '''
        these methods connect to QubicStudio and require pystudio
        '''
        from .acquisition import\
            connect_QubicStudio,\
            verify_QS_connection,\
            configure_PID,\
            compute_offsets,\
            feedback_offsets,\
            get_amplitude,\
            get_mean,\
            integrate_scientific_data,\
            get_NPIXELS,\
            get_nsamples,\
            get_chunksize,\
            get_RawMask,\
            get_bias,\
            get_PID,\
            set_Rfeedback,\
            set_VoffsetTES,\
            get_iv_data,\
            get_iv_timeline,\
            get_ASD

        from .squids import\
            squid_test

    else:
        def connect_QubicStudio(self,silent=True):
            return None

    def __init__(self):
        self.assign_defaults()
        return

    def debugmsg(self,msg,level=1):
        if self.debuglevel>=level:
            if self.logfile is None:
                print('DEBUG %s : %s' % (dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC'),msg))
            else:
                self.writelog('DEBUG: %s' % msg)
        return

    def printmsg(self,msg,verbosity=1):
        '''
        print a message to screen
        '''
        if self.verbosity>=verbosity:
            print(msg)
        return
    
    
    
