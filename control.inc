;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-08-27T17:52:26+0000
;; Last-Updated: 2013-12-05T21:10:24+0000
;;           By: Sebastian Luque
;; ------------------------------------------------------------------------
;;; Commentary: 
;;
;; This is a control file where variables for flux_capacitor.pro are
;; defined.
;;; -----------------------------------------------------------------------
;;; Code:

;; Temporarily setting these here; they will be set from main GUI program.
root_dir=expand_path('~/Data/ArcticNet/2011')
templates_dir=expand_path('~/Scripts/Projects/Flux_Capacitor/Templates')

;;;_ : RMC

rmc_ship_stamp='RMC_SHIP'
rmc_omg_stamp='RMC_OMG'
rmc_daily_rate=1
rmc_ship_raw_keep_fields=['latitude', 'longitude', 'sog', 'cog']
rmc_omg_raw_keep_fields=rmc_ship_raw_keep_fields

;;;_  . Ship

rmc_ship_cog_fld=15                  ; COG field in STD files
rmc_ship_sog_fld=14                  ; SOG field in STD files
;; Raw files (*.RAW files)
rmc_ship_idir=root_dir + path_sep() + 'RMC' + path_sep() + 'Ship'
;; STD files
rmc_ship_std_dir=rmc_ship_idir + path_sep() + 'STD'
;; Daily files
rmc_ship_daily_dir=rmc_ship_idir + path_sep() + 'Daily'
;; Averaged files
rmc_ship_avg_dir=rmc_ship_idir + path_sep() + '1min'
;; Period files
rmc_ship_period_dir=rmc_ship_daily_dir + path_sep() + 'Periods'

;;;_  . OMG

rmc_omg_cog_fld=9                  ; COG field in STD files
rmc_omg_sog_fld=8                  ; SOG field in STD files
;; Raw files (*.RAW files)
rmc_omg_idir=root_dir + path_sep() + 'RMC' + path_sep() + 'OMG'
;; STD files
rmc_omg_std_dir=rmc_omg_idir + path_sep() + 'STD'
;; Daily files
rmc_omg_daily_dir=rmc_omg_idir + path_sep() + 'Daily'
;; Averaged files
rmc_omg_avg_dir=rmc_omg_idir + path_sep() + '1min'
;; Period files
rmc_omg_period_dir=rmc_omg_daily_dir + path_sep() + 'Periods'

;;;_  . Templates

rmc_ship_raw_template=templates_dir + path_sep() + 'rmc_ship_raw.sav'
rmc_omg_raw_template=templates_dir + path_sep() + 'rmc_omg_raw.sav'
rmc_ship_std_template=templates_dir + path_sep() + 'rmc_ship_std.sav'
rmc_omg_std_template=templates_dir + path_sep() + 'rmc_omg_std.sav'
rmc_ship_avg_template=templates_dir + path_sep() + 'rmc_ship_avg.sav'
rmc_omg_avg_template=templates_dir + path_sep() + 'rmc_omg_avg.sav'


;;;_ : GYRO

gyro_ship_stamp='GYRO'
gyro_omg_stamp='GYRO_OMG'
gyro_hdg_fld=6                  ; heading field index in STD files
gyro_daily_rate=0.1
gyro_ship_raw_keep_fields=['heading']
gyro_omg_raw_keep_fields=gyro_ship_raw_keep_fields

;;;_  . Ship

;; Raw files
gyro_ship_idir=root_dir + path_sep() + 'GYRO' + path_sep() + 'Ship'
;; STD files
gyro_ship_std_dir=gyro_ship_idir + path_sep() + 'STD'
;; Filtered (10Hz) files
gyro_ship_std10Hz_dir=gyro_ship_idir + path_sep() + '10Hz'
;; Daily files
gyro_ship_daily_dir=gyro_ship_idir + path_sep() + 'Daily'
;; Averaged files
gyro_ship_avg_dir=gyro_ship_idir + path_sep() + '1min'
;; Period files
gyro_ship_period_dir=gyro_ship_daily_dir + path_sep() + 'Periods'

;;;_  . OMG

;; Raw files
gyro_omg_idir=root_dir + path_sep() + 'GYRO' + path_sep() + 'OMG'
;; STD files
gyro_omg_std_dir=gyro_omg_idir + path_sep() + 'STD'
;; Filtered (10Hz) files
gyro_omg_std10Hz_dir=gyro_omg_idir + path_sep() + '10Hz'
;; Daily files
gyro_omg_daily_dir=gyro_omg_idir + path_sep() + 'Daily'
;; Averaged files
gyro_omg_avg_dir=gyro_omg_idir + path_sep() + '1min'
;; Period files
gyro_omg_period_dir=gyro_omg_daily_dir + path_sep() + 'Periods'

;;;_  . Templates

gyro_ship_raw_template=templates_dir + path_sep() + 'gyro_ship_raw.sav'
gyro_omg_raw_template=templates_dir + path_sep() + 'gyro_omg_raw.sav'
gyro_ship_std_template=templates_dir + path_sep() + 'gyro_ship_std.sav'
gyro_omg_std_template=templates_dir + path_sep() + 'gyro_omg_std.sav'
;; gyro_ship_avg_template=templates_dir + path_sep() + 'gyro_avg_template.sav'
gyro_omg_avg_template=templates_dir + path_sep() + 'gyro_omg_avg.sav'


;;;_ : MET

met_stamp='MET'
met_daily_rate=60
met_raw_keep_fields=['prog_version', 'battery_voltage', 'panel_temperature', $
                     'pressure', 'air_temperature', 'rh_percent', $
                     'surface_temperature', 'wind_speed', $
                     'wind_direction', 'wind_direction_sd', 'par', $
                     'pitch', 'roll']
;; See PROCESS_MET for the meaning of fields below
;; Fields to pull from RMC STD files when calling PROCESS_MET (see docs)
rmc_pull_flds=[6, 7, 8, 9, 10, 11]
;; Fields to pull from Gyro STD files when calling PROCESS_MET (see docs)
gyro_pull_flds=gyro_hdg_fld
;; Location of angle fields in processed MET files
met_proc_angle_flds=[14, 17, 18, 21, 25, 26]
;; Location of magnitude fields in processed MET files
met_proc_mag_flds=[13, -1, -1, 22, -1, 27]
met_proc_diag_fld=33

;; Raw files
met_idir=root_dir + path_sep() + 'MET'
;; STD files
met_std_dir=met_idir + path_sep() + 'STD'
;; Daily files
met_daily_dir=met_idir + path_sep() + 'Daily'
;; Processed files
met_proc_dir=met_daily_dir + path_sep() + 'Processed'
;; Processed and averaged files
met_proc_avg_dir=met_proc_dir + path_sep() + '20min'

;;;_  . Templates
met_raw_template=templates_dir + path_sep() + 'met_raw.sav'
met_std_template=templates_dir + path_sep() + 'met_std.sav'
met_proc_template=templates_dir + path_sep() + 'met_processed.sav'
met_proc_avg_pre_template=templates_dir + path_sep() + $
                          'met_processed_avg_pre.sav'
met_proc_avg_template=templates_dir + path_sep() + 'met_processed_avg.sav'
met_proc_avg_diag_template=templates_dir + path_sep() + $
                           'met_processed_avg_diag.sav'

;;;_ : NAV

nav_stamp='NAV'
nav_daily_rate=1
;; These are the fields from raw NAV files that should be kept
nav_raw_keep_fields=['prog_version', 'latitude', 'longitude', 'sog', 'cog', $
                     'heading', 'pitch', 'roll', 'accel_x', 'accel_y', $
                     'accel_z']
;; nav_std_header='server_year, server_month, server_day, server_hour, ' + $
;;                'server_minute, server_second, prog_version, ' + $
;;                'utc_year, utc_month, utc_day, utc_hour, utc_min, utc_sec, ' + $
;;                'latitude, longitude, SOG, COG, heading, pitch, roll, ' + $
;;                'accel_x, accel_y, accel_z'
;; Corresponding field type codes
nav_raw_keep_types=[7, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
;; Angle and magnitude fields (args angle_fields and magnitude_fields in
;; AVERAGE_SERIES; see docs)
nav_angle_flds=[16, 17, 18, 19]
nav_mag_flds=[15, -1, -1, -1]

;; Raw files
nav_idir=root_dir + path_sep() + 'NAV'
;; STD files
nav_std_dir=nav_idir + path_sep() + 'STD'
;; Daily files
nav_daily_dir=nav_idir + path_sep() + 'Daily'
;; Averaged files
nav_avg_dir=nav_idir + path_sep() + '1min'

;;;_  . Templates

nav_raw_template=templates_dir + path_sep() + 'nav_raw.sav'
nav_std_template=templates_dir + path_sep() + 'nav_std.sav'


;;;_ : Radiation (RAD)

rad_raw_keep_fields=['battery_voltage', 'panel_temperature', $
                     'battery_voltage_sd', 'panel_temperature_sd', $
                     'k_down', 'thermopile', 'Tcase', 'Tdome', 'lw_in', $
                     'par', 'temperature_uv', 'uv_b', 'uv_a', 'uv_broad', $
                     'k_down_sd', 'thermopile_sd', 'Tcase_sd', $
                     'Tdome_sd', 'lw_in_sd', 'par_sd', 'temperature_uv_sd', $
                     'uv_b_sd', 'uv_a_sd', 'uv_broad_sd']
rad_stamp='RAD'
rad_daily_rate=60
;; See PROCESS_RAD docs for the two below.
rmc_latlon_flds=[6, 7]          ; lat and lon field location in RMC STD files
met_par_fld=16                  ; location of PAR field in MET STD files

;; Raw files
rad_idir=root_dir + path_sep() + 'RAD'
;; STD files
rad_std_dir=rad_idir + path_sep() + 'STD'
;; Daily files
rad_daily_dir=rad_idir + path_sep() + 'Daily'
;; Processed files (BURBA correction)
rad_proc_dir=rad_daily_dir + path_sep() + 'Processed'
;; Period files
rad_period_dir=rad_proc_dir + path_sep() + 'Periods'

;;;_  . Templates

rad_raw_template=templates_dir + path_sep() + 'rad_raw.sav'
rad_std_template=templates_dir + path_sep() + 'rad_std.sav'
rad_proc_template=templates_dir + path_sep() + 'rad_processed.sav'

;;;_ : Logs

logs_idir=root_dir + path_sep() + 'Logs'
met_log_file=logs_idir + path_sep() + 'met_log.csv'
closed_path_log_file=logs_idir + path_sep() + 'closed_path_log.csv'
cl_log_status_fld=10

;;;_  . Templates

met_log_template=templates_dir + path_sep() + 'met_log.sav'
closed_path_log_template=templates_dir + path_sep() + 'closed_path_log.sav'


;;;_ : Flux (eddy-covariance)

ec_daily_rate=0.1               ; sampling rate (s)
ec_stamp='EC'
ec_period=1200                  ; duration of flux period (s)
;; Fields in MET data (averaged over the flux period) used for filtering.
;; This is the filter_idx argument for FILTER_MET (see docs): mean wind
;; direction, mean SOG, mean heading, and mean true wind speed.
met_filter_flds=[15, 22, 19, 23]
;; And the corresponding thresholds
met_filter_thrs=[45.0, 5.0, 45.0, !VALUES.F_NAN]

ec_41_42_raw_keep_fields=['prog_version', 'accel_x', 'accel_y', 'accel_z', $
                          'rate_x', 'rate_y', 'rate_z', 'u_wind_analogue', $
                          'v_wind_analogue', 'w_wind_analogue', $
                          'sonic_temperature_analogue', 'CO2_op', $
                          'H2O_op', 'pressure_op', 'diag_op', 'CO2_cl', $
                          'H2O_cl', 'pressure_cl', 'temperature_cl', $
                          'diag_cl', 'u_wind_serial', 'v_wind_serial', $
                          'w_wind_serial', 'sonic_temperature_serial', $
                          'SoS_serial', 'gillstat_serial', $
                          'pressure_tower', 'temperature_base', $
                          'temperature_spar', 'temperature_bulb', $
                          'LI7000_mode', 'CO2_opsh', 'H2O_opsh', $
                          'pressure_opsh', 'diag_opsh', 'CO2_dry', $
                          'H2O_dry', 'CO2_LI7200', 'H2O_LI7200', $
                          'AGC_LI7200', 'pressure_LI7200', $
                          'temperature_in_LI7200', 'temperature_out_LI7200', $
                          'temperature_avg_LI7200', 'diag_LI7200']
ec_43_raw_keep_fields=['prog_version', 'accel_x', 'accel_y', 'accel_z', $
                       'rate_x', 'rate_y', 'rate_z', 'u_wind_analogue', $
                       'v_wind_analogue', 'w_wind_analogue', $
                       'sonic_temperature_analogue', 'CO2_op', 'H2O_op', $
                       'pressure_op', 'temperature_op', 'CO2_cl', 'H2O_cl', $
                       'pressure_cl', 'temperature_cl', 'diag_cl', $
                       'u_wind_serial', 'v_wind_serial', 'w_wind_serial', $
                       'sonic_temperature_serial', 'SoS_serial', $
                       'gillstat_serial', 'pressure_tower', 'LI7000_mode', $
                       'CO2_opsh', 'H2O_opsh', 'pressure_opsh', 'diag_opsh', $
                       'CO2_dry', 'H2O_dry', 'CO2_LI7200', 'H2O_LI7200', $
                       'AGC_LI7200', 'pressure_LI7200', $
                       'temperature_in_LI7200', 'temperature_out_LI7200', $
                       'temperature_avg_LI7200', 'diag_LI7200']
ec_44_45_raw_keep_fields=['prog_version', 'accel_x', 'accel_y', 'accel_z', $
                          'rate_x', 'rate_y', 'rate_z', 'u_wind_analogue', $
                          'v_wind_analogue', 'w_wind_analogue', $
                          'sonic_temperature_analogue', 'CO2_op', 'H2O_op', $
                          'pressure_op', 'temperature_op', 'diag_op', $
                          'CO2_cl', 'H2O_cl', 'pressure_cl', $
                          'temperature_cl', 'diag_cl', 'u_wind_serial', $
                          'v_wind_serial', 'w_wind_serial', $
                          'sonic_temperature_serial', 'SoS_serial', $
                          'gillstat_serial', 'pressure_tower', $
                          'LI7000_mode', 'CO2_opsh', 'H2O_opsh', 'pressure_opsh', $
                          'diag_opsh', 'CO2_dry', 'H2O_dry', 'CO2_LI7200', $
                          'H2O_LI7200', 'AGC_LI7200', 'pressure_LI7200', $
                          'temperature_in_LI7200', 'temperature_out_LI7200', $
                          'temperature_avg_LI7200', 'diag_LI7200']
ec_46_47_raw_keep_fields=['prog_version', 'accel_x', 'accel_y', 'accel_z', $
                          'rate_x', 'rate_y', 'rate_z', 'u_wind_analogue', $
                          'v_wind_analogue', 'w_wind_analogue', $
                          'sonic_temperature_analogue', 'CO2_op', 'H2O_op', $
                          'pressure_op', 'temperature_op', 'diag_op', $
                          'CO2_cl', 'H2O_cl', 'pressure_cl', 'temperature_cl', $
                          'diag_cl', 'u_wind_serial', 'v_wind_serial', $
                          'w_wind_serial', 'sonic_temperature_serial', $
                          'SoS_serial', 'gillstat_serial', 'pressure_tower', $
                          'LI7000_mode', 'CO2_opsh', 'H2O_opsh', 'pressure_opsh', $
                          'temperature_opsh', 'diag_opsh', 'CO2_dry', 'H2O_dry', $
                          'CO2_LI7200', 'H2O_LI7200', 'AGC_LI7200', $
                          'pressure_LI7200', 'temperature_in_LI7200', $
                          'temperature_out_LI7200', 'temperature_avg_LI7200', $
                          'diag_LI7200']
ec_std_fields=['prog_version', 'accel_x', 'accel_y', 'accel_z', $
               'rate_x', 'rate_y', 'rate_z', 'u_wind_analogue', $
               'v_wind_analogue', 'w_wind_analogue', $
               'sonic_temperature_analogue', 'CO2_op', 'H2O_op', $
               'pressure_op', 'temperature_op', 'diag_op', 'CO2_cl', $
               'H2O_cl', 'pressure_cl', 'temperature_cl', 'diag_cl', $
               'u_wind_serial', 'v_wind_serial', 'w_wind_serial', $
               'sonic_temperature_serial', 'SoS_serial', $
               'gillstat_serial', 'pressure_tower', 'LI7000_mode', $
               'CO2_opsh', 'H2O_opsh', 'pressure_opsh', $
               'temperature_opsh', 'diag_opsh', 'CO2_dry', 'H2O_dry', $
               'CO2_LI7200', 'H2O_LI7200', 'AGC_LI7200', $
               'pressure_LI7200', 'temperature_in_LI7200', $
               'temperature_out_LI7200', 'temperature_avg_LI7200', $
               'diag_LI7200', 'temperature_base', 'temperature_spar', $
               'temperature_bulb']
ec_std_types=[7, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 7, 4, 4, 4, 4, $
              7, 4, 4, 4, 4, 4, 4, 4, 7, 4, 4, 4, 4, 7, 4, 4, 4, 4, 4, 4, $
              4, 4, 4, 7, 4, 4, 4]

;; Raw files
ec_idir=root_dir + path_sep() + 'Flux'
ec_idir1=ec_idir + path_sep() + 'V4.1_4.2'
ec_idir2=ec_idir + path_sep() + 'V4.3'
ec_idir3=ec_idir + path_sep() + 'V4.4_4.5'
ec_idir4=ec_idir + path_sep() + 'V4.6_4.7'
;; STD files
ec_std_dir=ec_idir + path_sep() + 'STD'
;; Daily files
ec_daily_dir=ec_idir + path_sep() + 'Daily'
;; Period files
ec_period_dir=ec_daily_dir + path_sep() + 'Periods'
;; Output directory for motion corrected files
ec_motcorr_dir=ec_period_dir + path_sep() + 'Motion_Corrected'
;; Output directory for footprint plots
ec_footprint_dir=ec_period_dir + path_sep() + 'Footprint'
;; Output file for fluxes
ec_fluxes_file=ec_idir + path_sep() + 'fluxes.csv'

;; See FLUX docs for meaning of variables below
motpak_offset=[1.700, 0, 2.725]
sog_thr=0.25
lfreq_thr=0.02
hfreq_thr=1.0
xover_freq_thr=0.02

;;;_  . Templates

ec_raw_template1=templates_dir + path_sep() + 'ec_41_42.sav'
ec_raw_template2=templates_dir + path_sep() + 'ec_43.sav'
ec_raw_template3=templates_dir + path_sep() + 'ec_44_45.sav'
ec_raw_template4=templates_dir + path_sep() + 'ec_46_47.sav'
ec_std_template=templates_dir + path_sep() + 'ec_std.sav'



;;;_ * Emacs local variables
;; Local Variables:
;; mode: idlwave
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; control.inc ends here
