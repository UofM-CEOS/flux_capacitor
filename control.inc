;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-08-27T17:52:26+0000
;; Last-Updated: 2013-11-25T22:14:28+0000
;;           By: Sebastian P. Luque
;; ------------------------------------------------------------------------
;;; Commentary: 
;;
;; This is a control file where variables for flux_capacitor.pro are
;; defined.
;;; -----------------------------------------------------------------------
;;; Code:

;;;_ : Flux (eddy-covariance)

;; Time (in minutes) for flux run periods
fluxperiod=20
;; Allowable range of degrees (+/-) of the mean wind direction within the
;; flux period
t_wind=45
;; Allowable range of kts (+/-) of the mean ship velocity within the flux
;; period
t_SOG=5.0
;; Allowable range of degrees (+/-) of the mean heading direction within
;; the flux period
t_head=45.0

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
ec_daily_rate=0.1
ec_stamp='EC'
ec_period=1200

;; Raw files
ec_idir=rootdir + path_sep() + 'Flux'
ec_dir1=ec_idir + path_sep() + 'V4.1_4.2'
ec_dir2=ec_idir + path_sep() + 'V4.3'
ec_dir3=ec_idir + path_sep() + 'V4.4_4.5'
ec_dir4=ec_idir + path_sep() + 'V4.6_4.7'
;; STD files
ec_std_dir=ec_idir + path_sep() + 'STD'
;; Daily files
ec_daily_dir=ec_idir + path_sep() + 'Daily'
;; Period files
ec_period_dir=ec_idir + path_sep() + 'Periods'

;;;_  . Templates

ec_raw_template1=templates_dir + path_sep() + 'ec_41_42_template.sav'
ec_raw_template2=templates_dir + path_sep() + 'ec_43_template.sav'
ec_raw_template3=templates_dir + path_sep() + 'ec_44_45_template.sav'
ec_raw_template4=templates_dir + path_sep() + 'ec_46_47_template.sav'
ec_std_template=templates_dir + path_sep() + 'ec_std_template.sav'

;;;_ : Directory for output file and plotted data based on ddscode

val=widget_info((*pState).irga, /button_set)
IF val EQ 0 THEN BEGIN
   ddscode='Both'
   proc_flux_dir=proc_flux_dir_root + 'FullDigital'
   outname='ArcticNet11_FullDigital.out'
ENDIF ELSE BEGIN
   ddscode='IRGA'
   proc_flux_dir=proc_flux_dir_root + 'FullDigital'
   outname='ArcticNet11_IRGADigital.out'
ENDELSE   

;;;_ : MET

met_stamp='MET'
met_daily_rate=60
met_raw_keep_fields=['prog_version', 'battery_voltage', 'panel_temperature', $
                     'pressure', 'air_temperature', 'rh_percent', $
                     'surface_temperature', 'wind_speed', $
                     'wind_direction', 'wind_direction_sd', 'par', $
                     'pitch', 'roll']

;; Raw files
met_idir=rootdir + path_sep() + 'MET'
;; STD files
met_std_dir=met_idir + path_sep() + 'STD'
;; Daily files
met_daily_dir=met_idir + path_sep() + 'Daily'
;; Processed files
met_proc_dir=met_idir + path_sep() + 'Processed'
;; Processed and averaged files
met_proc_avg_dir=met_proc_dir + path_sep() + '20min'

;;;_  . Templates
met_raw_template=templates_dir + path_sep() + 'met_raw_template.sav'
met_std_template=templates_dir + path_sep() + 'met_std_template.sav'
met_log_template=templates_dir + path_sep() + 'met_log_template.sav'
met_proc_template=templates_dir + path_sep() + 'met_processed_template.sav'
met_proc_avg_pre_template=templates_dir + path_sep() + $
                          'met_processed_avg_pre_template.sav'
met_proc_avg_template=templates_dir + path_sep() + $
                      'met_processed_avg_template.sav'
met_proc_avg_diag_template=templates_dir + path_sep() + $
                           'met_processed_avg_diag_template.sav'

;;;_ : NAV

nav_stamp='NAV'
nav_daily_rate=1
;; These are the fields from raw NAV files that should be kept
nav_raw_keep_fields=['prog_version', 'latitude', 'longitude', 'sog', 'cog', $
                     'heading', 'pitch', 'roll', 'accel_x', 'accel_y', $
                     'accel_z']
;; Corresponding field type codes
nav_raw_keep_types=[7, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]

;; Raw files
nav_idir=rootdir + path_sep() + 'NAV'
;; STD files
nav_std_dir=nav_idir + path_sep() + 'STD'
;; Daily files
nav_daily_dir=nav_idir + path_sep() + 'Daily'
;; Averaged files
nav_avg_dir=nav_idir + path_sep() + '1min'

;;;_  . Templates

nav_raw_template=templates_dir + path_sep() + 'nav_raw_template.sav'
nav_std_template=templates_dir + path_sep() + 'nav_std_template.sav'


;;;_ : GYRO

gyro_ship_stamp='GYRO'
gyro_omg_stamp='GYRO_OMG'
gyro_daily_rate=0.1
gyro_ship_raw_keep_fields=['heading']
gyro_omg_raw_keep_fields=gyro_ship_raw_keep_fields

;;;_  . Ship

;; Raw files
gyro_ship_idir=rootdir + path_sep() + 'Gyro'
;; STD files
gyro_ship_std_dir=rmc_ship_idir + path_sep() + 'STD'
;; Filtered (10Hz) files
gyro_ship_std_dir=rmc_ship_idir + path_sep() + '10Hz'
;; Daily files
gyro_ship_daily_dir=rmc_ship_idir + path_sep() + 'Daily'
;; Averaged files
gyro_ship_avg_dir=rmc_ship_idir + path_sep() + '1min'

;;;_  . OMG

;; Raw files
gyro_omg_idir=rootdir + path_sep() + 'OMG' + path_sep() + 'HDG'
;; STD files
gyro_omg_std_dir=rmc_omg_idir + path_sep() + 'STD'
;; Filtered (10Hz) files
gyro_omg_std_dir=rmc_omg_idir + path_sep() + '10Hz'
;; Daily files
gyro_omg_daily_dir=rmc_omg_idir + path_sep() + 'Daily'
;; Averaged files
gyro_omg_avg_dir=rmc_omg_idir + path_sep() + '1min'

;;;_   : Templates

gyro_ship_raw_template=templates_dir + path_sep() + 'gyro_raw_template.sav'
gyro_omg_raw_template=templates_dir + path_sep() + 'omg_hdg_raw_template.sav'
gyro_ship_std_template=templates_dir + path_sep() + 'gyro_std_template.sav'
gyro_omg_std_template=templates_dir + path_sep() + 'omg_hdg_std_template.sav'
;; gyro_ship_avg_template=templates_dir + path_sep() + 'gyro_avg_template.sav'
gyro_omg_avg_template=templates_dir + path_sep() + $
                      'omg_hdg_std_avg_template.sav'


;;;_ : RMC

rmc_ship_stamp='RMC_SHIP'
rmc_omg_stamp='RMC_OMG'
rmc_daily_rate=1
rmc_ship_raw_keep_fields=['latitude', 'longitude', 'sog', 'cog']
rmc_omg_raw_keep_fields=rmc_ship_raw_keep_fields

;;;_  . Ship

;; Raw files (*.RAW files)
rmc_ship_idir=rootdir + path_sep() + 'RMC'
;; STD files
rmc_ship_std_dir=rmc_ship_idir + path_sep() + 'STD'
;; Daily files
rmc_ship_daily_dir=rmc_ship_idir + path_sep() + 'Daily'
;; Averaged files
rmc_ship_avg_dir=rmc_ship_idir + path_sep() + '1min'

;;;_  . OMG

;; Raw files (*.RAW files)
rmc_omg_idir=rootdir + path_sep() + 'OMG' + path_sep() + 'NVG'
;; STD files
rmc_omg_std_dir=rmc_omg_idir + path_sep() + 'STD'
;; Daily files
rmc_omg_daily_dir=rmc_omg_idir + path_sep() + 'Daily'
;; Averaged files
rmc_omg_avg_dir=rmc_omg_idir + path_sep() + '1min'

;;;_  . Templates

rmc_ship_raw_template=templates_dir + path_sep() + 'rmc_raw_template.sav'
rmc_omg_raw_template=templates_dir + path_sep() + 'omg_nav_raw_template.sav'
rmc_ship_std_template=templates_dir + path_sep() + 'rmc_std_template.sav'
rmc_omg_std_template=templates_dir + path_sep() + 'omg_nav_std_template.sav'
rmc_ship_avg_template=templates_dir + path_sep() + 'rmc_avg_template.sav'
rmc_omg_avg_template=templates_dir + path_sep() + $
                     'omg_nav_std_avg_template.sav'


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

;; Raw files
rad_idir=rootdir + path_sep() + 'RAD'
;; STD files
rad_std_dir=rad_idir + path_sep() + 'STD'
;; Daily files
rad_daily_dir=rad_dir + path_sep() + 'Daily'
;; Processed files (BURBA correction)
rad_proc_dir=rad_dir + path_sep() + 'Processed'

;;;_  . Templates

rad_raw_template=templates_dir + 'rad_raw_template.sav'
rad_std_template=templates_dir + 'rad_std_template.sav'


;;;_ : Log files and directories

met_log_dir=rootdir + 'TowerData' + path_sep() + 'Logs'
met_log_name='met_log.csv'
closed_path_log='closed_path_log.csv'
full_irga_log=rootdir + 'TowerData' + path_sep() + 'Logs' + path_sep() + $
              'IRGA_full_log.csv'
user_IRGA_log=rootdir + 'TowerData' + path_sep() + 'Logs' + path_sep() + $
              'closed_path_log.csv'



;;;_ * Emacs local variables
;; Local Variables:
;; mode: idlwave
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; control.inc ends here
