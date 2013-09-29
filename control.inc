;;; control.inc --- Variables for flux_capacitor.pro
;; Author: Sebastian Luque
;; Created: 2013-08-27T17:52:26+0000
;; Last-Updated: 2013-09-28T11:50:04+0000
;;           By: Sebastian P. Luque
;; ------------------------------------------------------------------------
;;; Commentary: 
;;
;; This is a control file where variables for flux_capacitor.pro are
;; defined.
;;; -----------------------------------------------------------------------
;;; Code:

;;;_ : Flux period parameters

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

;;;_ : Key directories for processed flux files

;; ArcticNet 2011 processed flux root directory
proc_flux_dir_root=rootdir + 'TowerData' + path_sep() + 'ProcessedFlux'
;; Directory for raw flux files
raw_ec_dir=rootdir + 'TowerData' + path_sep() + 'FLUX'
;; Directory for standardized raw flux files
std_ec_dir=rootdir + 'TowerData' + path_sep() + 'FLUX' + path_sep() + $
           'std-FLUX'
;; Directory for extracted EC flux runs
ec_run_dir=rootdir + 'TowerData' + path_sep() + 'ProcessedFlux' + $
           path_sep() + 'Runs'
;; Directory for extracted GYRO runs
gyro_run_dir=rootdir + 'TowerData' + path_sep() + 'ProcessedFlux' + $
             path_sep() + 'GYRORuns'
;; Directory for extracted RMC runs
rmc_run_dir=rootdir + 'TowerData' + path_sep() + 'ProcessedFlux' + $
            path_sep() + 'RMCRuns'

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

;;;_ : Directories to store processed flux plots

lag_plot_dir=proc_flux_dir + 'LagPlots'
spec_plot_dir=proc_flux_dir + 'SpecPlots'
ogvie_plot_dir=proc_flux_dir + 'OgviePlots'
sonic_timeseries_dir=proc_flux_dir + 'SonicTimeSeries'
co2_timeseries_dir=proc_flux_dir + 'CO2TimeSeries'
nav_timeseries_dir=proc_flux_dir + 'NavTimeSeries'
foot_plot_dir=proc_flux_dir + 'footprint_plots'
lag_TU_plot_dir=proc_flux_dir + 'LagTU_plots'

;;;_ : MET directories

;; Raw met files
met_idir=rootdir + 'TowerData' + path_sep() + 'MET' + path_sep() + $
         'raw-MET'
;; Daily met files
met_daily_dir=rootdir + 'TowerData' + path_sep() + 'MET' + path_sep() + $
             'daily-MET'
;; Processed met files
met_proc_dir=rootdir + 'TowerData' + path_sep() + 'MET' + path_sep() + $
            'proc-MET'

;;;_ : NAV directories

;; Raw nav files
nav_raw_dir=rootdir + 'TowerData' + path_sep() + 'NAV' + path_sep() + $
            'raw-NAV'
;; Std nav files
nav_std_dir=rootdir + 'TowerData' + path_sep() + 'NAV' + path_sep() + $
            'std-NAV'
;; Daily nav files
nav_daily_dir=rootdir + 'TowerData' + path_sep() + 'NAV' + path_sep() + $
              'daily-NAV'
;; 1 min nav files
nav_min_dir=rootdir + 'TowerData' + path_sep() + 'NAV/1min-NAV'

;;;_ : GYRO directories and templates

;; GYRO files directly from ship (*.RAW files)
raw_gyro_dir=rootdir + 'GYRO' + path_sep() + 'raw-GYRO'
;; Pre-processed GYRO files
pre_gyro_dir=rootdir + 'GYRO' + path_sep() + 'pre-GYRO'
;; GYRO files that have standardized datestamp and data format
std_gyro_dir=rootdir + 'GYRO' + path_sep() + 'std-GYRO'
;; Time-corrected high frequency GYRO files
cor_gyro_dir=rootdir + 'GYRO' + path_sep() + 'cor-GYRO'
;; 1 min GYRO files
min_gyro_dir=rootdir + 'GYRO/1min-GYRO'
;; Daily GYRO files
daily_gyro_dir=rootdir + 'GYRO' + path_sep() + 'daily-GYRO'
;; Template for raw GYRO files
raw_gyro_temp=rootdir + 'Templates' + path_sep() + 'rawgyrotemplate.sav'
;; Template for standardized GYRO files
gyro_template=rootdir + 'Templates' + path_sep() + 'gyro_template.sav'
;; Template for daily GYRO files
daily_gyro_temp=rootdir + 'Templates' + path_sep() + 'daily_gyro_template.sav'
;; Directory for time conversion files
conv_dir=rootdir + 'GYRO' + path_sep() + 'time_conversion'
;; Template for pre-processed GYRO files
pre_GYRO_tem=rootdir + 'Templates' + path_sep() + 'pre_GYRO_tem.sav'

;;;_ : RMC directories and templates

;; RMC files directly from ship (*.RAW files)
raw_rmc_dir=rootdir + 'RMC' + path_sep() + 'raw-RMC'
;; Pre-processed RMC files
pre_rmc_dir=rootdir + 'RMC' + path_sep() + 'pre-RMC'
;; RMC files that have standardized datestamp and data format
std_rmc_dir=rootdir + 'RMC' + path_sep() + 'std-RMC'
;; 1 min RMC files
min_rmc_dir=rootdir + 'RMC/1min-RMC'
;; daily RMC files
daily_rmc_dir=rootdir + 'RMC' + path_sep() + 'daily-RMC'
;; Template for raw RMC files
raw_rmc_temp=rootdir + 'Templates' + path_sep() + 'rawrmctemplate.sav'
;; Template for standardized RMC files
rmc_template=rootdir + 'Templates' + path_sep() + 'rmc_template.sav'
;; Template for daily RMC files
daily_rmc_temp=rootdir + 'Templates' + path_sep() + 'daily_rmc_template.sav'
;; Template for pre-processed RMC files
pre_RMC_tem=rootdir + 'Templates' + path_sep() + 'pre_RMC_template.sav'

;;;_ : Other template files

met_template=rootdir + 'Templates' + path_sep() + 'met_template.sav'
nav_std_template=rootdir + 'Templates' + path_sep() + 'nav_template.sav'
dailymet_temp=rootdir + 'Templates' + path_sep() + 'dailymet_temp.sav'
dailynav_temp=rootdir + 'Templates' + path_sep() + 'dailynav_temp.sav'
oneminnav_temp=rootdir + 'Templates' + path_sep() + 'oneminnav_temp.sav'
met_proctemp=rootdir + 'Templates' + path_sep() + 'met_proctemp.sav'
fixlagtemplate=rootdir + 'Templates' + path_sep() + 'fixlagtemplate.sav'
Rad_template=rootdir + 'Templates' + path_sep() + 'RAD_template.sav'
met_log_template=rootdir + 'Templates' + path_sep() + 'full_log_template.sav'

;;;_ : Log files and directories

met_log_dir=rootdir + 'TowerData' + path_sep() + 'Logs'
met_log_name='met_log.csv'
closed_path_log='closed_path_log.csv'
full_irga_log=rootdir + 'TowerData' + path_sep() + 'Logs' + path_sep() + $
              'IRGA_full_log.csv'
user_IRGA_log=rootdir + 'TowerData' + path_sep() + 'Logs' + path_sep() + $
              'closed_path_log.csv'

;;;_ : Output test files

run_file='AN11_fluxruns.txt'
run_file_full=rootdir + 'TowerData' + path_sep() + 'ProcessedFlux' + $
              path_sep() + 'AN11_fluxruns.txt'
fixlagfile=rootdir + 'TowerData' + path_sep() + 'ProcessedFlux' + $
           path_sep() + 'fixedlags.txt'
gill_digi_analysis_file=rootdir + 'TowerData' + path_sep() + $
                        'ProcessedFlux' + path_sep() + 'gill_digi_analysis.txt'
IRGA_digi_analysis_file=rootdir + 'TowerData' + path_sep() + $
                        'ProcessedFlux' + path_sep() + 'LI7000_digi_analysis.txt'

;;;_ : Directories for radiation data

;; Directory of raw radiation files
raw_rad_dir=rootdir + 'RadData' + path_sep() + 'raw-RAD'
;; Directory of daily radiation files
daily_rad_dir=rootdir + 'RadData' + path_sep() + 'daily-RAD'
;; processed radiation data (BURBA correction)
proc_rad_dir=rootdir + 'RadData' + path_sep() + 'proc-RAD'

;;;_ : IRGA directories

IRGA_dir=rootdir + 'TowerData' + path_sep() + 'IRGA'

;;_ : Set file parameters and text for file headers

;; header     - text string that defines the header
;; timing     - sample rate of collected data (e.g. 1 or 60 sec)
;; outcolumns - number of columns in the output file
;; stamp      - identification stamp for type of data
;; prefix     - number of characters BEFORE the date of the input file name
;; bear_field - field number which has bearing (compass) data (set to 0 for no data)
;; vec_field  - field number which has vector magnitude data associated
;;              with the bearing field (set to 0 for no data)

met_std_header='year, month, day, hour, minute, second, ' + $
               'prog_version, battery_volt, panel_temperature, ' + $
               'pressure, air_temperature, RH_percent, ' + $
               'surface_temperature, wind_speed, wind_direction, ' + $
               'wind_sd, PAR, pitch, roll'
met_header='year, month, day, hour, minute, second, prog_version, ' + $
           'battery_voltage, panel_temperature, pressure, ' + $
           'air_temperature, RH_percent, surface_temperature, ' + $
           'wind_speed, wind_direction, wind_stdev, ' + $
           'PAR, pitch, roll, battery_stddev, ' + $
           'panel_temperature_stddev, pressure_stdev, ' + $
           'air_temperature_stdev, RH_stdev, ' + $
           'surface_temperature_stdev, PAR_stdev'
met_timing=60
met_outcolumns=26
met_stamp='MET'
met_prefix=8

nav_std_header='year, month, day, hour, minute, second, prog_version, ' + $
               'GPSyear, GPSmonth, GPSday, GPShour, GPSmin, GPSsec, ' + $
               'latitude, longitude, SOG, COG, heading, pitch, roll, ' + $
               'accel_x, accel_y, accel_z'
nav_timing=1
nav_outcolumns=23
nav_stamp='NAV'
nav_prefix=8

rad_header='year, month, day, hour, minute, second, battery_voltage, ' + $
           'panel_temperature, battery_stdev, ' + $
           'panel_temperature_stdev, Kdown, thermopile, Tcase, ' + $
           'Tdome_avg, LWin, PAR, Tuv, UV_B, UV_A, UV_broad, ' + $
           'Kdown_stdev, thermopile_stdev, Tcase_stdev, Tdome_stdev, ' + $
           'LWin_stdev, PAR_stdev, Tuv_stdev, UV_B_stdev, UV_A_stdev, ' + $
           'UV_broad_stdev' 
rad_timing=60
rad_outcolumns=30
rad_stamp='RAD'
rad_prefix=8

gyro_header='year, month, day, hour, min, sec, heading'
gyro_timing=1
gyro_outcolumns=7
gyro_stamp='GYRO'
gyro_prefix=7
gyro_bear_field=6
gyro_vec_field=0

rmc_header='year, month, day, hour, min, sec, latitude, longitude, SOG, COG'
rmc_timing=1
rmc_outcolumns=10
rmc_stamp='RMC'
rmc_prefix=6
rmc_bear_field=9
rmc_vec_field=8


;;; Emacs local variables
;; Local Variables:
;; mode: idlwave
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; control.inc ends here
