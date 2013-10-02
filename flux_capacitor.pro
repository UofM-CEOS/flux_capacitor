;; $Id$
;;; flux_capacitor.pro --- 
;; Author: Bruce Johnson, Ryan Smith, Sebastian Luque
;; Created: 2013-08-23T22:24:25+0000
;; Last-Updated: 2013-10-02T18:20:21+0000
;;           By: Sebastian Luque
;; ------------------------------------------------------------------------
;;; Commentary:
;; 
;; The purpose of this program is to simplify the PROCESS_AN11.pro program
;; written by Bruce and make it easier to modify elements of the EC
;; program.
;;
;; Things to clarify with Ryan:
;;
;; o Where are the header descriptions?
;;
;; o Radiation section doesn't match the flow chart.
;;
;; o Argument 'out_columns' for std_nav_AN11 seems useless since the input
;;   columns are hard-coded with a specific constant meaning, so the number
;;   of columns will always be 23.  It is also determined by the
;;   'nav_header' variable, so it's redundant.  I have demoted this
;;   argument to a variable inside the procedure, but are there any gotchas
;;   for this?
;;
;; o I have replaced the function 'jd_to_calendar' with a more general and
;;   simpler solution 'doy2calendar'.
;;
;; o Should an '/overwrite' keyword be added to procedures/functions to
;;   overwrite output files?
;;
;; ------------------------------------------------------------------------
;;; Code:

;;;_ : Main widget setting initial variables

;;;_  . Utility procedure for copying dates in main widget

PRO copydates, event

  widget_control, event.top, get_uvalue=pState
  widget_control, (*pState).gyro_sdate, get_value=gyro_sdate
  widget_control, (*pState).gyro_edate, get_value=gyro_edate
  gyro_sdate=strtrim(gyro_sdate, 2)
  gyro_edate=strtrim(gyro_edate, 2)
  widget_control, (*pState).rmc_sdate, set_value=gyro_sdate
  widget_control, (*pState).rmc_edate, set_value=gyro_edate
  widget_control, (*pState).met_sdate, set_value=gyro_sdate
  widget_control, (*pState).met_edate, set_value=gyro_edate
  widget_control, (*pState).nav_sdate, set_value=gyro_sdate
  widget_control, (*pState).nav_edate, set_value=gyro_edate
  widget_control, (*pState).rad_sdate, set_value=gyro_sdate
  widget_control, (*pState).rad_edate, set_value=gyro_edate
  widget_control, (*pState).flux_sdate, set_value=gyro_sdate
  widget_control, (*pState).flux_edate, set_value=gyro_edate
  widget_control, (*pState).calc_sdate, set_value=gyro_sdate
  widget_control, (*pState).calc_edate, set_value=gyro_edate

END

;;;_  . Utility procedure for choosing directories in main widget

PRO chooserootdirectory, event

  widget_control, event.top, get_uvalue=pState
  filedir=dialog_pickfile(/directory, get_path=rootdir_select, /must_exist)
  draw_msg, (*pState).viewer, 'Selected Dir: ' + rootdir_select
  widget_control, (*pState).rootdir, set_value=rootdir_select
  
END

;;;_  . Clean up

PRO flux_capacitor_cleanup, base

  widget_control, base, get_uvalue=pState ;free the pointers
  ptr_free, pState

END

PRO flux_capacitor_event, event

END

;;;_  . Main widget

PRO flux_capacitor

  ;; First, check to make sure all routines are being called by the program
  ;; correctly
  
  resolve_all, /continue_on_error, unresolved=unresolved
  
  ;; result=file_search('day_splitter_GPS.pro',/fully_qualify_path)
  ;; un_file=filepath(root_dir=workingdir,unresolved[0])
  ;; cd,current=workingdir
  
  ;; Define the base widget, from which all other widgets will be added to...

  base=widget_base(title='Flux Capacitor 2012', scr_xsize=700, scr_ysize=800, $
                   col=3, /align_left)
  col1=widget_base(base, col=1, /align_left)
  top=widget_base(col1, col=2, /align_left)
  topleft=widget_base(top, col=1)
  label=widget_label(topleft, value=' ')
  label=widget_label(topleft, value='Welcome to Flux Capacitor 2012', $
                     /align_left)
  label=widget_label(topleft, value=' ')
  label=widget_label(topleft, $
                     value='Select root directory for ArcticNet data:', $
                     /align_left)
  rootdirbase=widget_base(topleft, col=2, /align_left)
  rootdir=widget_text(rootdirbase, /editable, scr_xsize=275)
  ;; bm=filepath('open.bmp',root_dir='resources',subdir='bitmaps')
  select_rootdir=widget_button(rootdirbase, scr_xsize=25, value='...', $
                               event_PRO='chooserootdirectory')

  label=widget_label(topleft, value='Select digital data stream type:', $
                     /align_left)
  ddcodebase=widget_base(topleft, /exclusive, col=1, /align_left, frame=1, $
                         scr_xsize=300)
  irga=widget_button(ddcodebase, value='IRGA')
  both=widget_button(ddcodebase, value='BOTH')
  widget_control, irga, /set_button
      
  topright=widget_base(top, col=1, /align_center)    
  viewer=widget_draw(topright, scr_xsize=370, scr_ysize=200)
  widget_control, viewer, set_uvalue=0    

  label=widget_label(col1, value='Set flux calculation options:', $
                     /align_left)
  moreoptions=widget_base(col1, col=3, /align_left, frame=1)
  ynbase=widget_base(moreoptions, /nonexclusive, col=1, /align_left, $
                     scr_xsize=325)
  rmc_pp=widget_button(ynbase, value='Pre-process the RMC and GYRO files?', $
                       scr_ysize=35)
  rmc_std=widget_button(ynbase, $
                        value='Standardize the RMC files?' + $
                        '                    -->', $
                        scr_ysize=35)
  gyro_std=widget_button(ynbase, $
                         value='Standardize the GYRO files?' + $
                         '                   -->', $
                         scr_ysize=35)
  gyro_split=widget_button(ynbase, $
                           value='Run the day splitter on the GYRO files?' + $
                           '       -->', $
                           scr_ysize=35) ; ---> more options if YES
  rmc_split=widget_button(ynbase, $
                          value='Run the day splitter on the RMC files?' + $
                          '        -->', $
                          scr_ysize=35) ; ---> more options if YES
  nav_std=widget_button(ynbase, value='Standardize the NAV files?', $
                        scr_ysize=35)
  met_split=widget_button(ynbase, $
                          value='Run the day splitter on the MET files?' + $
                          '        -->', $
                          scr_ysize=34) ; ---> more options if YES
  nav_split=widget_button(ynbase, $
                          value='Run the day splitter on the NAV files?' + $
                          '        -->', $
                          scr_ysize=34) ; ---> more options if YES
  met_proc=widget_button(ynbase, $
                         value='Process the MET files?', scr_ysize=35)
  rad_split=widget_button(ynbase, $
                          value='Run the day splitter on the Radiation files?' + $
                          ' -->', $
                          scr_ysize=34) ; ---> more options if YES
  flux_std=widget_button(ynbase, value='Standardize the flux files?', $
                         scr_ysize=35)
  flux_ext=widget_button(ynbase, $
                         value='Find and extract flux runs?' + $
                         '                   -->', $
                         scr_ysize=34) ; ---> more options if YES
  flux_calc=widget_button(ynbase,value= 'Process the flux calculations?' + $
                          '                -->', $
                          scr_ysize=34)  
  enterdates=widget_base(moreoptions, col=1)
  blank=widget_label(enterdates, value=' ', scr_ysize=35)
  omg_nvg_base=widget_base(enterdates, /nonexclusive)
  omg_nvg_std=widget_button(omg_nvg_base, value='Use OMG files?')
  omg_hdg_base=widget_base(enterdates, /nonexclusive)
  omg_hdg_std=widget_button(omg_hdg_base, value='Use OMG files?')
  gyro_split_dates_base=widget_base(enterdates, col=2, scr_ysize=28)
  gyro_sdate=widget_text(gyro_split_dates_base, scr_xsize=150, $
                         value='20110725', /editable)
  gyro_edate=widget_text(gyro_split_dates_base, scr_xsize=150, $
                         value='20110725', /editable)
  rmc_split_dates_base=widget_base(enterdates, col=2, scr_ysize=28)
  rmc_sdate=widget_text(rmc_split_dates_base, scr_xsize=150, $
                        value='Start Date: (yyyymmdd)', /editable)
  rmc_edate=widget_text(rmc_split_dates_base, scr_xsize=150, $
                        value='End Date: (yyyymmdd)', /editable)
  blank=widget_label(enterdates,value=' ', scr_ysize=35)
  met_split_dates_base=widget_base(enterdates, col=2, scr_ysize=28)
  met_sdate=widget_text(met_split_dates_base, scr_xsize=150, $
                        value='Start Date: (yyyymmdd)', /editable)
  met_edate=widget_text(met_split_dates_base, scr_xsize=150, $
                        value='End Date: (yyyymmdd)', /editable)
  nav_split_dates_base=widget_base(enterdates,col=2, scr_ysize=28)
  nav_sdate=widget_text(nav_split_dates_base, scr_xsize=150, $
                        value='Start Date: (yyyymmdd)', /editable)
  nav_edate=widget_text(nav_split_dates_base, scr_xsize=150, $
                        value='End Date: (yyyymmdd)', /editable)           
  blank=widget_label(enterdates,value=' ',scr_ysize=35)
  rad_split_dates_base=widget_base(enterdates, col=2, scr_ysize=28)
  rad_sdate=widget_text(rad_split_dates_base, scr_xsize=150, $
                        value='Start Date: (yyyymmdd)', /editable)
  rad_edate=widget_text(rad_split_dates_base, scr_xsize=150, $
                        value='End Date: (yyyymmdd)', /editable)              
  blank=widget_label(enterdates,value=' ', scr_ysize=35)
  flux_split_dates_base=widget_base(enterdates, col=2, scr_ysize=28)
  flux_sdate=widget_text(flux_split_dates_base, scr_xsize=150, $
                         value='Start Date: (yyyymmdd)', /editable)
  flux_edate=widget_text(flux_split_dates_base, scr_xsize=150, $
                         value='End Date: (yyyymmdd)', /editable) 
  flux_calc_dates_base=widget_base(enterdates,col=2, scr_ysize=28)
  calc_sdate=widget_text(flux_calc_dates_base, scr_xsize=150, $
                         value='Start Date: (yyyymmdd)', /editable)
  calc_edate=widget_text(flux_calc_dates_base, scr_xsize=150, $
                         value='End Date: (yyyymmdd)', /editable)
               
  copydatesbase=widget_base(moreoptions, col=1)
  blank=widget_label(copydatesbase, value=' ', scr_ysize=35)
  blank=widget_label(copydatesbase, value=' ', scr_ysize=35)
  blank=widget_label(copydatesbase, value=' ', scr_ysize=35)
  ;; bm=filepath('shift_down.bmp',root_dir='resources',subdir='bitmaps')
  copydates=widget_button(copydatesbase, scr_ysize=35, value=' V ', $
                          tooltip='Copy these start/end values to ' + $
                          'all boxes below ...', event_PRO='copydates')
                
  label=widget_label(col1, value=' ')
  start=widget_button(col1, value='Process ...', $
                      event_PRO='process', scr_xsize=300)
   
  col2=widget_base(base, col=1, /align_left)  
  
  state={rootdir:rootdir, irga:irga, rmc_pp:rmc_pp, rmc_std:rmc_std, $
         omg_nvg_std:omg_nvg_std, gyro_std:gyro_std, gyro_split:gyro_split, $
         rmc_split:rmc_split, nav_std:nav_std, met_split:met_split, $
         nav_split:nav_split, met_proc:met_proc, rad_split:rad_split, $
         flux_std:flux_std, flux_ext:flux_ext, flux_calc:flux_calc, $
         gyro_sdate:gyro_sdate, gyro_edate:gyro_edate, rmc_sdate:rmc_sdate, $
         rmc_edate:rmc_edate, met_sdate:met_sdate, met_edate:met_edate, $
         nav_sdate:nav_sdate, nav_edate:nav_edate, rad_sdate:rad_sdate, $
         rad_edate:rad_edate, flux_sdate:flux_sdate, flux_edate:flux_edate, $
         calc_sdate:calc_sdate, calc_edate:calc_edate, $
         omg_hdg_std:omg_hdg_std, viewer:viewer}

  pState = ptr_new(state, /no_copy)
  widget_control, base, set_uvalue=pState, /no_copy
  widget_control, base, /realize

  widget_control, viewer, get_value=drawid
  wset, drawid
  device, decomposed=0
  device, retain=2

  xmanager, 'flux_capacitor', base, cleanup='flux_capacitor_cleanup'

END

;;;_ : Process all options

;;;_  . Drawing function used throughout the main procedure

PRO draw_msg, ptr_viewer, msg

  widget_control, ptr_viewer, get_uvalue=xval
  x=float(xval)
  IF x GT 0.99 THEN BEGIN
     x=0.1
     erase
  ENDIF ELSE BEGIN
     x=x + 0.1
  ENDELSE
  widget_control, ptr_viewer, set_uvalue=x
  xyouts, 0.05, 1.0 - x, msg, /normal

END

;;;_  . Check value of a given variable selected in main widget

FUNCTION check_ctrl_var, strc_var

  val=widget_info(strc_var, /button_set)

  RETURN, val

END

;;;_  . Main procedure for processing all options

;; Make sure to set proper values for variables in control.inc 

PRO process, event

  widget_control, event.top, get_uvalue=pState
  widget_control, (*pState).rootdir, get_value=rootdir
  rootdir=strtrim(rootdir, 2)
  rootdir=rootdir[0]
  widget_control, (*pState).gyro_sdate, get_value=GYRO_SDATE
  widget_control, (*pState).gyro_edate, get_value=GYRO_EDATE
  widget_control, (*pState).rmc_sdate, get_value=rmc_sdate
  widget_control, (*pState).rmc_edate, get_value=rmc_edate
  widget_control, (*pState).met_sdate, get_value=met_sdate
  widget_control, (*pState).met_edate, get_value=met_edate
  widget_control, (*pState).nav_sdate, get_value=nav_sdate
  widget_control, (*pState).nav_edate, get_value=nav_edate
  widget_control, (*pState).rad_sdate, get_value=rad_sdate
  widget_control, (*pState).rad_edate, get_value=rad_edate
  widget_control, (*pState).flux_sdate, get_value=flux_sdate
  widget_control, (*pState).flux_edate, get_value=flux_edate
  widget_control, (*pState).calc_sdate, get_value=calc_sdate
  widget_control, (*pState).calc_edate, get_value=calc_edate
  
  ;; cd, current=workingdir

  @control.inc

  ;; ===========================================================
  ;; Call process if NAV standardization selector was set to Yes

  IF check_ctrl_var((*pState).nav_std) EQ 1 THEN BEGIN
     draw_msg, (*pState).viewer, 'Standardizing NAV Files ...'
     ;; std_nav, nav_raw_dir, nav_std_dir, nav_std_header, /overwrite
  ENDIF ELSE BEGIN
     print, 'Skipping standardization of NAV files ...'
  ENDELSE

  ;; ================================================================
  ;; Call process if NAV day splitter selector was set to Yes

  IF check_ctrl_var((*pState).nav_split) EQ 1 THEN BEGIN
     draw_msg, (*pState).viewer, 'Splitting NAV files ...'
     ;; day_splitter, nav_SDATE[0], nav_EDATE[0], nav_std_dir, $
     ;;               nav_daily_dir, nav_std_template, 0, nav_timing, $
     ;;               nav_stamp    ;, /overwrite
     ;; draw_msg, (*pState).viewer, 'Creating 1-min averages of NAV files ...'
     ;; nav_avg, nav_daily_dir, nav_min_dir, nav_timing, 60, 15, 16, 17, $
     ;;          nav_stamp, nav_std_template, 0
  ENDIF ELSE BEGIN
     print, 'Skipping splitting of NAV files ...'
  ENDELSE

  ;; ================================================================
  ;; Call process if MET day splitter selector was set to Yes

  IF check_ctrl_var((*pState).met_split) EQ 1 THEN BEGIN
     draw_msg, (*pState).viewer, 'Splitting MET files ...'
     day_splitter, met_SDATE[0], met_EDATE[0], met_std_dir, met_daily_dir, $
                   met_std_template, 0, met_timing, met_stamp, /overwrite
  ENDIF ELSE BEGIN
     print, 'Skipping splitting of MET files ...'
  ENDELSE

  ;; ==========================================================
  ;; Call process if RMC standardization selector was set to Yes

  rmc_std=check_ctrl_var((*pState).rmc_std)
  IF rmc_std EQ 1 THEN BEGIN
     ;; Check whether we want to standardize the RMC or OMG files
     omg_nvg_std=check_ctrl_var((*pState).omg_nvg_std)
     IF omg_nvg_std EQ 1 THEN BEGIN
        draw_msg, (*pState).viewer, 'Creating standardized NVG file ...'
        ;; omg_navigation, rootdir
     ENDIF ELSE BEGIN
        draw_msg, (*pState).viewer, 'Creating standardized RMC Files ...'
        ;; stand_r, pre_rmc_dir, std_rmc_dir, pre_RMC_tem, conv_dir
        draw_msg, (*pState).viewer, '[VERIFY THAT PRE-PROCESSED GYRO ' + $
                  'FILES AND TIME CONVERSION FILES HAVE SAME DATE/TIME STAMP]'
     ENDELSE
  ENDIF ELSE BEGIN
     print, 'Skipping standardization of RMC files ...'
  ENDELSE

  ;; ==========================================================
  ;; Call process if GYRO standardization selector was set to Yes

  gyro_std=check_ctrl_var((*pState).gyro_std)
  IF gyro_std EQ 1 THEN BEGIN 
     ;; Check whether we want to standardize the RMC or OMG files
     omg_hdg_std=check_ctrl_var((*pState).omg_hdg_std)
     IF omg_hdg_std EQ 1 THEN BEGIN
        draw_msg, (*pState).viewer, 'Creating standardized HDG file ...'
        ;; o_head, rootdir
     ENDIF ELSE BEGIN
        draw_msg, (*pState).viewer, 'Creating standardized GYRO files ...'
        ;; standardize_gyro_v3, pre_gyro_dir, std_gyro_dir, pre_GYRO_tem, $
        ;;                      conv_dir, cor_gyro_dir
     ENDELSE
  ENDIF ELSE BEGIN
     print, 'Skipping standardization of GYRO files ...'
  ENDELSE

  ;; ===========================================================
  ;; Call processes if GYRO day splitter selector was set to Yes

  IF check_ctrl_var((*pState).gyro_split) EQ 1 THEN BEGIN
     draw_msg, (*pState).viewer, 'Splitting GYRO into daily files ...'
     ;; day_splitter_GPS, gyro_SDATE, gyro_EDATE, std_gyro_dir, $
     ;;                   daily_gyro_dir, gyro_template, gyro_header, gyro_timing, $
     ;;                   gyro_outcolumns, gyro_stamp, gyro_prefix, rootdir
     draw_msg, (*pState).viewer, 'Creating 1-min averages of GYRO files ...'
     ;; one_min_vAN11, daily_gyro_dir, min_gyro_dir, daily_gyro_temp, $
     ;;                gyro_outcolumns, gyro_bear_field, gyro_vec_field, $
     ;;                gyro_stamp, gyro_header, rootdir
  ENDIF ELSE BEGIN
     print, 'Skipping splitting of GYRO files ...'
  ENDELSE
  
  ;; ========================================================
  ;; Call process if RMC day splitter selector was set to Yes

  IF check_ctrl_var((*pState).rmc_split) EQ 1 THEN BEGIN
     draw_msg, (*pState).viewer, 'Splitting RMC into daily files ...'
     ;; day_splitter_GPS, rmc_SDATE, rmc_EDATE, std_rmc_dir, daily_rmc_dir, $
     ;;                   rmc_template, rmc_header, rmc_timing, $
     ;;                   rmc_outcolumns, rmc_stamp, rmc_prefix
     draw_msg, (*pState).viewer, 'Creating 1-min averages of RMC files ...'
     ;; one_min_vRMC_AN11, daily_rmc_dir, min_rmc_dir, daily_rmc_temp, $
     ;;                    rmc_outcolumns, rmc_bear_field, rmc_vec_field, $
     ;;                    rmc_stamp, rmc_header
  ENDIF ELSE BEGIN
     print, 'Skipping splitting of RMC files ...'
  ENDELSE
  
  ;; ==============================================================
  ;; Call process if radiation day splitter selector was set to Yes

  ;; NOTE [SPL]: It's interesting that this section doesn't match the flow
  ;; chart; process_RAD is called *only* if we choose to do the daily
  ;; splitter and, further, there's no decision about whether to use
  ;; OMG files, as shown in the flow chart.

  IF check_ctrl_var((*pState).rad_split) EQ 1 THEN BEGIN
     draw_msg, (*pState).viewer, 'Creating daily radiation files ...'
     ;; day_splitter_MET, rad_SDATE, rad_EDATE, raw_rad_dir, daily_rad_dir, $
     ;;                   Rad_template, rad_header, rad_timing, rad_outcolumns, $
     ;;                   rad_stamp, rad_prefix
     draw_msg, (*pState).viewer, 'Processing radiation data ...'
     ;; process_RAD, daily_rad_dir, proc_rad_dir, min_rmc_dir, $
     ;;              met_daily_dir, rootdir
  ENDIF ELSE BEGIN
     print, 'Skipping processing of radiation data ...'
  ENDELSE
  
  ;; ======================================================
  ;; Call process if MET processing selector was set to Yes

  IF check_ctrl_var((*pState).met_proc) EQ 1 THEN BEGIN
     draw_msg, (*pState).viewer, 'Processing MET Files ...'
     ;; process_met_an11, met_daily_dir, min_gyro_dir, min_rmc_dir, $
     ;;                   met_proc_dir, dailymet_temp, met_log_name, $
     ;;                   met_log_dir, met_log_template,rootdir
  ENDIF ELSE BEGIN
     print, 'Skipping processing of MET files ...'
  ENDELSE
  
  ;; =================================================================
  ;; Call process if flux file standardization selector was set to Yes

  IF check_ctrl_var((*pState).flux_std) EQ 1 THEN BEGIN
     draw_msg, (*pState).viewer, 'Standardizing flux files ...'
     ;; ec_tables2011, raw_ec_dir, std_ec_dir, rootdir
  ENDIF ELSE BEGIN
     print, 'Skipping standardization of flux files ...'
  ENDELSE

  ;; ===========================================================
  ;; Call process if flux run extraction selector was set to Yes

  IF check_ctrl_var((*pState).flux_ext) EQ 1 THEN BEGIN
     draw_msg, (*pState).viewer, 'Finding flux runs ...'
     ;; filter_flux_runs_AN11, fluxperiod, flux_SDATE, flux_EDATE, t_wind, $
     ;;                        t_SOG, t_head, run_file, met_proc_dir, $
     ;;                        proc_flux_dir_root, met_proctemp
     ;; Run the script for extracting suitable flux runs
     draw_msg, (*pState).viewer, 'Extracting flux runs ...'
     ;; Extract_Flux_Runs_AN11, std_ec_dir, ec_run_dir, run_file_full
     ;; Run the script to extract GYRO runs
     draw_msg, (*pState).viewer, 'Extracting GYRO Runs ...'
     ;; Extract_GYRO_Runs_AN11, cor_gyro_dir, gyro_run_dir, run_file_full, $
     ;;                         rootdir
     ;; Run the script to extract RMC runs
     draw_msg, (*pState).viewer, 'Extracting RMC Runs ...'
     ;; Extract_RMC_Runs_AN11, std_rmc_dir, rmc_run_dir, run_file_full, $
     ;;                        rootdir
  ENDIF ELSE BEGIN
     print, 'Skipping flux run extraction ... '
  ENDELSE

  ;; ========================================================
  ;; Call process if flux calculation selector was set to Yes

  IF check_ctrl_var((*pState).flux_calc) EQ 1 THEN BEGIN
     draw_msg, (*pState).viewer, 'Calculating fluxes ...'
     ;; resolve_routine, 'ec_main_an11v2'
     ;; ec_main_AN11v2, DigitalDataStream=ddscode, outname, proc_flux_dir, $
     ;;                 met_log_dir, closed_path_log, met_log_template, $
     ;;                 lag_plot_dir, spec_plot_dir, ogvie_plot_dir, $
     ;;                 sonic_timeseries_dir, co2_timeseries_dir, $
     ;;                 nav_timeseries_dir, ec_run_dir, gyro_run_dir, $
     ;;                 rmc_run_dir, run_file_full, dailynav_temp, $
     ;;                 proc_rad_dir, foot_plot_dir, fixlagtemplate, $
     ;;                 fixlagfile, lag_TU_plot_dir, rootdir, calc_sdate, $
     ;;                 calc_edate, pState  
  ENDIF ELSE BEGIN
     print, 'Skipping calculation of fluxes ...'
     draw_msg, (*pState).viewer, '*FINISHED PROCESSING*'
  ENDELSE
  
  RETURN
  
END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; flux_capacitor.pro ends here
