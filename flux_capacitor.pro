;; Author: Sebastian P. Luque
;; Created: 2013-11-26T23:41:42+0000
;; Last-Updated: 2015-06-30T19:50:23+0000
;;           By: Sebastian P. Luque
;; ------------------------------------------------------------------------
;;; Code:

@control.inc

;; NAV

;; STD_NAV, expand_path('~/tmp/ArcticNet2011/NAV'), $
;;          expand_path('~/tmp/ArcticNet2011/NAV/STD'), $
;;          nav_std_header, /overwrite
;; Below is what should be used, but NAV files are messed up.  Fortunately,
;; we can use the server time from them and provide that to day_splitter.
;; This is because, luckily, server and UTC time seem to be almost the same
;; in these files (there's a difference of a few seconds), and although UTC
;; times are often missing, the server times are not.  So note that we are
;; using the server time *as if* it was UTC (this may introduce some time
;; mismatches if the differences are large).
STD2, nav_idir, nav_std_dir, nav_raw_template, 6, 1, nav_raw_keep_fields, $
      keep_types=nav_raw_keep_types, file_type='NAV', /overwrite
DAY_SPLITTER, '20110719', '20111022', $
              nav_std_dir, nav_daily_dir, nav_std_template, 0, $
              nav_daily_rate, nav_stamp, /overwrite
AVERAGE_SERIES, nav_daily_dir, nav_avg_dir, nav_std_template, 0, $
                nav_daily_rate, 60, $
                angle_fields=nav_angle_flds, $
                magnitude_fields=nav_mag_flds, /overwrite

;; MET
STD, met_idir, met_std_dir, met_raw_template, 1, met_raw_keep_fields, $
     file_type='MET', /overwrite
DAY_SPLITTER, '20110719', '20111022', $
              met_std_dir, met_daily_dir, met_std_template, 0, $
              met_daily_rate, met_stamp, /overwrite

;; RMC

;; Ship
STD2, rmc_ship_idir, rmc_ship_std_dir, rmc_ship_raw_template, 3, 0, $
      rmc_ship_raw_keep_fields, file_type='RMC', /overwrite
;; Ocean Mapping Group
STD, rmc_omg_idir, rmc_omg_std_dir, rmc_omg_raw_template, 0, $
     rmc_omg_raw_keep_fields, file_type='OMG', /overwrite


;; GYRO

;; Ship
STD_GYRO, gyro_ship_idir, gyro_ship_std_dir, gyro_ship_raw_template, 0, $
          rmc_ship_std_dir, rmc_ship_std_template, 0, 6, $
          gyro_ship_raw_keep_fields, /overwrite
;; Here we need to generate exactly 10 Hz
FILTER_SERIES, gyro_ship_std_dir, gyro_ship_std10Hz_dir, $
               gyro_ship_std_template, 0, gyro_hdg_fld, $
               gyro_daily_rate, /overwrite
;; Ocean Mapping Group
STD, gyro_omg_idir, gyro_omg_std_dir, gyro_omg_raw_template, 0, $
     gyro_omg_raw_keep_fields, file_type='OMG', /overwrite
;; Here we need to generate exactly 10 Hz
FILTER_SERIES, gyro_omg_std_dir, gyro_omg_std10Hz_dir, $
               gyro_omg_std_template, 0, gyro_hdg_fld, $
               gyro_daily_rate, /overwrite


;; Split RMC/GYRO and OMG NVG/HDG
;; Note that we're creating RMC daily files at 1-s intervals, even though
;; the input is at 2 s.  This is required because we can't tell whether the
;; series follows odd or even seconds.  Perhaps data should be interpolated
;; at 1 s via filter_series (just as was done for gyro files) prior to
;; this.
DAY_SPLITTER, '20110719', '20111022', $
              rmc_ship_std_dir, rmc_ship_daily_dir, rmc_ship_std_template, $
              0, rmc_daily_rate, rmc_ship_stamp, /overwrite
;; For splitting gyro files, we'll have to loop through each day or a small
;; number of days to avoid memory problems
d0=julday(07, 19, 2011, 0)      ; start date (full request)
d1=julday(10, 22, 2011, 0)      ; end date (full request)
subn=10                         ; n days to process at a time
days=timegen(start=d0, final=d1 + 1, step_size=1)
idx=lindgen(n_elements(days))
blks=histogram(idx, binsize=subn) ; how to allocate blocks of days
idx_ends=long(total(blks, /cumulative)) - 1
idx_begs=idx_ends - (blks - 1)
FOREACH i, lindgen(n_elements(idx_begs)) DO BEGIN
   caldat, days[idx_begs[i]], mo0, dd0, yyyy0
   caldat, days[idx_ends[i]], mo1, dd1, yyyy1
   begi=strcompress(string(yyyy0, format='(i04)') + $
                    string(mo0, format='(i02)') + $
                    string(dd0, format='(i02)'))
   endi=strcompress(string(yyyy1, format='(i04)') + $
                    string(mo1, format='(i02)') + $
                    string(dd1, format='(i02)'))
   ;; print, begi, endi
   DAY_SPLITTER, begi, endi, $
                 gyro_ship_std10Hz_dir, gyro_ship_daily_dir, $
                 gyro_ship_std_template, 0, gyro_daily_rate, $
                 gyro_ship_stamp, /overwrite
ENDFOREACH
;; DAY_SPLITTER, '20110719', '20111022', $
;;               expand_path('~/tmp/ArcticNet2011/GYRO/10Hz'), $
;;               expand_path('~/tmp/ArcticNet2011/GYRO/Daily'), $
;;               'Templates/gyro_std_template.sav', 0, gyro_daily_rate, $
;;               gyro_stamp, /overwrite
DAY_SPLITTER, '20110719', '20111022', $
              rmc_omg_std_dir, rmc_omg_daily_dir, rmc_omg_std_template, $
              0, rmc_daily_rate, rmc_omg_stamp, /overwrite
FOREACH i, lindgen(n_elements(idx_begs)) DO BEGIN
   caldat, days[idx_begs[i]], mo0, dd0, yyyy0
   caldat, days[idx_ends[i]], mo1, dd1, yyyy1
   begi=strcompress(string(yyyy0, format='(i04)') + $
                    string(mo0, format='(i02)') + $
                    string(dd0, format='(i02)'))
   endi=strcompress(string(yyyy1, format='(i04)') + $
                    string(mo1, format='(i02)') + $
                    string(dd1, format='(i02)'))
   ;; print, begi, endi
   DAY_SPLITTER, begi, endi, $
                 gyro_omg_std10Hz_dir, gyro_omg_daily_dir, $
                 gyro_omg_std_template, 0, gyro_daily_rate, $
                 gyro_omg_stamp, /overwrite
ENDFOREACH
;; DAY_SPLITTER, '20110719', '20111022', $
;;               expand_path('~/tmp/ArcticNet2011/OMG/HDG/10Hz'), $
;;               expand_path('~/tmp/ArcticNet2011/OMG/HDG/Daily'), $
;;               'Templates/omg_hdg_std_template.sav', 0, gyro_daily_rate, $
;;               omg_hdg_stamp, /overwrite
;; 1-min averages
AVERAGE_SERIES, rmc_ship_daily_dir, rmc_ship_avg_dir, $
                rmc_ship_std_template, 0, rmc_daily_rate, 60, $
                angle_fields=rmc_ship_cog_fld, $
                magnitude_fields=rmc_ship_sog_fld, /overwrite
AVERAGE_SERIES, gyro_ship_daily_dir, gyro_ship_avg_dir, $
                gyro_ship_std_template, 0, gyro_daily_rate, 60, $
                angle_fields=gyro_hdg_fld, magnitude_fields=[-1], /overwrite
AVERAGE_SERIES, rmc_omg_daily_dir, rmc_omg_avg_dir, $
                rmc_omg_std_template, 0, rmc_daily_rate, 60, $
                angle_fields=rmc_omg_cog_fld, $
                magnitude_fields=rmc_omg_sog_fld, /overwrite
AVERAGE_SERIES, gyro_omg_daily_dir, gyro_omg_avg_dir, $
                gyro_omg_std_template, 0, gyro_daily_rate, 60, $
                angle_fields=gyro_hdg_fld, magnitude_fields=[-1], /overwrite

;; RAD
STD, rad_idir, rad_std_dir, rad_raw_template, 1, rad_raw_keep_fields, $
     file_type='RAD', /overwrite
DAY_SPLITTER, '20110719', '20111022', $
              rad_std_dir, rad_daily_dir, rad_std_template, 0, $
              rad_daily_rate, rad_stamp, /overwrite
;; Process RAD.  Here one can choose to use the ship or OMG files for RMC
;; data; simply specify the chosen directory for the RMC_DIR argument for
;; the PROCESS_RAD procedure.
PROCESS_RAD, rad_daily_dir, rad_proc_dir, rad_std_template, 0, $
             rmc_omg_avg_dir, rmc_omg_avg_template, 0, rmc_latlon_flds, $
             met_daily_dir, met_std_template, 0, met_par_fld, /overwrite

;; Processing MET
PROCESS_MET, met_daily_dir, met_proc_dir, met_std_template, 0, $
             rmc_omg_avg_dir, rmc_omg_avg_template, 0, rmc_pull_flds, $
             gyro_omg_avg_dir, gyro_omg_avg_template, 0, gyro_pull_flds, $
             met_log_file, met_log_template, 0, 5, 10, /overwrite

;; Flux (eddy-covariance) files

;; Note that we re-organized input files, placing them according to common
;; structure (e.g. V4.1_4.2 holds input files that can be read with the
;; same template, and so on)
STD_EC, ec_idir1, ec_std_dir, ec_raw_template1, 1, $
        ec_41_42_raw_keep_fields, ec_std_fields, ec_std_types, /overwrite
STD_EC, ec_idir2, ec_std_dir, ec_raw_template2, 1, $
        ec_43_raw_keep_fields, ec_std_fields, ec_std_types, /overwrite
STD_EC, ec_idir3, ec_std_dir, ec_raw_template3, 1, $
        ec_44_45_raw_keep_fields, ec_std_fields, ec_std_types, /overwrite
STD_EC, ec_idir4, ec_std_dir, ec_raw_template4, 1, $
        ec_46_47_raw_keep_fields, ec_std_fields, ec_std_types, /overwrite
;; For splitting the flux files, we'll have to loop through each day
;; or a small number of days to avoid memory problems
;; d0=julday(07, 19, 2011, 0)      ; start date (full request)
d0=julday(10, 02, 2011, 0)      ; start date (full request)
d1=julday(10, 22, 2011, 0)      ; end date (full request)
subn=5                         ; n days to process at a time
days=timegen(start=d0, final=d1 + 1, step_size=1)
idx=lindgen(n_elements(days))
blks=histogram(idx, binsize=subn) ; how to allocate blocks of days
idx_ends=long(total(blks, /cumulative)) - 1
idx_begs=idx_ends - (blks - 1)
FOREACH i, lindgen(n_elements(idx_begs)) DO BEGIN
   caldat, days[idx_begs[i]], mo0, dd0, yyyy0
   caldat, days[idx_ends[i]], mo1, dd1, yyyy1
   begi=strcompress(string(yyyy0, format='(i04)') + $
                    string(mo0, format='(i02)') + $
                    string(dd0, format='(i02)'))
   endi=strcompress(string(yyyy1, format='(i04)') + $
                    string(mo1, format='(i02)') + $
                    string(dd1, format='(i02)'))
   ;; print, begi, endi
   DAY_SPLITTER, begi, endi, $
                 ec_std_dir, ec_daily_dir, ec_std_template, 0, $
                 ec_daily_rate, ec_stamp, /overwrite
ENDFOREACH

;; Average processed MET files for flux periods (this is for checking
;; whether flux periods have valid MET data)
AVERAGE_SERIES, met_proc_dir, met_proc_avg_dir, $
                met_proc_template, 0, met_daily_rate, ec_period, $
                angle_fields=met_proc_angle_flds, $
                magnitude_fields=met_proc_mag_flds, /overwrite
;; Small utility to remove some fields so that FILTER_MET interprets fields
;; correctly
REMOVE_FIELDS_ASCII, met_proc_avg_dir, met_proc_avg_pre_template, $
                     ['wind_direction_sd', 'cog_sd', 'sog_sd']
FILTER_MET, met_proc_avg_dir, met_proc_dir, met_proc_avg_template, 0, $
            ec_period, met_proc_template, 0, met_daily_rate, $
            met_filter_flds, met_filter_thrs

;; Extract flux periods

;; We will need to filter the standard GYRO to 10 Hz
SUBSET_FLUX, ec_daily_dir, ec_std_template, 0, ec_daily_rate, $
             met_proc_avg_dir, met_proc_avg_diag_template, 0, $
             met_proc_diag_fld, ec_period, rmc_omg_daily_dir, $
             rmc_omg_std_template, 0, gyro_omg_daily_dir, $
             gyro_omg_std_template, 0, rad_proc_dir, $
             rad_std_template, 0, /overwrite
;; ;; Run a subset.
;; SUBSET_FLUX, expand_path('~/tmp/ArcticNet2011_Oct/Flux'), $
;;              'Templates/ec_std_template.sav', 0, ec_daily_rate, $
;;              expand_path('~/tmp/ArcticNet2011_Oct/MET'), $
;;              'Templates/met_processed_avg_diag_template.sav', 0, 33, $
;;              ec_period, $
;;              expand_path('~/tmp/ArcticNet2011/OMG/NVG/Daily'), $
;;              'Templates/omg_nav_std_template.sav', 0, $
;;              expand_path('~/tmp/ArcticNet2011/OMG/HDG/Daily'), $
;;              'Templates/omg_hdg_std_template.sav', 0, $
;;              expand_path('~/tmp/ArcticNet2011/RAD/Processed'), $
;;              'Templates/rad_std_template.sav', 0, /overwrite


;; Calculate fluxes

FLUX, ec_period_dir, ec_std_template, 0, ec_daily_rate, $
      met_proc_avg_dir, met_proc_avg_diag_template, 0, met_proc_diag_fld, $
      ec_period, rmc_omg_period_dir, rmc_omg_std_template, 0, $
      gyro_omg_period_dir, gyro_omg_std_template, 0, $
      rad_period_dir, rad_proc_template, 0, $
      motpak_offset, sog_thr, lfreq_thr, hfreq_thr, xover_freq_thr, $
      closed_path_log_file, closed_path_log_template, 0, 5, $
      cl_log_status_fld, ec_motcorr_dir, ec_footprint_dir, ec_fluxes_file, $
      /overwrite
;; Debugging
FLUX, expand_path('~/tmp/EC'), ec_std_template, 0, ec_daily_rate, $
      met_proc_avg_dir, met_proc_avg_diag_template, 0, met_proc_diag_fld, $
      ec_period, rmc_omg_period_dir, rmc_omg_std_template, 0, $
      gyro_omg_period_dir, gyro_omg_std_template, 0, $
      rad_period_dir, rad_proc_template, 0, $
      motpak_offset, sog_thr, lfreq_thr, hfreq_thr, xover_freq_thr, $
      closed_path_log_file, closed_path_log_template, 0, 5, $
      cl_log_status_fld, expand_path('~/tmp/EC/Motion_Corrected'), $
      ec_footprint_dir, expand_path('~/tmp/fluxes_junk.csv'), /overwrite



;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
