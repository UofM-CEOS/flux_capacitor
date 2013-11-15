;; Author: Sebastian Luque
;; Created: 2013-11-12T17:07:28+0000
;; Last-Updated: 2013-11-15T20:51:46+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;; 
;; 
;; PURPOSE:
;; 
;; 
;; 
;; CATEGORY:
;; 
;; 
;; 
;; CALLING SEQUENCE:
;; 
;; 
;; 
;; INPUTS:
;; 
;; 
;; 
;; OPTIONAL INPUTS:
;; 
;; 
;; 
;; KEYWORD PARAMETERS:
;; 
;; 
;; 
;; OUTPUTS:
;; 
;; 
;; 
;; OPTIONAL OUTPUTS:
;; 
;; 
;; 
;; COMMON BLOCKS:
;; 
;; 
;; 
;; SIDE EFFECTS:
;; 
;; 
;; 
;; RESTRICTIONS:
;; 
;; 
;; 
;; PROCEDURE:
;; 
;; 
;; 
;; EXAMPLE:
;; 
;; 
;; 
;;- ----------------------------------------------------------------------
;;; Code:

PRO FLUX, IDIR, ITEMPLATE_SAV, TIME_IDX, ISAMPLE_RATE, $
          DIAG_DIR, DIAG_ITEMPLATE_SAV, DIAG_TIME_IDX, DIAG_IDX, $
          EC_PERIOD, RMC_DIR, RMC_ITEMPLATE_SAV, RMC_TIME_IDX, $
          GYRO_DIR, GYRO_ITEMPLATE_SAV, GYRO_TIME_IDX, $
          RAD_DIR, RAD_ITEMPLATE_SAV, RAD_TIME_IDX, $
          MOTPAK_OFFSET, SOG_THR, LFREQ_THR, HFREQ_THR, XOVER_FREQ_THR, $
          LOG_FILE, LOG_ITEMPLATE_SAV, LOG_TIME_BEG_IDX, $
          LOG_TIME_END_IDX, LOG_STATUS_IDX, MOT_CORR_ODIR, OFILE, $
          SERIAL=SERIAL, OVERWRITE=OVERWRITE

  log_file_info=file_info(log_file)
  IF log_file_info.regular NE 1 THEN $
     message, 'Log file is not a regular file.  Exiting'
  IF ((n_elements(log_time_beg_idx) NE 1) OR (log_time_beg_idx LT 0)) THEN $
     message, 'LOG_TIME_BEG_IDX must be a scalar >= zero'
  IF ((n_elements(log_time_end_idx) NE 1) OR (log_time_end_idx LT 0)) THEN $
     message, 'LOG_TIME_END_IDX must be a scalar >= zero'
  IF ((n_elements(log_status_idx) NE 1) OR (log_status_idx LT 0)) THEN $
     message, 'LOG_STATUS_IDX must be a scalar >= zero'
  idir_files=file_search(idir + path_sep() + '*', count=nidir_files, $
                         /nosort, /fold_case, /test_regular)
  diag_files=file_search(diag_dir + path_sep() + '*', count=ndiag_files, $
                         /nosort, /fold_case, /test_regular)
  rmc_files=file_search(rmc_dir + path_sep() + '*', count=nrmc_files, $
                        /nosort, /fold_case, /test_regular)
  gyro_files=file_search(gyro_dir + path_sep() + '*', count=ngyro_files, $
                         /nosort, /fold_case, /test_regular)
  rad_files=file_search(rad_dir + path_sep() + '*', count=nrad_files, $
                        /nosort, /fold_case, /test_regular)
  IF (nidir_files EQ 0) OR (ndiag_files EQ 0) OR $
     (nrmc_files EQ 0) OR (ngyro_files EQ 0) OR $
     (nrad_files EQ 0) THEN $
        message, 'No files in at least one of IDIR, DIAG_DIR, ' + $
                 'RMC_DIR, GYRO_DIR, or RAD_DIR'

  ;; Parse MET Diag files
  restore, diag_itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  diag_template=itemplate
  diag_field_names=strlowcase(diag_template.FIELDNAMES)
  ;; Break file names and extract the piece to match
  diag_filesl=strsplit(diag_files, '_.', /extract)
  diag_files_a=diag_filesl.toArray(/transpose) ; array with pieces in cols
  diag_files_mstr_dims=size(diag_files_a, /dimensions)
  ;; We want to match against the concatenated 4th to last piece (period
  ;; start is obtained from each line later)
  diag_files_mstr=diag_files_a[diag_files_mstr_dims[0] - 4, *]

  ;; Parse RMC files
  restore, rmc_itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  rmc_template=itemplate
  rmc_field_names=strlowcase(rmc_template.FIELDNAMES)
  ;; Break file names and extract the piece to match
  rmc_filesl=strsplit(rmc_files, '_.', /extract)
  rmc_files_a=rmc_filesl.toArray(/transpose) ; array with pieces in cols
  rmc_files_mstr_dims=size(rmc_files_a, /dimensions)
  ;; We want to match against the concatenated 3rd and 2nd to last pieces
  rmc_files_mstr=rmc_files_a[rmc_files_mstr_dims[0] - 3, *] + $
                 rmc_files_a[rmc_files_mstr_dims[0] - 2, *]

  ;; Parse Gyro files
  restore, gyro_itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  gyro_template=itemplate
  gyro_field_names=strlowcase(gyro_template.FIELDNAMES)
  ;; Break file names and extract the piece to match
  gyro_filesl=strsplit(gyro_files, '_.', /extract)
  gyro_files_a=gyro_filesl.toArray(/transpose) ; array with pieces in cols
  gyro_files_mstr_dims=size(gyro_files_a, /dimensions)
  ;; We want to match against the concatenated 3rd and 2nd to last pieces
  gyro_files_mstr=gyro_files_a[gyro_files_mstr_dims[0] - 3, *] + $
                  gyro_files_a[gyro_files_mstr_dims[0] - 2, *]

  ;; Parse RAD files
  restore, rad_itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  rad_template=itemplate
  rad_field_names=strlowcase(rad_template.FIELDNAMES)
  ;; Break file names and extract the piece to match
  rad_filesl=strsplit(rad_files, '_.', /extract)
  rad_files_a=rad_filesl.toArray(/transpose) ; array with pieces in cols
  rad_files_mstr_dims=size(rad_files_a, /dimensions)
  ;; We want to match against the concatenated 4th and 2nd to last pieces
  rad_files_mstr=rad_files_a[rad_files_mstr_dims[0] - 4, *] + $
                 rad_files_a[rad_files_mstr_dims[0] - 2, *]

  ;; Parse log template
  restore, log_itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  log_template=itemplate
  log_field_names=strlowcase(log_template.FIELDNAMES)
  log_is_begtime_fld=log_template.FIELDGROUPS EQ log_time_beg_idx
  log_is_endtime_fld=log_template.FIELDGROUPS EQ log_time_end_idx

  ;; Parse input flux files (last because we don't make copy of template)
  restore, itemplate_sav
  field_names=strlowcase(itemplate.FIELDNAMES)
  time_flds=where(itemplate.FIELDGROUPS EQ time_idx)
  ;; Break file names and extract the piece to match
  flux_filesl=strsplit(idir_files, '_.', /extract)
  flux_files_a=flux_filesl.toArray(/transpose) ; array with pieces in cols
  flux_files_mstr_dims=size(flux_files_a, /dimensions)
  ;; We want to match against the concatenated 3rd and 2nd to last pieces
  flux_files_mstr=flux_files_a[flux_files_mstr_dims[0] - 3, *] + $
                  flux_files_a[flux_files_mstr_dims[0] - 2, *]
  ;; [Original comment: Temporarily set g for filtering purposes... true g
  ;; will be calculated later... we'll set this low to make sure we filter
  ;; properly]
  g=8.0

  ;; Read log file and obtain beginning and end times
  log=read_std2_file(log_file, log_template, $
                     log_time_beg_idx, log_time_end_idx)
  log_names=strlowcase(tag_names(log))
  log_tbeg_loc=where(log_names EQ log_field_names[log_time_beg_idx])
  log_tend_loc=where(log_names EQ log_field_names[log_time_end_idx])
  log_beg_times=log.(log_tbeg_loc)
  log_end_times=log.(log_tend_loc)
  ;; Convert to julian
  log_tbeg_jd=reform(julday(long(log_beg_times[1, *]), $
                            long(log_beg_times[2, *]), $
                            long(log_beg_times[0, *]), $
                            long(log_beg_times[3, *]), $
                            long(log_beg_times[4, *]), $
                            fix(log_beg_times[5, *])))
  log_tend_jd=reform(julday(long(log_end_times[1, *]), $
                            long(log_end_times[2, *]), $
                            long(log_end_times[0, *]), $
                            long(log_end_times[3, *]), $
                            long(log_end_times[4, *]), $
                            fix(log_end_times[5, *])))
  log_status=log.(where(log_names EQ log_field_names[log_status_idx]))
  logn=(size(log_beg_times, /dimensions))[1]

  ;; Perhaps set up hashes to hold output data from all valid flux runs at
  ;; this point, before starting outermost loop.

  FOR k=0, ndiag_files - 1 DO BEGIN
     dfile=diag_files[k]
     diag=read_std_file(dfile, diag_template, diag_time_idx)
     diag_names=strlowcase(tag_names(diag))
     time_loc=where(diag_names EQ diag_field_names[diag_time_idx])
     diag_times=diag.(time_loc)
     diag_flag=diag.(where(diag_names EQ diag_field_names[diag_idx]))
     fluxable=where(diag_flag EQ 0, nfluxable)
     IF nfluxable EQ 0 THEN CONTINUE
     dfile_mstr=diag_files_mstr[k]

     FOREACH fperiod, fluxable DO BEGIN
        ;; Time stamp.  Note we are flooring the seconds (as in SUBSET_FLUX)
        tstamp=string(diag_times[3, fperiod] + diag_times[4, fperiod] + $
                      diag_times[5, fperiod], format='(i06)')

        ;; Read matching flux file
        flux_pair=where(flux_files_mstr EQ (dfile_mstr + tstamp), mcount)
        IF mcount LT 1 THEN BEGIN
           message, 'No matching flux file found. Skipping.', /CONTINUE
           CONTINUE
        ENDIF
        flux=read_std_file(idir_files[flux_pair], itemplate, time_idx)
        flux_names=strlowcase(tag_names(flux))
        flux_time_loc=where(flux_names EQ field_names[time_idx])
        flux_times=flux.(flux_time_loc)
        flux_times_dims=size(flux_times, /dimensions)
        ftimes_s=flux_times[5, *]
        ftimes_s=fix(ftimes_s) + $
                 (round((double(ftimes_s) - fix(ftimes_s)) * 10) / 10.0)
        flux_jd=reform(julday(long(flux_times[1, *]), $
                              long(flux_times[2, *]), $
                              long(flux_times[0, *]), $
                              long(flux_times[3, *]), $
                              long(flux_times[4, *]), $
                              double(ftimes_s)))
        ;; Get a file name prefix to be shared by the output files from
        ;; this period
        iname=strsplit(file_basename(idir_files[flux_pair]), '.', /extract)
        iname_prefix=iname[0]

        ;; [Original comment: no need recalibrate motion sensor for this
        ;; experiment.  Read in accelerations in RH coordinate system,
        ;; convert to m/s2].  All these values could just be placed back on
        ;; the intput structure, and avoid cluttering the workspace and
        ;; memory so much
        accel_x=flux.accel_z * 9.81  ; accel_z on the tower
        accel_y=-flux.accel_x * 9.81 ; accel_x on the tower
        accel_z=-flux.accel_y * 9.81 ; accel_y on the tower
        ;; Put acceleration components in 3-column array and make copy to
        ;; keep uncorrected data
        accel=transpose([[accel_x], [accel_y], [accel_z]])
        accel_raw=accel
        ;; Original comment: read in angular rates in RH coordinate system,
        ;; convert to rad/s
        rate_phi=flux.rate_z * !DTOR    ; rate_z on the tower
        rate_theta=-flux.rate_x * !DTOR ; rate_x on the tower
        rate_shi=-flux.rate_y * !DTOR   ; rate_y on the tower
        ;; Put rate components in 3-column array and make copy to keep
        ;; uncorrected data
        rate=transpose([[rate_phi], [rate_theta], [rate_shi]])
        rate_raw=rate
        ;; Extract the wind components based on whether we want serial or
        ;; analogue data
        wind_u=keyword_set(serial) ? $
               flux.u_wind_serial : $
               flux.u_wind_analogue
        wind_v=keyword_set(serial) ? $
               flux.v_wind_serial : $
               flux.v_wind_analogue
        wind_w=keyword_set(serial) ? $
               flux.w_wind_serial : $
               flux.w_wind_analogue
        ;; Put wind components together in 3-column array, and make a copy
        ;; to keep uncorrected data
        wind=transpose([[wind_u], [wind_v], [wind_w]])
        wind_raw=wind
        sonic_temperature=keyword_set(serial) ? $
                          flux.sonic_temperature_serial : $
                          flux.sonic_temperature_analogue

        ;; [Original comment: use the closed path log to identify any data
        ;; that should be removed... Set to 'NAN'].  Check the log status
        ;; flag (initial value=0 -> not used in log, so we use to OK data)
        status_flag=intarr(flux_times_dims[1])
        FOREACH logi, indgen(logn) DO BEGIN
           flagi=where((flux_jd GE log_tbeg_jd[logi]) AND $
                       (flux_jd LT log_tend_jd[logi]), nflagi)
           IF nflagi GT 0 THEN status_flag[flagi]=log_status[logi]
        ENDFOREACH
        bad_status=where(status_flag NE 0, nbad_status)
        IF nbad_status GT 0 THEN BEGIN
           flux.co2_cl[bad_status]=!VALUES.D_NAN
           flux.h20_cl[bad_status]=!VALUES.D_NAN
           flux.pressure_cl[bad_status]=!VALUES.D_NAN
           flux.temperature_cl[bad_status]=!VALUES.D_NAN
        ENDIF

        ;; [Original comment: create flags for the 4 possible sources of
        ;; "bad" data, flag=0 means data good]
        open_flag=0
        closed_flag=0
        sonic_flag=0
        motion_flag=0

        ;; [Original comment: check for any significant number of 'NAN's
        ;; (not worried about the odd one scattered here and there)]
        bad_co2_op=where(~finite(flux.co2_op), nbad_co2_op)
        bad_h2o_op=where(~finite(flux.h2o_op), nbad_h2o_op)
        bad_diag_op=where(~finite(fix(flux.diag_op)), nbad_diag_op)
        ;; [Original comment: set open flag if gt 2% of records are 'NAN']
        IF (float(nbad_co2_op) / flux_times_dims[1] * 100.0 GT 2) OR $
           (float(nbad_h2o_op) / flux_times_dims[1] * 100.0 GT 2) OR $
           (float(nbad_diag_op) / flux_times_dims[1] * 100.0 GT 2) THEN $
              open_flag=1

        bad_u=where(~finite(wind_u), nbad_u)
        bad_v=where(~finite(wind_v), nbad_v)
        bad_w=where(~finite(wind_w), nbad_w)
        bad_sonic_temperature=where(~finite(sonic_temperature), $
                                    nbad_sonic_temperature)
        sonicNANs=float(nbad_u / flux_times_dims[1] * 100.0)
        ;; [Original comment: set wind flag if gt 2% of records are 'NAN']
        IF (float(nbad_u) / flux_times_dims[1] * 100.0 GT 2) OR $
           (float(nbad_v) / flux_times_dims[1] * 100.0 GT 2) OR $
           (float(nbad_w) / flux_times_dims[1] * 100.0 GT 2) OR $
           (float(nbad_sonic_temperature) / $
            flux_times_dims[1] * 100.0 GT 2) THEN sonic_flag=1

        bad_accel_x=where(~finite(accel_x), nbad_accel_x)
        bad_accel_y=where(~finite(accel_y), nbad_accel_y)
        bad_accel_z=where(~finite(accel_z), nbad_accel_z)
        bad_rate_phi=where(~finite(rate_phi), nbad_rate_phi)
        bad_rate_theta=where(~finite(rate_theta), nbad_rate_theta)
        bad_rate_shi=where(~finite(rate_shi), nbad_rate_shi)
        ;; [Original comment: set motion flag if gt 2% of records are 'NAN']
        IF (float(nbad_accel_x) / flux_times_dims[1] * 100.0 GT 2) OR $
           (float(nbad_accel_y) / flux_times_dims[1] * 100.0 GT 2) OR $
           (float(nbad_accel_z) / flux_times_dims[1] * 100.0 GT 2) OR $
           (float(nbad_rate_phi) / flux_times_dims[1] * 100.0 GT 2) OR $
           (float(nbad_rate_theta) / flux_times_dims[1] * 100.0 GT 2) OR $
           (float(nbad_rate_shi) / flux_times_dims[1] * 100.0 GT 2) THEN $
              motion_flag=1

        bad_co2_cl=where(~finite(flux.co2_cl), nbad_co2_cl)
        bad_h2o_cl=where(~finite(flux.h2o_cl), nbad_h2o_cl)
        bad_pressure_cl=where(~finite(flux.pressure_cl), nbad_pressure_cl)
        irgaNANs=float(nbad_co2_cl / flux_times_dims[1] * 100.0)
        ;; [Original comment: set closed flag if gt 2% of records are 'NAN']
        IF (float(nbad_co2_cl) / flux_times_dims[1] * 100.0 GT 2) OR $
           (float(nbad_h2o_cl) / flux_times_dims[1] * 100.0 GT 2) OR $
           (float(nbad_pressure_cl) / flux_times_dims[1] * 100.0 GT 2) THEN $
              closed_flag=1

        ;; [Original comment: now that we have looked for NANs, we may as
        ;; well fill in the NANs and any spikes using the shot filter]
        IF sonic_flag NE 1 THEN BEGIN
           wind_u=shot_filter(wind_u)
           wind_v=shot_filter(wind_v)
           wind_w=shot_filter(wind_w)
           sonic_temperature=shot_filter(sonic_temperature)
        ENDIF
        IF open_flag NE 1 THEN BEGIN
           flux.co2_op=shot_filter(flux.co2_op)
           flux.h2o_op=shot_filter(flux.h2o_op)
           flux.pressure_op=shot_filter(flux.pressure_op)
           ;; [Original comment: this is necessary to check if there is
           ;; still ugly shot noise... if there is, we need to skip this]
           bad_co2_op=where(abs(flux.co2_op - mean(flux.co2_op, /NAN)) GT $
                            (6.0 * stddev(flux.co2_op, /NAN)), nbad_co2_op)
           bad_h2o_op=where(abs(flux.h2o_op - mean(flux.h2o_op, /NAN)) GT $
                            (6.0 * stddev(flux.h2o_op, /NAN)), nbad_h2o_op)
           IF (nbad_co2_op GT 0) OR (nbad_h2o_op GT 0) THEN $
              open_flag=1
        ENDIF
        IF closed_flag NE 1 THEN BEGIN
           flux.co2_cl=shot_filter(flux.co2_cl)
           flux.h2o_cl=shot_filter(flux.h2o_cl)
           flux.pressure_cl=shot_filter(flux.pressure_cl)
           flux.temperature_cl=shot_filter(flux.temperature_cl)
           fco2_cl=flux.co2_cl
           bad_co2_cl=where(abs(fco2_cl - mean(fco2_cl, /NAN)) GT $
                            (6.0 * stddev(fco2_cl, /NAN)), nbad_co2_cl)
           fh2o_cl=flux.h2o_cl
           bad_h2o_cl=where(abs(fh2o_cl - mean(fh2o_cl, /NAN)) GT $
                            (6.0 * stddev(fh2o_cl, /NAN)), nbad_h2o_cl)
           fp_cl=flux.pressure_cl
           bad_pressure_cl=where(abs(fp_cl - mean(fp_cl, /NAN)) GT $
                                 (6.0 * stddev(fp_cl, /NAN)), $
                                 nbad_pressure_cl)
           ft_cl=flux.temperature_cl
           bad_temperature_cl=where(abs(ft_cl - mean(ft_cl, /NAN)) GT $
                                    (6.0 * stddev(ft_cl, /NAN)), $
                                    nbad_temperature_cl)
           IF (nbad_co2_cl GT 0) OR (nbad_h2o_cl GT 0) OR $
              (nbad_pressure_cl GT 0) THEN closed_flag=1
        ENDIF

        ;; Check for high or low diagnostic flags (set open flag if more
        ;; than 2% of records have them)
        bad_diags=where((flux.diag_op GT 249) OR (flux.diag_op LT 240), $
                        nbad_diags)
        IF (float(nbad_diags) / flux_times_dims[1] * 100.0 GT 2) THEN $
           open_flag=1

        ;; [Original comment: check for bad wind data: bad wind data can usually be diagnosed
        ;; by unusually high wind velocities.  this is most obvious in the
        ;; vertical wind where we wouldn't expect high values bad sonic
        ;; data can also turn up in the Tsonic before the wind, check the
        ;; deviation between Tsonic and mean air T]
        bad_vert_wind=where(abs(wind_w) GT 7, nbad_vert_wind)
        t_avg=diag.air_temperature[fperiod]
        bad_sonic_temperature=where(abs(sonic_temperature - t_avg) GT 7, $
                                    nbad_sonic_temperature)
        ;; ;; sonic_count appears to not be working... I temporarily disabled
        ;; ;; this check... RS Feb 2013
        ;; sonic_count=0  
        ;; Set wind flag high if gt 0.5% of records are frost contaminated
        IF (float(nbad_vert_wind) / flux_times_dims[1] * 100.0 GT 0.5) OR $
           (float(nbad_sonic_temperature) / $
            flux_times_dims[1] * 100.0 GT 0.5) THEN sonic_flag = 1

        ;; Check for Motion Pak data that are out of range
        bad_motionpak=where((accel_x GT g) OR (accel_y GT g), nbad_motionpak)
        IF nbad_motionpak GT 0 THEN motion_flag=1

        ;; if sonic data are bad, skip this period
        IF sonic_flag GT 0 THEN BEGIN
           message, 'Bad sonic anemometer data.  Skipping period', $
                    /CONTINUE
           CONTINUE
        ENDIF

        ;; [Original comment: check critical low f variabiles]
        IF (~finite(diag.rh_percent[fperiod])) OR $
           (~finite(t_avg)) THEN BEGIN
           message, 'RH or mean air temperature data unavailable.  ' + $
                    'Skipping period', /CONTINUE
           CONTINUE
        ENDIF

        ;; Read matching RAD file
        rad_pair=where(rad_files_mstr EQ (dfile_mstr + tstamp), mcount)
        IF mcount LT 1 THEN BEGIN
           message, 'No matching RAD file found. Skipping.', /CONTINUE
           CONTINUE
        ENDIF
        rad=read_std_file(rad_files[rad_pair], rad_template, $
                          rad_time_idx)
        ;; Mean radiation from RAD structure
        sw_avg=mean(rad.k_down, /nan)
        lw_avg=mean(rad.lw_in, /nan)

        ;; Read matching RMC file
        rmc_pair=where(rmc_files_mstr EQ (dfile_mstr + tstamp), mcount)
        IF mcount LT 1 THEN BEGIN
           message, 'No matching RMC file found. Skipping.', /CONTINUE
           ;; [Original comment: if we can't find the rmc file, we'll set
           ;; the motion flag to 1 and skip loading the GPS].  [SPL:
           ;; Skipping the whole period, for now, as I don't understand
           ;; what should happen in this case.]
           motion_flag=1
           CONTINUE
        ENDIF
        rmc=read_std_file(rmc_files[rmc_pair], rmc_template, rmc_time_idx)
        rmc_names=strlowcase(tag_names(rmc))
        rmc_time_loc=where(rmc_names EQ rmc_field_names[rmc_time_idx])
        rmc_times=rmc.(rmc_time_loc)
        rmc_times_dims=size(rmc_times, /dimensions)
        rtimes_s=rmc_times[5, *]
        rtimes_s=fix(rtimes_s) + $
                 (round((double(rtimes_s) - fix(rtimes_s)) * 10) / 10.0)
        rmc_jd=reform(julday(long(rmc_times[1, *]), $
                             long(rmc_times[2, *]), $
                             long(rmc_times[0, *]), $
                             long(rmc_times[3, *]), $
                             long(rmc_times[4, *]), $
                              double(rtimes_s)))
        ;; Find matching times
        match2, flux_jd, rmc_jd, flux_in_rmc
        flux_matches=where(flux_in_rmc GE 0, mcount, /null)
        IF mcount LT 1 THEN BEGIN
           message, 'No matching RMC records found.  Skipping file.', $
                    /informational
           ;; [Original comment: if we can't find the rmc file, we'll set
           ;; the motion flag to 1 and skip loading the GPS].  [SPL:
           ;; Skipping the whole period, for now, as I don't understand
           ;; what should happen in this case.]
           motion_flag=1
           CONTINUE
        ENDIF
        latitude=make_array(flux_times_dims[1], type=4, $
                            value=!VALUES.F_NAN)
        longitude=latitude
        sog=latitude
        cog=latitude
        latitude[flux_matches]=rmc.latitude[flux_matches]
        longitude[flux_matches]=rmc.longitude[flux_matches]
        sog[flux_matches]=rmc.sog[flux_matches]
        cog[flux_matches]=rmc.cog[flux_matches]

        ;; Read matching Gyro file
        gyro_pair=where(gyro_files_mstr EQ (dfile_mstr + tstamp), mcount)
        IF mcount LT 1 THEN BEGIN
           message, 'No matching Gyro file found. Skipping.', /CONTINUE
           ;; [Original comment: if we can't find the rmc file, we'll set
           ;; the motion flag to 1 and skip loading the GPS].  [SPL:
           ;; Skipping the whole period, for now, as I don't understand
           ;; what should happen in this case.]
           motion_flag=1
           CONTINUE
        ENDIF
        gyro=read_std_file(gyro_files[gyro_pair], gyro_template, $
                           gyro_time_idx)
        gyro_names=strlowcase(tag_names(gyro))
        gyro_time_loc=where(gyro_names EQ gyro_field_names[gyro_time_idx])
        gyro_times=gyro.(gyro_time_loc)
        gyro_times_dims=size(gyro_times, /dimensions)
        gtimes_s=gyro_times[5, *]
        gtimes_s=fix(gtimes_s) + $
                 (round((double(gtimes_s) - fix(gtimes_s)) * 10) / 10.0)
        gyro_jd=reform(julday(long(gyro_times[1, *]), $
                              long(gyro_times[2, *]), $
                              long(gyro_times[0, *]), $
                              long(gyro_times[3, *]), $
                              long(gyro_times[4, *]), $
                              double(gtimes_s)))
        ;; Find matching times
        match2, flux_jd, gyro_jd, flux_in_gyro
        flux_matches=where(flux_in_gyro GE 0, mcount, /null)
        IF mcount LT 1 THEN BEGIN
           message, 'No matching Gyro records found.  Skipping file.', $
                    /informational
           ;; [Original comment: if we can't find the rmc file, we'll set
           ;; the motion flag to 1 and skip loading the GPS].  [SPL:
           ;; Skipping the whole period, for now, as I don't understand
           ;; what should happen in this case.]
           motion_flag=1
           CONTINUE
        ENDIF
        heading=make_array(flux_times_dims[1], type=4, $
                           value=!VALUES.F_NAN)
        heading[flux_matches]=gyro.heading[flux_matches]

        ;; [Original comment: now fill in the gaps by applying a moving
        ;; average... In this case, we use a 100 sample window (10 sec)
        ;; moving average... may need to tweak this value].  [SPL: I THINK
        ;; THIS STEP SHOULD HAVE BEEN DONE EARLY DURING RMC PROCESSING.
        ;; Also, perhaps a simple linear interpolation is better; I don't
        ;; know why this moving average is used, where a window must be
        ;; specified and seems to be introducing bias.  Why aren't latitude
        ;; and longitude not similarly interpolated?]
        cog_xy=decompose(cog, sog)
        cog_x=smooth(cog_xy[0, *], 100, /nan, /edge_truncate)
        cog_y=smooth(cog_xy[1, *], 100, /nan, /edge_truncate)
        cog_sm=recompose(cog_x, cog_y)
        cog=reform(cog_sm[0, *])
        sog=reform(cog_sm[1, *])
        ;; Dummy heading magnitude for decomposition purposes only
        vheading=rebin([1], flux_times_dims[1], /sample)
        ;; Guard if missing heading
        noheading=where(~finite(heading), nnoheading)
        IF nnoheading GT 0 THEN vheading[noheading]=!VALUES.D_NAN
        heading_xy=decompose(heading, vheading)
        heading_x=smooth(heading_xy[0, *], 100, /nan, /edge_truncate)
        heading_y=smooth(heading_xy[1, *], 100, /nan, /edge_truncate)
        heading_sm=recompose(heading_x, heading_y)
        heading=reform(heading_sm[0, *])

        ;; [Original comment: Check to make sure that no 'NaNs' dropped
        ;; through... if they did, we'll have to skip this one]
        no_cog=where(~finite(cog), nno_cog)
        no_sog=where(~finite(sog), nno_sog)
        no_heading=where(~finite(heading), nno_heading)
        IF (nno_cog GT 0) OR (nno_sog GT 0) OR (nno_heading GT 0) THEN $
           motion_flag=1

        ;; Level sonic anemometer

        ;; [Original comment: here we level the sonic anemometer, as long
        ;; as we know the MEAN roll/pitch angles from an inclinometer, and
        ;; as long as those are in a L.H.S.] [SPL: WATCH LEVEL_SONIC
        ;; FUNCTION]
        IF finite(diag.roll[fperiod]) AND $
           finite(diag.pitch[fperiod]) THEN $
              wind=level_sonic(wind, diag.roll[fperiod] * !DTOR, $
                               -diag.pitch[fperiod] * !DTOR)

        ;; [Original comment: now let's calculate the raw mean w
        ;; value... this is going to be important for sort of tracking flow
        ;; distortion... We could do something more in depth, but we'll
        ;; keep it like this for now]
        wind_v_mean=mean(wind[2, *], /nan)

        ;; High frequency motion correction

        IF diag.sog[fperiod] GT sog_thr THEN BEGIN
           ;; [Original comment: shot filter the motion channels... this
           ;; helps with a problem where unreasonably high accelerations
           ;; cause a 'NaN' calculation
           accel[0, *]=shot_filter(accel[0, *])
           accel[1, *]=shot_filter(accel[1, *])
           accel[2, *]=shot_filter(accel[2, *])
           rate[0, *]=shot_filter(rate[0, *])
           rate[1, *]=shot_filter(rate[1, *])
           rate[2, *]=shot_filter(rate[2, *])

           ;; [Original comment: synthetically level the Motion Pak] [SPL:
           ;; WATCH LEVEL_MOTIONPAK FUNCTION]
           level=level_motionpak(accel, rate, 4, 4)
           IF ~finite(level[0]) THEN BEGIN
              motion_flag=1
              message, 'Motion-flagged', /informational
           ENDIF
           accel[0, *]=level[0, *]
           accel[1, *]=level[1, *]
           accel[2, *]=level[2, *]
           rate[0, *]=level[3, *]
           rate[1, *]=level[4, *]
           rate[2, *]=level[5, *]

           ;; [Original comment: demean the rates by extracting the
           ;; fluctuating component and setting the mean to 0]
           rate[0, *]=(rate[0, *] - mean(rate[0, *], /NAN, /DOUBLE)) + 0
           rate[1, *]=(rate[1, *] - mean(rate[1, *], /NAN, /DOUBLE)) + 0
           rate[2, *]=(rate[2, *] - mean(rate[2, *], /NAN, /DOUBLE)) + 0
  
           ;; Gravity (what was the purpose of the g variables at the
           ;; beginning of the procedure?)
           g=(mean(accel[2, *], /NAN, /DOUBLE))

           ;; Check to make sure the x/y accelerations are not greater than
           ;; g... this is indicative of a problem w/ the Motion Pak, and
           ;; the run needs to be skipped
           bad_mp=where((accel[0, *] GT g) OR $
                        (accel[1, *] GT g), nbad_mp)
           IF nbad_mp GT 0 THEN BEGIN
              motion_flag=1
              message, 'Invalid Motion Pak accelerations', /continue
              CONTINUE
           ENDIF

           ;; Integrate angular rates to give angles.  This also performs a
           ;; high pass filter, cutting off all frequencies below the
           ;; cutoff period (in this case, the cutoff period is set to 20s
           ;; (or 0.05Hz)). [SPL: WATCH INT1BYF FUNCTION]
           hf_pitch=int1byf(rate[0, *], 10, xover_freq_thr)
           hf_roll=int1byf(rate[1, *], 10, xover_freq_thr)
           hf_yaw=int1byf(rate[2, *], 10, xover_freq_thr)

           ;; Use the accelerometer data to calculate low frequency angle
           ;; information, then add that to the high frequency angles Low
           ;; pass filter cuts off frequencies abve the cutoff period (in
           ;; this case 20s or 0.05Hz).
           lf_pitch=lowpass_filter(asin(accel[0, *] / g), 10, $
                                   xover_freq_thr)
           lf_roll=lowpass_filter(asin(accel[1, *] / g), 10, $
                                  xover_freq_thr)
           pitch=detrend(hf_pitch) - lf_pitch
           roll=detrend(hf_roll) + lf_roll

           ;; Here, we are getting the low frequency part of the compass
           ;; heading out.  Convert compass heading of ship gyro to radians
           ;; and redefine coordinate system as: north=0, east=PI/2,
           ;; south=PI, west=-PI/2
           heading_rad=atan(sin(heading * !DTOR), cos(heading * !DTOR))
    
           ;; Demean low pass filter that
           gyro_sf=flux_times_dims[1] / double(ec_period) ; sampling freq
           head_dmean=heading_rad - mean(heading_rad, /NAN)
           flux_lf_yaw=lowpass_filter([head_dmean], gyro_sf, $
                                      xover_freq_thr)
  
           ;; Add the lowfreq and highfreq components --> BUT multiply
           ;; lf_yaw by -1 to convert to L.H. system
           yaw=reform(-flux_lf_yaw + hf_yaw)
           angle=transpose([[pitch], [roll], [yaw]])
  
           ;; Call motion correction routine
           u_true=motcorr(wind, accel, angle, motpak_offset, $
                          double(isample_rate) / 10, lfreq_thr, $
                          hfreq_thr, g)
           wind=u_true

        ENDIF

        ;; Low frequency motion correction

        ;; [SPL: WATCH LOW_FREQ_CORR FUNCTION, which redundantly uses
        ;; truewind.pro code, and also works with a loop that is no longer
        ;; needed, since the new TRUEWIND can process the entire array.]
        IF diag.sog[fperiod] GT 0.25 THEN BEGIN
           U_TRUE_2=low_freq_corr(wind[0, *], wind[1, *], cog, $
                                  sog / 1.9438449, heading, 337.0)
           wind[0, *]=U_TRUE_2[0, *]
           wind[1, *]=U_TRUE_2[1, *]
           true_son=bearing_avg(U_TRUE_2[3, *], U_TRUE_2[2, *])
           true_sonic_vel=true_son[0, 1]
           true_sonic_dir=true_son[0, 0]
           raw_son=bearing_avg(U_TRUE_2[5, *], U_TRUE_2[4, *])
           raw_sonic_vel=raw_son[0, 1]
           raw_sonic_dir=raw_son[0, 0]
        ENDIF

        ;; Output motion corrected data
        omot_corr=create_struct(field_names[time_flds[0]], $
                                reform(flux_times[0, *]))
        FOREACH tfld, time_flds[1:*] DO BEGIN ; include all time data
           omot_corr=create_struct(omot_corr, field_names[tfld], $
                                   reform(flux_times[tfld, *]))
        ENDFOREACH
        omot_corr=create_struct(omot_corr, $
                                'wind_speed_u', $
                                reform(wind_raw[0, *]), $
                                'wind_speed_v', $
                                reform(wind_raw[1, *]), $
                                'wind_speed_w', $
                                reform(wind_raw[2, *]), $
                                'wind_speed_u_corr', $
                                reform(wind[0, *]), $
                                'wind_speed_v_corr', $
                                reform(wind[1, *]), $
                                'wind_speed_w_corr', $
                                reform(wind[2, *]), $
                                'sonic_temperature', $
                                sonic_temperature, $
                                'co2_op', $
                                flux.co2_op, $
                                'h2o_op', $
                                flux.h2o_op, $
                                'pressure_op', $
                                flux.pressure_op, $
                                'diag_op', $
                                flux.diag_op, $
                                'accel_x', $
                                reform(accel_raw[0, *]), $
                                'accel_y', $
                                reform(accel_raw[1, *]), $
                                'accel_z', $
                                reform(accel_raw[2, *]), $
                                'accel_x_corr', $
                                reform(accel[0, *]), $
                                'accel_y_corr', $
                                reform(accel[1, *]), $
                                'accel_z_corr', $
                                reform(accel[2, *]), $
                                'rate_phi', $
                                reform(rate_raw[0, *]), $
                                'rate_theta', $
                                reform(rate_raw[1, *]), $
                                'rate_shi', $
                                reform(rate_raw[2, *]), $
                                'rate_phi_corr', $
                                reform(rate[0, *]), $
                                'rate_theta_corr', $
                                reform(rate[1, *]), $
                                'rate_shi_corr', $
                                reform(rate[2, *]), $
                                'co2_cl', $
                                flux.co2_cl, $
                                'h2o_cl', $
                                flux.h2o_cl, $
                                'pressure_cl', $
                                flux.pressure_cl, $
                                'temperature_cl', $
                                flux.temperature_cl, $
                                'latitude', $
                                latitude, $
                                'longitude', $
                                longitude, $
                                'sog', $
                                sog, $
                                'cog', $
                                cog, $
                                'heading', $
                                heading)

        mc_ofile_name=strcompress(mot_corr_odir + path_sep() + $
                                  iname_prefix + '_' + 'mc.' + $
                                  iname[1], /remove_all)
        mc_ofile_stamp=file_basename(mc_ofile_name)
        mc_out_list=file_search(mot_corr_odir + path_sep() + '*.' + $
                                iname[1], /nosort, /fold_case, /test_regular)
        matchfiles=where(mc_ofile_stamp EQ file_basename(mc_out_list), $
                         matchfilecount)
        IF matchfilecount GT 0 THEN BEGIN
           msg1='Motion corrected file '
           IF keyword_set(overwrite) THEN BEGIN
              message, msg1 + mc_ofile_stamp + $
                       ' already exists.  Overwriting', /informational
              write_csv, mc_ofile_name, omot_corr, $
                         header=strlowcase(tag_names(omot_corr))
           ENDIF ELSE BEGIN
              message, msg1 + mc_ofile_stamp + $
                       ' already exists.  Not overwriting', /informational
           ENDELSE
        ENDIF ELSE BEGIN
           write_csv, mc_ofile_name, omot_corr, $
                      header=strlowcase(tag_names(omot_corr))
        ENDELSE
        delvar, omot_corr

        ;; Eddy covariance calculations
        sf_hz=float(isample_rate) * 10 ; sampling freq (Hz)
        mom=ec_momentum(wind, sonic_temperature, $
                        diag.air_temperature[fperiod], $
                        diag.rh_percent[fperiod], $
                        diag.pressure[fperiod], $
                        float(ec_period) / 60, $ ; in minutes
                        sf_hz, CORR_MASSMAN=[isample_rate, 0.145])

     ENDFOREACH

  ENDFOR

END



;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; flux.pro ends here
