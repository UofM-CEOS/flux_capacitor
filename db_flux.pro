;; Author: Sebastian Luque
;; Created: 2013-11-12T17:07:28+0000
;; Last-Updated: 2016-02-23T21:39:06+0000
;;           By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     FLUX
;; 
;; PURPOSE:
;; 
;;     This procedure controls all flux processing.
;; 
;; CALLING SEQUENCE:
;; 
;;     DB_FLUX, Idir, Itemplate_Sav, Time_Idx, Isample_Rate, Ec_Period, $
;;              Motpak_Offset, SOG_Thr, Lfreq_Thr, Hfreq_Thr, $
;;              Xover_freq_Thr, Mot_Corr_Odir, Footprint_Odir, Ofile
;; 
;; INPUTS:
;; 
;;     Idir:                 Input directory (no trailing separator), with
;;                           flux files.
;;     Itemplate_Sav:        Ascii template to read Idir files.
;;     Time_Idx:             Index (in template) where time in Idir files.
;;     Isample_Rate:         Scalar indicating the frequency (s) with which
;;                           input Idir data were sampled.
;;     Ec_Period:            Scalar indicating the duration (s) of flux
;;                           study periods.
;;     Motpak_Offset:        Vector showing displacement between motionPak
;;                           and anemometer (in m), x, y, z.
;;     SOG_Thr:              Threshold SOG (knots) below which we do not do
;;                           motion correction (i.e. we assume ship is
;;                           drifting due to ice, etc.)
;;     Lfreq_Thr:            Low frequency threshold (Hz) for high pass
;;                           filter operation during accelerometer
;;                           integration.
;;     Hfreq_Thr:            High frequency threshold (Hz) for low pass
;;                           filter operation on derived
;;                           angles/translational velocities (filtering
;;                           motionpak noise.)
;;     Xover_freq_Thr:       Crossover period (Hz) which serves as a
;;                           threshold frequency for low pass filtering the
;;                           linear accelerations and for high pass
;;                           filtering the angular rates.
;;     Mot_Corr_Odir:        Output directory for motion-corrected files.
;;     Ofile:                Output file path for fluxes.
;;     Footprint_Odir:       Output directory for footprint plots.
;; 
;; KEYWORD PARAMETERS:
;; 
;;     OVERWRITE:            Whether to overwrite files in output
;;                           directories.
;; 
;; SIDE EFFECTS:
;; 
;;     Writes a file Ofile with the eddy covariance calculations.  Writes
;;     motion-corrected files in Mot_Corr_Odir, and footprint plots in
;;     Footprint_Odir.
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

PRO DB_FLUX, IDIR, ITEMPLATE_SAV, TIME_IDX, ISAMPLE_RATE, $
             EC_PERIOD, MOTPAK_OFFSET, SOG_THR, LFREQ_THR, HFREQ_THR, $
             XOVER_FREQ_THR, MOT_CORR_ODIR, FOOTPRINT_ODIR, $
             OFILE, SERIAL=SERIAL, OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 13) THEN $
     message, 'Usage: FLUX, IDIR, ITEMPLATE_SAV, TIME_IDX, ' + $
              'ISAMPLE_RATE, EC_PERIOD, MOTPAK_OFFSET, ' + $
              'SOG_THR, LFREQ_THR, HFREQ_THR, XOVER_FREQ_THR, ' + $
              'MOT_CORR_ODIR, FOOTPRINT_ODIR, OFILE'
  idir_info=file_info(idir)
  itpl_info=file_info(itemplate_sav)
  IF (~idir_info.directory) THEN $
     message, 'IDIR must be a string pointing to an existing directory'
  IF (~itpl_info.read) THEN $
     message, 'ITEMPLATE_SAV must be a string pointing to a readable file'
  IF ((n_elements(time_idx) NE 1) OR $
      ((size(time_idx, /type) NE 2) || time_idx LT 0)) THEN $
         message, 'TIME_IDX must be an integer scalar >= zero'
  IF ((n_elements(isample_rate) NE 1) OR (isample_rate LT 0)) THEN $
     message, 'ISAMPLE_RATE must be a scalar >= zero'
  sf_hz=float(isample_rate) * 100 ; sampling frequency (Hz)
  IF ((n_elements(ec_period) NE 1) OR (ec_period LT 0)) THEN $
     message, 'EC_PERIOD must be a scalar >= zero'
  IF (n_elements(motpak_offset) NE 3) THEN $
         message, 'MOTPAK_OFFSET must be a 3-element numerical array'
  IF ((n_elements(sog_thr) NE 1) OR (sog_thr LT 0)) THEN $
     message, 'SOG_THR must be a scalar >= zero'
  IF ((n_elements(lfreq_thr) NE 1) OR (lfreq_thr LT 0)) THEN $
     message, 'LFREQ_THR must be a scalar >= zero'
  IF ((n_elements(hfreq_thr) NE 1) OR (hfreq_thr LT 0)) THEN $
     message, 'HFREQ_THR must be a scalar >= zero'
  IF ((n_elements(xover_freq_thr) NE 1) OR (xover_freq_thr LT 0)) THEN $
     message, 'XOVER_FREQ_THR must be a scalar >= zero'
  IF ((n_elements(mot_corr_odir) NE 1) OR $
      (size(mot_corr_odir, /type) NE 7)) THEN $
         message, 'MOT_CORR_ODIR must be a string scalar'
  IF ((n_elements(footprint_odir) NE 1) OR $
      (size(footprint_odir, /type) NE 7)) THEN $
         message, 'FOOTPRINT_ODIR must be a string scalar'
  IF ((n_elements(ofile) NE 1) OR (size(ofile, /type) NE 7)) THEN $
         message, 'OFILE must be a string scalar'
  idir_files=file_search(idir + path_sep() + '*', count=nidir_files, $
                         /nosort, /fold_case, /test_regular)
  IF (nidir_files EQ 0) THEN $
        message, 'No files in at least one of IDIR or DIAG_DIR, '

  ;; Standard time names (output from read_iso_file())
  stdtime_names=['year', 'month', 'day', 'hour', 'minute', 'second']
  ;; Parse input flux files
  restore, itemplate_sav
  field_names=strlowcase(itemplate.FIELDNAMES)
  time_flds=where(itemplate.FIELDGROUPS EQ time_idx)
  ;; Break file names and extract the piece to match
  ;; [Original comment: Temporarily set g for filtering purposes... true g
  ;; will be calculated later... we'll set this low to make sure we filter
  ;; properly].  CHECK. [SPL: I've renamed this to avoid confusion.]
  g_thr=8.0

  ;; Set up a hash (or structure) to hold output data from all valid flux
  ;; runs at this point, before starting outermost loop.
  ;; BE VERY CAREFUL THESE EXIST IN THE TEMPLATE!  I'm only checking that
  ;; "pressure" gets interpreted as "atmospheric_pressure", and
  ;; "rh_percent" as "relative_humidity" below. Building names in steps, as
  ;; we will need to use each bit at the end.
  okeys_diag=['latitude', 'longitude', 'speed_over_ground', 'wind_speed', $
              'wind_direction', 'true_wind_speed', $
              'true_wind_direction', 'air_temperature', $
              'relative_humidity', 'surface_temperature', $
              'atmospheric_pressure']
  okeys_mom=['cf_wu', 'Ustar', 'Tau', 'MO_L', 'cov_w_psAirT', 'psH', $
             'cov_w_u']
  okeys_op=['CO2_op', 'H2O_op', 'cov_w_airT', 'cf_wAirT', 'H', $
            'cov_w_CO2_op', 'cf_CO2_op', 'FCO2_op', 'WPLcont', $
            'CO2Burba_mt', 'CO2Burba_ln', 'H2OBurba_mt', 'H2OBurba_ln', $
            'cov_w_H2O_op', 'cf_wH2O_op', 'E_op', 'Qe_op', 'lag_op', $
            'pkt_FCO2_op', 'pkt_loop', 'cflux_det', 'dRH_by_dq', $
            'dc_by_dRH_initial', 'dc_by_dRH_final', 'dc_by_dq_initial', $
            'dc_by_dq_final']
  okeys_calc=['sonic_speed', 'sonic_direction', 'true_sonic_speed', $
                'true_sonic_direction', 'vertical', 'open_flag', $
                'sonic_flag', 'motion_flag', 'sonic_NAN_pct', $
                'SW', 'CDm']
  okeys_micro=['Um', 'U10', 'CD10', 'z0', 'CHm', 'CH10', 'zT', 'CEm', $
               'CE10', 'zQ', 'peakF', 'dist90', 'psim', 'psih', 'U10N', $
               'U10Nocean']
  okeys=['time', 'DOY', okeys_diag, okeys_mom, okeys_op, 'diag_op', $
         okeys_calc, okeys_micro]
  fluxes=hash(okeys)            ; just empty keys; we'll be appending data

  FOREACH ifile, idir_files DO BEGIN

     ;; Read flux file
     message, 'Processing ' + ifile, /informational
     flux=read_iso_file(ifile, itemplate, time_idx)
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
     ;; Temporary approach (GENERALIZE THIS LATER)
     diag_jd=reform(julday(long(strmid(flux.time_20min[0], 5, 2)), $
                             long(strmid(flux.time_20min[0], 8, 2)), $
                             long(strmid(flux.time_20min[0], 0, 4)), $
                             long(strmid(flux.time_20min[0], 11, 2)), $
                             long(strmid(flux.time_20min[0], 14, 2)), $
                             double(strmid(flux.time_20min[0], 17, 1))))
     diag_tstamp=jul2timestamp(diag_jd)
     diag_doy=calendar2doy(long(strmid(flux.time_20min[0], 0, 4)), $
                             long(strmid(flux.time_20min[0], 5, 2)), $
                             long(strmid(flux.time_20min[0], 8, 2)))
     ;; Get a file name prefix to be shared by the output files from
     ;; this period
     iname=strsplit(file_basename(ifile), '.', /extract)
     iname_prefix=iname[0]

     ;; [Original comment: no need recalibrate motion sensor for this
     ;; experiment.  Read in accelerations in RH coordinate system,
     ;; convert to m/s2].  All these values could just be placed back on
     ;; the intput structure, and avoid cluttering the workspace and
     ;; memory so much
     accel_x=flux.acceleration_z * 9.81     ; accel_z on the tower
     accel_y=-flux.acceleration_x * 9.81    ; accel_x on the tower
     accel_z=-flux.acceleration_y * 9.81    ; accel_y on the tower
     ;; Put acceleration components in 3-column array and make copy to
     ;; keep uncorrected data
     accel=transpose([[accel_x], [accel_y], [accel_z]])
     accel_raw=accel
     ;; Original comment: read in angular rates in RH coordinate system,
     ;; convert to rad/s
     rate_phi=flux.rate_z * !DTOR       ; rate_z on the tower
     rate_theta=-flux.rate_x * !DTOR    ; rate_x on the tower
     rate_shi=-flux.rate_y * !DTOR      ; rate_y on the tower
     ;; Put rate components in 3-column array and make copy to keep
     ;; uncorrected data
     rate=transpose([[rate_phi], [rate_theta], [rate_shi]])
     rate_raw=rate
     ;; Extract the wind components based on whether we want serial or
     ;; analogue data. KEEPING THEM THE SAME FOR NOW.
     wind_u=keyword_set(serial) ? $
            flux.wind_speed_u : $
            flux.wind_speed_u
     wind_v=keyword_set(serial) ? $
            flux.wind_speed_v : $
            flux.wind_speed_v
     wind_w=keyword_set(serial) ? $
            flux.wind_speed_w : $
            flux.wind_speed_w
     ;; Put wind components together in 3-column array, and make a copy
     ;; to keep uncorrected data.  KEEPING SERIAL AND ANALOGUE THE SAME
     ;; FOR NOW, until we determine a stable workflow.
     wind=transpose([[wind_u], [wind_v], [wind_w]])
     wind_raw=wind
     sonic_temperature=keyword_set(serial) ? $
                       flux.air_temperature_sonic : $
                       flux.air_temperature_sonic

     ;; [Original comment: create flags for the 4 possible sources of
     ;; "bad" data, flag=0 means data good]
     open_flag=0
     sonic_flag=0
     motion_flag=0

     ;; [Original comment: check for any significant number of 'NAN's
     ;; (not worried about the odd one scattered here and there)]
     bad_co2_op=where(~finite(flux.op_co2_density), nbad_co2_op)
     bad_h2o_op=where(~finite(flux.op_h2o_density), nbad_h2o_op)
     ;; bad_diag_op=where(~finite(fix(flux.diag_op)), nbad_diag_op)
     ;; ;; [Original comment: set open flag if gt 2% of records are 'NAN']
     ;; IF (float(nbad_co2_op) / flux_times_dims[1] * 100.0 GT 2) OR $
     ;;    (float(nbad_h2o_op) / flux_times_dims[1] * 100.0 GT 2) OR $
     ;;    (float(nbad_diag_op) / flux_times_dims[1] * 100.0 GT 2) THEN $
     ;;       open_flag=1
     ;; [Original comment: set open flag if gt 2% of records are 'NAN']
     IF (float(nbad_co2_op) / flux_times_dims[1] * 100.0 GT 2) OR $
        (float(nbad_h2o_op) / flux_times_dims[1] * 100.0 GT 2) THEN $
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

     ;; [Original comment: now that we have looked for NANs, we may as
     ;; well fill in the NANs and any spikes using the shot filter].
     ;; [SPL: these changes are done outside the WIND array, which is
     ;; the one that is used later for motion correction, etc., so they
     ;; are lost.]
     IF sonic_flag NE 1 THEN BEGIN
        wind_u=shot_filter(wind_u)
        wind_v=shot_filter(wind_v)
        wind_w=shot_filter(wind_w)
        sonic_temperature=shot_filter(sonic_temperature)
     ENDIF
     IF open_flag NE 1 THEN BEGIN
        flux.op_co2_density=shot_filter(flux.op_co2_density)
        flux.op_h2o_density=shot_filter(flux.op_h2o_density)
        flux.op_pressure=shot_filter(flux.op_pressure)
        ;; [Original comment: this is necessary to check if there is
        ;; still ugly shot noise... if there is, we need to skip this]
        bad_co2_op=where(abs(flux.op_co2_density - $
                               mean(flux.op_co2_density, /NAN)) GT $
                           (6.0 * stddev(flux.op_co2_density, /NAN)), $
                           nbad_co2_op)
        bad_h2o_op=where(abs(flux.op_h2o_density - $
                               mean(flux.op_h2o_density, /NAN)) GT $
                           (6.0 * stddev(flux.op_h2o_density, /NAN)), $
                           nbad_h2o_op)
        IF (nbad_co2_op GT 0) OR (nbad_h2o_op GT 0) THEN open_flag=1
     ENDIF

     ;;       ;; Check for high or low diagnostic flags (set open flag if more
     ;;       ;; than 2% of records have them)
     ;;       bad_diags=where((flux.diag_op GT 249) OR (flux.diag_op LT 240), $
     ;;                       nbad_diags)
     ;;       IF (float(nbad_diags) / flux_times_dims[1] * 100.0 GT 2) THEN $
     ;;          open_flag=1

     ;; [Original comment: check for bad wind data: bad wind data can
     ;; usually be diagnosed by unusually high wind speeds.  this is
     ;; most obvious in the vertical wind where we wouldn't expect high
     ;; values bad sonic data can also turn up in the Tsonic before the
     ;; wind, check the deviation between Tsonic and mean air T]
     bad_vert_wind=where(abs(wind_w) GT 7, nbad_vert_wind)
     t_avg=flux.air_temperature[0]
     bad_sonic_temperature=where(abs(sonic_temperature - t_avg) GT 7, $
                                   nbad_sonic_temperature)
     ;; ;; sonic_count appears to not be working... I temporarily disabled
     ;; ;; this check... RS Feb 2013
     ;; sonic_count=0  
     ;; Set wind flag high if gt 0.5% of records are frost contaminated
     IF (float(nbad_vert_wind) / flux_times_dims[1] * 100.0 GT 0.5) OR $
        (float(nbad_sonic_temperature) / $
           flux_times_dims[1] * 100.0 GT 0.5) THEN sonic_flag=1

     ;; Check for Motion Pak data that are out of range
     bad_motionpak=where((accel_x GT g_thr) OR (accel_y GT g_thr), $
                           nbad_motionpak)
     IF nbad_motionpak GT 0 THEN motion_flag=1

     ;; if sonic data are bad, skip this period
     IF sonic_flag GT 0 THEN BEGIN
        message, 'Bad sonic anemometer data. Skipping.', $
                 /CONTINUE
        CONTINUE
     ENDIF

     ;; [Original comment: check critical low f variabiles]
     IF (~finite(flux.relative_humidity[0])) OR $
        (~finite(t_avg)) THEN BEGIN
        message, 'RH or mean air temperature unavailable. ' + $
                 'Skipping.', /CONTINUE
        CONTINUE
     ENDIF

     ;; Mean radiation
     sw_avg=flux.k_down[0]
     lw_avg=flux.lw_in[0]

     ;; Get RMC
     latitude=flux.latitude
     longitude=flux.longitude
     sog=flux.speed_over_ground
     cog=flux.course_over_ground
     sog_avg=mean(sog, /nan)    ; for comparing against threshold later
     ;; Get GYRO
     heading=flux.heading

     ;; [Original comment: now fill in the gaps by applying a moving
     ;; average... In this case, we use a 100 sample window (10 sec)
     ;; moving average... may need to tweak this value].  [SPL: I THINK
     ;; THIS STEP SHOULD HAVE BEEN DONE EARLIER DURING RMC PROCESSING.
     ;; Also, perhaps a simple linear interpolation is better; I don't
     ;; know why this moving average is used, where a window must be
     ;; specified and may be introducing bias.  Why aren't latitude and
     ;; longitude not similarly interpolated?]
     cog_xy=smooth_angle(cog, sog, 100)
     cog=reform(cog_xy[0, *])
     sog=reform(cog_xy[1, *])
     ;; Dummy heading magnitude for decomposition purposes only
     vheading=float(rebin([1], flux_times_dims[1], /sample))
     ;; Guard if missing heading
     noheading=where(~finite(heading), nnoheading)
     IF nnoheading GT 0 THEN vheading[noheading]=!VALUES.F_NAN
     heading_xy=smooth_angle(heading, vheading, 100)
     heading=reform(heading_xy[0, *])

     ;; [Original comment: Check to make sure that no 'NaNs' dropped
     ;; through... if they did, we'll have to skip this one]
     no_cog=where(~finite(cog), nno_cog, ncomplement=nok_cog)
     no_sog=where(~finite(sog), nno_sog, ncomplement=nok_sog)
     no_heading=where(~finite(heading), nno_heading, $
                        ncomplement=nok_heading)
     IF (nno_cog GT 0) OR (nno_sog GT 0) OR (nno_heading GT 0) THEN $
        motion_flag=1
     ;; If we have no good COG, SOG, or heading, then we should skip
     ;; processing entirely
     IF (nok_cog LT 1) OR (nok_sog LT 1) OR $
        (nok_heading LT 1) THEN BEGIN
        message, 'Unusable COG, SOG, or heading records. Skipping.', $
                 /CONTINUE
        CONTINUE
     ENDIF

     IF sog_avg GT sog_thr THEN BEGIN
        ;; [Original comment: shot filter the motion channels... this helps
        ;; with a problem where unreasonably high accelerations cause a
        ;; 'NaN' calculation]
        accel[0, *]=shot_filter(accel[0, *])
        accel[1, *]=shot_filter(accel[1, *])
        accel[2, *]=shot_filter(accel[2, *])
        rate[0, *]=shot_filter(rate[0, *])
        rate[1, *]=shot_filter(rate[1, *])
        rate[2, *]=shot_filter(rate[2, *])

        ;; [Original comment: synthetically level the Motion Pak] [SPL:
        ;; WATCH LEVEL_MOTIONPAK FUNCTION.  Note how this calls the
        ;; function with the fixed value of 4 for R_RANGE and P_RANGE
        ;; arguments.  Why?]
        level=level_motionpak(accel, rate, 4, 4, status=status)
        IF status NE 0 THEN BEGIN
           message, 'Levelling Motion Pak data failed. Skipping.', $
                    /CONTINUE
           CONTINUE
        ENDIF
        ;; What is the condition below trying to test?  BE had this to
        ;; catch if something went wrong during LEVEL_MOTIONPAK, beyond the
        ;; status variable above.  I guess this simply flags the data, as
        ;; opposed to simply rejecting the entire period.
        IF ~finite(level[0]) THEN BEGIN
           motion_flag=1
           message, 'Motion-flagged', /informational
        ENDIF
        accel_lev=level[0:2, *]
        rate_lev=level[3:5, *]
     
        ;; [Original comment: demean the rates by extracting the
        ;; fluctuating component and setting the mean to 0]
        rate_lev[0, *]=(rate_lev[0, *] - $
                          mean(rate_lev[0, *], /NAN, /DOUBLE)) + 0
        rate_lev[1, *]=(rate_lev[1, *] - $
                          mean(rate_lev[1, *], /NAN, /DOUBLE)) + 0
        rate_lev[2, *]=(rate_lev[2, *] - $
                          mean(rate_lev[2, *], /NAN, /DOUBLE)) + 0
     
        ;; Gravity calculation
        g=(mean(accel_lev[2, *], /NAN, /DOUBLE))

        ;; Check to make sure the x/y accelerations are not greater than
        ;; g... this is indicative of a problem w/ the Motion Pak, and the
        ;; run needs to be skipped
        bad_mp=where((accel_lev[0, *] GT g) OR $
                       (accel_lev[1, *] GT g), nbad_mp)
        IF nbad_mp GT 0 THEN BEGIN
           motion_flag=1
           message, 'Invalid Motion Pak accelerations. Skipping.', $
                    /CONTINUE
           CONTINUE
        ENDIF

        ;; Integrate angular rates to give angles.  This also performs a
        ;; high pass filter, cutting off all frequencies below the cutoff
        ;; period (in this case, the cutoff period is set to 20s (or
        ;; 0.05Hz)). [SPL: WATCH INT1BYF FUNCTION]
        hf_pitch=int1byf(rate_lev[0, *], sf_hz, xover_freq_thr)
        hf_roll=int1byf(rate_lev[1, *], sf_hz, xover_freq_thr)
        hf_yaw=int1byf(rate_lev[2, *], sf_hz, xover_freq_thr)

        ;; Use the accelerometer data to calculate low frequency angle
        ;; information, then add that to the high frequency angles Low pass
        ;; filter cuts off frequencies abve the cutoff period (in this case
        ;; 20s or 0.05Hz).
        lf_pitch=-lowpass_filter(reform(asin(accel_lev[0, *] / g)), sf_hz, $
                                   xover_freq_thr)
        lf_roll=lowpass_filter(reform(asin(accel_lev[1, *] / g)), sf_hz, $
                                 xover_freq_thr)
        pitch=detrend(hf_pitch) + lf_pitch
        roll=detrend(hf_roll) + lf_roll
        pitch_avg=(bearing_avg(pitch, 1))[0]
        roll_avg=(bearing_avg(roll, 1))[0]

        ;; Here, we are getting the low frequency part of the compass
        ;; heading out.  Convert compass heading of ship gyro to radians
        ;; and redefine coordinate system as: north=0, east=PI/2, south=PI,
        ;; west=-PI/2
        heading_rad=atan(sin(heading * !DTOR), cos(heading * !DTOR))
        
        ;; Demean and low pass filter
        gyro_sf=flux_times_dims[1] / double(ec_period) ; sampling freq
        head_dmean=heading_rad - mean(heading_rad, /NAN)
        flux_lf_yaw=-lowpass_filter(reform([head_dmean]), gyro_sf, $
                                      xover_freq_thr)
        ;; We need to stop here if low pass filter failed
        lpf_ok=where(finite(flux_lf_yaw), nlpf_ok)
        IF nlpf_ok EQ 0 THEN BEGIN
           message, 'Low pass filter for heading failed. Skipping.', $
                    /CONTINUE
           CONTINUE
        ENDIF
        ;; Add the lowfreq and highfreq components --> BUT multiply lf_yaw
        ;; by -1 to convert to L.H. system
        yaw=reform(hf_yaw + flux_lf_yaw)
        angle=transpose([[pitch], [roll], [yaw]])

        ;; Level sonic anemometer

        ;; [Original comment: here we level the sonic anemometer, as long
        ;; as we know the MEAN roll/pitch angles from an inclinometer, and
        ;; as long as those are in a L.H.S.] [SPL: WATCH LEVEL_SONIC
        ;; FUNCTION]
        IF finite(roll_avg) AND finite(pitch_avg) THEN $
              wind_lev=level_sonic(wind, roll_avg * !DTOR, $
                                     -pitch_avg * !DTOR)
        ;; [Original comment: now let's calculate the raw mean w
        ;; value... this is going to be important for sort of tracking flow
        ;; distortion... We could do something more in depth, but we'll
        ;; keep it like this for now]
        wind_v_mean=mean(wind_lev[2, *], /nan)

        ;; Call motion correction routine
        wind_corr=motcorr(wind_lev, accel_lev, angle, motpak_offset, $
                            sf_hz, lfreq_thr, hfreq_thr, g)
        ;; ;; Attempting to do it with Miller's code, ported to IDL. Debugging
        ;; ;; (these variables should be moved to the control file)
        ;; butter_num=transpose([0.9693, -3.8772, 5.8157, -3.8772, 0.9693])
        ;; butter_den=transpose([1, -3.9376, 5.8148, -3.8167, 0.9395])
        ;; butter_coeffs=[butter_num, butter_den]
        ;; cutoffs=float([80, 80, 40])
        ;; motion, wind_lev, accel_lev, rate_lev, heading, sog, $
        ;;         motpak_offset, sf_hz, 20, cutoffs, butter_coeffs, $
        ;;         [0, 0], [0, 0]
     
        ;; Low frequency motion correction
     
        ;; [SPL: Watch TRUEWIND_SONIC function, which redundantly uses
        ;; truewind.pro code.  See note above regarding these variables
        ;; needed later; what to do if the test in this block (e.g. SOG is
        ;; too low) fails?]
        u_true=truewind_sonic(wind_corr[0, *], wind_corr[1, *], cog, $
                                sog / 1.9438449, heading, 337.0)
        wind_corr[0, *]=u_true[0, *]
        wind_corr[1, *]=u_true[1, *]
        true_son=bearing_avg(u_true[3, *], u_true[2, *])
        true_sonic_spd=true_son[0, 1]
        true_sonic_dir=true_son[0, 0]
        raw_son=bearing_avg(u_true[5, *], u_true[4, *])
        raw_sonic_spd=raw_son[0, 1]
        raw_sonic_dir=raw_son[0, 0]

     ENDIF ELSE BEGIN
        wind_lev=wind
        wind_corr=wind
     ENDELSE

     ;; Output motion corrected (if it was needed) data
     omc_tags=[stdtime_names, $
                 'wind_speed_u', $
                 'wind_speed_v', $
                 'wind_speed_w', $
                 'wind_speed_lev_u', $
                 'wind_speed_lev_v', $
                 'wind_speed_lev_w', $
                 'wind_speed_u_corr', $
                 'wind_speed_v_corr', $
                 'wind_speed_w_corr', $
                 'sonic_temperature', $
                 'CO2_op', $
                 'H2O_op', $
                 'pressure_op', $
                 'diag_op', $
                 'accel_x', $
                 'accel_y', $
                 'accel_z', $
                 'accel_x_lev', $
                 'accel_y_lev', $
                 'accel_z_lev', $
                 'rate_phi', $
                 'rate_theta', $
                 'rate_shi', $
                 'rate_phi_lev', $
                 'rate_theta_lev', $
                 'rate_shi_lev', $
                 'latitude', $
                 'longitude', $
                 'SOG', $
                 'COG', $
                 'heading']
     omot_corr=create_struct(stdtime_names[0], reform(flux_times[0, *]))
     FOREACH tfld, stdtime_names[1:*] DO BEGIN ; include all time data
        tfld_idx=where(stdtime_names EQ tfld)
        omot_corr=create_struct(omot_corr, tfld, $
                                  reform(flux_times[tfld_idx, *]))
     ENDFOREACH
     omot_corr=create_struct(omot_corr, $
                               omc_tags[n_elements(stdtime_names):*], $
                               reform(wind_raw[0, *]), $
                               reform(wind_raw[1, *]), $
                               reform(wind_raw[2, *]), $
                               reform(wind_lev[0, *]), $
                               reform(wind_lev[1, *]), $
                               reform(wind_lev[2, *]), $
                               reform(wind_corr[0, *]), $
                               reform(wind_corr[1, *]), $
                               reform(wind_corr[2, *]), $
                               sonic_temperature, $
                               flux.op_co2_density, $
                               flux.op_h2o_density, $
                               flux.op_pressure, $
                               flux.op_analyzer_status, $
                               reform(accel_raw[0, *]), $
                               reform(accel_raw[1, *]), $
                               reform(accel_raw[2, *]), $
                               reform(accel_lev[0, *]), $
                               reform(accel_lev[1, *]), $
                               reform(accel_lev[2, *]), $
                               reform(rate_raw[0, *]), $
                               reform(rate_raw[1, *]), $
                               reform(rate_raw[2, *]), $
                               reform(rate_lev[0, *]), $
                               reform(rate_lev[1, *]), $
                               reform(rate_lev[2, *]), $
                               latitude, $
                               longitude, $
                               sog, $
                               cog, $
                               heading)
     file_mkdir, mot_corr_odir
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
           write_csv, mc_ofile_name, omot_corr, header=omc_tags
        ENDIF ELSE BEGIN
           message, msg1 + mc_ofile_stamp + $
                    ' already exists.  Not overwriting', /informational
        ENDELSE
     ENDIF ELSE BEGIN
        write_csv, mc_ofile_name, omot_corr, header=omc_tags
     ENDELSE
     delvar, omot_corr
  
     ;; Eddy covariance calculations
     mom=ec_momentum(wind_corr, sonic_temperature, $
                       flux.air_temperature[0], $
                       flux.relative_humidity[0], $
                       flux.atmospheric_pressure[0], $
                       float(ec_period) / 60, $ ; in minutes
                       sf_hz, CORR_MASSMAN=[isample_rate, 0.145])

     ;; Open path calculations.  IF TEST FAILS, THEN WE CREATE AN ARRAY
     ;; WITH THE SAME DIMENSIONS AS THAT RETURNED BY EC_OPEN, for
     ;; printing purposes.
     IF open_flag NE 1 THEN BEGIN
        open_path=ec_open(wind_corr, sonic_temperature, flux.op_co2_density, $
                            flux.op_h2o_density, flux.op_pressure, $
                            float(ec_period) / 60, $ ; in minutes
                            sf_hz, $
                            CORR_MASSMAN=[isample_rate, 0.145, -0.04, $
                                            0.38, !VALUES.D_NAN, 0.01, $
                                            0.125], $
                            BURBA=[sw_avg, lw_avg, raw_sonic_spd], $
                            pkt=mom[2])
     ENDIF ELSE open_path=make_array(27, value=!VALUES.D_NAN)

     ;; Closed-path flow rate for 2010 = 11.5 LPM.  Sample tube length
     ;; for 2010 = 8.0m.  Do closed path calculations if available.  IF
     ;; CLOSED FLAG IS UP, THEN WE CREATE AN ARRAY OF SAME DIMENSIONS AS
     ;; THAT RETURNED BY EC_CLOSED, for printing purposes.
     IF closed_flag EQ 0 THEN BEGIN
        cl_flow=11.5            ; NEED LUT FOR FLOW RATE
        closed_path=ec_closed(wind_corr, flux.co2_cl, flux.h2o_cl, $
                                flux.pressure_cl, flux.temperature_cl, $
                                flux.air_temperature[0], $
                                flux.relative_humidity[0], $
                                flux.atmospheric_pressure[0], [-10, 50], $
                                float(ec_period) / 60, $ ; in minutes
                                sf_hz, $
                                CORR_MASSMAN=[isample_rate, 0.10, -0.06, $
                                                0.44, !VALUES.D_NAN, $
                                                !VALUES.D_NAN, $
                                                !VALUES.D_NAN, cl_flow, $
                                                0.005, 8.0])
     ENDIF ELSE closed_path=make_array(11, value=!VALUES.D_NAN)

     ;; MICROMET CALCULATIONS

     ;; Here, we make a number of micromet calculations that calculate
     ;; the following things:
     ;;  - 10m wind speed
     ;;  - zo
     ;;  - footprint modeling
     ;;  - CD_10m, CH_10m, CV_10m
     
     ;; The following publications are heavily used, and I will try to
     ;; cite them and specific eq'ns: Andreas et al. 2005, BLM
     ;; 114:439-460 Jordan et al. 1999, JGR 104 No. C4:7785-7806
     
     ;; First, we need to extract the measurement height... we'll just
     ;; give it a value of 14.1m... come back to that
     zm=14.1 
     ;; couple of constants
     vonk=0.4                   ; von karman constant
     Rval=8.31451               ; j/mol/k universal gas constant
     mv=18.02                   ; g/mol molecular weight for water, Stull(1995)
     ma=28.96                   ; g/mol molecular weight for dry air, Stull(1995)
     ;; J/g/K specific heat for dry air at constant pressure Stull(1995)
     cpd=1004.67 / 1000.0
     
     ;; Now, we need to calculate the air density.  Sat'n vapour
     ;; pressure (Pa)
     es_met=[6.112 * exp((17.67 * flux.air_temperature[0]) / $
                           (flux.air_temperature[0] + 243.5))] * $
            100.0
     ;; Vapour pressure (Pa)
     ev_met=(flux.relative_humidity[0] / double(100)) * es_met
     ;;  mol 
     c_h2o_met=ev_met / (Rval * (flux.air_temperature[0] + 273.15))
     rho_v_met=c_h2o_met * mv   ; g/m3
     rho_d_met=((flux.atmospheric_pressure[0] * double(1000) - ev_met) / $
                  (Rval * (flux.air_temperature[0] + 273.15))) * ma
     rhoa=(rho_v_met + rho_d_met) / 1000.0 ; kg/m3

     ;; Calculate CD at measurement height.  CD at measurement height,
     ;; Andreas (2005), eqn, 10a)
     CDm=mom[3] / (rhoa * true_sonic_spd ^ double(2))
     
     ;; Drag coefficient should be positive... if not, do no more calcs
     IF CDm GT 0 THEN BEGIN
        ;; Calculate the wind profile modifiers.  Calculate the profile
        ;; modifier for stable conditions (Jordan, 1999, eq 33)
        IF zm / mom[4] GT 0.5 THEN BEGIN    
           psim=-[((0.70 * zm) / mom[4]) + 0.75 * $
                    ((zm / motpak_offset) - 14.3) * $
                    exp((-0.35 * zm) / mom[4] ) + 10.7]
           psih=psim
        ENDIF
        ;; Calculate the profile modifier for neutral/slightly stable
        ;; (Jordan, 1999, eq 32)
        IF zm / mom[4] GE 0 AND zm / mom[4] LE 0.5 THEN BEGIN
           psim=-6.0 * (zm / mom[4])
           psih=psim
        ENDIF
        ;; Calculate the profile modifier for unstable conditions
        ;; (Jordan, 1999, eq 30)
        IF zm / mom[4] LT 0 THEN BEGIN
           xval=(1.0 - 16.0 * (zm / mom[4])) ^ (1.0 / 4.0)
           psim=alog((1.0 + xval ^ 2.0) / 2.0) + 2.0 * $
                alog((1.0 + xval) / 2.0) - 2.0 * atan(xval) + $
                (!PI / 2.0)
           ;; Jordan (1999), eq 31
           psih=2.0 * alog((1.0 + xval ^ 2.0) / 2.0)
        ENDIF
        ;; Calculate the roughness length (Andreas, 2005, eq 12a)
        z0=zm * EXP( -[vonk * CDm ^ (-0.5) + psim * (zm / mom[4])])
        ;; Calculate the Drag coefficient at 10m (Andreas, 2005, eq 11a)
        CD10=vonk ^ 2.0 / $
             (alog(10.0 / z0) - psim * (10.0 / mom[4])) ^ 2.0
        ;; Calculate wind speed at 10m & for interest sake, we'll see
        ;; what our profile gives us for the measurement height wind
        ;; speed...  Stull (1988), eq'n 9.7.5g
        U10=[(1.0 / vonk) * $
               (alog(10.0 / z0) + psim * (10.0 / mom[4]))] * mom[2]
        Um=[(1.0 / vonk) * $
              (alog(zm / z0) + psim * (zm / mom[4]))] * mom[2]
        ;; Calculate the wind speed at 10m assuming neutral stability
        U10N=[(double(1) / vonk) * (alog(zm / z0))] * mom[2]
        ;; Calculate the wind speed at 10m assuming neutral stability
        ;; AND that we're over the ocean
        z0charnock=(0.016 * mom[2] ^ 2) / 9.81 ; Stull 9.7.2c
        U10Nocean=[(double(1) / vonk) * $
                     (alog(zm / z0charnock))] * mom[2]
        ;; ;; Make call to the footprint routine, which follows Hsieh et
        ;; ;; al. 2000, Advances in Water Resources 23 (2000) 765-777.
        ;; ;; [SPL: Turning this off; it is taking way too much, stuck on
        ;; ;; the while loop in hkt_footprint.pro.]
        ;; foot_ofile_name=strcompress(footprint_odir + path_sep() + $
        ;;                             iname_prefix + '_' + 'spec.ps', $
        ;;                             /remove_all)
        ;; IF finite(z0) EQ 1 THEN BEGIN
        ;;    foot=hkt_footprint(mom[4], z0, zm, 1.0, 100.0, 0.99, $
        ;;                       plot_file=foot_ofile_name)
        ;; ENDIF
        ;; peakF=foot[0]
        ;; dist90=foot[1]
        ;; Now we can also calculate the heat and water vapour transfer
        ;; coefficients...  BUT BEWARE... fundamentally, we need surface
        ;; T... we have to use our unreliable IR transducer for this...
        ;; NOTE - using a new type of IR transducer starting in 2010
        ;; which is more reliable
        ;;  
        ;; Calculate a few terms we'll need from the MET data.  Mean
        ;; specific humidity (g_h2o/g_moist air)
        mean_qh2o=rho_v_met / (rhoa * 1000.0)
        ;; Mean specific heat capacity of moist air (J/g/K)
        mean_cp=cpd * (1.0 + 0.84 * mean_qh2o)
        IF (flux.air_temperature[0] GE 0) THEN $ ; vaporization
           Lv=(2.50057 - 0.00245715 * $
                 flux.air_temperature[0]) * 1000.0
        IF (flux.air_temperature[0] LT 0) THEN $ ; sublimation
           lv=(2.83539 - 0.000135713 * $
                 flux.air_temperature[0]) * 1000.0
        ;; Sat'n vapour pressure (Pa)
        ev_surf=[6.112 * $
                   EXP((17.67 * flux.surface_temperature[0]) / $
                         (flux.surface_temperature[0] + 243.5))] * 100.0
        c_h2o_surf=ev_surf / $
                   (Rval * $
                      (flux.surface_temperature[0] + 273.15)) ; mol
        rho_v_surf=c_h2o_surf * mv                            ; g/m3
        surf_qh2o=rho_v_surf / (rhoa * 1000.0) 
        ;; Calculate the H coeficient and the T roughness length... from
        ;; that get CH at 10m. Andreas (2005) eq'n 10b
        CHm=open_path[2] / $
            ((rhoa * 1000.0) * mean_cp * true_sonic_spd * $
               (flux.surface_temperature[0] - $
                  flux.air_temperature[0]) )
        ;; Andreas (2005) eq'n 12b
        zT=zm * EXP( -(vonk * CDm ^ (1.0 / 2.0) * $
                         CHm ^ (-1.0) + psih * (zm / mom[4])))
        ;; Andreas (2005) eq'n 11b
        CH10=vonk ^ 2.0 / $
             [(alog(10.0 / z0) - psim * (10.0 / mom[4])) * $
                (alog(10.0 / zT) - psih * (10.0 / mom[4])) ]
        ;; Now calculate the E coefficient, and the Q roughness
        ;; length... from that get CE at 10m
        CEm=open_path[9] / [(rhoa * 1000.0) * Lv * raw_sonic_spd * $
                              (surf_qh2o - mean_qh2o)]
        zQ=zm * exp(-(vonk * CDm ^ (1.0 / 2.0) * CEm ^ (-1.0) + $
                        psih * (zm / mom[4]))) ;Andreas (2005) eq'n 12c
        ;; Andreas (2005) eq'n 11c
        CE10=vonk ^ 2.0 / $
             [(alog(10.0 / z0) - psim * (10.0 / mom[4])) * $
                (alog(10.0 / zQ) - psih * (10.0 / mom[4]))]
     ENDIF

     ;; Time to append data to output hash.  Start with MET summary
     ;; (diag) data.
     fluxes['time']=[fluxes['time'], diag_tstamp]
     fluxes['DOY']=[fluxes['DOY'], diag_doy]
     FOREACH diag_fld, okeys_diag DO BEGIN
        match_fld=where(strlowcase(tag_names(flux)) EQ $
                          strlowcase(diag_fld))
        fluxes[diag_fld]=[fluxes[diag_fld], $
                            reform(flux.(match_fld)[0])]
     ENDFOREACH
     ;; Next we have momentum output.  We need to find out which
     ;; elements from mom object correspond to field names in okeys_mom
     mom_idx=[1, 2, 3, 4, 6, 7, 0] ; same # elements as okeys_mom
     FOREACH mom_fld, indgen(n_elements(okeys_mom)) DO BEGIN
        fluxes[okeys_mom[mom_fld]]=[fluxes[okeys_mom[mom_fld]], $
                                      mom[mom_idx[mom_fld]]]
     ENDFOREACH
     ;; Next we have open path output.  We need to (be careful these
     ;; indices match) find out which elements from open_path object
     ;; correspond to field names in okeys_op.  Number of elements must
     ;; be the same.
     op_idx=[11, 12, 0, 1, 2, 3, 4, 5, 13, 14, 15, 16, 17, 6, 7, 8, 9, $
               10, 18, 19, 20, 21, 22, 23, 24, 25]
     FOREACH op_fld, indgen(n_elements(okeys_op)) DO BEGIN
        fluxes[okeys_op[op_fld]]=[fluxes[okeys_op[op_fld]], $
                                    open_path[op_idx[op_fld]]]
     ENDFOREACH
     ;; The mean diagnostic value for the open path flux data.  Not sure
     ;; what the meaning of such a value would be.  FIX THIS (below should
     ;; work, but 2013 we have all 'NAN' string and fails mean()/fix()
     ;; fluxes['diag_op']=[fluxes['diag_op'], $
     ;;                      mean(fix(flux.op_analyzer_status), /NAN)]
     fluxes['diag_op']=[fluxes['diag_op'], !VALUES.F_NAN]
     ;; For the calculated fields, we'll have to do it pseudo-manually,
     ;; since they are (so far) just strewn across the workspace
     calc_vals=[raw_sonic_spd, raw_sonic_dir, true_sonic_spd, $
                  true_sonic_dir, wind_v_mean, open_flag, $
                  sonic_flag, motion_flag, sonicNANs, sw_avg, CDm]
     FOREACH calc_fld, indgen(n_elements(okeys_calc)) DO BEGIN
        fluxes[okeys_calc[calc_fld]]=[fluxes[okeys_calc[calc_fld]], $
                                        calc_vals[calc_fld]]
     ENDFOREACH
     ;; What to do with these?
     ;; micro_vals=[Um, U10, CD10, z0, CHm, CH10, zT, CEm, CE10, $
     ;;             zQ, !VALUES.D_NAN, !VALUES.D_NAN, psim, psih, U10N, $
     ;;             U10Nocean]

  ENDFOREACH

  ;; Write the full hash, ordering fields in some way
  odata=create_struct('time', fluxes['time'], 'DOY', fluxes['DOY'])
  FOREACH fld, okeys_diag DO BEGIN ; MET summary data
     odata=create_struct(odata, $
                         okeys_diag[where(okeys_diag EQ fld)], $
                         fluxes[fld])
  ENDFOREACH
  FOREACH fld, okeys_mom DO BEGIN ; momentum data
     odata=create_struct(odata, $
                         okeys_mom[where(okeys_mom EQ fld)], fluxes[fld])
  ENDFOREACH
  FOREACH fld, okeys_op DO BEGIN ; open path data
     odata=create_struct(odata, $
                         okeys_op[where(okeys_op EQ fld)], fluxes[fld])
  ENDFOREACH
  odata=create_struct(odata, 'diag_op', fluxes['diag_op'])
  FOREACH fld, okeys_calc DO BEGIN ; calculated data
     odata=create_struct(odata, $
                         okeys_calc[where(okeys_calc EQ fld)], fluxes[fld])
  ENDFOREACH
  write_csv, ofile, odata, $
             header=['time', 'DOY', okeys_diag, okeys_mom, okeys_op, $
                     'diag_op', okeys_calc]

END



;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
