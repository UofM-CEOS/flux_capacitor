;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-10-15T15:43:56+0000
;; Last-Updated: 2013-10-28T22:00:20+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     PROCESS_MET
;; 
;; PURPOSE:
;; 
;;     Process MET data.
;; 
;; CALLING SEQUENCE:
;; 
;;     PROCESS_MET, Idir, Odir, Itemplate_Sav, Time_Beg_Idx, RMC_Dir, $
;;                  RMC_Itemplate_Sav, RMC_Time_Idx, RMC_Pull_Idx, $
;;                  GYRO_Dir, GYRO_Itemplate_Sav, GYRO_Time_Idx, $
;;                  GYRO_Pull_Idx, Log_File, Log_Itemplate_Sav, $
;;                  Log_Time_Beg_Idx, Log_Time_End_Idx
;; 
;; INPUTS:
;; 
;;     Idir:                Input directory (no trailing separator).
;;     Odir:                Output directory (no trailing separator).
;;     Itemplate_Sav:       Ascii template to read input files.
;;     Time_Beg_Idx:        Index (in template) where time is.
;;     RMC_Dir:             Directory where RMC files are found (no
;;                          trailing separator).
;;     RMC_Itemplate_Sav:   Ascii template to read RMC files.
;;     RMC_Time_Idx:        Index (in template) where time is in RMC files.
;;     RMC_Pull_Idx:        Integer array with indices (in template) of
;;                          fields to pull from RMC files.
;;     GYRO_Dir:            Directory where GYRO files are found (no
;;                          trailing separator).
;;     GYRO_Itemplate_Sav:  Ascii template to read GYRO files.
;;     GYRO_Time_Idx:       Index (in template) where time is in GYRO files.
;;     GYRO_Pull_Idx:       Integer array with indices (in template) of
;;                          fields to pull from GYRO files.
;;     Log_File:            Path to file with MET log.
;;     Log_Itemplate_Sav:   Ascii template to read MET log file.
;;     Log_Time_Beg_Idx:    Index (in template) where starting time is in
;;                          log file.
;;     Log_Time_End_Idx:    Index (in template) where ending time is in
;;                          log file.
;;     Log_Status_Idx:      Index (in template) where status flag is in
;;                          log file.
;; 
;; KEYWORD PARAMETERS:
;; 
;;     WIND_IDX_OFFSET:     2-element numeric array with the index of the
;;                          wind direction field and offset (due to the
;;                          position of the tower) in MET data.  Default:
;;                          [14, 23.0]
;;     OVERWRITE:           Whether to overwrite files in Odir.
;; 
;; RESTRICTIONS:
;; 
;;     This relies on the existence of the following field names being
;;     included via RMC_Pull_Idx and GYRO_Pull_Idx in the respective
;;     templates (case irrelevant): cog, sog, heading, wind_direction,
;;     wind_speed, true_wind_direction, true_wind_velocity, rh_percent,
;;     air_temperature, surface_temperature, pressure, pitch.
;; 
;; PROCEDURE:
;; 
;; 
;; 
;; EXAMPLE:
;; 
;;     process_MET, expand_path('~/tmp/ArcticNet2011/MET/Daily'), $
;;                  expand_path('~/tmp/ArcticNet2011/MET/Processed'), $
;;                  'met_std_template.sav', 0, $
;;                  expand_path('~/tmp/ArcticNet2011/OMG/NVG/1min'), $
;;                  'omg_nav_std_avg_template.sav', 0, [6, 7, 8, 9, 10, 11], $
;;                  expand_path('~/tmp/ArcticNet2011/OMG/HDG/1min'), $
;;                  'omg_hdg_std_avg_template.sav', 0, 6, $
;;                  expand_path('~/tmp/ArcticNet2011/Logs/met_log.csv'), $
;;                  'met_log_template.sav', 0, 5, 10
;; 
;;- -----------------------------------------------------------------------
;;; Code:

PRO PROCESS_MET, IDIR, ODIR, ITEMPLATE_SAV, TIME_BEG_IDX, RMC_DIR, $
                 RMC_ITEMPLATE_SAV, RMC_TIME_IDX, RMC_PULL_IDX, $
                 GYRO_DIR, GYRO_ITEMPLATE_SAV, GYRO_TIME_IDX, $
                 GYRO_PULL_IDX, LOG_FILE, LOG_ITEMPLATE_SAV, $
                 LOG_TIME_BEG_IDX, LOG_TIME_END_IDX, LOG_STATUS_IDX, $
                 WIND_IDX_OFFSET=WIND_IDX_OFFSET, OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 17) THEN $
     message, 'Usage: PROCESS_MET, IDIR, ODIR, ITEMPLATE_SAV, ' + $
              'TIME_BEG_IDX, RMC_DIR, RMC_ITEMPLATE_SAV, ' + $
              'RMC_TIME_IDX, RMC_PULL_IDX, GYRO_DIR, ' + $
              'GYRO_ITEMPLATE_SAV, GYRO_TIME_IDX, GYRO_PULL_IDX, ' + $
              'LOG_FILE, LOG_ITEMPLATE_SAV, LOG_TIME_BEG_IDX, ' + $
              'LOG_TIME_END_IDX, LOG_STATUS_IDX'
  IF ((n_elements(idir) EQ 0) OR (idir EQ '')) THEN $
     message, 'IDIR is undefined or is empty string'
  IF ((n_elements(odir) EQ 0) OR (odir EQ '')) THEN $
     message, 'ODIR is undefined or is empty string'
  IF ((n_elements(itemplate_sav) EQ 0) OR (itemplate_sav EQ '')) THEN $
     message, 'ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(time_beg_idx) NE 1) OR (time_beg_idx LT 0)) THEN $
     message, 'TIME_BEG_IDX must be a scalar >= zero'
  IF ((n_elements(rmc_dir) EQ 0) OR (rmc_dir EQ '')) THEN $
     message, 'RMC_DIR is undefined or is empty string'
  IF ((n_elements(rmc_itemplate_sav) EQ 0) OR $
      (rmc_itemplate_sav EQ '')) THEN $
     message, 'RMC_ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(rmc_time_idx) NE 1) OR (rmc_time_idx LT 0)) THEN $
     message, 'RMC_TIME_IDX must be a scalar >= zero'
  IF ((n_elements(rmc_pull_idx) LT 1) OR $
      (size(rmc_pull_idx, /type) NE 2)) THEN BEGIN
     message, 'RMC_PULL_IDX must be an integer array'
  ENDIF ELSE BEGIN
     rpi=where(rmc_pull_idx LT 0, nrpineg)
     IF nrpineg GT 0 THEN $
        message, 'RMC_PULL_IDX cannot have negative indices'
  ENDELSE
  IF ((n_elements(gyro_dir) EQ 0) OR (gyro_dir EQ '')) THEN $
     message, 'GYRO_DIR is undefined or is empty string'
  IF ((n_elements(gyro_itemplate_sav) EQ 0) OR $
      (gyro_itemplate_sav EQ '')) THEN $
     message, 'GYRO_ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(gyro_time_idx) NE 1) OR (gyro_time_idx LT 0)) THEN $
     message, 'GYRO_TIME_IDX must be a scalar >= zero'
  IF ((n_elements(gyro_pull_idx) LT 1) OR $
      (size(gyro_pull_idx, /type) NE 2)) THEN BEGIN
     message, 'GYRO_PULL_IDX must be an integer array'
  ENDIF ELSE BEGIN
     gpi=where(gyro_pull_idx LT 0, ngpineg)
     IF ngpineg GT 0 THEN $
        message, 'GYRO_PULL_IDX cannot have negative indices'
  ENDELSE
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
  IF nidir_files LT 1 THEN message, 'No input files found'
  rmc_files=file_search(rmc_dir + path_sep() + '*', count=nrmc_files, $
                        /nosort, /fold_case, /test_regular)
  IF nrmc_files LT 1 THEN message, 'No RMC files found.  Exiting'
  gyro_files=file_search(gyro_dir + path_sep() + '*', count=ngyro_files, $
                        /nosort, /fold_case, /test_regular)
  IF ngyro_files LT 1 THEN message, 'No GYRO files found.  Exiting'
  IF keyword_set(wind_idx_offset) THEN BEGIN
     IF n_elements(wind_idx_offset) NE 2 THEN $
        message, 'WIND_IDX_OFFSET must be a 2-element array'
     IF wind_idx_offset[0] LT 0 THEN $
        message, 'The index in WIND_IDX_OFFSET cannot be negative'
  ENDIF ELSE BEGIN
     wind_idx_offset=[14, 23.0] ; default
  ENDELSE

  ;; Parse input template
  restore, itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  met_template=itemplate
  field_names=strlowcase(met_template.FIELDNAMES)
  field_types=met_template.FIELDTYPES
  is_time_field=met_template.FIELDGROUPS EQ time_beg_idx
  ;; Ignore other groups when reading the data
  met_template.FIELDGROUPS=indgen(met_template.FIELDCOUNT)
  met_template.FIELDGROUPS[where(is_time_field)]=time_beg_idx
  non_time_fields=where(~is_time_field)
  non_time_field_names=field_names[non_time_fields]
  tags2remove=where(field_names EQ field_names[time_beg_idx])
  ;; Times
  tfields=where(is_time_field, /NULL)
  tnames=field_names[tfields]
  tnamesl=strsplit(tnames, '_', /extract)
  tnames_last=strarr(n_elements(tnamesl))
  tnames_id=strjoin((tnamesl[0])[0:n_elements(tnamesl[0]) - 2], '_')
  FOR i=0L, n_elements(tnames) - 1 DO $
     tnames_last[i]=tnamesl[i, n_elements(tnamesl[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  time_locs=locate_time_strings(tnames_last)

  ;; Parse log template
  restore, log_itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  log_template=itemplate
  log_field_names=strlowcase(log_template.FIELDNAMES)
  log_field_types=log_template.FIELDTYPES
  log_is_begtime_fld=log_template.FIELDGROUPS EQ log_time_beg_idx
  log_is_endtime_fld=log_template.FIELDGROUPS EQ log_time_end_idx
  ;; Beginning Times
  log_tbegfields=where(log_is_begtime_fld, /NULL)
  log_tbegnames=log_field_names[log_tbegfields]
  log_tbegnamesl=strsplit(log_tbegnames, '_', /extract)
  log_tbegnames_last=strarr(n_elements(log_tbegnamesl))
  log_tbegnames_id=(log_tbegnamesl[0])[0:n_elements(log_tbegnamesl[0]) - 2]
  log_tbegnames_id=strjoin(log_tbegnames_id, '_')
  FOR i=0L, n_elements(log_tbegnames) - 1 DO BEGIN
     lasti=log_tbegnamesl[i, n_elements(log_tbegnamesl[i]) - 1]
     log_tbegnames_last[i]=lasti
  ENDFOR
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  log_tbeg_locs=locate_time_strings(log_tbegnames_last)
  ;; Ending Times
  log_tendfields=where(log_is_endtime_fld, /NULL)
  log_tendnames=log_field_names[log_tendfields]
  log_tendnamesl=strsplit(log_tendnames, '_', /extract)
  log_tendnames_last=strarr(n_elements(log_tendnamesl))
  log_tendnames_id=(log_tendnamesl[0])[0:n_elements(log_tendnamesl[0]) - 2]
  log_tendnames_id=strjoin(log_tendnames_id, '_')
  FOR i=0L, n_elements(log_tendnames) - 1 DO BEGIN
     lasti=log_tendnamesl[i, n_elements(log_tendnamesl[i]) - 1]
     log_tendnames_last[i]=lasti
  ENDFOR
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  log_tend_locs=locate_time_strings(log_tendnames_last)

  ;; Parse RMC template
  restore, rmc_itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  rmc_template=itemplate
  rmc_field_names=strlowcase(rmc_template.FIELDNAMES)
  rmc_field_types=rmc_template.FIELDTYPES
  rmc_is_time_field=rmc_template.FIELDGROUPS EQ rmc_time_idx
  ;; Ignore other groups when reading the data
  rmc_template.FIELDGROUPS=indgen(rmc_template.FIELDCOUNT)
  rmc_template.FIELDGROUPS[where(rmc_is_time_field)]=rmc_time_idx
  rmc_non_time_fields=where(~rmc_is_time_field)
  rmc_non_time_field_names=rmc_field_names[rmc_non_time_fields]
  ;; Times
  rmc_tfields=where(rmc_is_time_field, /NULL)
  rmc_tnames=rmc_field_names[rmc_tfields]
  rmc_tnamesl=strsplit(rmc_tnames, '_', /extract)
  rmc_tnames_last=strarr(n_elements(rmc_tnamesl))
  rmc_tnames_id=strjoin((rmc_tnamesl[0])[0:n_elements(rmc_tnamesl[0]) - 2], '_')
  FOR i=0L, n_elements(rmc_tnames) - 1 DO $
     rmc_tnames_last[i]=rmc_tnamesl[i, n_elements(rmc_tnamesl[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  rmc_time_locs=locate_time_strings(rmc_tnames_last)
  ;; Break file names and extract the piece to match
  rmc_filesl=strsplit(rmc_files, '_.', /extract)
  rmc_files_a=rmc_filesl.toArray(/transpose) ; array with pieces in cols
  rmc_files_mstr_dims=size(rmc_files_a, /dimensions)
  ;; We want the 3rd to last piece
  rmc_files_mstr=rmc_files_a[rmc_files_mstr_dims[0] - 3, *]

  ;; Parse Gyro template
  restore, gyro_itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  gyro_template=itemplate
  gyro_field_names=strlowcase(gyro_template.FIELDNAMES)
  gyro_field_types=gyro_template.FIELDTYPES
  gyro_is_time_field=gyro_template.FIELDGROUPS EQ gyro_time_idx
  ;; Ignore other groups when reading the data
  gyro_template.FIELDGROUPS=indgen(gyro_template.FIELDCOUNT)
  gyro_template.FIELDGROUPS[where(gyro_is_time_field)]=gyro_time_idx
  gyro_non_time_fields=where(~gyro_is_time_field)
  gyro_non_time_field_names=gyro_field_names[gyro_non_time_fields]
  ;; Times
  gyro_tfields=where(gyro_is_time_field, /NULL)
  gyro_tnames=gyro_field_names[gyro_tfields]
  gyro_tnamesl=strsplit(gyro_tnames, '_', /extract)
  gyro_tnames_last=strarr(n_elements(gyro_tnamesl))
  gyro_tnames_id=strjoin((gyro_tnamesl[0])[0:n_elements(gyro_tnamesl[0]) - 2], '_')
  FOR i=0L, n_elements(gyro_tnames) - 1 DO $
     gyro_tnames_last[i]=gyro_tnamesl[i, n_elements(gyro_tnamesl[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  gyro_time_locs=locate_time_strings(gyro_tnames_last)
  ;; Break file names and extract the piece to match
  gyro_filesl=strsplit(gyro_files, '_.', /extract)
  gyro_files_a=gyro_filesl.toArray(/transpose) ; array with pieces in cols
  gyro_files_mstr_dims=size(gyro_files_a, /dimensions)
  ;; We want the 3rd to last piece
  gyro_files_mstr=gyro_files_a[gyro_files_mstr_dims[0] - 3, *]

  ;; Read log file
  log=read_ascii(log_file, template=log_template)
  log_names=strlowcase(tag_names(log))
  log_tbeg_loc=where(log_names EQ log_field_names[log_time_beg_idx])
  log_tend_loc=where(log_names EQ log_field_names[log_time_end_idx])
  log_beg_times=log.(log_tbeg_loc)
  log_end_times=log.(log_tend_loc)
  ;; Remove quotes and separators
  IF size(log_beg_times, /type) EQ 7 THEN BEGIN
     FOREACH fld, indgen((size(log_beg_times, $
                               /dimensions))[0]) DO BEGIN
        ok=strsplit(log_beg_times[fld, *], '" -/:', /extract)
        ok=(temporary(ok)).toArray()
        ok=strjoin(transpose(temporary(ok)))
        log_beg_times[fld, *]=ok
     ENDFOREACH
  ENDIF
  IF size(log_end_times, /type) EQ 7 THEN BEGIN
     FOREACH fld, indgen((size(log_end_times, $
                               /dimensions))[0]) DO BEGIN
        ok=strsplit(log_end_times[fld, *], '" -/:', /extract)
        ok=(temporary(ok)).toArray()
        ok=strjoin(transpose(temporary(ok)))
        log_end_times[fld, *]=ok
     ENDFOREACH
  ENDIF
  log_data_flds=where(~ (log_is_begtime_fld OR log_is_endtime_fld))
  FOREACH fld, (indgen(n_tags(log)))[log_data_flds] DO BEGIN
     IF size(log.(fld), /type) EQ 7 THEN BEGIN
        ok=strsplit(log.(fld), '" ', /extract)
        log.(fld)=ok.toArray()
     ENDIF
  ENDFOREACH
  ;; Obtain full log beginning time details and transform to julian
  log_beg_times_std=parse_times(log_beg_times, log_tbegnames_last, $
                                log_tbeg_locs)
  log_tbeg_jd=reform(julday(long(log_beg_times_std[1, *]), $
                            long(log_beg_times_std[2, *]), $
                            long(log_beg_times_std[0, *]), $
                            long(log_beg_times_std[3, *]), $
                            long(log_beg_times_std[4, *]), $
                            float(log_beg_times_std[5, *])))
  ;; Obtain full log ending time details and transform to julian
  log_end_times_std=parse_times(log_end_times, log_tendnames_last, $
                                log_tend_locs)
  log_tend_jd=reform(julday(long(log_end_times_std[1, *]), $
                            long(log_end_times_std[2, *]), $
                            long(log_end_times_std[0, *]), $
                            long(log_end_times_std[3, *]), $
                            long(log_end_times_std[4, *]), $
                            float(log_end_times_std[5, *])))
  logn=(size(log_beg_times, /dimensions))[1]

  ;; Loop through files in input directory
  FOR k=0, nidir_files - 1 DO BEGIN
     iname=strsplit(file_basename(idir_files[k]), '.', /extract)
     ;; Get a path for the file, check if it already exists
     ofile_name=strcompress(odir + path_sep() + iname[0] + '_proc.' + $
                            iname[1], /remove_all)
     ofile_stamp=file_basename(ofile_name)
     out_list=file_search(odir + path_sep() + '*.' + iname[1], $
                          /nosort, /fold_case, /test_regular)
     matchfiles=where(ofile_stamp EQ file_basename(out_list), $
                      matchfilecount)
     IF matchfilecount GT 0 THEN BEGIN
        IF keyword_set(overwrite) THEN BEGIN
           message, 'Processed file ' + ofile_stamp + $
                    ' already exists.  Overwriting', /informational
        ENDIF ELSE BEGIN
        message, 'Processed file ' + ofile_stamp + $
                 ' already exists.  Not overwriting', /informational
        CONTINUE
     ENDELSE
     ENDIF

     ifile=idir_files[k]
     message, 'Processing ' + ifile, /informational
     ;; Read input file
     idata=read_ascii(ifile, template=met_template)
     idata_names=strlowcase(tag_names(idata))
     time_loc=where(idata_names EQ field_names[time_beg_idx])
     idata_times=idata.(time_loc)
     ;; Number of lines in input
     lines=n_elements(idata_times[0, *])
     ;; Remove quotes and separators
     IF size(idata_times, /type) EQ 7 THEN BEGIN
        FOREACH fld, indgen((size(idata_times, $
                                  /dimensions))[0]) DO BEGIN
           ok=strsplit(idata_times[fld, *], '" -/:', /extract)
           ok=(temporary(ok)).toArray()
           ok=strjoin(transpose(temporary(ok)))
           idata_times[fld, *]=ok
        ENDFOREACH
     ENDIF
     match2, idata_names, field_names[tags2remove], is_time
     FOREACH fld, (indgen(n_tags(idata)))[where(is_time LT 0)] DO BEGIN
        IF size(idata.(fld), /type) EQ 7 THEN BEGIN
           ok=strsplit(idata.(fld), '" ', /extract)
           idata.(fld)=ok.toArray()
        ENDIF
     ENDFOREACH
     ;; Adjust wind direction
     wdfld=where(idata_names EQ field_names[wind_idx_offset[0]], /null)
     idata.(wdfld)=idata.(wdfld) - wind_idx_offset[1]
     x_wdir=where(idata.(wdfld) LT 1e-5, nx_wdir)
     IF nx_wdir GT 0 THEN $
        idata.(wdfld)[x_wdir]=360.0 + idata.(wdfld)[x_wdir]

     ;; Obtain full input time details
     itimes_std=parse_times(idata_times, tnames_last, time_locs)
     met_jd=reform(julday(long(itimes_std[1, *]), $
                          long(itimes_std[2, *]), $
                          long(itimes_std[0, *]), $
                          long(itimes_std[3, *]), $
                          long(itimes_std[4, *]), $
                          float(itimes_std[5, *])))

     ifile_strl=strsplit(ifile, '_.', /extract) ; break string
     ifile_mstr=ifile_strl[n_elements(ifile_strl) - 2]

     ;; Read matching RMC file
     rmc_pair=where(rmc_files_mstr EQ ifile_mstr, mcount)
     IF mcount LT 1 THEN BEGIN
        message, 'No matching RMC file found. Skipping.', /CONTINUE
        CONTINUE
     ENDIF
     rmc=read_ascii(rmc_files[rmc_pair], template=rmc_template)
     rmc_names=strlowcase(tag_names(rmc))
     ;; Obtain times and convert to Julian
     rmc_time_loc=where(rmc_names EQ rmc_field_names[rmc_time_idx])
     rmc_times=rmc.(rmc_time_loc)
     rmc_times_dims=size(rmc_times, /dimensions)
     ;; Remove quotes
     IF size(rmc_times, /type) EQ 7 THEN BEGIN
        FOREACH fld, indgen((size(rmc_times, $
                                  /dimensions))[0]) DO BEGIN
           ok=strsplit(rmc_times[fld, *], '" -/:', /extract)
           ok=(temporary(ok)).toArray()
           ok=strjoin(transpose(temporary(ok)))
           rmc_times[fld, *]=ok
        ENDFOREACH
     ENDIF
     rmc_jd=reform((rmc_times_dims[0] EQ 6) ? $
                   julday(long(rmc_times[1, *]), $
                          long(rmc_times[2, *]), $
                          long(rmc_times[0, *]), $
                          long(rmc_times[3, *]), $
                          long(rmc_times[4, *]), $
                          float(rmc_times[5, *])) : $
                   julday(long(rmc_times[1, *]), $
                          long(rmc_times[2, *]), $
                          long(rmc_times[0, *]), $
                          long(rmc_times[3, *]), $
                          long(rmc_times[4, *]), $
                          float(rmc_times[5, *] + '.' + $
                                rmc_times[6, *])))
     ;; Find matching times
     match2, met_jd, rmc_jd, met_in_rmc
     met_matches=where(met_in_rmc GE 0, mcount, /null)
     IF mcount LT 1 THEN BEGIN
        message, 'No matching RMC records found.  Skipping file.', $
                 /informational
        CONTINUE
     ENDIF
     FOREACH fld, rmc_pull_idx DO BEGIN
        match_fld=where(rmc_names EQ rmc_field_names[fld])
        idata=create_struct(idata, rmc_field_names[fld], $
                            rmc.(match_fld)[met_matches])
     ENDFOREACH

     ;; Read matching GYRO file
     gyro_pair=where(gyro_files_mstr EQ ifile_mstr, mcount)
     IF mcount LT 1 THEN BEGIN
        message, 'No matching GYRO file found. Skipping.', /CONTINUE
        CONTINUE
     ENDIF
     gyro=read_ascii(gyro_files[gyro_pair], template=gyro_template)
     gyro_names=strlowcase(tag_names(gyro))
     ;; Obtain times and convert to Julian
     gyro_time_loc=where(gyro_names EQ gyro_field_names[gyro_time_idx])
     gyro_times=gyro.(gyro_time_loc)
     gyro_times_dims=size(gyro_times, /dimensions)
     ;; Remove quotes
     IF size(gyro_times, /type) EQ 7 THEN BEGIN
        FOREACH fld, indgen((size(gyro_times, $
                                  /dimensions))[0]) DO BEGIN
           ok=strsplit(gyro_times[fld, *], '" -/:', /extract)
           ok=(temporary(ok)).toArray()
           ok=strjoin(transpose(temporary(ok)))
           gyro_times[fld, *]=ok
        ENDFOREACH
     ENDIF
     gyro_jd=reform((gyro_times_dims[0] EQ 6) ? $
                    julday(long(gyro_times[1, *]), $
                           long(gyro_times[2, *]), $
                           long(gyro_times[0, *]), $
                           long(gyro_times[3, *]), $
                           long(gyro_times[4, *]), $
                           float(gyro_times[5, *])) : $
                    julday(long(gyro_times[1, *]), $
                           long(gyro_times[2, *]), $
                           long(gyro_times[0, *]), $
                           long(gyro_times[3, *]), $
                           long(gyro_times[4, *]), $
                           float(gyro_times[5, *] + '.' + $
                                 gyro_times[6, *])))
     ;; Find matching times
     match2, met_jd, gyro_jd, met_in_gyro
     met_matches=where(met_in_gyro GE 0, mcount, /null) ; overwriting var
     IF mcount LT 1 THEN BEGIN
        message, 'No matching GYRO records found.  Skipping file.', $
                 /informational
        CONTINUE
     ENDIF
     FOREACH fld, gyro_pull_idx DO BEGIN
        match_fld=where(gyro_names EQ gyro_field_names[fld])
        idata=create_struct(idata, gyro_field_names[fld], $
                            gyro.(match_fld)[met_matches])
     ENDFOREACH

     ;; Calculate true wind direction and speed
     cog=idata.cog
     ;; Convert from nautical miles to m/s
     sog=idata.sog / 1.9438449
     ;; This is the bit about "beware of the COG when SOG is 0"
     cog_flag=finite(cog, /nan, sign=-1)
     cog0=where(cog_flag, ncog0)
     IF ncog0 GT 0 THEN $
        cog[cog0]=0
     twinds=truewind(0, cog, sog, idata.heading, idata.wind_direction, $
                     idata.wind_speed)
     tw_names=['true_wind_direction', 'true_wind_velocity']
     FOREACH fld, tw_names DO BEGIN
        idata=create_struct(idata, fld, twinds[*, where(tw_names EQ fld)])
     ENDFOREACH

     ;; Check the log status flag (initial value=0 -> not used in log, so
     ;; we use to OK data)
     status_flag=intarr(lines)
     log_status=log.(where(log_names EQ log_field_names[log_status_idx]))
     FOREACH logi, indgen(logn) DO BEGIN
        flagi=where((met_jd GE log_tbeg_jd[logi]) AND $
                    (met_jd LT log_tend_jd[logi]), nflagi)
        IF nflagi GT 0 THEN status_flag[flagi]=log_status[logi]
     ENDFOREACH

     ;; I don't understand why the flag is used to fiddle with variables
     ;; below here.  I think the flag should be used directly when
     ;; selecting appropriate records for calculation when needed.
     ;; Furthermore, the original code does this also assuming we have
     ;; standard deviations, which are only available if data correspond
     ;; from some averaging algorithm.  In 2011, input MET data are sampled
     ;; at 1-min, so no averaging was done, hence no standard deviations.
     ;; It seems overly manipulative.  Keeping all of it for now, until we
     ;; meet to sort all of this out.

     ;; Filter out bad wind direction
     baddir=where((idata.wind_direction GT 100) AND $
                  (idata.wind_direction LT 260), nbaddir)
     IF nbaddir GT 0 THEN BEGIN
        idata.true_wind_direction[baddir]=!VALUES.F_NAN
        idata.true_wind_velocity[baddir]=!VALUES.F_NAN
     ENDIF
     ;; For "very" bad wind direction, also remove T/RH data. NOTE: what we
     ;; are really filtering out here is from 170 to 190 deg (due to
     ;; anemometer offset for this year)
     vbaddir=where((idata.wind_direction GT 170) AND $
                   (idata.wind_direction LT 190), nvbaddir)
     IF nvbaddir GT 0 THEN BEGIN
        idata.rh_percent[vbaddir]=!VALUES.F_NAN
        idata.air_temperature[vbaddir]=!VALUES.F_NAN
        ;; We just don't have standard deviations in 2011, since we're not
        ;; averaging data and working straight from daily files.
        ;; idata.rh_percent_sd[vbaddir]=!VALUES.F_NAN
        ;; idata.surface_temperature_sd[vbaddir]=!VALUES.F_NAN
     ENDIF

     ;; Filter out ice breaking (where SOG stdev > 2, or COG stdev > 10)
     badice=where((idata.sog_sd GT 2) OR (idata.cog_sd GT 10), nbadice)
     IF nbadice GT 0 THEN BEGIN
        idata.true_wind_direction[badice]=!VALUES.F_NAN
        idata.true_wind_velocity[badice]=!VALUES.F_NAN
     ENDIF

     ;; Filter out of range Patm
     badPatm=where(idata.pressure LT 94, nbadpatm)
     IF nbadpatm GT 0 THEN BEGIN
        idata.pressure[badPatm]=!VALUES.F_NAN ;filter out P_avg
        ;; See note above about standard deviations
        ;; work_arr(21,badPatm) = 'NaN' ;filter out P_std
        status_flag[badPatm]=5
     ENDIF

     ;; Filter out of range pitch/roll; we can use the compass angle data
     ;; as a "second check" on tower down situations
     badroll=where(idata.pitch GT 50, nbadroll)
     IF nbadroll GT 0 THEN BEGIN
        idata.rh_percent[badroll]=!VALUES.F_NAN
        idata.air_temperature[badroll]=!VALUES.F_NAN
        idata.wind_speed[badroll]=!VALUES.F_NAN
        idata.wind_direction[badroll]=!VALUES.F_NAN
        ;; Not doing anything to standard deviations, as per note above
        idata.true_wind_direction[badroll]=!VALUES.F_NAN
        idata.true_wind_velocity[badroll]=!VALUES.F_NAN
        status_flag[badroll]=1
     endif

     flag1=where(status_flag EQ 1, nflag1)
     IF nflag1 GT 0 THEN BEGIN
        idata.air_temperature[flag1]=!VALUES.F_NAN
        idata.rh_percent[flag1]=!VALUES.F_NAN
        idata.wind_direction[flag1]=!VALUES.F_NAN
        idata.wind_speed[flag1]=!VALUES.F_NAN
        idata.true_wind_direction[flag1]=!VALUES.F_NAN
        idata.true_wind_velocity[flag1]=!VALUES.F_NAN
     ENDIF
     flag2=where(status_flag EQ 2, nflag2)
     IF nflag2 GT 0 THEN BEGIN
        idata.wind_speed[flag2]=!VALUES.F_NAN
        idata.true_wind_direction[flag1]=!VALUES.F_NAN
        idata.true_wind_velocity[flag1]=!VALUES.F_NAN
     ENDIF
     flag3=where(status_flag EQ 3, nflag3)
     IF nflag3 GT 0 THEN BEGIN
        idata.surface_temperature[flag3]=!VALUES.F_NAN
     ENDIF
     flag4=where(status_flag EQ 4, nflag4)
     IF nflag4 GT 0 THEN BEGIN
        idata.air_temperature[flag4]=!VALUES.F_NAN
        idata.rh_percent[flag4]=!VALUES.F_NAN
     ENDIF
     flag5=where(status_flag EQ 5, nflag5)
     IF nflag5 GT 0 THEN BEGIN
        idata.pressure[flag5]=!VALUES.F_NAN
     ENDIF

     idata=create_struct(idata, 'diag', status_flag)
     odata=remove_structure_tags(idata, field_names[tags2remove])
     delvar, idata
     ;; OK, how else to just extract the time info into a structure
     revtidx=reverse(indgen((size(idata_times, /dimensions))[0]))
     FOREACH fld, revtidx DO BEGIN
        fld_name=tnames[fld]
        odata=create_struct(fld_name, reform(idata_times[fld, *]), odata)
     ENDFOREACH
    
     write_csv, ofile_name, odata, header=strlowcase(tag_names(odata))

  ENDFOR


END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; process_MET.pro ends here
