;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-10-29T14:29:43+0000
;; Last-Updated: 2013-12-05T19:29:37+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     FILTER_MET
;; 
;; PURPOSE:
;; 
;;     This procedure checks whether processed and averaged MET data within
;;     a particular time period (e.g. 20 min) meets the conditions required
;;     for eddy covariance calculations.  It generates a file with
;;     requested fields from the processed MET file, consisting of the mean
;;     values within the requested period, starting from a requested start
;;     time to and end time.
;;
;;     A flag is added to the flux-period-averaged MET files where 0
;;     indicates the period is suitable for flux analyses, and integers > 0
;;     indicate the failure reason: 1 - bad quadrant, 2 - some wind
;;     directions are too far (given threshold) from the mean, 3 - some
;;     ship speeds (SOG) are too far (given threshold) from the mean, 4 -
;;     some headings are too far (given threshold) from the mean, 5 - some
;;     true wind speeds not available.
;; 
;; CALLING SEQUENCE:
;; 
;;     FILTER_MET, Avg_Dir, Full_Dir, Avg_Itemplate_Sav, Avg_Time_Idx, $
;;                 Avg_Period, Full_Itemplate_Sav, Full_Time_Idx, $
;;                 Full_Sample_Rate, Filter_Idx, Filter_Thr
;; 
;; INPUTS:
;; 
;;     Avg_Dir:               Directory (no trailing separator) containing
;;                            the files from Full_Dir averaged over
;;                            Avg_Period's.
;;     Full_Dir:              Directory (no trailing separator) containing
;;                            the processed MET files; typically 1-min
;;                            daily files.
;;     Avg_Itemplate_Sav:     Ascii template to read files in Avg_Dir.
;;     Avg_Time_Idx:          Scalar with index (in template) where time is
;;                            in files in Avg_Dir.
;;     Avg_Period:            Scalar with duration (s) of each period in
;;                            files under Avg_Dir.
;;     Full_Itemplate_Sav:    Ascii template to read files in Full_Dir.
;;     Full_Time_Idx:         Scalar with index (in template) where time is
;;                            in files in Full_Dir.
;;     Full_Sample_Rate:      Scalar with sample rate (s) of files under
;;                            Full_Dir.
;;     Filter_Idx:            Integer array with indices (in template) of
;;                            the following fields for files in Avg_Dir
;;                            (order is relevant): mean wind direction
;;                            (raw), mean SOG, mean heading, and mean true
;;                            wind speed.
;;     Filter_Thr:            Array with values to determine thresholds for
;;                            the fields in Filter_Idx.  Threshold is
;;                            computed as average +- value.
;; 
;; KEYWORD PARAMETERS:
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
;;- -----------------------------------------------------------------------
;;; Code:

PRO FILTER_MET, AVG_DIR, FULL_DIR, AVG_ITEMPLATE_SAV, AVG_TIME_IDX, $
                AVG_PERIOD, FULL_ITEMPLATE_SAV, FULL_TIME_IDX, $
                FULL_SAMPLE_RATE, FILTER_IDX, FILTER_THR

  ;; Check parameters
  IF (n_params() NE 10) THEN $
     message, 'Usage: FILTER_MET, AVG_DIR, FULL_DIR, ' + $
              'AVG_ITEMPLATE_SAV, AVG_TIME_IDX, AVG_PERIOD, ' + $
              'FULL_ITEMPLATE_SAV, FULL_TIME_IDX, FULL_SAMPLE_RATE, ' + $
              'FILTER_IDX, FILTER_THR'
  avg_dir_info=file_info(avg_dir)
  avg_tpl_info=file_info(avg_itemplate_sav)
  full_dir_info=file_info(full_dir)
  full_tpl_info=file_info(full_itemplate_sav)
  IF (~avg_dir_info.directory) THEN $
     message, 'AVG_DIR must be a string pointing to an existing directory'
  IF (~avg_tpl_info.read) THEN $
     message, 'AVG_ITEMPLATE_SAV must be a string pointing to a ' + $
              'readable file'
  IF ((n_elements(avg_time_idx) NE 1) OR $
      ((size(avg_time_idx, /type) NE 2) || avg_time_idx LT 0)) THEN $
         message, 'AVG_TIME_IDX must be an integer scalar >= zero'
  IF (~full_dir_info.directory) THEN $
     message, 'FULL_DIR must be a string pointing to an existing directory'
  IF (~full_tpl_info.read) THEN $
     message, 'FULL_ITEMPLATE_SAV must be a string pointing to a ' + $
              'readable file'
  IF ((n_elements(full_time_idx) NE 1) OR $
      ((size(full_time_idx, /type) NE 2) || full_time_idx LT 0)) THEN $
         message, 'FULL_TIME_IDX must be an integer scalar >= zero'
  IF ((n_elements(avg_period) NE 1) OR (avg_period LT 0)) THEN $
     message, 'AVG_PERIOD must be a scalar >= zero'
  IF ((n_elements(full_sample_rate) NE 1) OR $
      (full_sample_rate LT 0)) THEN $
         message, 'FULL_SAMPLE_RATE must be a scalar >= zero'
  n_fi=n_elements(filter_idx)
  n_ft=n_elements(filter_thr)
  IF n_fi NE 4 THEN message, 'FILTER_IDX must be a 4-element integer array'
  IF n_fi NE n_ft THEN $
     message, 'FILTER_THR must be supplied along with FILTER_IDX'
  IF n_ft GT 0 THEN BEGIN
     fi_size=size(filter_idx)
     IF (fi_size[0] GT 1) && (fi_size[2] NE 2) THEN $
        message, 'FILTER_IDX must be an integer scalar or vector'
     ft_size=size(filter_thr)
     IF (ft_size[0] GT 1) && (ft_size[2] NE 2) THEN $
        message, 'FILTER_THR must be an integer scalar or vector'
     IF ft_size[1] NE fi_size[1] THEN $
        message, 'FILTER_THR and FILTER_IDX must have the same ' + $
                 'number of elements'
  ENDIF
  full_files=file_search(full_dir + path_sep() + '*', $
                         count=nfull_files, $
                         /nosort, /fold_case, /test_regular)
  avg_files=file_search(avg_dir + path_sep() + '*', $
                        count=navg_files, $
                        /nosort, /fold_case, /test_regular)
  IF (nfull_files EQ 0) OR (navg_files EQ 0) THEN $
     message, 'No files in FULL_DIR or AVG_DIR'
  IF nfull_files NE navg_files THEN $
     message, 'Unequal number of files in FULL_DIR and AVG_DIR'

  ;; Parse average template
  restore, avg_itemplate_sav
  ;; We make a copy, since the RMC template is also called itemplate
  avg_template=itemplate
  avg_field_names=strlowcase(avg_template.FIELDNAMES)
  ;; Avg file times
  tfields_avg=where(avg_template.FIELDGROUPS EQ avg_time_idx, /NULL)
  tnames_avg=avg_field_names[tfields_avg]
  tnamesl_avg=strsplit(tnames_avg, '_', /extract)
  tnames_last_avg=strarr(n_elements(tnamesl_avg))
  ntn=n_elements(tnamesl_avg[0])
  tnames_id_avg=strjoin((tnamesl_avg[0])[0:ntn - 2], '_')
  FOR i=0L, n_elements(tnames_avg) - 1 DO $
     tnames_last_avg[i]=tnamesl_avg[i, n_elements(tnamesl_avg[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  time_locs_avg=locate_time_strings(tnames_last_avg)
  ;; Break Avg file names and extract the piece to match
  avg_filesl=strsplit(avg_files, '_.', /extract)
  avg_files_a=avg_filesl.toArray(/transpose) ; array with pieces in cols
  avg_files_mstr_dims=size(avg_files_a, /dimensions)
  ;; We want the 4th to last piece
  avg_files_mstr=avg_files_a[avg_files_mstr_dims[0] - 4, *]

  ;; Parse full template
  restore, full_itemplate_sav
  ;; We make a copy, since the RMC template is also called itemplate
  full_template=itemplate
  full_field_names=strlowcase(full_template.FIELDNAMES)
  ;; Full file times
  tfields_full=where(full_template.FIELDGROUPS EQ full_time_idx, /NULL)
  tnames_full=full_field_names[tfields_full]
  tnamesl_full=strsplit(tnames_full, '_', /extract)
  tnames_last_full=strarr(n_elements(tnamesl_full))
  ntn=n_elements(tnamesl_full[0])
  tnames_id_full=strjoin((tnamesl_full[0])[0:ntn - 2], '_')
  FOR i=0L, n_elements(tnames_full) - 1 DO $
     tnames_last_full[i]=tnamesl_full[i, n_elements(tnamesl_full[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  time_locs_full=locate_time_strings(tnames_last_full)
  ;; Break Full file names and extract the piece to match
  full_filesl=strsplit(full_files, '_.', /extract)
  full_files_a=full_filesl.toArray(/transpose) ; array with pieces in cols
  full_files_mstr_dims=size(full_files_a, /dimensions)
  ;; We want the 3rd to last piece
  full_files_mstr=full_files_a[full_files_mstr_dims[0] - 3, *]

  tags2remove=where(avg_field_names EQ avg_field_names[avg_time_idx])
  ;; Shape of matrix to analyse full time series (as in average_series)
  ncols_avgs=long(avg_period / full_sample_rate)
  nrows_avgs=86400 / long(full_sample_rate * ncols_avgs)

  ;; Loop through files in input directory
  FOR k=0, navg_files - 1 DO BEGIN
     ifile=avg_files[k]
     message, 'Processing ' + ifile, /informational
     ;; Read input file
     idata=read_ascii(ifile, template=avg_template)
     idata_names=strlowcase(tag_names(idata))
     avg_time_loc=where(idata_names EQ avg_field_names[avg_time_idx])
     idata_times_avg=idata.(avg_time_loc)
     ;; Number of lines in input
     lines=n_elements(idata_times_avg[0, *])
     ;; Remove quotes and separators
     IF size(idata_times_avg, /type) EQ 7 THEN BEGIN
        FOREACH fld, indgen((size(idata_times_avg, $
                                  /dimensions))[0]) DO BEGIN
           ok=strsplit(idata_times_avg[fld, *], '" -/:', /extract)
           ok=(temporary(ok)).toArray()
           ok=strjoin(transpose(temporary(ok)))
           idata_times_avg[fld, *]=ok
        ENDFOREACH
     ENDIF
     match2, idata_names, avg_field_names[tags2remove], is_time
     FOREACH fld, (indgen(n_tags(idata)))[where(is_time LT 0)] DO BEGIN
        IF size(idata.(fld), /type) EQ 7 THEN BEGIN
           ok=strsplit(idata.(fld), '" ', /extract)
           idata.(fld)=ok.toArray()
           ;; Replace NANs (bad) with empty string since they are
           ;; automatically turned into NaN elsewhere
           badnan=where(idata.(fld) EQ 'NAN', nbadnan)
           IF nbadnan GT 0 THEN idata.(fld)[badnan]=''
        ENDIF
     ENDFOREACH
     ;; Obtain full Avg time details
     itimes_avg_std=parse_times(idata_times_avg, tnames_last_avg, $
                                time_locs_avg)
     avg_jd=reform(julday(long(itimes_avg_std[1, *]), $
                          long(itimes_avg_std[2, *]), $
                          long(itimes_avg_std[0, *]), $
                          long(itimes_avg_std[3, *]), $
                          long(itimes_avg_std[4, *]), $
                          float(itimes_avg_std[5, *])))

     ;; Read matching Full file
     ifile_strl=strsplit(ifile, '_.', /extract) ; break string
     ifile_mstr=ifile_strl[n_elements(ifile_strl) - 4]
     full_pair=where(full_files_mstr EQ ifile_mstr, mcount)
     IF mcount LT 1 THEN BEGIN
        message, 'No matching file found. Skipping.', /CONTINUE
        CONTINUE
     ENDIF
     full=read_ascii(full_files[full_pair], template=full_template)
     full_names=strlowcase(tag_names(full))
     ;; Obtain Full times and convert to Julian
     full_time_loc=where(full_names EQ full_field_names[full_time_idx])
     full_times=full.(full_time_loc)
     full_times_dims=size(full_times, /dimensions)
     ;; Remove quotes
     IF size(full_times, /type) EQ 7 THEN BEGIN
        FOREACH fld, indgen((size(full_times, $
                                  /dimensions))[0]) DO BEGIN
           ok=strsplit(full_times[fld, *], '" -/:', /extract)
           ok=(temporary(ok)).toArray()
           ok=strjoin(transpose(temporary(ok)))
           full_times[fld, *]=ok
        ENDFOREACH
     ENDIF
     full_jd=reform(julday(long(full_times[1, *]), $
                           long(full_times[2, *]), $
                           long(full_times[0, *]), $
                           long(full_times[3, *]), $
                           long(full_times[4, *]), $
                           float(full_times[5, *])))

     ;; Checks
     
     avg_wd_loc=where(idata_names EQ avg_field_names[filter_idx[0]])
     avg_sog_loc=where(idata_names EQ avg_field_names[filter_idx[1]])
     avg_head_loc=where(idata_names EQ avg_field_names[filter_idx[2]])
     avg_ws_loc=where(idata_names EQ avg_field_names[filter_idx[3]])
     full_wd_loc=where(full_names EQ avg_field_names[filter_idx[0]])
     full_sog_loc=where(full_names EQ avg_field_names[filter_idx[1]])
     full_head_loc=where(full_names EQ avg_field_names[filter_idx[2]])
     full_ws_loc=where(full_names EQ avg_field_names[filter_idx[3]])
     ;; Subset check fields
     wd_avg=idata.(avg_wd_loc)
     sog_avg=idata.(avg_sog_loc)
     head_avg=idata.(avg_head_loc)
     ws_avg=idata.(avg_ws_loc)
     wd_full=full.(full_wd_loc)
     sog_full=full.(full_sog_loc)
     head_full=full.(full_head_loc)
     ws_full=full.(full_ws_loc)
     ;; Set the flag; 0s are OK, and the rest integers indicate which check
     ;; caused the period to fail
     fluxable=intarr(lines)

     ;; 1: Flag periods where raw wind direction is outside of right
     ;; quadrant
     badquadrant=where((wd_avg LT 270) AND (wd_avg GT 90), nbadquadrant)
     IF nbadquadrant GT 0 THEN fluxable[badquadrant]=1

     ;; 2: Flag periods where raw wind direction was too far from the mean
     ;; for the period, even for a single record
     full_wd2d=reform(wd_full, ncols_avgs, nrows_avgs)
     wd_lo=wd_avg - filter_thr[0] ; lower threshold
     ;; Check if we went below 360 for threshold and correct it
     wd_lo_neg=where(wd_lo LT 0, nwd_lo_neg)
     IF nwd_lo_neg GT 0 THEN $
        wd_lo[wd_lo_neg]=360 + wd_lo[wd_lo_neg]
     wd_hi=wd_avg + filter_thr[0] ; upper threshold
     ;; Check if we went above 360 for threshold and correct it
     wd_hi_over=where(wd_hi GT 360, nwd_hi_over)
     IF nwd_hi_over GT 0 THEN $
        wd_hi[wd_hi_over]=(wd_hi[wd_hi_over] * !DTOR) MOD (2 * !PI) / !DTOR
     ;; For each series of angles we check whether the lower bound is
     ;; larger than the upper bound.  In these cases, we reject periods
     ;; where any angle < lower bound *and* > upper bound.  In other cases,
     ;; i.e. when we are not crossing the 360 border, we proceed normally
     ;; by rejecting periods where any angle < lower bound *or* > upper
     ;; bound.
     FOR i=0L, nrows_avgs - 1 DO BEGIN
        wd_ok=where(finite(full_wd2d[*, i]) GT 0, nwd_ok, $
                    complement=wd_bad)
        IF nwd_ok EQ ncols_avgs THEN BEGIN ; period with complete data
           IF wd_hi[i] GE wd_lo[i] THEN BEGIN ; not crossing
              badwd=where((full_wd2d[*, i] LT wd_lo[i]) OR $
                          (full_wd2d[*, i] GT wd_hi[i]), nbadwd)
           ENDIF ELSE BEGIN     ; crossing 360
              badwd=where((full_wd2d[*, i] LT wd_lo[i]) AND $
                          (full_wd2d[*, i] GT wd_hi[i]), nbadwd)
           ENDELSE
           IF nbadwd GT 0 THEN fluxable[i]=2 ; got some outlier angles
        ENDIF ELSE fluxable[i]=5 ; same as flag 5
     ENDFOR

     ;; 3: Flag periods where ship speed (SOG) was too far from the mean
     ;; for the period
     full_sog2d=reform(sog_full, ncols_avgs, nrows_avgs)
     sog_lo=sog_avg - filter_thr[1] ; lower threshold
     sog_hi=sog_avg + filter_thr[1]   ; upper threshold
     FOR i=0L, nrows_avgs - 1 DO BEGIN
        badsog=where((full_sog2d[*, i] LT sog_lo[i]) OR $
                     (full_sog2d[*, i] GT sog_hi[i]), nbadsog)
        IF nbadsog GT 0 THEN fluxable[i]=3
     ENDFOR

     ;; 4: Flag periods where ship heading was too far from the mean for
     ;; the period, even for a single record.  Proceed as for wind
     ;; direction.
     full_head2d=reform(head_full, ncols_avgs, nrows_avgs)
     head_lo=head_avg - filter_thr[2] ; lower threshold
     head_lo_neg=where(head_lo LT 0, nhead_lo_neg)
     IF nhead_lo_neg GT 0 THEN $
        head_lo[head_lo_neg]=360 + head_lo[head_lo_neg]
     head_hi=head_avg + filter_thr[2] ; upper threshold
     head_hi_over=where(head_hi GT 360, nhead_hi_over)
     IF nhead_hi_over GT 0 THEN $
        head_hi[head_hi_over]=(head_hi[head_hi_over] * !DTOR) MOD $
                              (2 * !PI) / !DTOR
     FOR i=0L, nrows_avgs - 1 DO BEGIN
        head_ok=where(finite(full_head2d[*, i]) GT 0, nhead_ok, $
                      complement=head_bad)
        IF nhead_ok EQ ncols_avgs THEN BEGIN ; period with complete data
           IF head_hi[i] GE head_lo[i] THEN BEGIN ; not crossing
              badhead=where((full_head2d[*, i] LT head_lo[i]) OR $
                            (full_head2d[*, i] GT head_hi[i]), nbadhead)
           ENDIF ELSE BEGIN     ; crossing 360
              badhead=where((full_head2d[*, i] LT head_lo[i]) AND $
                            (full_head2d[*, i] GT head_hi[i]), nbadhead)
           ENDELSE
           IF nbadhead GT 0 THEN fluxable[i]=4 ; got some outlier angles
        ENDIF ELSE fluxable[i]=5               ; same as flag 5
     ENDFOR

     ;; 5: flag periods where any of the true wind speeds is unavailable
     full_ws2d=reform(full.(full_ws_loc), ncols_avgs, nrows_avgs)
     FOR i=0L, nrows_avgs - 1 DO BEGIN
        ws_bad=where(~ finite(full_ws2d[*, i]), nws_bad, complement=ws_ok)
        IF nws_bad GT 0 THEN fluxable[i]=5
     ENDFOR

     idata=create_struct(idata, 'diag_flux', fluxable)
     odata=remove_structure_tags(idata, avg_field_names[tags2remove])
     delvar, idata
     ;; OK, how else to just extract the time info into a structure
     revtidx=reverse(indgen((size(idata_times_avg, /dimensions))[0]))
     FOREACH fld, revtidx DO BEGIN
        fld_name=tnames_avg[fld]
        odata=create_struct(fld_name, $
                            reform(idata_times_avg[fld, *]), odata)
     ENDFOREACH

     write_csv, ifile, odata, header=strlowcase(tag_names(odata))

  ENDFOR

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; filter_MET.pro ends here
