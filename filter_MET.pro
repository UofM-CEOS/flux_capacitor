;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-10-29T14:29:43+0000
;; Last-Updated: 2013-11-02T00:54:13+0000
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
;;     ship velocities (SOG) are too far (given threshold) from the mean, 4
;;     - some headings are too far (given threshold) from the mean, 5 -
;;     some true wind speeds not available.
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
;;                            the following fields of files in Avg_Dir
;;                            (order is relevant): mean wind direction
;;                            (raw), mean SOG, and mean heading.
;;     Filter_Thr:            Array with thresholds for the fields in
;;                            Filter_Idx.
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
  IF ((n_elements(avg_dir) EQ 0) OR (avg_dir EQ '')) THEN $
     message, 'AVG_DIR is undefined or is empty string'
  IF ((n_elements(full_dir) EQ 0) OR (full_dir EQ '')) THEN $
     message, 'FULL_DIR is undefined or is empty string'
  IF ((n_elements(avg_itemplate_sav) EQ 0) OR $
      (avg_itemplate_sav EQ '')) THEN $
         message, 'AVG_ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(avg_time_idx) NE 1) OR (avg_time_idx LT 0)) THEN $
     message, 'AVG_TIME_IDX must be a scalar >= zero'
  IF ((n_elements(avg_period) NE 1) OR (avg_period LT 0)) THEN $
     message, 'AVG_PERIOD must be a scalar >= zero'
  IF ((n_elements(full_itemplate_sav) EQ 0) OR $
      (full_itemplate_sav EQ '')) THEN $
         message, 'FULL_ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(full_time_idx) NE 1) OR (full_time_idx LT 0)) THEN $
     message, 'FULL_TIME_IDX must be a scalar >= zero'
  IF ((n_elements(full_sample_rate) NE 1) OR $
      (full_sample_rate LT 0)) THEN $
         message, 'FULL_SAMPLE_RATE must be a scalar >= zero'
  n_fi=n_elements(filter_idx)
  n_ft=n_elements(filter_thr)
  IF n_fi NE 3 THEN message, 'FILTER_IDX must be a 3-element integer array'
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
     
     wind_direction_loc=where(idata_names EQ avg_field_names[filter_idx[0]])
     sog_loc=where(idata_names EQ avg_field_names[filter_idx[1]])
     heading_loc=where(idata_names EQ avg_field_names[filter_idx[2]])

     ;; Set the flag; 0s are OK, and the rest integers indicate which check
     ;; caused the period to fail
     fluxable=intarr(lines)

     ;; 1st flag periods outside of right quadrant
     badquadrant=where((idata.(wind_direction_loc) LT 270) AND $
                       (idata.(wind_direction_loc) GT 90), nbadquadrant)
     IF nbadquadrant GT 0 THEN $
        fluxable[badquadrant]=1

;;      if ((avg_wind_dir lt 270) and (avg_wind_dir gt 90)) then begin
;;   printf, 4, 'Flux Run   ' + strcompress(fix(filt_arr(0,x))) + strcompress(fix(filt_arr(1,x))) + strcompress(fix(filt_arr(2,x))) + $
;;               strcompress(fix(filt_arr(3,x))) + strcompress(fix(filt_arr(4,x))) + ' - ' + strcompress(fix(filt_arr(0,x+period-1))) + $
;;               strcompress(fix(filt_arr(1,x+period-1))) + strcompress(fix(filt_arr(2,x+period-1))) + strcompress(fix(filt_arr(3,x+period-1))) + $
;;               strcompress(fix(filt_arr(4,x+period-1))) + '   failed due to mean wind outside of quadrant'
;;   print, 'Flux Run   ' + strcompress(fix(filt_arr(0,x))) + strcompress(fix(filt_arr(1,x))) + strcompress(fix(filt_arr(2,x))) + $
;;               strcompress(fix(filt_arr(3,x))) + strcompress(fix(filt_arr(4,x))) + ' - ' + strcompress(fix(filt_arr(0,x+period-1))) + $
;;               strcompress(fix(filt_arr(1,x+period-1))) + strcompress(fix(filt_arr(2,x+period-1))) + strcompress(fix(filt_arr(3,x+period-1))) + $
;;               strcompress(fix(filt_arr(4,x+period-1))) + '   failed due to mean wind outside of quadrant'              
;;   goto, fail
;; endif


  ENDFOR

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; filter_MET.pro ends here
