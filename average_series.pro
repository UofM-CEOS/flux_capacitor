;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-09-17T14:59:07+0000
;; Last-Updated: 2013-11-29T04:19:35+0000
;;           By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;;
;;     AVERAGE_SERIES
;;
;; PURPOSE:
;;
;;     Calculate temporal averages for input file, taking into account
;;     angle and corresponding magnitude fields, if present.
;;
;; CALLING SEQUENCE:
;;
;;     AVERAGE_SERIES, Idir, Odir, Itemplate_Sav, Time_Beg_Idx, $
;;                     Isample_Rate, Osample_Rate
;;
;; INPUTS:
;;
;;     Idir:                  Input directory (no trailing separator).
;;     Odir:                  Output directory (no trailing separator).
;;     Itemplate_Sav:         Ascii template to read input files.
;;     Time_Idx:              Index (in template) where time is.
;;     Isample_Rate:          Scalar indicating the frequency (s) with
;;                            which input data were sampled.
;;     Osample_Rate:          Scalar indicating the frequency (s) with
;;                            which output data should be averaged over.
;;
;; KEYWORD PARAMETERS:
;;
;;     ANGLE_FIELDS:          An integer array indicating which fields in
;;                            the input template should be treated as
;;                            angular data.
;;     MAGNITUDE_FIELDS:      An integer array indicating which fields in
;;                            the input template should be treated as
;;                            magnitude for each field in ANGLE_FIELDS.  If
;;                            an element is -1, then magnitude is assumed
;;                            to be 1 for the corresponding angle.
;;     OVERWRITE:             Whether to overwrite files in Odir.
;;
;; RESTRICTIONS:
;;
;;     Input MUST have a full day's data; i.e. daily split.
;;
;; PROCEDURE:
;;
;;
;;
;; EXAMPLE:
;;
;;     AVERAGE_SERIES, expand_path('~/tmp/ArcticNet2011/NAV/Daily'), $
;;                     expand_path('~/tmp/ArcticNet2011/NAV/1min'), $
;;                     'nav_std_template.sav', 0, 1, 60, $
;;                     angle_fields=[16, 17], magnitude_fields=[15, -1], $
;;                     /overwrite
;;
;;
;;- -----------------------------------------------------------------------
;;; Code:

PRO AVERAGE_SERIES, IDIR, ODIR, ITEMPLATE_SAV, TIME_IDX, ISAMPLE_RATE, $
                    OSAMPLE_RATE, ANGLE_FIELDS=ANGLE_FIELDS, $
                    MAGNITUDE_FIELDS=MAGNITUDE_FIELDS, $
                    OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 6) THEN $
     message, 'Usage: AVERAGE_SERIES, IDIR, ODIR, ITEMPLATE_SAV, ' + $
              'TIME_IDX, ISAMPLE_RATE, OSAMPLE_RATE, ' + $
              'ANGLE_FIELDS, MAGNITUDE_FIELDS, STAMP'
  idir_info=file_info(idir)
  itpl_info=file_info(itemplate_sav)
  IF (~idir_info.directory) THEN $
     message, 'IDIR must be a string pointing to an existing directory'
  IF (~itpl_info.read) THEN $
     message, 'ITEMPLATE_SAV must be a string pointing to a readable file'
  IF ((n_elements(odir) NE 1) OR (size(odir, /type) NE 7)) THEN $
         message, 'ODIR must be a string scalar'
  IF ((n_elements(time_idx) NE 1) OR $
      ((size(time_idx, /type) NE 2) || time_idx LT 0)) THEN $
         message, 'TIME_IDX must be an integer scalar >= zero'
  IF ((n_elements(isample_rate) NE 1) OR (isample_rate LT 0)) THEN $
     message, 'ISAMPLE_RATE must be a scalar >= zero'
  IF ((n_elements(osample_rate) NE 1) OR (osample_rate LT 0)) THEN $
     message, 'OSAMPLE_RATE must be a scalar >= zero'
  IF ((osample_rate MOD isample_rate) NE 0) THEN $
     message, 'ISAMPLE_RATE must be an integer divisor of OSAMPLE_RATE'
  n_af=n_elements(angle_fields)
  n_mf=n_elements(magnitude_fields)
  IF n_af NE n_mf THEN $
     message, 'MAGNITUDE_FIELDS must be supplied along with ANGLE_FIELDS'
  IF n_af GT 0 THEN BEGIN
     af_size=size(angle_fields)
     IF (af_size[0] GT 1) && (af_size[2] NE 2) THEN $
        message, 'ANGLE_FIELDS must be an integer scalar or vector'
     mf_size=size(magnitude_fields)
     IF (mf_size[0] GT 1) && (mf_size[2] NE 2) THEN $
        message, 'MAGNITUDE_FIELDS must be an integer scalar or vector'
     IF mf_size[1] NE af_size[1] THEN $
        message, 'MAGNITUDE FIELDS and ANGLE_FIELDS must have the same ' + $
                 'number of elements'
  ENDIF                         ; else, we don't have angles/magnitudes
  idir_files=file_search(idir + path_sep() + '*', count=nidir_files, $
                         /nosort, /fold_case, /test_regular)
  IF nidir_files LT 1 THEN message, 'No input files found'

  restore, itemplate_sav
  field_names=strlowcase(itemplate.FIELDNAMES)
  field_types=itemplate.FIELDTYPES
  is_time_field=itemplate.FIELDGROUPS EQ time_idx
  ;; Ignore other groups when reading the data
  itemplate.FIELDGROUPS=indgen(itemplate.FIELDCOUNT)
  itemplate.FIELDGROUPS[where(is_time_field)]=time_idx
  non_time_fields=where(~is_time_field)
  non_time_field_names=field_names[non_time_fields]
  ;; Check which magnitude fields we have, since negative indices mean
  ;; magnitude is just scalar 1
  IF n_mf GT 0 THEN BEGIN       ; got angles/magnitudes?
     ;; Check if MAGNITUDE_FIELDS has any actual fields (mcount: how many)
     mag_yes=where(magnitude_fields GE 0, mcount)
     noang_cols_maybe=(mcount GT 0) ? $
                      cgSetDifference(non_time_fields, $
                                      [angle_fields, $
                                       magnitude_fields[mag_yes]], $
                                      count=nnoang) : $
                      cgSetDifference(non_time_fields, $
                                      angle_fields, count=nnoang)
     IF nnoang GT 0 THEN BEGIN
        ;; We got fields that are not angles/magnitudes and also not time
        ;; fields.  We only average those fields among these that are not
        ;; strings.  Becareful here.
        nostr=where(field_types[noang_cols_maybe] NE 7, navg)
        IF navg GT 0 THEN noang_cols=noang_cols_maybe[nostr]
        avg_cols=(navg GT 0) ? $
                 [noang_cols, angle_fields, magnitude_fields[mag_yes]] : $
                 [angle_fields, magnitude_fields[mag_yes]]
     ENDIF ELSE BEGIN
        ;; We now can say the non-time fields are all included in those
        ;; given as angles/magnitudes (Note cgSetDifference() returns the
        ;; first argument in this case)
        avg_cols=noang_cols_maybe
     ENDELSE
     ;; Now determine what the differences are
     noang=cgSetDifference(avg_cols, $
                           [angle_fields, magnitude_fields[mag_yes]], $
                           noresult=-1)
  ENDIF ELSE BEGIN              ; don't have any angles/magnitudes
     nostr=where(field_types[non_time_fields] NE 7, navg)
     IF navg GT 0 THEN BEGIN
        avg_cols=non_time_fields[nostr]
     ENDIF ELSE BEGIN
        message, 'Cannot average any non-time fields using' + itemplate_sav
     ENDELSE
  ENDELSE
  n_ifields=itemplate.FIELDCOUNT ; N fields in template
  tags2remove=where(field_names EQ field_names[time_idx])
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
  ncols_osample=long(osample_rate / isample_rate)
  nrows_osample=86400 / long(isample_rate * ncols_osample)

  FOR k=0, nidir_files - 1 DO BEGIN
     iname=strsplit(file_basename(idir_files[k]), '.', /extract)
     ;; Get a path for the file, check if it already exists
     ofile_name=strcompress(odir + path_sep() + iname[0] + '_' + $
                            strtrim(osample_rate, 2) + 's.' + iname[1], $
                            /remove_all)
     ofile_stamp=file_basename(ofile_name)
     out_list=file_search(odir + path_sep() + '*.' + iname[1], $
                          /nosort, /fold_case, /test_regular)
     matchfiles=where(ofile_stamp EQ file_basename(out_list), $
                      matchfilecount)
     IF matchfilecount GT 0 THEN BEGIN
        IF keyword_set(overwrite) THEN BEGIN
           message, 'Averaged file ' + ofile_stamp + $
                    ' already exists.  Overwriting', /informational
        ENDIF ELSE BEGIN
        message, 'Averaged file ' + ofile_stamp + $
                 ' already exists.  Not overwriting', /informational
        CONTINUE
     ENDELSE
     ENDIF

     ifile=idir_files[k]
     message, 'Producing ' + strtrim(osample_rate, 2) + 's average for ' + $
              ifile, /informational
     idata=read_ascii(ifile, template=itemplate)
     idata_names=strlowcase(tag_names(idata))
     time_loc=where(idata_names EQ field_names[time_idx])
     idata_times=idata.(time_loc)
     ;; Number of lines in input
     lines=n_elements(idata_times[0, *])
     ;; Extract times and remove quotes or spaces from strings
     IF size(idata_times, /type) EQ 7 THEN BEGIN
        FOREACH fld, indgen((size(idata_times, /dimensions))[0]) DO BEGIN
           ok=strsplit(idata_times[fld, *], '" ', /extract)
           idata_times[fld, *]=ok.toArray()
        ENDFOREACH
     ENDIF
     match2, idata_names, field_names[tags2remove], is_time
     FOREACH fld, (indgen(n_tags(idata)))[where(is_time LT 0)] DO BEGIN
        IF size(idata.(fld), /type) EQ 7 THEN BEGIN
           ok=strsplit(idata.(fld), '" ', /extract)
           idata.(fld)=ok.toArray()
        ENDIF
     ENDFOREACH

     ;; Obtain full time details array
     itimes_std=parse_times(idata_times, tnames_last, time_locs)
     beg_jd=julday(itimes_std[1, 0], $
                   itimes_std[2, 0], $
                   itimes_std[0, 0], 0)
     otimes=timegen(start=beg_jd, final=beg_jd + (86399.0 / 86400), $
                    step_size=osample_rate, units='seconds')
     otimes=jul2timestamp(otimes)
     ;; Set up hash that will contain output data
     ohash=hash(field_names)
     ;; Note we use our time_locs variable to locate each time info
     ohash[tnames[time_locs[0]]]=strmid(otimes, 0, 4)
     ohash[tnames[time_locs[1]]]=strmid(otimes, 5, 2)
     ohash[tnames[time_locs[2]]]=strmid(otimes, 8, 2)
     ohash[tnames[time_locs[3]]]=strmid(otimes, 11, 2)
     ohash[tnames[time_locs[4]]]=strmid(otimes, 14, 2)
     ohash[tnames[time_locs[5]]]=strmid(otimes, 17)
     FOREACH fld, non_time_field_names DO BEGIN
        fld_idx=where(field_names EQ fld)
        ifield_type=field_types[fld_idx]
        val=(ifield_type EQ 4) ? !VALUES.F_NAN : ''
        ohash[fld]=make_array(n_elements(otimes), type=ifield_type, $
                              value=val)
     ENDFOREACH
     ;; Calculate average for non-angular/magnitude fields
     IF (n_mf EQ 0) OR $
        ((n_elements(noang) GT 0) && noang[0] NE -1) THEN BEGIN
        flds=(n_elements(noang) EQ 0) ? avg_cols : noang
        FOREACH col, flds DO BEGIN
           match_fld=where(idata_names EQ field_names[col])
           col_idata=reform(idata.(match_fld), $
                            ncols_osample, nrows_osample)
           ;; This will throw "floating point illegal operand" arithmetic
           ;; error warning when all elements in each row are NaN, which
           ;; should should just be ignored.  In these cases, the result is
           ;; -NaN.
           avgs=mean(col_idata, dimension=1, /nan)
           ohash[field_names[col]]=avgs
        ENDFOREACH
     ENDIF
     ;; Calculate averages for angular/magnitude fields
     IF n_mf GT 0 THEN BEGIN
        FOREACH idx, indgen(n_elements(angle_fields)) DO BEGIN
           ang_name=field_names[angle_fields[idx]]
           ang=idata.(where(idata_names EQ ang_name))
           ang2d=reform(ang, ncols_osample, nrows_osample)
           ;; Add standard deviation field
           stdev_name=ang_name + '_stdev'
           ohash[stdev_name]=make_array(n_elements(otimes), $
                                        type=4, value=!VALUES.F_NAN)
           FOR i=0L, nrows_osample - 1 DO BEGIN
              ang_ok=where(finite(ang2d[*, i]) GT 0, n_ang_ok, $
                           complement=ang_bad)
              IF n_ang_ok GT 10 THEN $
                 ohash[stdev_name, i]=stddev(ang2d[ang_ok, i])
           ENDFOR
           IF (magnitude_fields[idx] GE 0) THEN BEGIN
              mag_name=field_names[magnitude_fields[idx]]
              mag_fld=where(idata_names EQ mag_name)
              mag=idata.(mag_fld)
              mag2d=reform(mag, ncols_osample, nrows_osample)
              avg=bearing_avg(ang2d, mag2d, dimension=1)
              ohash[ang_name]=avg[*, 0]
              ohash[mag_name]=avg[*, 1]
              ;; Add standard deviation field
              stdev_name=mag_name + '_stdev'
              ohash[stdev_name]=make_array(n_elements(otimes), $
                                           type=4, value=!VALUES.F_NAN)
              FOR i=0L, nrows_osample - 1 DO BEGIN
                 mag_ok=where(finite(mag2d[*, i]) GT 0, n_mag_ok, $
                              complement=mag_bad)
                 IF n_mag_ok GT 10 THEN $
                    ohash[stdev_name, i]=stddev(mag2d[mag_ok, i])
              ENDFOR
           ENDIF ELSE BEGIN
              avg=bearing_avg(ang2d, 1, dimension=1)
              ohash[ang_name]=avg[*, 0]
           ENDELSE
        ENDFOREACH
     ENDIF
     
     file_stamp=file_basename(ifile, '.dat') + '_' + $
                strtrim(osample_rate, 2) + 's.dat'

     IF n_mf EQ 0 THEN BEGIN    ; no angles/magnitudes
        okeys=[tnames, field_names[avg_cols]]
     ENDIF ELSE BEGIN
        okeys=(mcount GT 0) ? $ ; any magnitudes
              [tnames, field_names[avg_cols], $
               field_names[angle_fields] + '_stdev', $
               field_names[magnitude_fields[mag_yes]] + '_stdev'] : $
              [tnames, field_names[avg_cols], $
               field_names[angle_fields] + '_stdev']
     ENDELSE
     ts=create_struct(okeys[0], ohash[okeys[0]])
     FOREACH fld, okeys[1:*] DO BEGIN
        ts=create_struct(ts, okeys[where(okeys EQ fld)], ohash[fld])
     ENDFOREACH

     write_csv, odir + path_sep() + file_stamp, ts, $
                header=strlowcase(tag_names(ts))
  ENDFOR
        
END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; nav_avg.pro ends here
