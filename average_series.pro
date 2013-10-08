;; $Id$
;; Author: Bruce Johnson, Sebastian Luque
;; Created: 2013-09-17T14:59:07+0000
;; Last-Updated: 2013-10-08T18:29:30+0000
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
;;- -----------------------------------------------------------------------
;;; Code:

PRO AVERAGE_SERIES, IDIR, ODIR, ITEMPLATE_SAV, TIME_BEG_IDX, $
                    ISAMPLE_RATE, OSAMPLE_RATE, STAMP, $
                    ANGLE_FIELDS=ANGLE_FIELDS, $
                    MAGNITUDE_FIELDS=MAGNITUDE_FIELDS, $
                    OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 7) THEN $
     message, 'Usage: AVERAGE_SERIES, IDIR, ODIR, ITEMPLATE_SAV, ' + $
              'TIME_BEG_IDX, ISAMPLE_RATE, OSAMPLE_RATE, ' + $
              'ANGLE_FIELDS, MAGNITUDE_FIELDS, STAMP'
  IF ((n_elements(idir) EQ 0) OR (idir EQ '')) THEN $
     message, 'IDIR is undefined or is empty string'
  IF ((n_elements(odir) EQ 0) OR (odir EQ '')) THEN $
     message, 'ODIR is undefined or is empty string'
  IF ((n_elements(itemplate_sav) EQ 0) OR (itemplate_sav EQ '')) THEN $
     message, 'ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(time_beg_idx) NE 1) OR $
      ((size(time_beg_idx, /type) NE 2) || time_beg_idx LT 0)) THEN $
         message, 'TIME_BEG_IDX must be an integer scalar >= zero'
  IF ((n_elements(isample_rate) EQ 0) OR (isample_rate EQ '')) THEN $
     message, 'ISAMPLE_RATE is undefined or is empty string'
  IF ((n_elements(osample_rate) EQ 0) OR (osample_rate EQ '')) THEN $
     message, 'OSAMPLE_RATE is undefined or is empty string'
  IF ((osample_rate MOD isample_rate) NE 0) THEN $
     message, 'ISAMPLE_RATE must be an integer divisor of OSAMPLE_RATE'
  IF ((n_elements(stamp) EQ 0) OR (stamp EQ '')) THEN $
     message, 'STAMP is undefined or is empty string'
  n_af=n_elements(angle_fields)
  n_mf=n_elements(magnitude_fields)
  IF (n_af GT 0 AND n_mf EQ 0) || (n_mf GT 0 AND n_af EQ 0) THEN $
     message, 'MAGNITUDE_FIELDS must be supplied along with ANGLE_FIELDS'
  IF (n_af GT 0 AND n_mf GT 0) THEN BEGIN
     af_size=size(angle_fields)
     IF (af_size[0] GT 1) || (af_size[2] NE 2) THEN $
        message, 'ANGLE_FIELDS must be an integer scalar or vector'
     mf_size=size(magnitude_fields)
     IF (mf_size[0] GT 1) || (mf_size[2] NE 2) THEN $
        message, 'MAGNITUDE_FIELDS must be an integer scalar or vector'
     IF mf_size[1] NE af_size[1] THEN $
        message, 'MAGNITUDE FIELDS and ANGLE_FIELDS must have the same ' + $
                 'number of elements'
  ENDIF                         ; else, we don't have angles/magnitudes
  idir_files=file_search(idir + path_sep() + '*', count=nidir_files, $
                         /nosort, /fold_case, /test_regular)
  IF nidir_files LT 1 THEN BEGIN
     message, 'No input files found', /informational
     RETURN
  ENDIF

  restore, itemplate_sav
  field_names=itemplate.FIELDNAMES
  field_types=itemplate.FIELDTYPES
  is_time_field=itemplate.FIELDGROUPS EQ time_beg_idx
  non_time_fields=where(~is_time_field)
  non_time_field_names=strlowcase(field_names[non_time_fields])
  ;; Locate non-angular fields; check if any via nnoang
  noang_cols_maybe=n_af GT 0 || n_mf GT 0 ? $
                   cgSetDifference(non_time_fields, $
                                   [angle_fields, magnitude_fields], $
                                   count=nnoang) : $
                   non_time_fields
  ;; Only interpolate those fields that are not strings (BECAREFUL HERE)
  noang_cols=noang_cols_maybe[where(field_types[noang_cols_maybe] NE 7)]
  n_ifields=itemplate.FIELDCOUNT ; N fields in template
  tags2remove=where(field_names EQ field_names[time_beg_idx])
  ;; Times
  tfields=where(is_time_field, /NULL)
  tnames=strlowcase(field_names[tfields])
  tnamesl=strsplit(tnames, '_', /extract)
  tnames_last=strarr(n_elements(tnamesl))
  tnames_id=strjoin((tnamesl[0])[0:n_elements(tnamesl[0]) - 2], '_')
  FOR i=0L, n_elements(tnames) - 1 DO $
     tnames_last[i]=tnamesl[i, n_elements(tnamesl[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  time_locs=locate_time_strings(tnames_last)

  ;; ang_cols=[bmag_field, bear_field, head_field]
  ;; stdev_str=['sog_stdev', 'cog_stdev', 'heading_stdev']
  ;; avg_cols_maybe=(cgSetDifference(non_time_fields, ang_cols))
  ;; ;; Only average those fields that are not strings (BECAREFUL HERE)
  ;; avg_cols=avg_cols_maybe[where(itemplate.FIELDTYPES[avg_cols_maybe] NE 7)]
  ;; n_avg_cols=size(avg_cols, /n_elements)
  ;; ncols_osample=osample_rate / isample_rate
  ;; nrows_osample=86400 / osample_rate
  ;; oheader=[[time_field_names, header[avg_cols], stdev_str]]
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
     idata=read_ascii(ifile, count=n_inputfile, template=itemplate)
     idata_names=strlowcase(tag_names(idata))
     time_loc=where(idata_names EQ $
                    strlowcase(field_names[time_beg_idx]))
     idata_times=idata.(time_loc)

     ;; Extract times and remove quotes or spaces from strings
     idata_times=idata.(time_beg_idx)
     idata=remove_structure_tags(idata, (tag_names(idata))[time_beg_idx])
     idata_names=strlowcase(tag_names(idata))
     IF size(idata_times, /type) EQ 7 THEN BEGIN
        FOREACH fld, indgen((size(idata_times, /dimensions))[0]) DO BEGIN
           ok=strsplit(idata_times[fld, *], '" ', /extract)
           idata_times[fld, *]=ok.toArray()
        ENDFOREACH
     ENDIF
     FOREACH fld, indgen(n_tags(idata)) DO BEGIN
        IF size(idata.(fld), /type) EQ 7 THEN BEGIN
           ok=strsplit(idata.(fld), '" ', /extract)
           idata.(fld)=ok.toArray()
        ENDIF
     ENDFOREACH
     beg_jd=julday(idata_times[month_subfield, 0], $
                   idata_times[day_subfield, 0], $
                   idata_times[year_subfield, 0], 0)
     otimes=timegen(start=beg_jd, final=beg_jd + (86399.0 / 86400), $
                    step_size=osample_rate, units='seconds')
     otimes=jul2timestamp(otimes)
     ;; Set up hash that will contain output data
     ohash=hash(oheader)
     FOREACH value, ohash, fld DO BEGIN
        CASE fld OF
           'year': ohash[fld]=strmid(otimes, 0, 4)
           'month': ohash[fld]=strmid(otimes, 5, 2)
           'day': ohash[fld]=strmid(otimes, 8, 2)
           'hour': ohash[fld]=strmid(otimes, 11, 2)
           'minute': ohash[fld]=strmid(otimes, 14, 2)
           'second': ohash[fld]=strmid(otimes, 17, 2)
           ELSE: BEGIN
              fld_idx=where(itemplate.FIELDNAMES EQ fld)
              ifield_type=itemplate.FIELDTYPES[fld_idx]
              val=(ifield_type EQ 4) ? !VALUES.F_NAN : ''
              ohash[fld]=make_array(n_elements(otimes), type=ifield_type, $
                                    value=val)
           END
        ENDCASE
     ENDFOREACH

     FOREACH col, avg_cols DO BEGIN
        col_idata=reform(idata.(where(idata_names EQ header[col])), $
                         ncols_osample, nrows_osample)
        ;; This will throw "floating point illegal operand" arithmetic
        ;; error warning when all elements in each row are NaN, which
        ;; should should just be ignored.  In these cases, the result is
        ;; -NaN.
        avgs=mean(col_idata, dimension=1, /nan)
        ohash[header[col]]=avgs
     ENDFOREACH
     ;; Means for SOG, COG, and heading
     cog=idata.(where(idata_names EQ header[bear_field]))
     brg2d=reform(cog, ncols_osample, nrows_osample)
     sog=idata.(where(idata_names EQ header[bmag_field]))
     bmag2d=reform(sog, ncols_osample, nrows_osample)
     head=idata.(where(idata_names EQ header[head_field]))
     head2d=reform(head, ncols_osample, nrows_osample)
     brgmag_mean=bearing_avg(brg2d, bmag2d, dimension=1)
     ohash[header[bmag_field]]=brgmag_mean[*, 1]
     ohash[header[bear_field]]=brgmag_mean[*, 0]
     head_mean=bearing_avg(head2d, 1, dimension=1)
     ohash[header[head_field]]=head_mean[*, 0]
     ;; Standard deviations
     FOR i=0, nrows_osample - 1 DO BEGIN
        sog_ok=where(finite(bmag2d[*, i]) GT 0, n_sog_ok, complement=sog_bad)
        cog_ok=where(finite(brg2d[*, i]) GT 0, n_cog_ok, complement=cog_bad)
        head_ok=where(finite(head2d[*, i]) GT 0, n_head_ok, $
                      complement=head_bad)
        IF n_sog_ok GT 10 THEN $
           ohash[stdev_str[0], i]=stddev(bmag2d[sog_ok, i])
        IF n_cog_ok GT 10 THEN $
           ohash[stdev_str[1], i]=stddev_yamartino(brg2d[cog_ok, i])
        IF n_head_ok GT 10 THEN $
           ohash[stdev_str[2], i]=stddev_yamartino(head2d[head_ok, i])
     ENDFOR
     
     file_stamp=file_basename(ifile, '.dat') + '_' + $
                strtrim(osample_rate, 2) + 's.dat'
     ts=create_struct(oheader[0], ohash[oheader[0]])
     FOREACH fld, oheader[1:*] DO BEGIN
        ts=create_struct(ts, oheader[where(oheader EQ fld)], ohash[fld])
     ENDFOREACH
     write_csv, odir + path_sep() + file_stamp, ts, header=oheader
  ENDFOR
        
END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; nav_avg.pro ends here
