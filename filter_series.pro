;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-10-04T17:25:14+0000
;; Last-Updated: 2013-10-26T09:41:10+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     FILTER_SERIES
;; 
;; PURPOSE:
;; 
;;     This procedure regularizes a time series, i.e. generates a time
;;     series with a regular time step from input having irregular time
;;     steps.  For instance, gyro files are nominally sampling at 1 Hz, but
;;     not exactly.  If a regular time series is required for further
;;     calculations, then this procedure can be used to generate such a
;;     time series.  Interpolation is performed to estimate values of
;;     angular and non-angular fields in the input files.
;; 
;; CALLING SEQUENCE:
;; 
;;     FILTER_SERIES, Idir, Odir, Itemplate_Sav, Time_Beg_Idx, $
;;     Angle_Fields, Sample_Rate
;; 
;; INPUTS:
;; 
;;     Idir:               Input directory (no trailing separator).
;;     Odir:               Output directory (no trailing separator).
;;     Itemplate_Sav:      Ascii template to read input files.
;;     Time_Beg_Idx:       Index (in template) where time is.
;;     Angle_Fields:       Integer array with the index of angular fields.
;;     Sample_Rate:        Scalar indicating the frequency (s) that
;;                         output data should have.
;; 
;; KEYWORD PARAMETERS:
;; 
;;     OVERWRITE:             Whether to overwrite files in Odir.
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
;;     FILTER_SERIES, expand_path('~/tmp/ArcticNet2011/GYRO/STD'), $
;;                    expand_path('~/tmp/ArcticNet2011/GYRO/1s'), $
;;                    'gyro_std_template.sav', 0, 6, 1L, /overwrite
;; 
;;- -----------------------------------------------------------------------
;;; Code:

PRO FILTER_SERIES, IDIR, ODIR, ITEMPLATE_SAV, TIME_BEG_IDX, ANGLE_FIELDS, $
                   SAMPLE_RATE, OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 6) THEN $
     message, 'Usage: FILTER_SERIES, IDIR, ODIR, ITEMPLATE_SAV, ' + $
              'TIME_BEG_IDX, ANGLE_FIELDS, SAMPLE_RATE'
  IF ((n_elements(idir) EQ 0) OR (idir EQ '')) THEN $
     message, 'IDIR is undefined or is empty string'
  IF ((n_elements(odir) EQ 0) OR (odir EQ '')) THEN $
     message, 'ODIR is undefined or is empty string'
  IF ((n_elements(itemplate_sav) EQ 0) OR (itemplate_sav EQ '')) THEN $
     message, 'ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(time_beg_idx) NE 1) OR (time_beg_idx LT 0)) THEN $
     message, 'TIME_BEG_IDX must be a scalar >= zero'
  IF (n_elements(angle_fields) EQ 0) THEN $
     message, 'ANGLE_FIELDS is undefined'
  IF ((n_elements(sample_rate) NE 1) OR (sample_rate LT 0)) THEN $
     message, 'SAMPLE_RATE must be a a scalar >= zero'
  idir_files=file_search(idir + path_sep() + '*', count=nidir_files, $
                         /nosort, /fold_case, /test_regular)
  IF nidir_files LT 1 THEN message, 'No input files found'

  restore, itemplate_sav
  field_names=strlowcase(itemplate.FIELDNAMES)
  field_types=itemplate.FIELDTYPES
  is_time_field=itemplate.FIELDGROUPS EQ time_beg_idx
  non_time_fields=where(~is_time_field)
  non_time_field_names=field_names[non_time_fields]
  ;; Locate non-angular fields; check if any via nnoang
  noang_cols_maybe=cgSetDifference(non_time_fields, angle_fields, $
                                   count=nnoang)
  ;; Only interpolate those fields that are not strings (BECAREFUL HERE)
  noang_cols=noang_cols_maybe[where(field_types[noang_cols_maybe] NE 7)]
  n_ifields=itemplate.FIELDCOUNT ; N fields in template
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

  ;; Loop through files in input directory
  FOR k=0, nidir_files - 1 DO BEGIN
     iname=strsplit(file_basename(idir_files[k]), '.', /extract)
     ;; Get a path for the file, check if it already exists
     ofile_name=strcompress(odir + path_sep() + iname[0] + '_' + $
                            strtrim(sample_rate, 2) + 's.' + iname[1], $
                            /remove_all)
     ofile_stamp=file_basename(ofile_name)
     out_list=file_search(odir + path_sep() + '*.' + iname[1], $
                          /nosort, /fold_case, /test_regular)
     matchfiles=where(ofile_stamp EQ file_basename(out_list), $
                      matchfilecount)
     IF matchfilecount GT 0 THEN BEGIN
        IF keyword_set(overwrite) THEN BEGIN
           message, 'Standardized file ' + ofile_stamp + $
                    ' already exists.  Overwriting', /informational
        ENDIF ELSE BEGIN
        message, 'Standardized file ' + ofile_stamp + $
                 ' already exists.  Not overwriting', /informational
        CONTINUE
     ENDELSE
     ENDIF

     ifile=idir_files[k]
     message, 'Processing ' + ifile, /informational
     ;; Read input file
     idata=read_ascii(ifile, template=itemplate)
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

     ;; Obtain full time details array
     itimes_std=parse_times(idata_times, tnames_last, time_locs)
     itimes_jd=reform(julday(itimes_std[1, *], $
                             itimes_std[2, *], $
                             itimes_std[0, *], $
                             itimes_std[3, *], $
                             itimes_std[4, *], $
                             itimes_std[5, *]))
     beg_jd=julday(itimes_std[1, 0], itimes_std[2, 0], itimes_std[0, 0], $
                   itimes_std[3, 0], itimes_std[4, 0], $
                   floor(fix(itimes_std[5, 0])))
     end_jd=julday(itimes_std[1, lines - 1], $
                   itimes_std[2, lines - 1], $
                   itimes_std[0, lines - 1], $
                   itimes_std[3, lines - 1], $
                   itimes_std[4, lines - 1], $
                   floor(fix(itimes_std[5, lines - 1])))
     otimes_jd=timegen(start=beg_jd, final=end_jd, $
                       step_size=sample_rate, units='seconds')
     caldat, otimes_jd, mo, dd, yyyy, hh, mm, ss

     ;; Set up output hash
     ohash=hash(field_names)
     ;; Note we use our time_locs variable to locate each time info
     ohash[tnames[time_locs[0]]]=string(yyyy, format='(i04)')
     ohash[tnames[time_locs[1]]]=string(mo, format='(i02)')
     ohash[tnames[time_locs[2]]]=string(dd, format='(i02)')
     ohash[tnames[time_locs[3]]]=string(hh, format='(i02)')
     ohash[tnames[time_locs[4]]]=string(mm, format='(i02)')
     ;; Becareful here with the decimal places formatting
     ohash[tnames[time_locs[5]]]=string(ss, format='(f06.3)')
     FOREACH fld, field_names[angle_fields] DO BEGIN
        fld_idx=where(idata_names EQ fld)
        xs=sin(idata.(fld_idx) * !DTOR)
        ys=cos(idata.(fld_idx) * !DTOR)
        oxs=interpol(xs, itimes_jd, otimes_jd)
        oys=interpol(ys, itimes_jd, otimes_jd)
        oang=atan(oxs, oys) / !DTOR
        negs=where(oang LE 0, nneg)
        IF nneg GT 0 THEN $
           oang[negs]=oang[negs] + 360
        ohash[fld]=oang
     ENDFOREACH
     IF nnoang GT 0 THEN BEGIN
        FOREACH fld, field_names[noang_fields] DO BEGIN
           fld_idx=where(idata_names EQ fld)
           noang=interpol(idata.(fld_idx), itimes_jd, otimes_jd)
           ohash[fld]=noang
        ENDFOREACH
     ENDIF
     delvar, idata

     ts=create_struct(field_names[0], ohash[field_names[0]])
     FOREACH fld, field_names[1:*] DO BEGIN
        ts=create_struct(ts, $
                         field_names[where(field_names EQ fld)], $
                         ohash[fld])
     ENDFOREACH
     write_csv, ofile_name, ts, header=strlowcase(tag_names(ts))

  ENDFOR

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; filter_series.pro ends here
