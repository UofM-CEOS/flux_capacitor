;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-10-04T17:25:14+0000
;; Last-Updated: 2013-10-06T19:51:57+0000
;;           By: Sebastian P. Luque
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
                         /nosort)
  IF nidir_files LT 1 THEN $
     message, 'No input files found'

  restore, itemplate_sav
  field_names=itemplate.FIELDNAMES
  is_time_field=itemplate.FIELDGROUPS EQ time_beg_idx
  non_time_fields=where(~is_time_field)
  non_time_field_names=strlowcase(field_names[non_time_fields])
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

  ;; Loop through files in input directory
  FOR k=0, nidir_files - 1 DO BEGIN
     iname=strsplit(file_basename(idir_files[k]), '.', /extract)
     ;; Get a path for the file, check if it already exists
     ofile_name=strcompress(odir + path_sep() + iname[0] + '_std.' + $
                            iname[1], /remove_all)
     ofile_stamp=file_basename(ofile_name)
     out_list=file_search(odir + path_sep() + '*.' + iname[1], /nosort)
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
     time_loc=where(idata_names EQ $
                    strlowcase(field_names[time_beg_idx]))
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
     match2, idata_names, strlowcase(field_names[tags2remove]), is_time
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
     ;; Set up output hash
     
     FOREACH fld, field_names[angle_fields] DO BEGIN
        fld_idx=where(idata_names EQ fld)
        xs=sin(idata.(fld_idx) * !DTOR)
        ys=cos(idata.(fld_idx) * !DTOR)
        oxs=interpol(xs, itimes_jd, otimes_jd)
        oys=interpol(ys, itimes_jd, otimes_jd)
        oang=atan(oxs, oys) / !DTOR
        negs=where(oang LE 0, /null)
        oang[negs]=oang[negs] + 360
     ENDFOREACH
     odata=remove_structure_tags(idata, field_names[tags2remove])
     otimes=1
     ;; Find indices to keep
     match2, strlowcase(tag_names(odata)), keep_fields, toss, keep
     tags2remove_odata=where(toss LT 0, nremove)
     IF nremove GT 0 THEN $
        odata=remove_structure_tags(odata, $
                                    (tag_names(odata))[tags2remove_odata])
     odata=create_struct('year', reform(itimes_std[0, *]), $
                         'month', reform(itimes_std[1, *]), $
                         'day', reform(itimes_std[2, *]), $
                         'hour', reform(itimes_std[3, *]), $
                         'minute', reform(itimes_std[4, *]), $
                         'second', reform(itimes_std[5, *]), odata)
     delvar, idata

     ;; Fix things for some years

     ;; For 2011 RH data set to NaN when sensor was not working
     ;; from 0931 UTC on JD204 through 1435 on JD207. We're sure
     ;; we have DOY in these raw files, so no need to test.
     cal_badbeg2011=doy2calendar(2011, 204)
     cal_badend2011=doy2calendar(2011, 207)
     jd_badbeg2011=julday(strmid(cal_badbeg2011, 4, 2), $
                          strmid(cal_badbeg2011, 6, 2), 2011, 9, 31)
     jd_badend2011=julday(strmid(cal_badend2011, 4, 2), $
                          strmid(cal_badend2011, 6, 2), 2011, 14, 35)
     jd=julday(odata.month, odata.day, odata.year, odata.hour, odata.minute)
     bad2011=where((jd GE jd_badbeg2011) AND (jd LE jd_badend2011), nbad)
     IF nbad GT 0 THEN odata.(11)[bad2011]=!VALUES.F_NAN

     write_csv, ofile_name, odata, header=strlowcase(tag_names(odata))

  ENDFOR

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; filter_series.pro ends here
