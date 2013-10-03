;; $Id$
;; Author: Brent Else, Bruce Johnson, Sebastian Luque
;; Created: 2013-09-20T17:13:48+0000
;; Last-Updated: 2013-10-03T20:31:22+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;; std_met
;; 
;; PURPOSE:
;; 
;;  Standardize MET files.  It will likely standardize other files in the
;;  future.
;; 
;; CATEGORY:
;; 
;;  General Input/Output
;; 
;; CALLING SEQUENCE:
;; 
;;  STD_MET, Idir, Odir, Itemplate_Sav, Time_Beg_Idx, Oheader
;; 
;; INPUTS:
;; 
;;  Idir: Input directory path, without trailing path separator.
;;  
;;  Odir: Output directory path, without trailing path separator.
;;  
;;  Itemplate_Sav: ascii_template template. Must contain a template
;;   structure named 'itemplate'.
;;  
;;  Time_Beg_Idx: Index (0-based) of first field with time definition.
;;   This is typically the year field.
;;
;;  Oheader: Output header string with each field name, separated by
;;   commas.
;; 
;; KEYWORD PARAMETERS:
;; 
;;  Overwrite: Whether to overwrite files in Odir.
;; 
;; SIDE EFFECTS:
;; 
;; Writes files in Odir.
;; 
;; EXAMPLE:
;; 
;; met_raw_keep_fields=['prog_version', 'battery_voltage', 'panel_temperature', $
;;                      'pressure', 'air_temperature', 'rh_percent', $
;;                      'surface_temperature', 'wind_speed', $
;;                      'wind_direction', 'wind_sd', 'par', 'pitch', 'roll']
;; STD_MET, expand_path('~/tmp/ArcticNet2011/MET'), $
;;          expand_path('~/tmp/ArcticNet2011/MET/STD'), $
;;          'met_raw_template.sav', 1, met_raw_keep_fields, /overwrite
;; 
;;- -----------------------------------------------------------------------
;;; Code:

PRO STD_MET, IDIR, ODIR, ITEMPLATE_SAV, TIME_BEG_IDX, KEEP_FIELDS, $
             OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 5) THEN $
     message, 'Usage: STD_MET, IDIR, ODIR, ITEMPLATE_SAV, ' + $
              'TIME_BEG_IDX, KEEP_FIELDS'
  IF ((n_elements(idir) EQ 0) OR (idir EQ '')) THEN $
     message, 'IDIR is undefined or is empty string'
  IF ((n_elements(odir) EQ 0) OR (odir EQ '')) THEN $
     message, 'ODIR is undefined or is empty string'
  IF ((n_elements(itemplate_sav) EQ 0) OR (itemplate_sav EQ '')) THEN $
     message, 'ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(time_beg_idx) NE 1) OR (time_beg_idx LT 0)) THEN $
     message, 'TIME_BEG_IDX must be a scalar >= zero'
  IF (n_elements(keep_fields) EQ 0) THEN $
     message, 'KEEP_FIELDS is undefined'
  idir_files=file_search(idir + path_sep() + '*.dat', count=nidir_files, $
                         /nosort)
  IF nidir_files LT 1 THEN $
     message, 'No input files found'

  restore, itemplate_sav
  field_names=itemplate.FIELDNAMES
  n_ifields=itemplate.FIELDCOUNT ; N fields in template
  tags2remove=where(field_names EQ field_names[time_beg_idx])
  ;; Times
  tfields=where(itemplate.FIELDGROUPS EQ time_beg_idx, /NULL)
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
     
     odata=remove_structure_tags(idata, field_names[tags2remove])
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

;; Below was moved from std_MET.pro, so should be added to looping through
;; records section above.

;;     ;Some extra processing for the RAD data
;;     ;--------------------------------------
;;     ;converting time of 2400 to time of 0000, and adding 1 to Julian Day
;;     IF stamp EQ 'RAD' THEN BEGIN
;;       tfclock = where(long(data.(3)[*]) EQ 2400, tfclock_count)
;;       IF tfclock_count GT 0 THEN BEGIN
;;         data.(3)[tfclock] = '0'
;;         data.(2)[tfclock] = data.(2)[tfclock]+1
;;      ENDIF
      
;;       ;for 2011 substituting 'NaN' when the UVS-AB-T sensor was not installed
;;       ;no data for all days before Julian Day 213 (Aug 1)
;;       no_uv_sensor = where((data.(1)[*] EQ 2011) AND (data.(2)[*] LT 213), no_uv_count)
;;       IF no_uv_count GT 0 THEN BEGIN
;;         data.(14)[no_uv_sensor] = 'NaN'
;;         data.(15)[no_uv_sensor] = 'NaN'
;;         data.(16)[no_uv_sensor] = 'NaN'
;;         data.(24)[no_uv_sensor] = 'NaN'
;;         data.(25)[no_uv_sensor] = 'NaN'
;;         data.(26)[no_uv_sensor] = 'NaN'
;;      ENDIF 
      
;;       ;no data for all times before 1515 UTC on Julian Day 213 (Aug 1)
;;       no_uv_sensor = where((data.(1)[*] EQ 2011) AND (data.(2)[*] EQ 213) AND (long(data.(3)[*] LT 1515)), no_uv_count)
;;       IF no_uv_count GT 0 THEN BEGIN
;;         data.(14)[no_uv_sensor] = 'NaN'
;;         data.(15)[no_uv_sensor] = 'NaN'
;;         data.(16)[no_uv_sensor] = 'NaN'
;;         data.(24)[no_uv_sensor] = 'NaN'
;;         data.(25)[no_uv_sensor] = 'NaN'
;;         data.(26)[no_uv_sensor] = 'NaN'
;;      ENDIF
;;    ENDIF 

;;     ;decode our date/time information
;;     input_hour = lonarr(1,n_elements(data.(3)[*]))
;;     input_min  = input_hour
;;     input_0    = where(long(data.(3)[*]) LT 60, match0)
;;     input_1    = where((long(data.(3)[*]) GT 59) AND (long(data.(3)[*] LT 1000)),match1)
;;     input_2    = where(long(data.(3)[*]) GT 999,match2)

;;     IF match0 GT 0 THEN BEGIN
;;       input_min(input_0) = long(data.(3)[input_0])
;;    ENDIF

;;     IF match1 GT 0 THEN BEGIN
;;       input_hour(input_1) = long(strmid(data.(3)[input_1],0,1))
;;       input_min(input_1)  = long(strmid(data.(3)[input_1],1,2))
;;    ENDIF

;;     IF match2 GT 0 THEN BEGIN
;;       input_hour(input_2) = long(strmid(data.(3)[input_2],0,2))
;;       input_min(input_2)  = long(strmid(data.(3)[input_2],2,2))
;;    ENDIF

;;     JD_arr=fltarr(1,n_recs)
;;     JD_arr(0,*)=current_JD



;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; std_met.pro ends here
