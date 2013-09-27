;;; std_rmc.pro --- Standardize RMC files
;; Author: Sebastian Luque
;; Created: 2013-09-26T21:14:01+0000
;; Last-Updated: 2013-09-26T22:42:52+0000
;;           By: Sebastian Luque
;; ------------------------------------------------------------------------
;;; Commentary: 
;; 
;;
;; + ----------------------------------------------------------------------
;; NAME:
;; 
;; std_rmc
;; 
;; PURPOSE:
;; 
;;  Standardize RMC files.  It will likely standardize other files in the
;;  future.
;; 
;; CATEGORY:
;; 
;;  General Input/Output
;; 
;; CALLING SEQUENCE:
;; 
;;  STD_RMC, Idir, Odir, Itemplate_Sav, Time_Beg_Idx, Oheader
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
;; rmc_std_header='year, month, day, hour, minute, second, ' + $
;;                'prog_version, battery_volt, panel_temperature, ' + $
;;                'pressure, air_temperature, RH_percent, ' + $
;;                'surface_temperature, wind_speed, wind_direction, ' + $
;;                'wind_sd, PAR, pitch, roll'
;; STD_RMC, expand_path('~/tmp/RMC'), expand_path('~/tmp/RMC/STD'), $
;;          'rmc_raw_template.sav', 1, rmc_std_header
;; 
;; - ----------------------------------------------------------------------
;;; Code:

PRO STD_RMC, IDIR, ODIR, ITEMPLATE_SAV, TIME_BEG_IDX, OHEADER, $
             OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 5) THEN $
     message, 'Usage: STD_RMC, IDIR, ODIR, ITEMPLATE_SAV, ' + $
              'TIME_BEG_IDX, OHEADER'
  IF ((n_elements(idir) EQ 0) OR (idir EQ '')) THEN $
     message, 'IDIR is undefined or is empty string'
  IF ((n_elements(odir) EQ 0) OR (odir EQ '')) THEN $
     message, 'ODIR is undefined or is empty string'
  IF ((n_elements(itemplate_sav) EQ 0) OR (itemplate_sav EQ '')) THEN $
     message, 'ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(time_beg_idx) NE 1) OR (time_beg_idx LT 0)) THEN $
     message, 'TIME_BEG_IDX must be a scalar >= zero'
  IF (n_elements(oheader) EQ 0) THEN $
     message, 'OHEADER is undefined'

  idir_files=file_search(idir + path_sep() + '*.raw', count=nidir_files, $
                         /nosort, /fold_case)
  IF nidir_files LT 1 THEN BEGIN
     print, 'No input files found'
     RETURN
  ENDIF
  restore, itemplate_sav
  n_ifields=itemplate.FIELDCOUNT ; N fields in template
  header=strsplit(oheader, ', ', /extract) ; assuming ", "
  n_ofields=n_elements(header) ; N fields as in input header
  time_fields=where(itemplate.FIELDGROUPS EQ time_beg_idx, /NULL)
  time_field_names=strlowcase(itemplate.FIELDNAMES[time_fields])
  yyyymmdd_subfield=where(time_field_names EQ 'yyyymmdd') ; merged?
  mmddyyyy_subfield=where(time_field_names EQ 'mmddyyyy') ; merged?
  year_subfield=where(time_field_names EQ 'year') ; locate year
  month_subfield=where(time_field_names EQ 'month') ; locate month
  day_subfield=where(time_field_names EQ 'day') ; locate day
  hhmmss_subfield=where(time_field_names EQ 'hhmmss') ; merged?
  hhmm_subfield=where(time_field_names EQ 'hhmm') ; merged?
  hour_subfield=where(time_field_names EQ 'hour') ; locate hour
  minute_subfield=where(time_field_names EQ 'minute') ; locate minute
  second_subfield=where(time_field_names EQ 'second') ; locate second
  ;; Check if we have a name with "second" as substring somewhere after the
  ;; first character; matches e.g.: "decisecond", "millisecond"
  subsecond_exists=strpos(time_field_names, 'second', 1)
  subsecond_subfield=where(subsecond_exists GE 0)
  doy_subfield=where((time_field_names EQ 'doy') OR $
                     (time_field_names EQ 'julday')) ; locate DOY
  IF (year_subfield LT 0) AND $
     ((yyyymmdd_subfield LT 0) AND (mmddyyyy_subfield LT 0)) THEN BEGIN
     message, 'Cannot determine year using this ITEMPLATE, ' + $
              'where a merged date variables is also absent'
  ENDIF
  IF (month_subfield LT 0) AND $
     ((yyyymmdd_subfield LT 0) AND (mmddyyyy_subfield LT 0) AND $
      (doy_subfield LT 0)) THEN BEGIN
     message, 'Cannot determine month using this ITEMPLATE, ' + $
              'where merged date variables or DOY is also absent'
  ENDIF
  IF (day_subfield LT 0) AND $
     ((yyyymmdd_subfield LT 0) AND (mmddyyyy_subfield LT 0) AND $
      (doy_subfield LT 0)) THEN BEGIN
     message, 'Cannot determine day using this ITEMPLATE, ' + $
              'where merged date variables or DOY is also absent'
  ENDIF
  IF ((hour_subfield LT 0) AND $
      (hhmm_subfield LT 0) AND (hhmmss_subfield LT 0)) THEN BEGIN
     message, 'Cannot determine hour using this ITEMPLATE, ' + $
              'where merged hour-minute, or hour-minute-second ' + $
              'is also absent'
  ENDIF
  IF ((minute_subfield LT 0) AND $
      (hhmm_subfield LT 0) AND (hhmmss_subfield LT 0)) THEN BEGIN
     message, 'Cannot determine minute using this ITEMPLATE, ' + $
              'where merged hour-minute, or hour-minute-second ' + $
              'is also absent'
  ENDIF
  IF ((second_subfield LT 0) AND $
      (hhmm_subfield LT 0) AND (hhmmss_subfield LT 0)) THEN BEGIN
     message, 'Cannot determine second using this ITEMPLATE, ' + $
              'where merged hour-minute, or hour-minute-second ' + $
              'is also absent'
  ENDIF

  ;; Loop through files in input directory
  FOR k=0, nidir_files - 1 DO BEGIN
     iname=strsplit(file_basename(idir_files[k]), '.', /extract)
     ;; Get a path for the file, check if it already exists
     ofile_name=strcompress(odir + path_sep() + iname[0] + '_std.' + $
                            iname[1], /remove_all)
     ofile_stamp=file_basename(ofile_name)
     out_list=file_search(odir + path_sep() + '*.dat', /nosort)
     matchfiles=where(ofile_stamp EQ file_basename(out_list), $
                      matchfilecount)
     IF matchfilecount GT 0 THEN BEGIN
        IF keyword_set(overwrite) THEN BEGIN
           print, 'Standardized MET file ' + ofile_stamp + $
                  ' already exists.  Overwriting'
        ENDIF ELSE BEGIN
        print, 'Standardized MET file ' + ofile_stamp + $
               ' already exists.  Not overwriting'
        CONTINUE
     ENDELSE
     ENDIF

     ifile=idir_files[k]
     print, 'Processing file: ' + ifile

     ;; Read input file
     idata=read_ascii(ifile, template=itemplate)
     ;; Number of lines in input
     lines=n_elements(idata.(time_beg_idx)[0, *])
     ;; Obtain full time details
     yyyy=reform(idata.(time_beg_idx)[year_subfield, *])
     IF yyyymmdd_subfield GE 0 THEN BEGIN
        yyyymmdd_str=string(idata.(time_beg_idx)[yyyymmdd_subfield, *])
        yyyymmdd_str=strsplit(temporary(yyyymmdd_str), '- ', /extract)
        yyyymmdd_str=transpose(yyyymmdd_str.toArray())
        mo=string(yyyymmdd_str[1, *], format='(i02)')
        dd=string(yyyymmdd_str[2, *], format='(i02)')
     ENDIF ELSE BEGIN
        IF month_subfield GE 0 THEN BEGIN 
           mo=reform(idata.(time_beg_idx)[month_subfield, *])
        ENDIF ELSE BEGIN
           doy=idata.(time_beg_idx)[doy_subfield, *]
           calstr=doy2calendar(yyyy, doy)
           mo=reform(strmid(calstr, 4, 2))
           dd=reform(strmid(calstr, 6, 2))
        ENDELSE
        IF day_subfield GE 0 THEN BEGIN
           dd=reform(idata.(time_beg_idx)[month_subfield, *])
        ENDIF
     ENDELSE
     IF hhmmss_subfield GE 0 THEN BEGIN
        hhmmss_str=string(idata.(time_beg_idx)[hhmmss_subfield, *])
        hhmmss_str=strsplit(temporary(hhmmss_str), ':- ', /extract)
        hhmmss_str=transpose(hhmmss_str.toArray())
        hh=string(hhmmss_str[0, *], format='(i02)')
        mm=string(hhmmss_str[1, *], format='(i02)')
        ss=string(hhmmss_str[2, *], format='(i02)')
     ENDIF ELSE BEGIN
        IF hhmm_subfield GE 0 THEN BEGIN
           hhmm_str=string(idata.(time_beg_idx)[hhmm_subfield, *])
           hhmm_str=strsplit(temporary(hhmm_str), ':- ', /extract)
           hhmm_str=transpose(hhmm_str.toArray())
           hh=string(hhmm_str[0, *], format='(i02)')
           mm=string(hhmm_str[1, *], format='(i02)')
        ENDIF ELSE BEGIN
           IF hour_subfield GE 0 THEN BEGIN
              hh=reform(idata.(time_beg_idx)[hour_subfield, *])
           ENDIF ELSE BEGIN
              hh=reform('00', lines)
           ENDELSE
           IF minute_subfield GE 0 THEN BEGIN
              mm=reform(idata.(time_beg_idx)[minute_subfield, *])
           ENDIF ELSE BEGIN
              mm=reform('00', lines)
           ENDELSE
           IF second_subfield GE 0 THEN BEGIN
              ss=reform(idata.(time_beg_idx)[second_subfield, *])
           ENDIF ELSE BEGIN
              ss=reform('00', lines)
           ENDELSE
        ENDELSE
     ENDELSE
     
     odata=remove_structure_tag(idata, $
                                reform((tag_names(idata))[0:time_beg_idx]))
     IF subsecond_subfield LT 0 THEN BEGIN
        odata=create_struct('year', yyyy, 'month', mo, 'day', dd, $
                            'hour', hh, 'minute', mm, 'second', ss, odata)
     ENDIF ELSE BEGIN
        ds=idata.(time_beg_idx)[subsecond_subfield, *]
        odata=create_struct('year', yyyy, 'month', mo, 'day', dd, $
                            'hour', hh, 'minute', mm, 'second', ss, $
                            'second_fraction', ds, odata)
     ENDELSE
     delvar, idata

     ;; ;; Loop through records (lines)
     ;; FOR exd=0L, lines[0] - 1  DO BEGIN
     ;;    ;; Fix things depending on year
     ;;    SWITCH yyyy[exd] OF
     ;;       '2011': BEGIN
     ;;          ;; For 2011 RH data set to NaN when sensor was not working
     ;;          ;; from 0931 UTC on JD204 through 1435 on JD207. We're sure
     ;;          ;; we have DOY in these raw files, so no need to test.
     ;;          cal_bad_beg=doy2calendar(2011, 204)
     ;;          cal_bad_end=doy2calendar(2011, 207)
     ;;          jd_badbeg=julday(strmid(cal_bad_beg, 4, 2), $
     ;;                           strmid(cal_bad_beg, 6, 2), 2011, 9, 31)
     ;;          jd_badend=julday(strmid(cal_bad_end, 4, 2), $
     ;;                           strmid(cal_bad_end, 6, 2), 2011, 14, 35)
     ;;          jd_now=julday(mo[exd], dd[exd], 2011, hh[exd], mm[exd])
     ;;          IF (jd_now GE jd_badbeg) AND $
     ;;             (jd_now LE jd_badend) THEN BEGIN
     ;;             odata.(11)[exd]=!VALUES.F_NAN
     ;;          ENDIF
     ;;       END
     ;;    ENDSWITCH
     ;; ENDFOR

     write_csv, ofile_name, odata, header=header

  ENDFOR

END



;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; std_rmc.pro ends here
