;;; std_rmc.pro --- Standardize RMC files
;; Author: Sebastian Luque
;; Created: 2013-09-26T21:14:01+0000
;; Last-Updated: 2013-10-01T12:29:07+0000
;;           By: Sebastian P. Luque
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
;; STD_RMC, expand_path('~/tmp/ArcticNet2011/RMC'), $
;;          expand_path('~/tmp/ArcticNet2011/RMC/STD'), $
;;          'rmc_raw_template.sav', 2, rmc_std_header, /overwrite
;; 
;; - ----------------------------------------------------------------------
;;; Code:

PRO STD_RMC, IDIR, ODIR, ITEMPLATE_SAV, UTC_TIME_IDX, GPS_TIME_IDX, $
             OHEADER, OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 6) THEN $
     message, 'Usage: STD_MET, IDIR, ODIR, ITEMPLATE_SAV, ' + $
              'UTC_TIME_IDX, OHEADER'
  IF ((n_elements(idir) EQ 0) OR (idir EQ '')) THEN $
     message, 'IDIR is undefined or is empty string'
  IF ((n_elements(odir) EQ 0) OR (odir EQ '')) THEN $
     message, 'ODIR is undefined or is empty string'
  IF ((n_elements(itemplate_sav) EQ 0) OR (itemplate_sav EQ '')) THEN $
     message, 'ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(utc_time_idx) NE 1) OR (utc_time_idx LT 0)) THEN $
     message, 'UTC_TIME_IDX must be a scalar >= zero'
  IF ((n_elements(gps_time_idx) NE 1) OR (gps_time_idx LT 0)) THEN $
     message, 'GPS_TIME_IDX must be a scalar  >= zero'
  IF (n_elements(oheader) EQ 0) THEN $
     message, 'OHEADER is undefined'

  idir_files=file_search(idir + path_sep() + '*.raw', count=nidir_files, $
                         /nosort, /fold_case)
  IF nidir_files LT 1 THEN BEGIN
     print, 'No input files found'
     RETURN
  ENDIF
  restore, itemplate_sav
  field_names=itemplate.FIELDNAMES
  n_ifields=itemplate.FIELDCOUNT ; N fields in template
  header=strsplit(oheader, ', ', /extract) ; assuming ", "
  n_ofields=n_elements(header) ; N fields as in input header
  tfields=where(itemplate.FIELDGROUPS EQ utc_time_idx, /NULL) ; time fields
  tnames=strlowcase(field_names[tfields])            ; time names
  tnamesl=strsplit(time_names, '_', /extract)                 ; split list
  tnames_last=strarr(n_elements(tnamesl)) ; set up 
  FOR i=0L, n_elements(tnames) - 1 DO $
     tnames_last[i]=tnamesl[i, n_elements(tnamesl[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  time_locs=locate_time_strings(tnames_last)
  ;; GPS times
  tfields_gps=where(itemplate.FIELDGROUPS EQ gps_time_idx, /NULL)
  tnames_gps=strlowcase(field_names[tfields_gps])
  tnamesl_gps=strsplit(tnames_gps, '_', /extract)
  tnames_last_gps=strarr(n_elements(tnamesl_gps))
  FOR i=0L, n_elements(tnames_gps) - 1 DO $
     tnames_last_gps[i]=tnamesl_gps[i, n_elements(tnamesl_gps[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  time_locs_gps=locate_time_strings(tnames_last_gps)
     
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
     idata_names=strlowcase(tag_names(idata))
     idata_times=idata.(where(idata_names EQ $
                              strlowcase(field_names[utc_time_idx])))
     ;; Number of lines in input
     lines=n_elements(idata_times[0, *])
     ;; Remove quotes and separators
     IF size(idata_times, /type) EQ 7 THEN BEGIN
        FOREACH fld, indgen((size(idata_times, /dimensions))[0]) DO BEGIN
           ok=strsplit(idata_times[fld, *], '" -/:', /extract)
           ok=(temporary(ok)).toArray()
           ok=strjoin(transpose(temporary(ok)))
           idata_times[fld, *]=ok
        ENDFOREACH
     ENDIF
     ;; Obtain full UTC time details
     CASE tnames_last[time_locs[0]] OF
        'year': yyyy=string(idata_times[time_locs[0], *], $
                            format='(i04)')
        'yyyymmdd': yyyy=strmid(idata_times[time_locs[0], *], 0, 4)
        'mmddyyyy': yyyy=strmid(idata_times[time_locs[0], *], 4, 4)
        'ddmmyyyy': yyyy=strmid(idata_times[time_locs[0], *], 4, 4)
        'ddmmyy': BEGIN
           message, 'Assuming current century', /informational
           tstamp=jul2timestamp(systime(/julian))
           yyyy=strmid(tstamp, 0, 2) + $
                strmid(idata_times[time_locs[0], *], 4, 2)
        END
        ELSE: message, 'Do not know how to extract year from this field'
     ENDCASE
     CASE tnames_last[time_locs[1]] OF
        'month': mo=string(idata_times[time_locs[1], *], $
                           format='(i02)')
        'yyyymmdd': mo=strmid(idata_times[time_locs[1], *], 4, 2)
        'mmddyyyy': mo=strmid(idata_times[time_locs[1], *], 0, 2)
        'ddmmyyyy': mo=strmid(idata_times[time_locs[1], *], 2, 2)
        'ddmmyy': mo=strmid(idata_times[time_locs[1], *], 2, 2)
        'doy': BEGIN
           calendar=doy2calendar(yyyy, idata_times[time_locs[1], *])
           mo=strmid(calendar, 4, 2)
        END
        ELSE: message, 'Do not know how to extract month from this field'
     ENDCASE
     CASE tnames_last[time_locs[2]] OF
        'day': dd=string(idata_times[time_locs[2], *], $
                         format='(i02)')
        'yyyymmdd': dd=strmid(idata_times[time_locs[2], *], 6, 2)
        'mmddyyyy': dd=strmid(idata_times[time_locs[2], *], 2, 2)
        'ddmmyyyy': dd=strmid(idata_times[time_locs[2], *], 0, 2)
        'ddmmyy': dd=strmid(idata_times[time_locs[2], *], 0, 2)
        'doy': dd=strmid(calendar, 6, 2) ; we already have calendar
        ELSE: message, 'Do not know how to extract day from this field'
     ENDCASE
     CASE tnames_last[time_locs[3]] OF
        'hour': hh=string(idata_times[time_locs[3], *], $
                          format='(i02)')
        'hhmmss': hh=strmid(idata_times[time_locs[3], *], 0, 2)
        'hhmm': hh=strmid(idata_times[time_locs[3], *], 0, 2)
        ELSE: message, 'Do not know how to extract hour from this field'
     ENDCASE
     CASE tnames_last[time_locs[4]] OF
        'minute': mm=string(idata_times[time_locs[4], *], $
                            format='(i02)')
        'hhmmss': mm=strmid(idata_times[time_locs[4], *], 2, 2)
        'hhmm': mm=strmid(idata_times[time_locs[4], *], 2, 2)
        ELSE: message, 'Do not know how to extract minute from this field'
     ENDCASE
     CASE tnames_last[time_locs[5]] OF
        'second': ss=string(idata_times[time_locs[5], *], $
                            format='(i02)')
        'hhmmss': ss=strmid(idata_times[time_locs[5], *], 4, 2)
        ELSE: message, 'Do not know how to extract second from this field'
     ENDCASE

     ;; Obtain full GPS time details
     gps_time_loc=where(idata_names EQ $
                        strlowcase(field_names[gps_time_idx]))
     idata_times_gps=idata.(gps_time_loc)
     CASE tnames_last_gps[time_locs_gps[0]] OF
        'year': $
           yyyy_gps=string(idata_times_gps[time_locs_gps[0], *], $
                           format='(i04)')
        'yyyymmdd': $
           yyyy_gps=strmid(idata_times_gps[time_locs_gps[0], *], $
                           0, 4)
        'mmddyyyy': $
           yyyy_gps=strmid(idata_times_gps[time_locs_gps[0], *], $
                           4, 4)
        'ddmmyyyy': $
           yyyy_gps=strmid(idata_times_gps[time_locs_gps[0], *], $
                           4, 4)
        'ddmmyy': BEGIN
           message, 'Assuming current century', /informational
           tstamp=jul2timestamp(systime(/julian))
           yyyy_gps=strmid(tstamp, 0, 2) + $
                    strmid(idata_times_gps[time_locs_gps[0], *], 4, 2)
        END
        ELSE: message, 'Do not know how to extract year from this field'
     ENDCASE
     CASE tnames_last_gps[time_locs_gps[1]] OF
        'month': mo_gps=string(idata_times_gps[time_locs_gps[1], *], $
                               format='(i02)')
        'yyyymmdd': $
           mo_gps=strmid(idata_times_gps[time_locs_gps[1], *], 4, 2)
        'mmddyyyy': $
           mo_gps=strmid(idata_times_gps[time_locs_gps[1], *], 0, 2)
        'ddmmyyyy': $
           mo_gps=strmid(idata_times_gps[time_locs_gps[1], *], 2, 2)
        'ddmmyy': $
           mo_gps=strmid(idata_times_gps[time_locs_gps[1], *], 2, 2)
        'doy': BEGIN
           calendar_gps=doy2calendar(yyyy, $
                                     idata_times_gps[time_locs_gps[1], *])
           mo_gps=strmid(calendar_gps, 4, 2)
        END
        ELSE: message, 'Do not know how to extract month from this field'
     ENDCASE
     CASE tnames_last_gps[time_locs[2]] OF
        'day': $
           dd_gps=string(idata_times_gps[time_locs_gps[2], *], $
                         format='(i02)')
        'yyyymmdd': $
           dd_gps=strmid(idata_times_gps[time_locs_gps[2], *], $
                         6, 2)
        'mmddyyyy': $
           dd_gps=strmid(idata_times_gps[time_locs_gps[2], *], $
                         2, 2)
        'ddmmyyyy': $
           dd_gps=strmid(idata_times_gps[time_locs_gps[2], *], $
                         0, 2)
        'ddmmyy': $
           dd_gps=strmid(idata_times_gps[time_locs_gps[2], *], $
                         0, 2)
        'doy': $
           dd_gps=strmid(calendar_gps, 6, 2) ; we already have calendar_gps
        ELSE: message, 'Do not know how to extract day from this field'
     ENDCASE
     CASE tnames_last_gps[time_locs_gps[3]] OF
        'hour': hh_gps=string(idata_times_gps[time_locs_gps[3], *], $
                              format='(i02)')
        'hhmmss': $
           hh_gps=strmid(idata_times_gps[time_locs_gps[3], *], 0, 2)
        'hhmm': $
           hh_gps=strmid(idata_times_gps[time_locs_gps[3], *], 0, 2)
        ELSE: message, 'Do not know how to extract hour from this field'
     ENDCASE
     CASE tnames_last_gps[time_locs_gps[4]] OF
        'minute': mm_gps=string(idata_times_gps[time_locs_gps[4], *], $
                                format='(i02)')
        'hhmmss': $
           mm_gps=strmid(idata_times_gps[time_locs_gps[4], *], 2, 2)
        'hhmm': $
           mm_gps=strmid(idata_times_gps[time_locs_gps[4], *], 2, 2)
        ELSE: message, 'Do not know how to extract minute from this field'
     ENDCASE
     CASE tnames_last_gps[time_locs_gps[5]] OF
        'second': ss_gps=string(idata_times_gps[time_locs_gps[5], *], $
                                format='(i02)')
        'hhmmss': $
           ss_gps=strmid(idata_times_gps[time_locs_gps[5], *], 4, 2)
        ELSE: message, 'Do not know how to extract second from this field'
     ENDCASE
     
     ;; keep=
     odata=remove_structure_tag(idata, (tag_names(idata))[utc_time_idx])
     IF time_locs[6] LT 0 THEN BEGIN
        odata=create_struct('year', reform(yyyy), 'month', reform(mo), $
                            'day', reform(dd), 'hour', reform(hh), $
                            'minute', reform(mm), 'second', reform(ss), $
                            odata)
     ENDIF ELSE BEGIN
        ds=idata_times[time_locs[6], *]
        odata=create_struct('year', reform(yyyy), 'month', reform(mo), $
                            'day', reform(dd), 'hour', reform(hh), $
                            'minute', reform(mm), 'second', reform(ss), $
                            'second_fraction', reform(ds), odata)
     ENDELSE
     delvar, idata

     ;; Loop through records (lines)
     FOR exd=0L, lines[0] - 1  DO BEGIN
        ;; Fix things depending on year
        SWITCH yyyy[exd] OF
           '2011': BEGIN
              ;; For 2011 RH data set to NaN when sensor was not working
              ;; from 0931 UTC on JD204 through 1435 on JD207. We're sure
              ;; we have DOY in these raw files, so no need to test.
              cal_bad_beg=doy2calendar(2011, 204)
              cal_bad_end=doy2calendar(2011, 207)
              jd_badbeg=julday(strmid(cal_bad_beg, 4, 2), $
                               strmid(cal_bad_beg, 6, 2), 2011, 9, 31)
              jd_badend=julday(strmid(cal_bad_end, 4, 2), $
                               strmid(cal_bad_end, 6, 2), 2011, 14, 35)
              jd_now=julday(mo[exd], dd[exd], 2011, hh[exd], mm[exd])
              IF (jd_now GE jd_badbeg) AND $
                 (jd_now LE jd_badend) THEN BEGIN
                 odata.(11)[exd]=!VALUES.F_NAN
              ENDIF
           END
        ENDSWITCH
     ENDFOR

     write_csv, ofile_name, odata, header=header

  ENDFOR

END



;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; std_rmc.pro ends here
