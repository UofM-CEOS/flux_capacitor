;; $Id$
;; Author: Sebastian Luque (from original by Brent Else, Bruce Johnson)
;; Created: 2013-09-20T17:13:48+0000
;; Last-Updated: 2013-12-02T02:24:49+0000
;;           By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     STD
;; 
;; PURPOSE:
;; 
;;     Standardize files with a single set of time stamps (e.g. UTC).
;; 
;; CALLING SEQUENCE:
;; 
;;     STD, Idir, Odir, Itemplate_Sav, Time_Beg_Idx, Keep_Fields, $
;;          Keep_types=Keep_Types, File_type=File_Type
;; 
;; INPUTS:
;; 
;;     Idir:                  Input directory (no trailing separator).
;;     Odir:                  Output directory (no trailing separator).
;;     Itemplate_Sav:         Ascii template to read input files.
;;     Time_Beg_Idx:          Index (in template) where time is.
;;     Keep_Fields:           String array or scalar with names of fields
;;                            to keep.
;; 
;; KEYWORD PARAMETERS:
;; 
;;     KEEP_TYPES:            Integer array with IDL data type codes to
;;                            convert the Keep_Fields to.
;;     FILE_TYPE:             String scalar specifying input file type.
;;     OVERWRITE:             Whether to overwrite files in Odir.
;; 
;; SIDE EFFECTS:
;; 
;; Writes files in Odir.
;; 
;; PROCEDURE:
;;
;;     The File_Type argument is used to perform corrections and
;;     manipulations depending on the file type.  This is performed after
;;     the standardization process.  Check the last case statement, where
;;     we may add clauses for different years/subsets.
;; 
;; EXAMPLE:
;; 
;;     met_raw_keep_fields=['prog_version', 'battery_voltage', $
;;                          'panel_temperature', 'pressure', $
;;                          'air_temperature', 'rh_percent', $
;;                          'surface_temperature', 'wind_speed', $
;;                          'wind_direction', 'wind_sd', 'par', 'pitch',
;;                          'roll']
;;     STD, expand_path('~/tmp/ArcticNet2011/MET'), $
;;          expand_path('~/tmp/ArcticNet2011/MET/STD'), $
;;          'met_raw_template.sav', 1, met_raw_keep_fields,
;;          file_type='MET', /OVERWRITE
;; 
;;- -----------------------------------------------------------------------
;;; Code:

PRO STD, IDIR, ODIR, ITEMPLATE_SAV, TIME_BEG_IDX, KEEP_FIELDS, $
         KEEP_TYPES=KEEP_TYPES, FILE_TYPE=FILE_TYPE, OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 5) THEN $
     message, 'Usage: STD, IDIR, ODIR, ITEMPLATE_SAV, ' + $
              'TIME_BEG_IDX, KEEP_FIELDS'
  idir_info=file_info(idir)
  itpl_info=file_info(itemplate_sav)
  IF (~idir_info.directory) THEN $
     message, 'IDIR must be a string pointing to an existing directory'
  IF (~itpl_info.read) THEN $
     message, 'ITEMPLATE_SAV must be a string pointing to a readable file'
  IF ((n_elements(time_beg_idx) NE 1) OR $
      ((size(time_beg_idx, /type) NE 2) || time_beg_idx LT 0)) THEN $
         message, 'TIME_BEG_IDX must be an integer scalar >= zero'
  IF ((n_elements(odir) EQ 0) OR (odir EQ '')) THEN $
     message, 'ODIR is undefined or is empty string'
  n_kf=n_elements(keep_fields)
  IF (n_kf EQ 0) THEN message, 'KEEP_FIELDS is undefined'
  n_kt=n_elements(keep_types)
  IF (n_kt GT 0) AND (n_kt NE n_kf) THEN $
     message, 'KEEP_TYPES and KEEP_FIELDS must have the same number of elements'
  IF ((n_elements(file_type) EQ 0) OR (file_type EQ '')) THEN $
     message, 'FILE_TYPE is undefined or is empty string'
  idir_files=file_search(idir + path_sep() + '*', count=nidir_files, $
                         /nosort, /fold_case, /test_regular)
  IF nidir_files LT 1 THEN message, 'No input files found'

  restore, itemplate_sav
  field_names=strlowcase(itemplate.FIELDNAMES)
  n_ifields=itemplate.FIELDCOUNT ; N fields in template
  tags2remove=where(field_names EQ field_names[time_beg_idx])
  ;; Times
  tfields=where(itemplate.FIELDGROUPS EQ time_beg_idx, /NULL)
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
     ofile_name=strcompress(odir + path_sep() + iname[0] + '_std.' + $
                            iname[1], /remove_all)
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
     times_dims=size(idata_times, /dimensions)
     valid_flag=make_array(times_dims[1], type=2, value=1)
     ;; Number of lines in input
     lines=n_elements(idata_times[0, *])
     ;; Remove quotes and separators. Also set invalid numbers to empty
     ;; string to avoid errors in parse_times()
     IF size(idata_times, /type) EQ 7 THEN BEGIN
        FOREACH fld, indgen((size(idata_times, /dimensions))[0]) DO BEGIN
           ok=strsplit(idata_times[fld, *], '" -/:', /extract)
           ok=(temporary(ok)).toArray()
           ok=strjoin(transpose(temporary(ok)))
           idata_times[fld, *]=ok
           is_valid=valid_num(idata_times[fld, *])
           ibad=where(~is_valid, bcount)
           IF bcount GT 0 THEN BEGIN
              idata_times[fld, ibad]=''
              valid_flag[ibad]=0 ; bad times -> invalid record
           ENDIF
        ENDFOREACH
     ENDIF
     match2, idata_names, field_names[tags2remove], is_time
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
     odata=remove_structure_tags(idata, field_names[tags2remove])
     ;; Find indices to keep
     match2, strlowcase(tag_names(odata)), strlowcase(keep_fields), $
             toss, keep
     tags2remove_odata=where(toss LT 0, nremove)
     IF nremove GT 0 THEN $
        odata=remove_structure_tags(odata, $
                                    (tag_names(odata))[tags2remove_odata])
     onames=tag_names(odata)
     ;; Re-check validity and subset
     ibad=where(~valid_flag, bcount, complement=ok, ncomplement=nok)
     IF nok LT 1 THEN BEGIN
        message, 'All time stamps are invalid.  Skipping this file', $
                 /continue
        CONTINUE
     ENDIF
     IF bcount GT 0 THEN BEGIN
        ohash=hash(odata)
        FOREACH value, ohash, fld DO BEGIN
           ohash[fld]=ohash[fld, ok]
        ENDFOREACH
        odata=create_struct(onames[0], ohash[onames[0]])
        FOREACH fld, onames[1:*] DO BEGIN
           odata=create_struct(odata, onames[where(onames EQ fld)], $
                               ohash[fld])
        ENDFOREACH
        idata_times=idata_times[*, ok]
     ENDIF
     ;; Type conversion, if requested
     IF n_kt GT 0 THEN BEGIN
        ohash=hash(odata)
        FOREACH value, ohash, fld DO BEGIN
           match_type=keep_types[where(onames EQ fld)]
           ohash[fld]=fix(ohash[fld], type=match_type)
        ENDFOREACH
        odata=create_struct(onames[0], ohash[onames[0]])
        FOREACH fld, onames[1:*] DO BEGIN
           odata=create_struct(odata, onames[where(onames EQ fld)], $
                               ohash[fld])
        ENDFOREACH
     ENDIF
     ;; Obtain full time details array
     itimes_std=parse_times(idata_times, tnames_last, time_locs)

     odata=create_struct(tnames_id + '_year', reform(itimes_std[0, *]), $
                         tnames_id + '_month', reform(itimes_std[1, *]), $
                         tnames_id + '_day', reform(itimes_std[2, *]), $
                         tnames_id + '_hour', reform(itimes_std[3, *]), $
                         tnames_id + '_minute', reform(itimes_std[4, *]), $
                         tnames_id + '_second', reform(itimes_std[5, *]), $
                         odata)
     delvar, idata

     ;; Fix things for certain file types
     CASE file_type OF
        'MET': BEGIN
           ;; Fix things for some years

           ;; For 2011 RH data set to NaN when sensor was not working
           ;; from 0931 UTC on JD204 through 1435 on JD207. We're sure
           ;; we have DOY in these raw files, so no need to test.
           cal_badbeg2011=doy2calendar(2011, 204)
           cal_badend2011=doy2calendar(2011, 207)
           jd_badbeg2011=julday(strmid(cal_badbeg2011, 4, 2), $
                                strmid(cal_badbeg2011, 6, 2), $
                                2011, 9, 31)
           jd_badend2011=julday(strmid(cal_badend2011, 4, 2), $
                                strmid(cal_badend2011, 6, 2), $
                                2011, 14, 35)
           jd=julday(odata.(1), odata.(2), odata.(0), odata.(3), odata.(4))
           bad2011=where((jd GE jd_badbeg2011) AND (jd LE jd_badend2011), $
                         nbad)
           IF nbad GT 0 THEN odata.(11)[bad2011]=!VALUES.F_NAN
        END
        'RAD': BEGIN
           ;; Silly time stamp
           badhour=where(odata.(3) EQ '24', nbad)
           IF nbad GT 0 THEN BEGIN
              odata.(3)[badhour]='00'
              jd_new=julday(odata.(1)[badhour], odata.(2)[badhour], $
                            odata.(0)[badhour], odata.(3)[badhour], $
                            odata.(4)[badhour], odata.(5)[badhour]) + 1
              caldat, jd_new, mo, dd
              odata.(2)[badhour]=string(dd, format='(i02)')
           ENDIF
           ;; Original comment: no data for all days before Julian Day 213
           ;; (Aug 1), UVS-AB-T sensor not installed. [SPL: I am adding a
           ;; minimum date, as too careless otherwise.  Assuming DOY 91
           ;; (2011-04-01).  Also, the code used two separate tests: one
           ;; for all dates prior to "Julian Day" - really DOY - 213 and
           ;; another for all dates prior to DOY 213 at 15:15, which should
           ;; really be a single one: i.e. the last one.  The variables
           ;; modified correspond to: temperature_uv, uv_b, uv_a,
           ;; temperature_uv_sd, uv_b_sd, and uv_a_sd.
           cal_badbeg2011=doy2calendar(2011, 91)
           cal_badend2011=doy2calendar(2011, 213)
           jd_badbeg2011=julday(strmid(cal_badbeg2011, 4, 2), $
                                strmid(cal_badbeg2011, 6, 2), $
                                strmid(cal_badbeg2011, 0, 4))
           jd_badend2011=julday(strmid(cal_badend2011, 4, 2), $
                                strmid(cal_badend2011, 6, 2), $
                                2011, 15, 15)
           jd=julday(odata.(1), odata.(2), odata.(0), odata.(3), odata.(4))
           bad2011=where((jd GE jd_badbeg2011) AND (jd LE jd_badend2011), $
                         nbad)
           IF nbad GT 0 THEN BEGIN
              FOREACH fld, [14, 15, 16, 24, 25, 26] DO BEGIN
                 match_fld=where(tag_names(odata) EQ $
                                 strupcase(field_names[fld]))
                 odata.(match_fld)[bad2011]=!VALUES.F_NAN
              ENDFOREACH
           ENDIF
        END
        ELSE: message, 'No further processing for ' + $
                       file_type, /informational
     ENDCASE

     write_csv, ofile_name, odata, header=strlowcase(tag_names(odata))

  ENDFOR

END



;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; std.pro ends here
