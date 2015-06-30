;; Author: Sebastian Luque
;; Created: 2013-09-20T17:13:48+0000
;; Last-Updated: 2015-06-30T19:53:28+0000
;;           By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     STD_EC
;; 
;; PURPOSE:
;; 
;;     Standardize flux files from various instruments and generate a
;;     standardized file with a fixed structure (number of fields and field
;;     order).  It is very similar to std.pro, so not sure how useful this
;;     will be in the long run.
;; 
;; CALLING SEQUENCE:
;; 
;;     STD_EC, Idir, Odir, Time_Beg_Idx, Keep_Fields, Std_Fields, $
;;             Keep_types=Keep_Types
;; 
;; INPUTS:
;; 
;;     Idir:                  Input directory (no trailing separator).
;;     Odir:                  Output directory (no trailing separator).
;;     Itemplate_Sav:         Ascii template to read input files.
;;     Time_Beg_Idx:          Index (in template) where time is.
;;     Std_Fields:            String array or scalar with names of fields
;;                            to output into standard file.
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
;; ec_41_42_raw_keep_fields=['prog_version', 'accel_x', 'accel_y', 'accel_z', $
;;                           'rate_x', 'rate_y', 'rate_z', 'u_wind_analogue', $
;;                           'v_wind_analogue', 'w_wind_analogue', $
;;                           'sonic_temperature_analogue', 'CO2_op', $
;;                           'H2O_op', 'pressure_op', 'diag_op', 'CO2_cl', $
;;                           'H2O_cl', 'pressure_cl', 'temperature_cl', $
;;                           'diag_cl', 'u_wind_serial', 'v_wind_serial', $
;;                           'w_wind_serial', 'sonic_temperature_serial', $
;;                           'SoS_serial', 'gillstat_serial', $
;;                           'pressure_tower', 'temperature_base', $
;;                           'temperature_spar', 'temperature_bulb', $
;;                           'LI7000_mode', 'CO2_opsh', 'H2O_opsh', $
;;                           'pressure_opsh', 'diag_opsh', 'CO2_dry', $
;;                           'H2O_dry', 'CO2_LI7200', 'H2O_LI7200', $
;;                           'AGC_LI7200', 'pressure_LI7200', $
;;                           'temperature_in_LI7200', 'temperature_out_LI7200', $
;;                           'temperature_avg_LI7200', 'diag_LI7200']
;; ec_std_fields=['prog_version', 'accel_x', 'accel_y', 'accel_z', $
;;                'rate_x', 'rate_y', 'rate_z', 'u_wind_analogue', $
;;                'v_wind_analogue', 'w_wind_analogue', $
;;                'sonic_temperature_analogue', 'CO2_op', 'H2O_op', $
;;                'pressure_op', 'temperature_op', 'diag_op', 'CO2_cl', $
;;                'H2O_cl', 'pressure_cl', 'temperature_cl', 'diag_cl', $
;;                'u_wind_serial', 'v_wind_serial', 'w_wind_serial', $
;;                'sonic_temperature_serial', 'SoS_serial', $
;;                'gillstat_serial', 'pressure_tower', 'LI7000_mode', $
;;                'CO2_opsh', 'H2O_opsh', 'pressure_opsh', $
;;                'temperature_opsh', 'diag_opsh', 'CO2_dry', 'H2O_dry', $
;;                'CO2_LI7200', 'H2O_LI7200', 'AGC_LI7200', $
;;                'pressure_LI7200', 'temperature_in_LI7200', $
;;                'temperature_out_LI7200', 'temperature_avg_LI7200', $
;;                'diag_LI7200', 'temperature_base', 'temperature_spar', $
;;                'temperature_bulb']
;; STD_EC, expand_path('~/tmp/ArcticNet2011/Flux/V4.1_4.2'), $
;;         expand_path('~/tmp/ArcticNet2011/Flux/STD'), $
;;         'ec_41_42_template.sav', 1, ec_41_42_raw_keep_fields, $
;;         ec_std_fields, /overwrite
;; 
;;- -----------------------------------------------------------------------
;;; Code:

PRO STD_EC, IDIR, ODIR, ITEMPLATE_SAV, TIME_BEG_IDX, KEEP_FIELDS, $
            STD_FIELDS, STD_TYPES, KEEP_TYPES=KEEP_TYPES, OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 7) THEN $
     message, 'Usage: STD_EC, IDIR, ODIR, ITEMPLATE_SAV, ' + $
              'TIME_BEG_IDX, KEEP_FIELDS, STD_FIELDS, STD_TYPES'
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
  IF n_kf EQ 0 THEN message, 'KEEP_FIELDS is undefined'
  n_kt=n_elements(keep_types)
  IF (n_kt GT 0) AND (n_kt NE n_kf) THEN $
     message, 'KEEP_TYPES and KEEP_FIELDS must have the same number of elements'
  n_sf=n_elements(std_fields)
  IF n_sf EQ 0 THEN message, 'STD_FIELDS is undefined'
  n_st=n_elements(std_types)
  IF (n_st GT 0) AND (n_st NE n_sf) THEN $
     message, 'STD_TYPES and STD_FIELDS must have the same number of elements'
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
     onames=strlowcase(tag_names(odata))
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

     delvar, idata
     ;; Set up output hash
     ohash=hash(strlowcase(std_fields))
     FOREACH value, ohash, fld DO BEGIN
        match_fld=where(onames EQ fld, nmatch_fld)
        IF nmatch_fld GT 0 THEN BEGIN
           ohash[fld]=odata.(where(onames EQ fld))
        ENDIF ELSE BEGIN 
           match_type=std_types[where(std_fields EQ fld)]
           ohash[fld]=make_array(lines, type=match_type, $
                                 value=match_type EQ 7 ? $
                                 '' : !VALUES.F_NAN)
        ENDELSE
     ENDFOREACH
     ofield_names=[tnames_id + '_year', tnames_id + '_month', $
                   tnames_id + '_day', tnames_id + '_hour', $
                   tnames_id + '_minute', tnames_id + '_second', $
                   strlowcase(std_fields)]
     ;; Add standard times
     ohash[ofield_names[0]]=reform(itimes_std[0, *])
     ohash[ofield_names[1]]=reform(itimes_std[1, *])
     ohash[ofield_names[2]]=reform(itimes_std[2, *])
     ohash[ofield_names[3]]=reform(itimes_std[3, *])
     ohash[ofield_names[4]]=reform(itimes_std[4, *])
     ohash[ofield_names[5]]=reform(itimes_std[5, *])

     ;; Order output
     odata=create_struct(ofield_names[0], ohash[ofield_names[0]])
     FOREACH fld, ofield_names[1:*] DO BEGIN
        odata=create_struct(odata, $
                            ofield_names[where(ofield_names EQ fld)], $
                            ohash[fld])
     ENDFOREACH
     delvar, ohash

     write_csv, ofile_name, odata, header=ofield_names

  ENDFOR

END



;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
