;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-10-01T20:08:28+0000
;; Last-Updated: 2013-12-02T02:46:42+0000
;;           By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;;
;;     STD_GYRO
;;
;; PURPOSE:
;;
;;     Standardize Gyro (ship) files, correcting time based on offset from
;;     corresponding RMC file.
;;
;; CALLING SEQUENCE:
;;
;;     STD_GYRO, Idir, Odir, Itemplate_Sav, ' + $
;;               Server_Time_Idx, Rmc_Std_Dir, Rmc_Std_Itemplate_Sav, ' + $
;;               Rmc_Utc_Time_Idx, Rmc_Server_Time_Idx, Keep_Fields
;;
;; INPUTS:
;;
;;     Idir:                  Input directory (no trailing separator).
;;     Odir:                  Output directory (no trailing separator).
;;     Itemplate_Sav:         Ascii template to read input files.
;;     Server_Time_Idx:       Index (in template) where server time is.
;;     Rmc_Std_Dir:           Directory of matching RMC files.
;;     Rmc_Std_Itemplate_Sav: Ascii template for reading RMC files.
;;     Rmc_Utc_Time_Idx:      Index (in template) where UTC time is.
;;     Rmc_Server_Time_Idx:   Index (in template) where server time is.
;;     Keep_Fields:           String array or scalar with names of fields
;;                            to keep.
;;
;; KEYWORD PARAMETERS:
;;
;;     OVERWRITE:             Whether to overwrite files in Odir.
;;
;; SIDE EFFECTS:
;;
;;     Writes files to Odir.
;;
;; EXAMPLE:
;;
;;     STD_GYRO, expand_path('~/tmp/ArcticNet2011/GYRO'), $
;;               expand_path('~/tmp/ArcticNet2011/GYRO/STD'), $
;;               'gyro_raw_template.sav', 0, $
;;               expand_path('~/tmp/ArcticNet2011/RMC/STD'), $
;;               'rmc_std_template.sav', 0, 6, gyro_raw_keep_fields,
;;               /OVERWRITE 
;;
;;- -----------------------------------------------------------------------
;;; Code:

PRO STD_GYRO, IDIR, ODIR, ITEMPLATE_SAV, SERVER_TIME_IDX, RMC_STD_DIR, $
              RMC_STD_ITEMPLATE_SAV, RMC_UTC_TIME_IDX, RMC_SERVER_TIME_IDX, $
              KEEP_FIELDS, OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 9) THEN $
     message, 'Usage: STD_GYRO, IDIR, ODIR, ITEMPLATE_SAV, ' + $
              'SERVER_TIME_IDX, RMC_STD_DIR, RMC_STD_ITEMPLATE_SAV, ' + $
              'RMC_UTC_TIME_IDX, RMC_SERVER_TIME_IDX, KEEP_FIELDS'
  idir_info=file_info(idir)
  itpl_info=file_info(itemplate_sav)
  rmc_dir_info=file_info(rmc_std_dir)
  rmc_tpl_info=file_info(rmc_std_itemplate_sav)
  IF (~idir_info.directory) THEN $
     message, 'IDIR must be a string pointing to an existing directory'
  IF (~itpl_info.read) THEN $
     message, 'ITEMPLATE_SAV must be a string pointing to a readable file'
  IF ((n_elements(server_time_idx) NE 1) OR $
      ((size(server_time_idx, /type) NE 2) || server_time_idx LT 0)) THEN $
         message, 'SERVER_TIME_IDX must be an integer scalar >= zero'
  IF ((n_elements(odir) EQ 0) OR (odir EQ '')) THEN $
     message, 'ODIR is undefined or is empty string'
  IF (~rmc_dir_info.directory) THEN $
     message, 'RMC_STD_DIR must be a string pointing to an existing directory'
  IF (~rmc_tpl_info.read) THEN $
     message, 'RMC_STD_ITEMPLATE_SAV must be a string pointing to a ' + $
              'readable file'
  IF ((n_elements(rmc_utc_time_idx) NE 1) OR $
      ((size(rmc_utc_time_idx, /type) NE 2) || $
       rmc_utc_time_idx LT 0)) THEN $
         message, 'RMC_UTC_TIME_IDX must be an integer scalar >= zero'
  IF ((n_elements(rmc_server_time_idx) NE 1) OR $
      ((size(rmc_server_time_idx, /type) NE 2) || $
       rmc_server_time_idx LT 0)) THEN $
         message, 'RMC_SERVER_TIME_IDX must be an integer scalar >= zero'
  IF (n_elements(keep_fields) EQ 0) THEN $
     message, 'KEEP_FIELDS is undefined'
  idir_files=file_search(idir + path_sep() + '*', count=nidir_files, $
                         /nosort, /fold_case, /test_regular)
  IF nidir_files LT 1 THEN message, 'No input files found'

  rmc_std_files=file_search(rmc_std_dir + path_sep() + '*', $
                            count=nrmc_files, /nosort, /fold_case, $
                            /test_regular)
  IF nrmc_files LT 1 THEN BEGIN
     message, 'No standard RMC files found, so ' + $
              'cannot correct GYRO files.  Exiting'
  ENDIF

  restore, itemplate_sav
  ;; We make a copy, since the RMC template is also called itemplate
  gyro_template=itemplate
  field_names=strlowcase(gyro_template.FIELDNAMES)
  n_ifields=gyro_template.FIELDCOUNT ; N fields in template
  tags2remove=where(field_names EQ field_names[server_time_idx])
  ;; Server times
  tfields_srv=where(gyro_template.FIELDGROUPS EQ server_time_idx, /NULL)
  tnames_srv=field_names[tfields_srv]
  tnamesl_srv=strsplit(tnames_srv, '_', /extract)
  tnames_last_srv=strarr(n_elements(tnamesl_srv))
  tnames_id_srv=strjoin((tnamesl_srv[0])[0:n_elements(tnamesl_srv[0]) - 2], '_')
  FOR i=0L, n_elements(tnames_srv) - 1 DO $
     tnames_last_srv[i]=tnamesl_srv[i, n_elements(tnamesl_srv[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  time_locs_srv=locate_time_strings(tnames_last_srv)
  ;; Break RMC file names and extract the piece to match
  rmc_filesl=strsplit(rmc_std_files, '_.', /extract)
  rmc_files_a=rmc_filesl.toArray(/transpose) ; array with pieces in cols
  rmc_files_mstr_dims=size(rmc_files_a, /dimensions)
  ;; We want the 2nd to last piece
  rmc_files_mstr=rmc_files_a[rmc_files_mstr_dims[0] - 3, *]
  
  ;; RMC template
  restore, rmc_std_itemplate_sav
  ;; We make a copy, since the RMC template is also called itemplate
  rmc_template=itemplate
  rmc_field_names=strlowcase(rmc_template.FIELDNAMES)
     
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
     idata=read_ascii(ifile, template=gyro_template)
     idata_names=strlowcase(tag_names(idata))
     srv_time_loc=where(idata_names EQ field_names[server_time_idx])
     idata_times_srv=idata.(srv_time_loc)
     ;; Number of lines in input
     lines=n_elements(idata_times_srv[0, *])
     ;; Remove quotes and separators
     IF size(idata_times_srv, /type) EQ 7 THEN BEGIN
        FOREACH fld, indgen((size(idata_times_srv, $
                                  /dimensions))[0]) DO BEGIN
           ok=strsplit(idata_times_srv[fld, *], '" -/:', /extract)
           ok=(temporary(ok)).toArray()
           ok=strjoin(transpose(temporary(ok)))
           idata_times_srv[fld, *]=ok
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

     ;; Obtain full server time details
     itimes_srv_std=parse_times(idata_times_srv, tnames_last_srv, $
                                time_locs_srv)
     gyro_srv_jd=reform(julday(long(itimes_srv_std[1, *]), $
                               long(itimes_srv_std[2, *]), $
                               long(itimes_srv_std[0, *]), $
                               long(itimes_srv_std[3, *]), $
                               long(itimes_srv_std[4, *]), $
                               float(itimes_srv_std[5, *])))

     ;; Read matching RMC file
     ifile_strl=strsplit(ifile, '_.', /extract) ; break string
     ifile_mstr=ifile_strl[n_elements(ifile_strl) - 2]
     rmc_pair=where(rmc_files_mstr EQ ifile_mstr, mcount)
     IF mcount LT 1 THEN BEGIN
        message, 'No matching file found. Skipping.', /CONTINUE
        CONTINUE
     ENDIF
     rmc=read_ascii(rmc_std_files[rmc_pair], template=rmc_template)
     rmc_names=strlowcase(tag_names(rmc))
     ;; Obtain server times and convert to Julian
     rmc_srv_time_loc=where(rmc_names EQ $
                            rmc_field_names[rmc_server_time_idx])
     rmc_times_srv=rmc.(rmc_srv_time_loc)
     rmc_times_srv_dims=size(rmc_times_srv, /dimensions)
     ;; Remove quotes
     IF size(rmc_times_srv, /type) EQ 7 THEN BEGIN
        FOREACH fld, indgen((size(rmc_times_srv, $
                                  /dimensions))[0]) DO BEGIN
           ok=strsplit(rmc_times_srv[fld, *], '" -/:', /extract)
           ok=(temporary(ok)).toArray()
           ok=strjoin(transpose(temporary(ok)))
           rmc_times_srv[fld, *]=ok
        ENDFOREACH
     ENDIF
     rmc_srv_jd=reform((rmc_times_srv_dims[0] EQ 6) ? $
                       julday(long(rmc_times_srv[1, *]), $
                              long(rmc_times_srv[2, *]), $
                              long(rmc_times_srv[0, *]), $
                              long(rmc_times_srv[3, *]), $
                              long(rmc_times_srv[4, *]), $
                              float(rmc_times_srv[5, *])) : $
                       julday(long(rmc_times_srv[1, *]), $
                              long(rmc_times_srv[2, *]), $
                              long(rmc_times_srv[0, *]), $
                              long(rmc_times_srv[3, *]), $
                              long(rmc_times_srv[4, *]), $
                              float(rmc_times_srv[5, *] + '.' + $
                                    rmc_times_srv[6, *])))
     ;; Obtain UTC times and convert to Julian
     rmc_utc_time_loc=where(rmc_names EQ $
                            rmc_field_names[rmc_utc_time_idx])
     rmc_times_utc=rmc.(rmc_utc_time_loc)
     rmc_times_utc_dims=size(rmc_times_utc, /dimensions)
     ;; Remove quotes
     IF size(rmc_times_utc, /type) EQ 7 THEN BEGIN
        FOREACH fld, indgen((size(rmc_times_utc, $
                                  /dimensions))[0]) DO BEGIN
           ok=strsplit(rmc_times_utc[fld, *], '" -/:', /extract)
           ok=(temporary(ok)).toArray()
           ok=strjoin(transpose(temporary(ok)))
           rmc_times_utc[fld, *]=ok
        ENDFOREACH
     ENDIF
     rmc_utc_jd=reform((rmc_times_utc_dims[0] EQ 6) ? $
                       julday(long(rmc_times_utc[1, *]), $
                              long(rmc_times_utc[2, *]), $
                              long(rmc_times_utc[0, *]), $
                              long(rmc_times_utc[3, *]), $
                              long(rmc_times_utc[4, *]), $
                              float(rmc_times_utc[5, *])) : $
                       julday(long(rmc_times_utc[1, *]), $
                              long(rmc_times_utc[2, *]), $
                              long(rmc_times_utc[0, *]), $
                              long(rmc_times_utc[3, *]), $
                              long(rmc_times_utc[4, *]), $
                              float(rmc_times_utc[5, *] + '.' + $
                                    rmc_times_utc[6, *])))

     ;; Find interpolate UTC time at the server time
     gyro_utc_jd=interpol(rmc_utc_jd, rmc_srv_jd, gyro_srv_jd)
     caldat, gyro_utc_jd, utc_month, utc_day, utc_year, $
             utc_hour, utc_minute, utc_second
     ;; Format these for output
     utc_year=string(utc_year, format='(i04)')
     utc_month=string(utc_month, format='(i02)')
     utc_day=string(utc_day, format='(i02)')
     utc_hour=string(utc_hour, format='(i02)')
     utc_minute=string(utc_minute, format='(i02)')
     ;; Becareful here with the format, since we may need more decimal
     ;; places, depending on the interpolation
     utc_second=string(utc_second, format='(f06.3)')

     odata=remove_structure_tags(idata, field_names[tags2remove])
     ;; Find indices to keep
     match2, strlowcase(tag_names(odata)), keep_fields, toss, keep
     tags2remove_odata=where(toss LT 0, nremove)
     IF nremove GT 0 THEN $
        odata=remove_structure_tags(odata, $
                                    (tag_names(odata))[tags2remove_odata])
     odata=create_struct('UTC_year', utc_year, $
                         'UTC_month', utc_month, $
                         'UTC_day', utc_day, $
                         'UTC_hour', utc_hour, $
                         'UTC_minute', utc_minute, $
                         'UTC_second', utc_second, $
                         odata)
     delvar, idata

     write_csv, ofile_name, odata, header=strlowcase(tag_names(odata))

  ENDFOR

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; std_gyro.pro ends here
