;; $Id: $
;; Author: Sebastian Luque
;; Created: 2013-10-15T15:43:56+0000
;; Last-Updated: 2013-10-24T21:10:57+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     PROCESS_MET
;; 
;; PURPOSE:
;; 
;;     Process MET data.
;; 
;; CALLING SEQUENCE:
;; 
;;     PROCESS_MET, Idir, Odir, Itemplate_Sav, Time_Beg_Idx, RMC_Dir, $
;;                  RMC_Itemplate_Sav, RMC_Time_Idx, RMC_Pull_Idx, $
;;                  GYRO_Dir, GYRO_Itemplate_Sav, GYRO_Time_Idx, $
;;                  GYRO_Pull_Idx, Log_File, Log_Itemplate_Sav, $
;;                  Log_TimeBeg_Idx, Log_TimeEnd_Idx
;; 
;; INPUTS:
;; 
;;     Idir:                Input directory (no trailing separator).
;;     Odir:                Output directory (no trailing separator).
;;     Itemplate_Sav:       Ascii template to read input files.
;;     Time_Beg_Idx:        Index (in template) where time is.
;;     RMC_Dir:             Directory where RMC files are found (no
;;                          trailing separator).
;;     RMC_Itemplate_Sav:   Ascii template to read RMC files.
;;     RMC_Time_Idx:        Index (in template) where time is in RMC files.
;;     RMC_Pull_Idx:        Integer array with indices (in template) of
;;                          fields to pull from RMC files.
;;     GYRO_Dir:            Directory where GYRO files are found (no
;;                          trailing separator).
;;     GYRO_Itemplate_Sav:  Ascii template to read GYRO files.
;;     GYRO_Time_Idx:       Index (in template) where time is in GYRO files.
;;     GYRO_Pull_Idx:       Integer array with indices (in template) of
;;                          fields to pull from GYRO files.
;;     Log_File:            Path to file with MET log.
;;     Log_Itemplate_Sav:   Ascii template to read MET log file.
;;     Log_TimeBeg_Idx:     Index (in template) where starting time is in
;;                          log file.
;;     Log_TimeEnd_Idx:     Index (in template) where ending time is in
;;                          log file.
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
;; 
;; 
;;- -----------------------------------------------------------------------
;;; Code:

PRO PROCESS_MET, IDIR, ODIR, ITEMPLATE_SAV, TIME_BEG_IDX, $
                 RMC_DIR, RMC_ITEMPLATE_SAV, RMC_TIME_IDX, $
                 GYRO_DIR, GYRO_ITEMPLATE_SAV, GYRO_TIME_IDX, $
                 LOG_FILE, LOG_ITEMPLATE_SAV, LOG_TIME_BEG_IDX, $
                 LOG_TIME_END_IDX, OVERWRITE=OVERWRITE

  idir_files=file_search(idir + path_sep() + '*', count=nidir_files, $
                         /nosort, /fold_case, /test_regular)
  IF nidir_files LT 1 THEN message, 'No input files found'
  rmc_files=file_search(rmc_dir + path_sep() + '*', count=nrmc_files, $
                        /nosort, /fold_case, /test_regular)
  IF nrmc_files LT 1 THEN message, 'No RMC files found.  Exiting'
  gyro_files=file_search(gyro_dir + path_sep() + '*', count=ngyro_files, $
                        /nosort, /fold_case, /test_regular)
  IF ngyro_files LT 1 THEN message, 'No GYRO files found.  Exiting'
  log_file_info=file_info(log_file)
  IF log_file_info.regular NE 1 THEN $
     message, 'Log file is not a regular file.  Exiting'

  ;; Parse input template
  restore, itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  met_template=itemplate
  field_names=strlowcase(met_template.FIELDNAMES)
  field_types=met_template.FIELDTYPES
  is_time_field=met_template.FIELDGROUPS EQ time_beg_idx
  ;; Ignore other groups when reading the data
  met_template.FIELDGROUPS=indgen(met_template.FIELDCOUNT)
  met_template.FIELDGROUPS[where(is_time_field)]=time_beg_idx
  non_time_fields=where(~is_time_field)
  non_time_field_names=field_names[non_time_fields]
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

  ;; Parse RMC template
  restore, rmc_itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  rmc_template=itemplate
  rmc_field_names=strlowcase(rmc_template.FIELDNAMES)
  rmc_field_types=rmc_template.FIELDTYPES
  rmc_is_time_field=rmc_template.FIELDGROUPS EQ rmc_time_idx
  ;; Ignore other groups when reading the data
  rmc_template.FIELDGROUPS=indgen(rmc_template.FIELDCOUNT)
  rmc_template.FIELDGROUPS[where(rmc_is_time_field)]=rmc_time_idx
  rmc_non_time_fields=where(~rmc_is_time_field)
  rmc_non_time_field_names=rmc_field_names[rmc_non_time_fields]
  ;; Times
  rmc_tfields=where(rmc_is_time_field, /NULL)
  rmc_tnames=rmc_field_names[rmc_tfields]
  rmc_tnamesl=strsplit(rmc_tnames, '_', /extract)
  rmc_tnames_last=strarr(n_elements(rmc_tnamesl))
  rmc_tnames_id=strjoin((rmc_tnamesl[0])[0:n_elements(rmc_tnamesl[0]) - 2], '_')
  FOR i=0L, n_elements(tnames) - 1 DO $
     rmc_tnames_last[i]=rmc_tnamesl[i, n_elements(rmc_tnamesl[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  rmc_time_locs=locate_time_strings(rmc_tnames_last)
  ;; Break file names and extract the piece to match
  rmc_filesl=strsplit(rmc_files, '_.', /extract)
  rmc_files_a=rmc_filesl.toArray(/transpose) ; array with pieces in cols
  rmc_files_mstr_dims=size(rmc_files_a, /dimensions)
  ;; We want the 3rd to last piece
  rmc_files_mstr=rmc_files_a[rmc_files_mstr_dims[0] - 3, *]

  ;; Parse Gyro template
  restore, gyro_itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  gyro_template=itemplate
  gyro_field_names=strlowcase(gyro_template.FIELDNAMES)
  gyro_field_types=gyro_template.FIELDTYPES
  gyro_is_time_field=gyro_template.FIELDGROUPS EQ gyro_time_idx
  ;; Ignore other groups when reading the data
  gyro_template.FIELDGROUPS=indgen(gyro_template.FIELDCOUNT)
  gyro_template.FIELDGROUPS[where(gyro_is_time_field)]=gyro_time_idx
  gyro_non_time_fields=where(~gyro_is_time_field)
  gyro_non_time_field_names=gyro_field_names[gyro_non_time_fields]
  ;; Times
  gyro_tfields=where(gyro_is_time_field, /NULL)
  gyro_tnames=rmc_field_names[gyro_tfields]
  gyro_tnamesl=strsplit(gyro_tnames, '_', /extract)
  gyro_tnames_last=strarr(n_elements(gyro_tnamesl))
  gyro_tnames_id=strjoin((gyro_tnamesl[0])[0:n_elements(gyro_tnamesl[0]) - 2], '_')
  FOR i=0L, n_elements(tnames) - 1 DO $
     gyro_tnames_last[i]=gyro_tnamesl[i, n_elements(gyro_tnamesl[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  gyro_time_locs=locate_time_strings(gyro_tnames_last)
  ;; Break file names and extract the piece to match
  gyro_filesl=strsplit(gyro_files, '_.', /extract)
  gyro_files_a=gyro_filesl.toArray(/transpose) ; array with pieces in cols
  gyro_files_mstr_dims=size(gyro_files_a, /dimensions)
  ;; We want the 3rd to last piece
  gyro_files_mstr=gyro_files_a[gyro_files_mstr_dims[0] - 3, *]

  ;; Loop through files in input directory
  FOR k=0, nidir_files - 1 DO BEGIN
     iname=strsplit(file_basename(idir_files[k]), '.', /extract)
     ;; Get a path for the file, check if it already exists
     ofile_name=strcompress(odir + path_sep() + iname[0] + '_proc.' + $
                            iname[1], /remove_all)
     ofile_stamp=file_basename(ofile_name)
     out_list=file_search(odir + path_sep() + '*.' + iname[1], $
                          /nosort, /fold_case, /test_regular)
     matchfiles=where(ofile_stamp EQ file_basename(out_list), $
                      matchfilecount)
     IF matchfilecount GT 0 THEN BEGIN
        IF keyword_set(overwrite) THEN BEGIN
           message, 'Processed file ' + ofile_stamp + $
                    ' already exists.  Overwriting', /informational
        ENDIF ELSE BEGIN
        message, 'Processed file ' + ofile_stamp + $
                 ' already exists.  Not overwriting', /informational
        CONTINUE
     ENDELSE
     ENDIF

     ifile=idir_files[k]
     message, 'Processing ' + ifile, /informational
     ;; Read input file
     idata=read_ascii(ifile, template=met_template)
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

     ;; Obtain full input time details
     itimes_std=parse_times(idata_times, tnames_last, time_locs)
     met_jd=reform(julday(long(itimes_std[1, *]), $
                          long(itimes_std[2, *]), $
                          long(itimes_std[0, *]), $
                          long(itimes_std[3, *]), $
                          long(itimes_std[4, *]), $
                          float(itimes_std[5, *])))

     ifile_strl=strsplit(ifile, '_.', /extract) ; break string
     ifile_mstr=ifile_strl[n_elements(ifile_strl) - 2]
     ;; Read matching RMC file
     rmc_pair=where(rmc_files_mstr EQ ifile_mstr, mcount)
     IF mcount LT 1 THEN BEGIN
        message, 'No matching RMC file found. Skipping.', /CONTINUE
        CONTINUE
     ENDIF
     rmc=read_ascii(rmc_files[rmc_pair], template=rmc_template)
     rmc_names=strlowcase(tag_names(rmc))
     ;; Obtain times and convert to Julian
     rmc_time_loc=where(rmc_names EQ rmc_field_names[rmc_time_idx])
     rmc_times=rmc.(rmc_time_loc)
     rmc_times_dims=size(rmc_times, /dimensions)
     ;; Remove quotes
     IF size(rmc_times, /type) EQ 7 THEN BEGIN
        FOREACH fld, indgen((size(rmc_times, $
                                  /dimensions))[0]) DO BEGIN
           ok=strsplit(rmc_times[fld, *], '" -/:', /extract)
           ok=(temporary(ok)).toArray()
           ok=strjoin(transpose(temporary(ok)))
           rmc_times[fld, *]=ok
        ENDFOREACH
     ENDIF
     rmc_jd=reform((rmc_times_dims[0] EQ 6) ? $
                   julday(long(rmc_times[1, *]), $
                          long(rmc_times[2, *]), $
                          long(rmc_times[0, *]), $
                          long(rmc_times[3, *]), $
                          long(rmc_times[4, *]), $
                          float(rmc_times[5, *])) : $
                   julday(long(rmc_times[1, *]), $
                          long(rmc_times[2, *]), $
                          long(rmc_times[0, *]), $
                          long(rmc_times[3, *]), $
                          long(rmc_times[4, *]), $
                          float(rmc_times[5, *] + '.' + $
                                rmc_times[6, *])))
     ;; Find matching times
     match2, met_jd, rmc_jd, met_in_rmc, rmc_in_met
     met_matches=where(met_in_rmc GE 0, mcount, /null)
     IF mcount LT 1 THEN BEGIN
        message, 'No matching RMC records found.  Skipping file.', $
                 /informational
        CONTINUE
     ENDIF
     rmc_matches=where(rmc_in_met GE 0)
     FOREACH fld, rmc_lonlat_idx DO BEGIN
        match_fld=where(rmc_names EQ rmc_field_names[fld])
        idata=create_struct(idata, rmc_field_names[fld], $
                            rmc.(match_fld)[rad_matches])
     ENDFOREACH

     ;; Read matching MET file
     met_pair=where(met_files_mstr EQ ifile_mstr, mcount)
     IF mcount LT 1 THEN BEGIN
        message, 'No matching MET file found. Skipping.', /CONTINUE
        CONTINUE
     ENDIF
     met=read_ascii(met_files[met_pair], template=met_template)
     met_names=strlowcase(tag_names(met))
     ;; Obtain times and convert to Julian
     met_time_loc=where(met_names EQ met_field_names[met_time_idx])
     met_times=met.(met_time_loc)
     met_times_dims=size(met_times, /dimensions)
     ;; Remove quotes
     IF size(met_times, /type) EQ 7 THEN BEGIN
        FOREACH fld, indgen((size(met_times, $
                                  /dimensions))[0]) DO BEGIN
           ok=strsplit(met_times[fld, *], '" -/:', /extract)
           ok=(temporary(ok)).toArray()
           ok=strjoin(transpose(temporary(ok)))
           met_times[fld, *]=ok
        ENDFOREACH
     ENDIF
     met_jd=reform((met_times_dims[0] EQ 6) ? $
                   julday(long(met_times[1, *]), $
                          long(met_times[2, *]), $
                          long(met_times[0, *]), $
                          long(met_times[3, *]), $
                          long(met_times[4, *]), $
                          float(met_times[5, *])) : $
                   julday(long(met_times[1, *]), $
                          long(met_times[2, *]), $
                          long(met_times[0, *]), $
                          long(met_times[3, *]), $
                          long(met_times[4, *]), $
                          float(met_times[5, *] + '.' + $
                                met_times[6, *])))
     ;; Find matching times
     match2, rad_jd, met_jd, rad_in_met, met_in_rad
     rad_matches=where(rad_in_met GE 0, mcount, /null)
     IF mcount LT 1 THEN BEGIN
        message, 'No matching MET records found.  Skipping file.', $
                 /informational
        CONTINUE
     ENDIF
     met_matches=where(met_in_rad GE 0)
     ;; Below may not be necessary now, but in the future we may want to
     ;; pull other data from MET.  This extra PAR field is from another
     ;; sensor in the deck.  Need explanation about why this is being
     ;; pulled in.  Also the daily files have a number of blank fields
     ;; (look like just holders) for standard deviations, which are ignored
     ;; in the current daily MET files.
     FOREACH fld, met_par_idx DO BEGIN
        match_fld=where(met_names EQ met_field_names[fld])
        idata=create_struct(idata, met_field_names[fld] + '_met', $
                            met.(match_fld)[rad_matches])
     ENDFOREACH

     odata=remove_structure_tags(idata, field_names[tags2remove])
     delvar, idata
     ;; OK, how else to just extract the time info into a structure
     revtidx=reverse(indgen((size(idata_times, /dimensions))[0]))
     FOREACH fld, revtidx DO BEGIN
        fld_name=tnames[fld]
        odata=create_struct(fld_name, reform(idata_times[fld, *]), odata)
     ENDFOREACH
    
     write_csv, ofile_name, odata, header=strlowcase(tag_names(odata))

  ENDFOR


END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; process_MET.pro ends here
