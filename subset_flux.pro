;; Author: Sebastian Luque
;; Created: 2013-11-03T18:49:19+0000
;; Last-Updated: 2015-06-30T19:53:47+0000
;;           By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     SUBSET_STD_FILE
;; 
;; PURPOSE:
;; 
;;     This function subsets any standardized file, extracting data within
;;     the given time bounds (bounds[0] <= data < bounds[1]).
;; 
;; CALLING SEQUENCE:
;; 
;;     subdata=subset_std_file(Ifile, Itemplate, Time_Idx, Offset, $
;;                             Time_Bounds, STATUS=STATUS)
;; 
;; INPUTS:
;; 
;;     Ifile:         Path of the CSV file to read.
;;     Itemplate:     ASCII template, as returned by ASCII_TEMPLATE.
;;     Time_Idx:      Index (in template) where time matrix is located.
;;     Offset:        A scalar value, usually corresponding to the sample
;;                    rate.  It is used to avoid numerical representation
;;                    issues when comparing against the time bounds.
;;     Time_Bounds:   A 2-element array with the time (julian) bounds.
;; 
;; KEYWORD PARAMETERS:
;; 
;;     STATUS:        A variable that, on output, will indicate whether the
;;                    subset operation was successful (0), or no records
;;                    where found within the time bounds.
;; 
;; OUTPUTS:
;;
;;     A structure like the ones returned by READ_ASCII, with the requested
;;     subset.
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

FUNCTION SUBSET_STD_FILE, IFILE, ITEMPLATE, TIME_IDX, OFFSET, $
                          TIME_BOUNDS, STATUS=STATUS

  ;; Parse input template
  field_names=itemplate.FIELDNAMES
  is_time_field=itemplate.FIELDGROUPS EQ time_idx
  ;; Ignore other groups when reading the data
  itemplate.FIELDGROUPS=indgen(itemplate.FIELDCOUNT)
  itemplate.FIELDGROUPS[where(is_time_field)]=time_idx
  tags2remove=where(field_names EQ field_names[time_idx])
  time_field_names=field_names[where(is_time_field)]
  ;; Offset fraction of the day
  offset_dfrac=double(offset) / double(86400)
  
  idata=read_std_file(ifile, itemplate, time_idx)
  idata_names=strlowcase(tag_names(idata))
  times=idata.(where(idata_names EQ field_names[time_idx]))
  idata=remove_structure_tags(idata, field_names[tags2remove])
  idata_names=strlowcase(tag_names(idata)) ; redoing after removing times
  times_s=times[5, *]
  times_s=fix(times_s) + $
          (round((double(times_s) - fix(times_s)) * 10) / 10.0)
  times_jd=reform(julday(long(times[1, *]), $
                         long(times[2, *]), $
                         long(times[0, *]), $
                         long(times[3, *]), $
                         long(times[4, *]), $
                         double(times_s)))
  ;; We have to subtract 0.5 the input sample rate to protect against
  ;; numerical representation issues in IDL...
  matches=where(times_jd GE time_bounds[0] AND $
                times_jd LT (time_bounds[1] - $
                             (offset_dfrac / 2)), $
                mcount)
  ;; Set the status error check
  status=1                      ; failed
  IF mcount LT 1 THEN RETURN, !VALUES.D_NAN
  odata=create_struct(time_field_names[0], reform(times[0, matches]))
  ;; Subset the rest of the data
  FOREACH tfld, (indgen(n_elements(time_field_names)))[1:*] DO BEGIN
     odata=create_struct(odata, time_field_names[tfld], $
                         reform(times[tfld, matches]))
  ENDFOREACH
  FOREACH fld, indgen(n_elements(idata_names)) DO BEGIN
     odata=create_struct(odata, idata_names[fld], $
                         idata.(fld)[matches])
  ENDFOREACH

  status=0                      ; success
  RETURN, odata

END


;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     SUBSET_FLUX
;; 
;; PURPOSE:
;; 
;;     This procedure subsets the (standardized) daily EC flux files,
;;     extracting the periods with suitable data, based on the diagnostic
;;     flag from MET data.  It also subsets the corresponding RMC and GYRO
;;     data for these periods.
;; 
;; CALLING SEQUENCE:
;; 
;;     SUBSET_FLUX, Idir, Itemplate_Sav, Time_Beg_Idx, Isample_Rate, $
;;                  Diag_Dir, Diag_Itemplate_Sav, Diag_Time_Idx, Diag_Idx, $
;;                  Ec_Period, RMC_Dir, RMC_Itemplate_Sav, RMC_Time_Idx,
;;                  Gyro_Dir, Gyro_Itemplate_Sav, Gyro_Time_Idx, RAD_Dir,
;;                  RAD_Itemplate_Sav, RAD_Time_Idx
;; 
;; INPUTS:
;; 
;;     Idir:                  Input directory (no trailing separator) for
;;                            (standardized) daily EC flux files.
;;     Itemplate_Sav:         Ascii template to read input files.
;;     Time_Beg_Idx:          Index (in template) where time is.
;;     Isample_Rate:          Scalar indicating the frequency (s) with
;;                            which input data were sampled.
;;     Diag_Dir:              Directory (no trailing separator) containing
;;                            processed MET files, averaged across flux
;;                            periods, and with a diagnostic field
;;                            indicating whether the period is suitable for
;;                            flux analyses (flag=0).
;;     Diag_Itemplate_Sav:    Ascii template to read Diag_Dir files.
;;     Diag_Time_Idx:         Index (in template) where time is in files in
;;                            Diag_Dir.
;;     Diag_Idx:              Index (in template) where the diagnostic
;;                            field is in files in Diag_Dir.
;;     Flux_Period:           Scalar indicating the duration (s) of flux
;;                            study periods.
;;     Rmc_Dir:               Directory (no trailing separator) containing
;;                            standardized daily RMC files.
;;     Rmc_Itemplate_Sav:     Ascii template to read Rmc_Dir files.
;;     Rmc_Time_Idx:          Index (in template) where time is in files in
;;                            Rmc_Dir.
;;     Gyro_Dir:              Directory (no trailing separator) containing
;;                            standardized daily GYRO files.
;;     Gyro_Itemplate_Sav:    Ascii template to read Gyro_Dir files.
;;     Gyro_Time_Idx:         Index (in template) where time is in files in
;;                            Gyro_Dir.
;;     RAD_Dir:               Directory (no trailing separator) containing
;;                            standardized daily processed RAD files.
;;     RAD_Itemplate_Sav:     Ascii template to read RAD_Dir files.
;;     RAD_Time_Idx:          Index (in template) where time is in files in
;;                            RAD_Dir.
;; 
;; KEYWORD PARAMETERS:
;; 
;;     OVERWRITE:             Whether to overwrite files in output
;;                            subdirectories.
;; 
;; SIDE EFFECTS:
;; 
;;     Files are written in Odir.
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

PRO SUBSET_FLUX, IDIR, ITEMPLATE_SAV, TIME_IDX, ISAMPLE_RATE, $
                 DIAG_DIR, DIAG_ITEMPLATE_SAV, DIAG_TIME_IDX, DIAG_IDX, $
                 EC_PERIOD, RMC_DIR, RMC_ITEMPLATE_SAV, RMC_TIME_IDX, $
                 GYRO_DIR, GYRO_ITEMPLATE_SAV, GYRO_TIME_IDX, $
                 RAD_DIR, RAD_ITEMPLATE_SAV, RAD_TIME_IDX, $
                 OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 18) THEN $
     message, 'Usage: SUBSET_FLUX, IDIR, ITEMPLATE_SAV, TIME_IDX, ' + $
              'ISAMPLE_RATE, DIAG_DIR, DIAG_ITEMPLATE_SAV, ' + $
              'DIAG_TIME_IDX, DIAG_IDX, EC_PERIOD, RMC_DIR, ' + $
              'RMC_ITEMPLATE_SAV, RMC_TIME_IDX, GYRO_DIR, ' + $
              'GYRO_ITEMPLATE_SAV, GYRO_TIME_IDX, ' + $
              'RAD_DIR, RAD_ITEMPLATE_SAV, RAD_TIME_IDX'
  idir_info=file_info(idir)
  itpl_info=file_info(itemplate_sav)
  diag_dir_info=file_info(diag_dir)
  diag_tpl_info=file_info(diag_itemplate_sav)
  rmc_dir_info=file_info(rmc_dir)
  rmc_tpl_info=file_info(rmc_itemplate_sav)
  gyro_dir_info=file_info(gyro_dir)
  gyro_tpl_info=file_info(gyro_itemplate_sav)
  rad_dir_info=file_info(rad_dir)
  rad_tpl_info=file_info(rad_itemplate_sav)
  IF (~idir_info.directory) THEN $
     message, 'IDIR must be a string pointing to an existing directory'
  IF (~itpl_info.read) THEN $
     message, 'ITEMPLATE_SAV must be a string pointing to a readable file'
  IF ((n_elements(time_idx) NE 1) OR $
      ((size(time_idx, /type) NE 2) || time_idx LT 0)) THEN $
         message, 'TIME_IDX must be an integer scalar >= zero'
  IF ((n_elements(isample_rate) NE 1) OR (isample_rate LT 0)) THEN $
     message, 'ISAMPLE_RATE must be a scalar >= zero'
  IF (~diag_dir_info.directory) THEN $
     message, 'DIAG_DIR must be a string pointing to an existing directory'
  IF (~diag_tpl_info.read) THEN $
     message, 'DIAG_ITEMPLATE_SAV must be a string pointing to a ' + $
              'readable file'
  IF ((n_elements(diag_time_idx) NE 1) OR $
      ((size(diag_time_idx, /type) NE 2) || diag_time_idx LT 0)) THEN $
         message, 'DIAG_TIME_IDX must be an integer scalar >= zero'
  IF ((n_elements(diag_idx) NE 1) OR $
      ((size(diag_idx, /type) NE 2) || diag_idx LT 0)) THEN $
         message, 'DIAG_IDX must be an integer scalar >= zero'
  IF ((n_elements(ec_period) NE 1) OR (ec_period LT 0)) THEN $
     message, 'EC_PERIOD must be a scalar >= zero'
  IF ((86400 MOD float(ec_period)) NE 0) THEN $
     message, 'EC_PERIOD must be an integer divisor of 86400 s'
  IF (~rmc_dir_info.directory) THEN $
     message, 'RMC_DIR must be a string pointing to an existing directory'
  IF (~rmc_tpl_info.read) THEN $
     message, 'RMC_ITEMPLATE_SAV must be a string pointing to a ' + $
              'readable file'
  IF ((n_elements(rmc_time_idx) NE 1) OR $
      ((size(rmc_time_idx, /type) NE 2) || rmc_time_idx LT 0)) THEN $
         message, 'RMC_TIME_IDX must be an integer scalar >= zero'
  IF (~gyro_dir_info.directory) THEN $
     message, 'GYRO_DIR must be a string pointing to an existing directory'
  IF (~gyro_tpl_info.read) THEN $
     message, 'GYRO_ITEMPLATE_SAV must be a string pointing to a ' + $
              'readable file'
  IF ((n_elements(gyro_time_idx) NE 1) OR $
      ((size(gyro_time_idx, /type) NE 2) || gyro_time_idx LT 0)) THEN $
         message, 'GYRO_TIME_IDX must be an integer scalar >= zero'
  IF (~rad_dir_info.directory) THEN $
     message, 'RAD_DIR must be a string pointing to an existing directory'
  IF (~rad_tpl_info.read) THEN $
     message, 'RAD_ITEMPLATE_SAV must be a string pointing to a ' + $
              'readable file'
  IF ((n_elements(rad_time_idx) NE 1) OR $
      ((size(rad_time_idx, /type) NE 2) || rad_time_idx LT 0)) THEN $
         message, 'RAD_TIME_IDX must be an integer scalar >= zero'
  idir_files=file_search(idir + path_sep() + '*', count=nidir_files, $
                         /nosort, /fold_case, /test_regular)
  diag_files=file_search(diag_dir + path_sep() + '*', count=ndiag_files, $
                         /nosort, /fold_case, /test_regular)
  rmc_files=file_search(rmc_dir + path_sep() + '*', count=nrmc_files, $
                        /nosort, /fold_case, /test_regular)
  gyro_files=file_search(gyro_dir + path_sep() + '*', count=ngyro_files, $
                         /nosort, /fold_case, /test_regular)
  rad_files=file_search(rad_dir + path_sep() + '*', count=nrad_files, $
                        /nosort, /fold_case, /test_regular)
  IF (nidir_files EQ 0) OR (ndiag_files EQ 0) OR $
     (nrmc_files EQ 0) OR (ngyro_files EQ 0) OR $
     (nrad_files EQ 0) THEN $
        message, 'No files in at least one of IDIR, DIAG_DIR, ' + $
                 'RMC_DIR, GYRO_DIR, or RAD_DIR'

  ;; Parse DIAG MET template
  restore, diag_itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  diag_template=itemplate
  diag_field_names=strlowcase(diag_template.FIELDNAMES)
  diag_field_types=diag_template.FIELDTYPES
  diag_is_time_field=diag_template.FIELDGROUPS EQ diag_time_idx
  ;; Ignore other groups when reading the data
  diag_template.FIELDGROUPS=indgen(diag_template.FIELDCOUNT)
  diag_template.FIELDGROUPS[where(diag_is_time_field)]=diag_time_idx
  diag_non_time_fields=where(~diag_is_time_field)
  diag_non_time_field_names=diag_field_names[diag_non_time_fields]
  diag_tags2remove=where(diag_field_names EQ diag_field_names[diag_time_idx])
  ;; Times
  diag_tfields=where(diag_is_time_field, /NULL)
  diag_tnames=diag_field_names[diag_tfields]
  diag_tnamesl=strsplit(diag_tnames, '_', /extract)
  diag_tnames_last=strarr(n_elements(diag_tnamesl))
  ntbits=n_elements(diag_tnamesl[0])
  diag_tnames_id=strjoin((diag_tnamesl[0])[0:ntbits - 2], '_')
  FOR i=0L, n_elements(diag_tnames) - 1 DO $
     diag_tnames_last[i]=diag_tnamesl[i, n_elements(diag_tnamesl[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  diag_time_locs=locate_time_strings(diag_tnames_last)

  ;; Parse RMC files
  restore, rmc_itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  rmc_template=itemplate
  ;; Break file names and extract the piece to match
  rmc_filesl=strsplit(rmc_files, '_.', /extract)
  rmc_files_a=rmc_filesl.toArray(/transpose) ; array with pieces in cols
  rmc_files_mstr_dims=size(rmc_files_a, /dimensions)
  ;; We want the 2nd to last piece
  rmc_files_mstr=rmc_files_a[rmc_files_mstr_dims[0] - 2, *]

  ;; Parse Gyro files
  restore, gyro_itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  gyro_template=itemplate
  ;; Break file names and extract the piece to match
  gyro_filesl=strsplit(gyro_files, '_.', /extract)
  gyro_files_a=gyro_filesl.toArray(/transpose) ; array with pieces in cols
  gyro_files_mstr_dims=size(gyro_files_a, /dimensions)
  ;; We want the 2nd to last piece
  gyro_files_mstr=gyro_files_a[gyro_files_mstr_dims[0] - 2, *]

  ;; Parse RAD files
  restore, rad_itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  rad_template=itemplate
  ;; Break file names and extract the piece to match
  rad_filesl=strsplit(rad_files, '_.', /extract)
  rad_files_a=rad_filesl.toArray(/transpose) ; array with pieces in cols
  rad_files_mstr_dims=size(rad_files_a, /dimensions)
  ;; We want the 3rd to last piece
  rad_files_mstr=rad_files_a[rad_files_mstr_dims[0] - 3, *]

  ;; Parse input flux files (last because we don't make a copy of template)
  restore, itemplate_sav
  ;; Break file names and extract the piece to match
  flux_filesl=strsplit(idir_files, '_.', /extract)
  flux_files_a=flux_filesl.toArray(/transpose) ; array with pieces in cols
  flux_files_mstr_dims=size(flux_files_a, /dimensions)
  ;; We want the 2nd to last piece
  flux_files_mstr=flux_files_a[flux_files_mstr_dims[0] - 2, *]

  ;; Loop through files in DIAG MET directory
  FOR k=0, ndiag_files - 1 DO BEGIN
     dfile=diag_files[k]
     message, 'Processing ' + dfile, /informational
     ;; Read input file
     diag=read_ascii(dfile, template=diag_template)
     diag_names=strlowcase(tag_names(diag))
     time_loc=where(diag_names EQ diag_field_names[diag_time_idx])
     diag_times=diag.(time_loc)
     ;; Number of lines in input
     lines=n_elements(diag_times[0, *])
     ;; Remove quotes and separators
     IF size(diag_times, /type) EQ 7 THEN BEGIN
        FOREACH fld, indgen((size(diag_times, $
                                  /dimensions))[0]) DO BEGIN
           ok=strsplit(diag_times[fld, *], '" -/:', /extract)
           ok=(temporary(ok)).toArray()
           ok=strjoin(transpose(temporary(ok)))
           diag_times[fld, *]=ok
        ENDFOREACH
     ENDIF
     match2, diag_names, diag_field_names[diag_tags2remove], is_time
     FOREACH fld, (indgen(n_tags(diag)))[where(is_time LT 0)] DO BEGIN
        IF size(diag.(fld), /type) EQ 7 THEN BEGIN
           ok=strsplit(diag.(fld), '" ', /extract)
           diag.(fld)=ok.toArray()
        ENDIF
     ENDFOREACH
     ;; At this stage we can skip records that indicate non-fluxable
     ;; periods (diag flag is 0)
     diag_flag=diag.(where(diag_names EQ diag_field_names[diag_idx]))
     fluxable=where(diag_flag EQ 0, nfluxable)
     IF nfluxable EQ 0 THEN CONTINUE
     ;; Obtain full input time details, rounding off seconds to nearest 0.1
     ;; s, to avoid numerical representation problems
     dtimes_std=parse_times(diag_times, diag_tnames_last, diag_time_locs)
     dtimes_s=dtimes_std[5, *]
     dtimes_s=fix(dtimes_s) + $
              (round((double(dtimes_s) - fix(dtimes_s)) * 10) / 10.0)
     flux_jd=reform(julday(long(dtimes_std[1, *]), $
                           long(dtimes_std[2, *]), $
                           long(dtimes_std[0, *]), $
                           long(dtimes_std[3, *]), $
                           long(dtimes_std[4, *]), $
                           double(dtimes_s)))
     dfile_strl=strsplit(dfile, '_.', /extract) ; break string
     dfile_mstr=dfile_strl[n_elements(dfile_strl) - 4]

     ;; Loop through fluxable periods
     FOREACH fperiod, fluxable DO BEGIN
        ;; Determine flux period time bounds
        ec_period_dfrac=double(ec_period) / double(86400)
        isample_rate_dfrac=double(isample_rate) / double(86400)
        bounds=[flux_jd[fperiod], flux_jd[fperiod] + ec_period_dfrac]
        ;; Time stamp.  Note we are flooring the seconds
        tstamp=string(dtimes_std[3, fperiod] + dtimes_std[4, fperiod] + $
                      dtimes_std[5, fperiod], format='(i06)')

        ;; Read matching RMC file
        rmc_pair=where(rmc_files_mstr EQ dfile_mstr, mcount)
        IF mcount LT 1 THEN BEGIN
           message, 'No matching RMC file found. Skipping.', /CONTINUE
           CONTINUE
        ENDIF
        ;; Build a name for output file and check existence
        iname=strsplit(file_basename(rmc_files[rmc_pair]), '.', /extract)
        ;; Get a path for the file, check if it already exists.  We are
        ;; placing this in a sub-directory in RMC_DIR.
        rmc_subdir=rmc_dir + path_sep() + 'Periods'
        file_mkdir, rmc_subdir
        oname_prefix=strcompress(rmc_subdir + path_sep() + iname[0])
        ofile_name=strcompress(oname_prefix + '_' + tstamp + '.' + $
                               iname[1], /remove_all)
        ofile_stamp=file_basename(ofile_name)
        out_list=file_search(rmc_subdir + path_sep() + '*.' + iname[1], $
                             /nosort, /fold_case, /test_regular)
        matchfiles=where(ofile_stamp EQ file_basename(out_list), $
                         matchfilecount)
        IF matchfilecount GT 0 THEN BEGIN
           IF keyword_set(overwrite) THEN BEGIN
              message, 'Flux period file ' + ofile_stamp + $
                       ' already exists.  Overwriting', /informational
           ENDIF ELSE BEGIN
              message, 'Flux period file ' + ofile_stamp + $
                       ' already exists.  Not overwriting', /informational
              CONTINUE
           ENDELSE
        ENDIF
        rmc_period=subset_std_file(rmc_files[rmc_pair], rmc_template, $
                                   rmc_time_idx, isample_rate, bounds, $
                                   status=status)
        IF status NE 0 THEN BEGIN
           message, 'No records found within bounds. Skipping file.', $
                    /informational
           CONTINUE
        ENDIF
        write_csv, ofile_name, rmc_period, $
                   header=strlowcase(tag_names(rmc_period))
        delvar, rmc_period

        ;; Read matching GYRO file
        gyro_pair=where(gyro_files_mstr EQ dfile_mstr, mcount)
        IF mcount LT 1 THEN BEGIN
           message, 'No matching GYRO file found. Skipping.', /CONTINUE
           CONTINUE
        ENDIF
        ;; Build a name for output file and check existence
        iname=strsplit(file_basename(gyro_files[gyro_pair]), '.', /extract)
        ;; Get a path for the file, check if it already exists.  We are
        ;; placing this in a sub-directory in GYRO_DIR.
        gyro_subdir=gyro_dir + path_sep() + 'Periods'
        file_mkdir, gyro_subdir
        oname_prefix=strcompress(gyro_subdir + path_sep() + iname[0])
        ofile_name=strcompress(oname_prefix + '_' + tstamp + '.' + $
                               iname[1], /remove_all)
        ofile_stamp=file_basename(ofile_name)
        out_list=file_search(gyro_subdir + path_sep() + '*.' + iname[1], $
                             /nosort, /fold_case, /test_regular)
        matchfiles=where(ofile_stamp EQ file_basename(out_list), $
                         matchfilecount)
        IF matchfilecount GT 0 THEN BEGIN
           IF keyword_set(overwrite) THEN BEGIN
              message, 'Flux period file ' + ofile_stamp + $
                       ' already exists.  Overwriting', /informational
           ENDIF ELSE BEGIN
              message, 'Flux period file ' + ofile_stamp + $
                       ' already exists.  Not overwriting', /informational
              CONTINUE
           ENDELSE
        ENDIF
        gyro_period=subset_std_file(gyro_files[gyro_pair], gyro_template, $
                                    gyro_time_idx, isample_rate, bounds, $
                                    status=status)
        IF status NE 0 THEN BEGIN
           message, 'No records found within bounds. Skipping file.', $
                    /informational
           CONTINUE
        ENDIF
        write_csv, ofile_name, gyro_period, $
                   header=strlowcase(tag_names(gyro_period))
        delvar, gyro_period

        ;; Read matching RAD file
        rad_pair=where(rad_files_mstr EQ dfile_mstr, mcount)
        IF mcount LT 1 THEN BEGIN
           message, 'No matching RAD file found. Skipping.', /CONTINUE
           CONTINUE
        ENDIF
        ;; Build a name for output file and check existence
        iname=strsplit(file_basename(rad_files[rad_pair]), '.', /extract)
        ;; Get a path for the file, check if it already exists.  We are
        ;; placing this in a sub-directory in RAD_DIR.
        rad_subdir=rad_dir + path_sep() + 'Periods'
        file_mkdir, rad_subdir
        oname_prefix=strcompress(rad_subdir + path_sep() + iname[0])
        ofile_name=strcompress(oname_prefix + '_' + tstamp + '.' + $
                               iname[1], /remove_all)
        ofile_stamp=file_basename(ofile_name)
        out_list=file_search(rad_subdir + path_sep() + '*.' + iname[1], $
                             /nosort, /fold_case, /test_regular)
        matchfiles=where(ofile_stamp EQ file_basename(out_list), $
                         matchfilecount)
        IF matchfilecount GT 0 THEN BEGIN
           IF keyword_set(overwrite) THEN BEGIN
              message, 'Flux period file ' + ofile_stamp + $
                       ' already exists.  Overwriting', /informational
           ENDIF ELSE BEGIN
              message, 'Flux period file ' + ofile_stamp + $
                       ' already exists.  Not overwriting', /informational
              CONTINUE
           ENDELSE
        ENDIF
        rad_period=subset_std_file(rad_files[rad_pair], rad_template, $
                                   rad_time_idx, isample_rate, bounds, $
                                   status=status)
        IF status NE 0 THEN BEGIN
           message, 'No records found within bounds. Skipping file.', $
                    /informational
           CONTINUE
        ENDIF
        write_csv, ofile_name, rad_period, $
                   header=strlowcase(tag_names(rad_period))
        delvar, rad_period

        ;; Read matching flux EC file
        flux_pair=where(flux_files_mstr EQ dfile_mstr, mcount)
        IF mcount LT 1 THEN BEGIN
           message, 'No matching EC flux file found. Skipping.', /CONTINUE
           CONTINUE
        ENDIF
        ;; Build a name for output file and check existence
        iname=strsplit(file_basename(idir_files[flux_pair]), '.', /extract)
        ;; Get a path for the file, check if it already exists.  We are
        ;; placing this in a sub-directory in IDIR
        odir=idir + path_sep() + 'Periods'
        file_mkdir, odir
        oname_prefix=strcompress(odir + path_sep() + iname[0])
        ofile_name=strcompress(oname_prefix + '_' + tstamp + '.' + $
                               iname[1], /remove_all)
        ofile_stamp=file_basename(ofile_name)
        out_list=file_search(odir + path_sep() + '*.' + iname[1], $
                             /nosort, /fold_case, /test_regular)
        matchfiles=where(ofile_stamp EQ file_basename(out_list), $
                         matchfilecount)
        IF matchfilecount GT 0 THEN BEGIN
           IF keyword_set(overwrite) THEN BEGIN
              message, 'Flux period file ' + ofile_stamp + $
                       ' already exists.  Overwriting', /informational
           ENDIF ELSE BEGIN
              message, 'Flux period file ' + ofile_stamp + $
                       ' already exists.  Not overwriting', /informational
              CONTINUE
           ENDELSE
        ENDIF
        flux_period=subset_std_file(idir_files[flux_pair], itemplate, $
                                    time_idx, isample_rate, bounds, $
                                    status=status)
        IF status NE 0 THEN BEGIN
           message, 'No records found within bounds. Skipping file.', $
                    /informational
           CONTINUE
        ENDIF
        write_csv, ofile_name, flux_period, $
                   header=strlowcase(tag_names(flux_period))
        delvar, flux_period

     ENDFOREACH

  ENDFOR

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
