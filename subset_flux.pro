;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-11-03T18:49:19+0000
;; Last-Updated: 2013-11-11T22:19:15+0000
;;           By: Sebastian Luque
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
;;     SUBSET_FLUX, Idir, Odir, Itemplate_Sav, Time_Beg_Idx, Isample_Rate, $
;;                  DIAG_DIR, Diag_Itemplate_Save, Diag_Time_Idx, Diag_Idx, $
;;                  Flux_Period
;; 
;; INPUTS:
;; 
;;     Idir:                  Input directory (no trailing separator) for
;;                            (standardized) daily EC flux files.
;;     Odir:                  Output directory (no trailing separator).
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
;; 
;; KEYWORD PARAMETERS:
;; 
;; 
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
;;; Code:

FUNCTION SUBSET_FILE, IFILE, ITEMPLATE_SAV, TIME_IDX, ISAMPLE_RATE, $
                      TIME_BOUNDS

  ;; Parse input flux template
  restore, itemplate_sav
  field_names=strlowcase(itemplate.FIELDNAMES)
  field_types=itemplate.FIELDTYPES
  is_time_field=itemplate.FIELDGROUPS EQ time_idx
  ;; Ignore other groups when reading the data
  itemplate.FIELDGROUPS=indgen(itemplate.FIELDCOUNT)
  itemplate.FIELDGROUPS[where(is_time_field)]=time_idx
  non_time_fields=where(~is_time_field)
  non_time_field_names=field_names[non_time_fields]
  tags2remove=where(field_names EQ field_names[time_idx])
  ;; Times
  tfields=where(is_time_field, /NULL)
  tnames=field_names[tfields]
  tnamesl=strsplit(tnames, '_', /extract)
  tnames_last=strarr(n_elements(tnamesl))
  ntbits=n_elements(tnamesl[0])
  tnames_id=strjoin((tnamesl[0])[0:ntbits - 2], '_')
  FOR i=0L, n_elements(tnames) - 1 DO $
     tnames_last[i]=tnamesl[i, n_elements(tnamesl[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  time_locs=locate_time_strings(tnames_last)
  isample_rate_dfrac=double(isample_rate) / double(86400)
  
  idata=read_ascii(ifile, template=itemplate)
  idata_names=strlowcase(tag_names(idata))
  ;; Obtain times and convert to Julian
  idata_time_loc=where(idata_names EQ field_names[time_idx])
  times=idata.(idata_time_loc)
  times_dims=size(times, /dimensions)
  ;; Remove quotes
  IF size(times, /type) EQ 7 THEN BEGIN
     FOREACH fld, indgen(times_dims[0]) DO BEGIN
        ok=strsplit(times[fld, *], '" -/:', /extract)
        ok=(temporary(ok)).toArray()
        ok=strjoin(transpose(temporary(ok)))
        times[fld, *]=ok
     ENDFOREACH
  ENDIF
  match2, idata_names, idata_field_names[tags2remove], is_time
  FOREACH fld, (indgen(n_tags(idata)))[where(is_time LT 0)] DO BEGIN
     IF size(idata.(fld), /type) EQ 7 THEN BEGIN
        ok=strsplit(idata.(fld), '" ', /extract)
        idata.(fld)=ok.toArray()
     ENDIF
  ENDFOREACH
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
                             (isample_rate_dfrac / 2)), $
                mcount)
  IF mcount LT 1 THEN BEGIN
     message, 'No matching records found.  Skipping file.', $
              /informational
     CONTINUE
  ENDIF
  odata=create_struct(tnames[0], reform(times[0, matches]))
  ;; Subset the rest of the time data
  FOREACH fld, (indgen(times_dims[0]))[1:*]  DO BEGIN
     odata=create_struct(odata, tnames[fld], $
                         reform(times[fld, matches]))
  ENDFOREACH
  ;; Add and subset the rest of the data
  FOREACH fld, (indgen(n_tags(idata)))[where(is_time LT 0)] DO BEGIN
     odata=create_struct(odata, idata_names[fld], $
                         idata.(fld)[matches])
  ENDFOREACH
  
  RETURN, odata
END

PRO SUBSET_FLUX, IDIR, ODIR, ITEMPLATE_SAV, TIME_BEG_IDX, ISAMPLE_RATE, $
                 DIAG_DIR, DIAG_ITEMPLATE_SAV, DIAG_TIME_IDX, DIAG_IDX, $
                 EC_PERIOD, RMC_DIR, RMC_ITEMPLATE_SAV, RMC_TIME_IDX, $
                 GYRO_DIR, GYRO_ITEMPLATE_SAV, GYRO_TIME_IDX, $
                 RAD_DIR, RAD_ITEMPLATE_SAV, RAD_TIME_IDX, $
                 OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 19) THEN $
     message, 'Usage: SUBSET_FLUX, IDIR, ODIR, ITEMPLATE_SAV, ' + $
              'TIME_BEG_IDX, ISAMPLE_RATE, DIAG_DIR, DIAG_ITEMPLATE_SAV, ' + $
              'DIAG_TIME_IDX, DIAG_IDX, EC_PERIOD, RMC_DIR, ' + $
              'RMC_ITEMPLATE_SAV, RMC_TIME_IDX, GYRO_DIR, ' + $
              'GYRO_ITEMPLATE_SAV, GYRO_TIME_IDX, ' + $
              'RAD_DIR, RAD_ITEMPLATE_SAV, RAD_TIME_IDX'
  IF ((n_elements(idir) EQ 0) OR (idir EQ '')) THEN $
     message, 'IDIR is undefined or is empty string'
  IF ((n_elements(odir) EQ 0) OR (odir EQ '')) THEN $
     message, 'ODIR is undefined or is empty string'
  IF ((n_elements(itemplate_sav) EQ 0) OR (itemplate_sav EQ '')) THEN $
         message, 'ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(time_beg_idx) NE 1) OR (time_beg_idx LT 0)) THEN $
     message, 'TIME_BEG_IDX must be a scalar >= zero'
  IF ((n_elements(isample_rate) NE 1) OR (isample_rate LT 0)) THEN $
     message, 'ISAMPLE_RATE must be a scalar >= zero'
  IF ((n_elements(diag_dir) EQ 0) OR (diag_dir EQ '')) THEN $
     message, 'DIAG_DIR is undefined or is empty string'
  IF ((n_elements(diag_itemplate_sav) EQ 0) OR $
      (diag_itemplate_sav EQ '')) THEN $
         message, 'DIAG_ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(diag_time_idx) NE 1) OR (diag_time_idx LT 0)) THEN $
     message, 'DIAG_TIME_IDX must be a scalar >= zero'
  IF ((n_elements(diag_idx) NE 1) OR (diag_idx LT 0)) THEN $
     message, 'DIAG_IDX must be a scalar >= zero'
  IF ((n_elements(ec_period) NE 1) OR (ec_period LT 0)) THEN $
     message, 'EC_PERIOD must be a scalar >= zero'
  IF ((86400 MOD float(ec_period)) NE 0) THEN $
     message, 'EC_PERIOD must be an integer divisor of 86400 s'
  IF ((n_elements(rmc_dir) EQ 0) OR (rmc_dir EQ '')) THEN $
     message, 'RMC_DIR is undefined or is empty string'
  IF ((n_elements(rmc_itemplate_sav) EQ 0) OR $
      (rmc_itemplate_sav EQ '')) THEN $
         message, 'RMC_ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(rmc_time_idx) NE 1) OR (rmc_time_idx LT 0)) THEN $
     message, 'RMC_TIME_IDX must be a scalar >= zero'
  IF ((n_elements(gyro_dir) EQ 0) OR (gyro_dir EQ '')) THEN $
     message, 'GYRO_DIR is undefined or is empty string'
  IF ((n_elements(gyro_itemplate_sav) EQ 0) OR $
      (gyro_itemplate_sav EQ '')) THEN $
         message, 'GYRO_ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(gyro_time_idx) NE 1) OR (gyro_time_idx LT 0)) THEN $
     message, 'GYRO_TIME_IDX must be a scalar >= zero'
  IF ((n_elements(rad_dir) EQ 0) OR (rad_dir EQ '')) THEN $
     message, 'RAD_DIR is undefined or is empty string'
  IF ((n_elements(rad_itemplate_sav) EQ 0) OR $
      (rad_itemplate_sav EQ '')) THEN $
         message, 'RAD_ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(rad_time_idx) NE 1) OR (rad_time_idx LT 0)) THEN $
     message, 'RAD_TIME_IDX must be a scalar >= zero'
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

  ;; Parse input flux template
  restore, itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  flux_template=itemplate
  flux_field_names=strlowcase(flux_template.FIELDNAMES)
  flux_field_types=flux_template.FIELDTYPES
  flux_is_time_field=flux_template.FIELDGROUPS EQ time_beg_idx
  ;; Ignore other groups when reading the data
  flux_template.FIELDGROUPS=indgen(flux_template.FIELDCOUNT)
  flux_template.FIELDGROUPS[where(flux_is_time_field)]=time_beg_idx
  flux_non_time_fields=where(~flux_is_time_field)
  flux_non_time_field_names=flux_field_names[flux_non_time_fields]
  flux_tags2remove=where(flux_field_names EQ flux_field_names[time_beg_idx])
  ;; Times
  flux_tfields=where(flux_is_time_field, /NULL)
  flux_tnames=flux_field_names[flux_tfields]
  flux_tnamesl=strsplit(flux_tnames, '_', /extract)
  flux_tnames_last=strarr(n_elements(flux_tnamesl))
  ntbits=n_elements(flux_tnamesl[0])
  flux_tnames_id=strjoin((flux_tnamesl[0])[0:ntbits - 2], '_')
  FOR i=0L, n_elements(flux_tnames) - 1 DO $
     flux_tnames_last[i]=flux_tnamesl[i, n_elements(flux_tnamesl[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  flux_time_locs=locate_time_strings(flux_tnames_last)
  ;; Break file names and extract the piece to match
  flux_filesl=strsplit(idir_files, '_.', /extract)
  flux_files_a=flux_filesl.toArray(/transpose) ; array with pieces in cols
  flux_files_mstr_dims=size(flux_files_a, /dimensions)
  ;; We want the 2nd to last piece
  flux_files_mstr=flux_files_a[flux_files_mstr_dims[0] - 2, *]

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
  rmc_tags2remove=where(rmc_field_names EQ rmc_field_names[rmc_time_idx])
  ;; Times
  rmc_tfields=where(rmc_is_time_field, /NULL)
  rmc_tnames=rmc_field_names[rmc_tfields]
  rmc_tnamesl=strsplit(rmc_tnames, '_', /extract)
  rmc_tnames_last=strarr(n_elements(rmc_tnamesl))
  ntbits=n_elements(rmc_tnamesl[0])
  rmc_tnames_id=strjoin((rmc_tnamesl[0])[0:ntbits - 2], '_')
  FOR i=0L, n_elements(rmc_tnames) - 1 DO $
     rmc_tnames_last[i]=rmc_tnamesl[i, n_elements(rmc_tnamesl[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  rmc_time_locs=locate_time_strings(rmc_tnames_last)
  ;; Break file names and extract the piece to match
  rmc_filesl=strsplit(rmc_files, '_.', /extract)
  rmc_files_a=rmc_filesl.toArray(/transpose) ; array with pieces in cols
  rmc_files_mstr_dims=size(rmc_files_a, /dimensions)
  ;; We want the 2nd to last piece
  rmc_files_mstr=rmc_files_a[rmc_files_mstr_dims[0] - 2, *]

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
  gyro_tags2remove=where(gyro_field_names EQ gyro_field_names[gyro_time_idx])
  ;; Times
  gyro_tfields=where(gyro_is_time_field, /NULL)
  gyro_tnames=gyro_field_names[gyro_tfields]
  gyro_tnamesl=strsplit(gyro_tnames, '_', /extract)
  gyro_tnames_last=strarr(n_elements(gyro_tnamesl))
  ntbits=n_elements(gyro_tnamesl[0])
  gyro_tnames_id=strjoin((gyro_tnamesl[0])[0:ntbits - 2], '_')
  FOR i=0L, n_elements(gyro_tnames) - 1 DO $
     gyro_tnames_last[i]=gyro_tnamesl[i, n_elements(gyro_tnamesl[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  gyro_time_locs=locate_time_strings(gyro_tnames_last)
  ;; Break file names and extract the piece to match
  gyro_filesl=strsplit(gyro_files, '_.', /extract)
  gyro_files_a=gyro_filesl.toArray(/transpose) ; array with pieces in cols
  gyro_files_mstr_dims=size(gyro_files_a, /dimensions)
  ;; We want the 2nd to last piece
  gyro_files_mstr=gyro_files_a[gyro_files_mstr_dims[0] - 2, *]

  ;; Parse RAD template
  restore, rad_itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  rad_template=itemplate
  rad_field_names=strlowcase(rad_template.FIELDNAMES)
  rad_field_types=rad_template.FIELDTYPES
  rad_is_time_field=rad_template.FIELDGROUPS EQ rad_time_idx
  ;; Ignore other groups when reading the data
  rad_template.FIELDGROUPS=indgen(rad_template.FIELDCOUNT)
  rad_template.FIELDGROUPS[where(rad_is_time_field)]=rad_time_idx
  rad_non_time_fields=where(~rad_is_time_field)
  rad_non_time_field_names=rad_field_names[rad_non_time_fields]
  rad_tags2remove=where(rad_field_names EQ rad_field_names[rad_time_idx])
  ;; Times
  rad_tfields=where(rad_is_time_field, /NULL)
  rad_tnames=rad_field_names[rad_tfields]
  rad_tnamesl=strsplit(rad_tnames, '_', /extract)
  rad_tnames_last=strarr(n_elements(rad_tnamesl))
  ntbits=n_elements(rad_tnamesl[0])
  rad_tnames_id=strjoin((rad_tnamesl[0])[0:ntbits - 2], '_')
  FOR i=0L, n_elements(rad_tnames) - 1 DO $
     rad_tnames_last[i]=rad_tnamesl[i, n_elements(rad_tnamesl[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  rad_time_locs=locate_time_strings(rad_tnames_last)
  ;; Break file names and extract the piece to match
  rad_filesl=strsplit(rad_files, '_.', /extract)
  rad_files_a=rad_filesl.toArray(/transpose) ; array with pieces in cols
  rad_files_mstr_dims=size(rad_files_a, /dimensions)
  ;; We want the 3rd to last piece
  rad_files_mstr=rad_files_a[rad_files_mstr_dims[0] - 3, *]

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
        rmc=read_ascii(rmc_files[rmc_pair], template=rmc_template)
        rmc_names=strlowcase(tag_names(rmc))
        ;; Obtain times and convert to Julian
        rmc_time_loc=where(rmc_names EQ rmc_field_names[rmc_time_idx])
        rmc_times=rmc.(rmc_time_loc)
        rmc_times_dims=size(rmc_times, /dimensions)
        ;; Remove quotes
        IF size(rmc_times, /type) EQ 7 THEN BEGIN
           FOREACH fld, indgen(rmc_times_dims[0]) DO BEGIN
              ok=strsplit(rmc_times[fld, *], '" -/:', /extract)
              ok=(temporary(ok)).toArray()
              ok=strjoin(transpose(temporary(ok)))
              rmc_times[fld, *]=ok
           ENDFOREACH
        ENDIF
        match2, rmc_names, rmc_field_names[rmc_tags2remove], is_time
        FOREACH fld, (indgen(n_tags(rmc)))[where(is_time LT 0)] DO BEGIN
           IF size(rmc.(fld), /type) EQ 7 THEN BEGIN
              ok=strsplit(rmc.(fld), '" ', /extract)
              rmc.(fld)=ok.toArray()
           ENDIF
        ENDFOREACH
        rtimes_s=rmc_times[5, *]
        rtimes_s=fix(rtimes_s) + $
                 (round((double(rtimes_s) - fix(rtimes_s)) * 10) / 10.0)
        rmc_jd=reform(julday(long(rmc_times[1, *]), $
                             long(rmc_times[2, *]), $
                             long(rmc_times[0, *]), $
                             long(rmc_times[3, *]), $
                             long(rmc_times[4, *]), $
                             double(rtimes_s)))
        ;; We have to subtract 0.5 the input sample rate to protect against
        ;; numerical representation issues in IDL.
        rmc_matches=where(rmc_jd GE bounds[0] AND $
                          rmc_jd LT (bounds[1] - $
                                     (isample_rate_dfrac / 2)), $
                          mcount)
        IF mcount LT 1 THEN BEGIN
           message, 'No matching RMC records found.  Skipping file.', $
                    /informational
           CONTINUE
        ENDIF
        ;; Build a name for output file and check existence
        iname=strsplit(file_basename(rmc_files[rmc_pair]), '.', /extract)
        ;; Get a path for the file, check if it already exists
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
        rmc_period=create_struct(rmc_tnames[0], $
                                 reform(rmc_times[0, rmc_matches]))
        ;; Subset the rest of the time data
        FOREACH fld, (indgen(rmc_times_dims[0]))[1:*]  DO BEGIN
           rmc_period=create_struct(rmc_period, rmc_tnames[fld], $
                                    reform(rmc_times[fld, rmc_matches]))
        ENDFOREACH
        ;; Add and subset the rest of the data
        FOREACH fld, (indgen(n_tags(rmc)))[where(is_time LT 0)] DO BEGIN
           rmc_period=create_struct(rmc_period, rmc_names[fld], $
                                    rmc.(fld)[rmc_matches])
        ENDFOREACH
        write_csv, ofile_name, rmc_period, header=rmc_field_names
        delvar, rmc, rmc_times

        ;; Read matching GYRO file
        gyro_pair=where(gyro_files_mstr EQ dfile_mstr, mcount)
        IF mcount LT 1 THEN BEGIN
           message, 'No matching GYRO file found. Skipping.', /CONTINUE
           CONTINUE
        ENDIF
        gyro=read_ascii(gyro_files[gyro_pair], template=gyro_template)
        gyro_names=strlowcase(tag_names(gyro))
        ;; Obtain times and convert to Julian
        gyro_time_loc=where(gyro_names EQ gyro_field_names[gyro_time_idx])
        gyro_times=gyro.(gyro_time_loc)
        gyro_times_dims=size(gyro_times, /dimensions)
        ;; Remove quotes
        IF size(gyro_times, /type) EQ 7 THEN BEGIN
           FOREACH fld, indgen((size(gyro_times, $
                                     /dimensions))[0]) DO BEGIN
              ok=strsplit(gyro_times[fld, *], '" -/:', /extract)
              ok=(temporary(ok)).toArray()
              ok=strjoin(transpose(temporary(ok)))
              gyro_times[fld, *]=ok
           ENDFOREACH
        ENDIF
        match2, gyro_names, gyro_field_names[gyro_tags2remove], is_time
        FOREACH fld, (indgen(n_tags(gyro)))[where(is_time LT 0)] DO BEGIN
           IF size(gyro.(fld), /type) EQ 7 THEN BEGIN
              ok=strsplit(gyro.(fld), '" ', /extract)
              gyro.(fld)=ok.toArray()
           ENDIF
        ENDFOREACH
        gtimes_s=gyro_times[5, *]
        gtimes_s=fix(gtimes_s) + $
                 (round((double(gtimes_s) - fix(gtimes_s)) * 10) / 10.0)
        gyro_jd=reform(julday(long(gyro_times[1, *]), $
                              long(gyro_times[2, *]), $
                              long(gyro_times[0, *]), $
                              long(gyro_times[3, *]), $
                              long(gyro_times[4, *]), $
                              double(gtimes_s)))
        gyro_matches=where(gyro_jd GE bounds[0] AND $
                           gyro_jd LT (bounds[1] - $
                                       (isample_rate_dfrac / 2)), $
                           mcount)
        IF mcount LT 1 THEN BEGIN
           message, 'No matching GYRO records found.  Skipping file.', $
                    /informational
           CONTINUE
        ENDIF
        ;; Build a name for output file and check existence
        iname=strsplit(file_basename(gyro_files[gyro_pair]), '.', /extract)
        ;; Get a path for the file, check if it already exists
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
        gyro_period=create_struct(gyro_tnames[0], $
                                  reform(gyro_times[0, gyro_matches]))
        ;; Subset the rest of the time data
        FOREACH fld, (indgen(gyro_times_dims[0]))[1:*]  DO BEGIN
           gyro_period=create_struct(gyro_period, gyro_tnames[fld], $
                                     reform(gyro_times[fld, gyro_matches]))
        ENDFOREACH
        ;; Add and subset the rest of the data
        FOREACH fld, (indgen(n_tags(gyro)))[where(is_time LT 0)] DO BEGIN
           gyro_period=create_struct(gyro_period, gyro_names[fld], $
                                     gyro.(fld)[gyro_matches])
        ENDFOREACH
        write_csv, ofile_name, gyro_period, header=gyro_field_names
        delvar, gyro, gyro_times

        ;; Read matching RAD file
        rad_pair=where(rad_files_mstr EQ dfile_mstr, mcount)
        IF mcount LT 1 THEN BEGIN
           message, 'No matching RAD file found. Skipping.', /CONTINUE
           CONTINUE
        ENDIF
        rad=read_ascii(rad_files[rad_pair], template=rad_template)
        rad_names=strlowcase(tag_names(rad))
        ;; Obtain times and convert to Julian
        rad_time_loc=where(rad_names EQ rad_field_names[rad_time_idx])
        rad_times=rad.(rad_time_loc)
        rad_times_dims=size(rad_times, /dimensions)
        ;; Remove quotes
        IF size(rad_times, /type) EQ 7 THEN BEGIN
           FOREACH fld, indgen((size(rad_times, $
                                     /dimensions))[0]) DO BEGIN
              ok=strsplit(rad_times[fld, *], '" -/:', /extract)
              ok=(temporary(ok)).toArray()
              ok=strjoin(transpose(temporary(ok)))
              rad_times[fld, *]=ok
           ENDFOREACH
        ENDIF
        match2, rad_names, rad_field_names[rad_tags2remove], is_time
        FOREACH fld, (indgen(n_tags(rad)))[where(is_time LT 0)] DO BEGIN
           IF size(rad.(fld), /type) EQ 7 THEN BEGIN
              ok=strsplit(rad.(fld), '" ', /extract)
              rad.(fld)=ok.toArray()
           ENDIF
        ENDFOREACH
        rdtimes_s=rad_times[5, *]
        rdtimes_s=fix(rdtimes_s) + $
                  (round((double(rdtimes_s) - fix(rdtimes_s)) * 10) / 10.0)
        rad_jd=reform(julday(long(rad_times[1, *]), $
                             long(rad_times[2, *]), $
                             long(rad_times[0, *]), $
                             long(rad_times[3, *]), $
                             long(rad_times[4, *]), $
                             double(rdtimes_s)))
        rad_matches=where(rad_jd GE bounds[0] AND $
                          rad_jd LT (bounds[1] - $
                                     (isample_rate_dfrac / 2)), $
                           mcount)
        IF mcount LT 1 THEN BEGIN
           message, 'No matching RAD records found.  Skipping file.', $
                    /informational
           CONTINUE
        ENDIF
        ;; Build a name for output file and check existence
        iname=strsplit(file_basename(rad_files[rad_pair]), '.', /extract)
        ;; Get a path for the file, check if it already exists
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
        rad_period=create_struct(rad_tnames[0], $
                                 reform(rad_times[0, rad_matches]))
        ;; Subset the rest of the time data
        FOREACH fld, (indgen(rad_times_dims[0]))[1:*]  DO BEGIN
           rad_period=create_struct(rad_period, rad_tnames[fld], $
                                    reform(rad_times[fld, rad_matches]))
        ENDFOREACH
        ;; Add and subset the rest of the data
        FOREACH fld, (indgen(n_tags(rad)))[where(is_time LT 0)] DO BEGIN
           rad_period=create_struct(rad_period, rad_names[fld], $
                                    rad.(fld)[rad_matches])
        ENDFOREACH
        write_csv, ofile_name, rad_period, header=rad_field_names
        delvar, rad, rad_times

        ;; Read matching flux EC file
        flux_pair=where(flux_files_mstr EQ dfile_mstr, mcount)
        IF mcount LT 1 THEN BEGIN
           message, 'No matching EC flux file found. Skipping.', /CONTINUE
           CONTINUE
        ENDIF
        flux=read_ascii(idir_files[flux_pair], template=flux_template)
        flux_names=strlowcase(tag_names(flux))
        ;; Obtain times and convert to Julian
        flux_time_loc=where(flux_names EQ flux_field_names[time_beg_idx])
        flux_times=flux.(flux_time_loc)
        flux_times_dims=size(flux_times, /dimensions)
        ;; Remove quotes
        IF size(flux_times, /type) EQ 7 THEN BEGIN
           FOREACH fld, indgen((size(flux_times, $
                                     /dimensions))[0]) DO BEGIN
              ok=strsplit(flux_times[fld, *], '" -/:', /extract)
              ok=(temporary(ok)).toArray()
              ok=strjoin(transpose(temporary(ok)))
              flux_times[fld, *]=ok
           ENDFOREACH
        ENDIF
        match2, flux_names, flux_field_names[flux_tags2remove], is_time
        FOREACH fld, (indgen(n_tags(flux)))[where(is_time LT 0)] DO BEGIN
           IF size(flux.(fld), /type) EQ 7 THEN BEGIN
              ok=strsplit(flux.(fld), '" ', /extract)
              flux.(fld)=ok.toArray()
           ENDIF
        ENDFOREACH
        ftimes_s=flux_times[5, *]
        ftimes_s=fix(ftimes_s) + $
                 (round((double(ftimes_s) - fix(ftimes_s)) * 10) / 10.0)
        flux_jd=reform(julday(long(flux_times[1, *]), $
                              long(flux_times[2, *]), $
                              long(flux_times[0, *]), $
                              long(flux_times[3, *]), $
                              long(flux_times[4, *]), $
                              double(ftimes_s)))
        flux_matches=where(flux_jd GE bounds[0] AND $
                           flux_jd LT (bounds[1] - $
                                       (isample_rate_dfrac / 2)), $
                           mcount)
        IF mcount LT 1 THEN BEGIN
           message, 'No matching EC flux records found.  Skipping file.', $
                    /informational
           CONTINUE
        ENDIF
        ;; Build a name for output file and check existence
        iname=strsplit(file_basename(idir_files[flux_pair]), '.', /extract)
        ;; Get a path for the file, check if it already exists
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
        flux_period=create_struct(flux_tnames[0], $
                                  reform(flux_times[0, flux_matches]))
        ;; Subset the rest of the time data
        FOREACH fld, (indgen(flux_times_dims[0]))[1:*]  DO BEGIN
           flux_period=create_struct(flux_period, flux_tnames[fld], $
                                     reform(flux_times[fld, flux_matches]))
        ENDFOREACH
        ;; Add and subset the rest of the data
        FOREACH fld, (indgen(n_tags(flux)))[where(is_time LT 0)] DO BEGIN
           flux_period=create_struct(flux_period, flux_names[fld], $
                                     flux.(fld)[flux_matches])
        ENDFOREACH
        write_csv, ofile_name, flux_period, header=flux_field_names
        delvar, flux, flux_times

     ENDFOREACH

  ENDFOR


END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; subset_flux.pro ends here
