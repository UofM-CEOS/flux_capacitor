;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-10-11T14:48:45+0000
;; Last-Updated: 2013-10-24T18:19:37+0000
;;           By: Sebastian Luque
;; 
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     PROCESS_RAD
;; 
;; PURPOSE:
;; 
;;     Process RAD data.
;; 
;; CALLING SEQUENCE:
;; 
;;     PROCESS_RAD, Idir, Odir, Itemplate_Sav, Time_Beg_Idx, RMC_Dir, $
;;                  RMC_Itemplate_Sav, RMC_Time_Idx, RMC_LonLat_Idx, $
;;                  MET_Dir, MET_Itemplate_Sav, MET_Time_Idx, MET_PAR_Idx
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
;;     RMC_LonLat_Idx:      Integer array with indices (in template) where
;;                          latitude and longitude are located.
;;     MET_Dir:             Directory where MET files are found (no
;;                          trailing separator).
;;     MET_Itemplate_Sav:   Ascii template to read MET files.
;;     MET_Time_Idx:        Index (in template) where time is in MET files.
;;     MET_PAR_Idx:         Integer with index (in template) where PAR
;;                          field is located.
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
;;     PROCESS_RAD, expand_path('~/tmp/ArcticNet2011/RAD/Daily'), $
;;                  expand_path('~/tmp/ArcticNet2011/RAD/Processed'), $
;;                  'rad_std_template.sav', 0, $
;;                  expand_path('~/tmp/ArcticNet2011/RMC/1min'), $
;;                  'rmc_avg_template.sav', 0, [6, 7], $
;;                  expand_path('~/tmp/ArcticNet2011/MET/Daily'), $
;;                  'met_std_template.sav', 0, 16, /overwrite
;; 
;;- -----------------------------------------------------------------------
;;; Code:

PRO PROCESS_RAD, IDIR, ODIR, ITEMPLATE_SAV, TIME_BEG_IDX, $
                 RMC_DIR, RMC_ITEMPLATE_SAV, RMC_TIME_IDX, $
                 RMC_LONLAT_IDX, MET_DIR, MET_ITEMPLATE_SAV, MET_TIME_IDX, $
                 MET_PAR_IDX, OVERWRITE=OVERWRITE


  idir_files=file_search(idir + path_sep() + '*', count=nidir_files, $
                         /nosort, /fold_case, /test_regular)
  IF nidir_files LT 1 THEN message, 'No input files found'
  rmc_files=file_search(rmc_dir + path_sep() + '*', count=nrmc_files, $
                        /nosort, /fold_case, /test_regular)
  IF nrmc_files LT 1 THEN message, 'No RMC files found.  Exiting'
  met_files=file_search(met_dir + path_sep() + '*', count=nmet_files, $
                        /nosort, /fold_case, /test_regular)
  IF nmet_files LT 1 THEN message, 'No MET files found.  Exiting'

  ;; Parse input template
  restore, itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  rad_template=itemplate
  field_names=strlowcase(rad_template.FIELDNAMES)
  field_types=rad_template.FIELDTYPES
  is_time_field=rad_template.FIELDGROUPS EQ time_beg_idx
  ;; Ignore other groups when reading the data
  rad_template.FIELDGROUPS=indgen(rad_template.FIELDCOUNT)
  rad_template.FIELDGROUPS[where(is_time_field)]=time_beg_idx
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

  ;; Parse MET template
  restore, met_itemplate_sav
  ;; Make a copy to avoid name collisions with the other templates
  met_template=itemplate
  met_field_names=strlowcase(met_template.FIELDNAMES)
  met_field_types=met_template.FIELDTYPES
  met_is_time_field=met_template.FIELDGROUPS EQ met_time_idx
  ;; Ignore other groups when reading the data
  met_template.FIELDGROUPS=indgen(met_template.FIELDCOUNT)
  met_template.FIELDGROUPS[where(met_is_time_field)]=met_time_idx
  met_non_time_fields=where(~met_is_time_field)
  met_non_time_field_names=met_field_names[met_non_time_fields]
  ;; Times
  met_tfields=where(met_is_time_field, /NULL)
  met_tnames=met_field_names[met_tfields]
  met_tnamesl=strsplit(met_tnames, '_', /extract)
  met_tnames_last=strarr(n_elements(met_tnamesl))
  met_tnames_id=strjoin((met_tnamesl[0])[0:n_elements(met_tnamesl[0]) - 2], '_')
  FOR i=0L, n_elements(tnames) - 1 DO $
     met_tnames_last[i]=met_tnamesl[i, n_elements(met_tnamesl[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  met_time_locs=locate_time_strings(met_tnames_last)
  ;; Break file names and extract the piece to match
  met_filesl=strsplit(met_files, '_.', /extract)
  met_files_a=met_filesl.toArray(/transpose) ; array with pieces in cols
  met_files_mstr_dims=size(met_files_a, /dimensions)
  ;; We want the 2nd to last piece
  met_files_mstr=met_files_a[met_files_mstr_dims[0] - 2, *]

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
     idata=read_ascii(ifile, template=rad_template)
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
     rad_jd=reform(julday(long(itimes_std[1, *]), $
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
     match2, rad_jd, rmc_jd, rad_in_rmc, rmc_in_rad
     rad_matches=where(rad_in_rmc GE 0, mcount, /null)
     IF mcount LT 1 THEN BEGIN
        message, 'No matching RMC records found.  Skipping file.', $
                 /continue
        CONTINUE
     ENDIF
     rmc_matches=where(rmc_in_rad GE 0)
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
                 /continue
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

     ;; Deal with out of range values.  [SPL: Note that we rely on field
     ;; names here, so be WATCH the names in the template.]
     ;; For shortwave -- < 0 or > 800
     ;; For longwave  -- < 0 or > 550
     ;; For PAR       -- < 0 or > 2000
     bad_sw=where(idata.k_down LT 0.0 OR idata.k_down GT 800.0, nbad_sw)
     bad_lw=where(idata.lw_in LT 0.0 OR idata.lw_in GT 550.0, nbad_lw)
     bad_par1=where(idata.par LT 0.0 OR idata.par GT 2000.0, nbad_par1)
     bad_par2=where(idata.par_met LT 0.0 OR idata.par_met GT 2000.0, $
                    nbad_par2)
     ;; If shortwave value is between 0 and -2, set value to zero,
     ;; otherwise set matching bad_sw records to 'NaN'
     IF nbad_sw GT 0 THEN BEGIN
        setzero=where((idata.k_down[bad_sw] LT 0) AND $
                      (idata.k_down[bad_sw] GT -2.0), nsetzero, $
                      complement=notzero)
        IF nsetzero GT 0 THEN BEGIN
           idata.k_down[bad_sw[setzero]]=0.0
           IF notzero[0] GT -1 THEN $
              idata.k_down[bad_sw[notzero]]=!VALUES.F_NAN
        ENDIF ELSE BEGIN
           idata.k_down[badsw]=!VALUES.F_NAN
        ENDELSE
     ENDIF
     ;; Set all matching bad_lw records to 'NaN'
     IF nbad_lw GT 0 THEN $
        idata.lw_in[badlw]=!VALUES.F_NAN
     ;; If PAR value is between 0 and -2, set value to zero, otherwise set
     ;; matching bad_par1 records to 'NaN'
     IF nbad_par1 GT 0 THEN BEGIN
        setzero=where((idata.par[bad_par1] LT 0) AND $
                      (idata.par[bad_par1] GT -2.0), nsetzero, $
                      complement=notzero)
        IF nsetzero GT 0 THEN BEGIN
           idata.par[bad_par1[setzero]]=0.0
           IF notzero[0] GT -1 THEN $
              idata.par[bad_par1[notzero]]=!VALUES.F_NAN
        ENDIF ELSE BEGIN
           idata.par[bad_par1]=!VALUES.F_NAN
        ENDELSE
     ENDIF
     ;; If PAR_MET value is between 0 and -2, set value to zero, otherwise
     ;; set matching bad_par2 records to 'NaN'
     IF nbad_par2 GT 0 THEN BEGIN
        setzero=where((idata.par_met[bad_par2] LT 0) AND $
                      (idata.par_met[bad_par2] GT -2.0), nsetzero, $
                      complement=notzero)
        IF nsetzero GT 0 THEN BEGIN
           idata.par_met[bad_par2[setzero]]=0.0
           IF notzero[0] GT -1 THEN $
              idata.par_met[bad_par2[notzero]]=!VALUES.F_NAN
        ENDIF ELSE BEGIN
           idata.par_met[bad_par2]=!VALUES.F_NAN
        ENDELSE
     ENDIF

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
;;; process_RAD.pro ends here
