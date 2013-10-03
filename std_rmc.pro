;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-09-26T21:14:01+0000
;; Last-Updated: 2013-10-03T20:54:21+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
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
;;  UTC_Time_Idx: Index (0-based) of first field with UTC time definition.
;;   This is typically the year field.
;;
;;  GPS_Time_Idx: Index (0-based) of first field with GPS time definition.
;;   This is typically the year field.  NOTE: We assume both UTC and GPS
;;   have the same resolution (i.e. if there's sub-seconds in one, there'd
;;   better be sub-seconds in the other!).
;;
;;  Keep_Fields: String listing the fields (excluding time fields above) to
;;   keeep in output file.  These must match names in Itemplate_Sav.
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
;; rmc_raw_keep_fields=['latitude', 'longitude', 'sog', 'cog']
;; STD_RMC, expand_path('~/tmp/ArcticNet2011/RMC'), $
;;          expand_path('~/tmp/ArcticNet2011/RMC/STD'), $
;;          'rmc_raw_template.sav', 3, 0, rmc_raw_keep_fields, /overwrite
;; 
;;- -----------------------------------------------------------------------
;;; Code:

PRO STD_RMC, IDIR, ODIR, ITEMPLATE_SAV, UTC_TIME_IDX, SERVER_TIME_IDX, $
             KEEP_FIELDS, OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 6) THEN $
     message, 'Usage: STD_MET, IDIR, ODIR, ITEMPLATE_SAV, ' + $
              'UTC_TIME_IDX, SERVER_TIME_IDX, KEEP_FIELDS'
  IF ((n_elements(idir) EQ 0) OR (idir EQ '')) THEN $
     message, 'IDIR is undefined or is empty string'
  IF ((n_elements(odir) EQ 0) OR (odir EQ '')) THEN $
     message, 'ODIR is undefined or is empty string'
  IF ((n_elements(itemplate_sav) EQ 0) OR (itemplate_sav EQ '')) THEN $
     message, 'ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(utc_time_idx) NE 1) OR (utc_time_idx LT 0)) THEN $
     message, 'UTC_TIME_IDX must be a scalar >= zero'
  IF ((n_elements(server_time_idx) NE 1) OR (server_time_idx LT 0)) THEN $
     message, 'SERVER_TIME_IDX must be a scalar  >= zero'
  IF (n_elements(keep_fields) EQ 0) THEN $
     message, 'KEEP_FIELDS is undefined'

  idir_files=file_search(idir + path_sep() + '*.raw', count=nidir_files, $
                         /nosort, /fold_case)
  IF nidir_files LT 1 THEN BEGIN
     message, 'No input files found', /informational
     RETURN
  ENDIF
  restore, itemplate_sav
  field_names=itemplate.FIELDNAMES
  n_ifields=itemplate.FIELDCOUNT ; N fields in template
  tags2remove=where(field_names EQ field_names[utc_time_idx] OR $
                    field_names EQ field_names[server_time_idx])
  tfields=where(itemplate.FIELDGROUPS EQ utc_time_idx, /NULL) ; time fields
  tnames=strlowcase(field_names[tfields])            ; time names
  tnamesl=strsplit(tnames, '_', /extract)                 ; split list
  tnames_last=strarr(n_elements(tnamesl)) ; set up
  tnames_id=strjoin((tnamesl[0])[0:n_elements(tnamesl[0]) - 2], '_')
  FOR i=0L, n_elements(tnames) - 1 DO $
     tnames_last[i]=tnamesl[i, n_elements(tnamesl[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  time_locs=locate_time_strings(tnames_last)
  ;; Server times
  tfields_srv=where(itemplate.FIELDGROUPS EQ server_time_idx, /NULL)
  tnames_srv=strlowcase(field_names[tfields_srv])
  tnamesl_srv=strsplit(tnames_srv, '_', /extract)
  tnames_last_srv=strarr(n_elements(tnamesl_srv))
  tnames_id_srv=strjoin((tnamesl_srv[0])[0:n_elements(tnamesl_srv[0]) - 2], '_')
  FOR i=0L, n_elements(tnames_srv) - 1 DO $
     tnames_last_srv[i]=tnamesl_srv[i, n_elements(tnamesl_srv[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  time_locs_srv=locate_time_strings(tnames_last_srv)
     
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
     idata_times=idata.(where(idata_names EQ $
                              strlowcase(field_names[utc_time_idx])))
     srv_time_loc=where(idata_names EQ $
                        strlowcase(field_names[server_time_idx]))
     idata_times_srv=idata.(srv_time_loc)
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
     IF size(idata_times_srv, /type) EQ 7 THEN BEGIN
        FOREACH fld, indgen((size(idata_times_srv, $
                                  /dimensions))[0]) DO BEGIN
           ok=strsplit(idata_times_srv[fld, *], '" -/:', /extract)
           ok=(temporary(ok)).toArray()
           ok=strjoin(transpose(temporary(ok)))
           idata_times_srv[fld, *]=ok
        ENDFOREACH
     ENDIF
     match2, idata_names, strlowcase(field_names[tags2remove]), is_time
     FOREACH fld, (indgen(n_tags(idata)))[where(is_time LT 0)] DO BEGIN
        IF size(idata.(fld), /type) EQ 7 THEN BEGIN
           ok=strsplit(idata.(fld), '" ', /extract)
           idata.(fld)=ok.toArray()
        ENDIF
     ENDFOREACH
     ;; Obtain full UTC time details
     itimes_utc_std=parse_times(idata_times, tnames_last, time_locs)
     ;; Obtain full server time details
     itimes_srv_std=parse_times(idata_times_srv, tnames_last_srv, $
                                time_locs_srv)
     
     odata=remove_structure_tags(idata, field_names[tags2remove])
     ;; Find indices to keep
     match2, strlowcase(tag_names(odata)), keep_fields, toss, keep
     tags2remove_odata=where(toss LT 0, nremove)
     IF nremove GT 0 THEN $
        odata=remove_structure_tags(odata, $
                                    (tag_names(odata))[tags2remove_odata])
     odata=create_struct(tnames_id + '_year', $
                         reform(itimes_utc_std[0, *]), $
                         tnames_id + '_month', $
                         reform(itimes_utc_std[1, *]), $
                         tnames_id + '_day', $
                         reform(itimes_utc_std[2, *]), $
                         tnames_id + '_hour', $
                         reform(itimes_utc_std[3, *]), $
                         tnames_id + '_minute', $
                         reform(itimes_utc_std[4, *]), $
                         tnames_id + '_second', $
                         reform(itimes_utc_std[5, *]), $
                         tnames_id_srv + '_year', $
                         reform(itimes_srv_std[0, *]), $
                         tnames_id_srv + '_month', $
                         reform(itimes_srv_std[1, *]), $
                         tnames_id_srv + '_day', $
                         reform(itimes_srv_std[2, *]), $
                         tnames_id_srv + '_hour', $
                         reform(itimes_srv_std[3, *]), $
                         tnames_id_srv + '_minute', $
                         reform(itimes_srv_std[4, *]), $
                         tnames_id_srv + '_second', $
                         reform(itimes_srv_std[5, *]), $
                         odata)
     delvar, idata

     ;; Fix things for certain years (first odata field)
     bad2011=where(odata.(0) EQ 2011, nbad, /null)
     IF nbad GT 0 THEN BEGIN
        ;; Fix silly format (DDMM.mm) for encoding degrees
        lat_str=string(odata.latitude[bad2011], format='(f010.5)')
        lon_str=string(odata.longitude[bad2011], format='(f011.5)')
        lat_deg=float(strmid(lat_str, 0, 2))
        lat_min=float(strmid(lat_str, 2, 8)) / 60
        odata.latitude[bad2011]=lat_deg + lat_min
        lon_deg=float(strmid(lon_str, 0, 3))
        lon_min=float(strmid(lon_str, 3, 8)) / 60
        odata.longitude[bad2011]=-1 * (lon_deg + lon_min)
     ENDIF

     write_csv, ofile_name, odata, header=strlowcase(tag_names(odata))

  ENDFOR

END



;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; std_rmc.pro ends here
