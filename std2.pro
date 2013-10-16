;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-09-26T21:14:01+0000
;; Last-Updated: 2013-10-16T05:36:03+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     STD2
;; 
;; PURPOSE:
;; 
;;     Standardize files with two sets of time stamps (usually UTC and
;;     server) by parsing times and selecting columns to keep from
;;     un-processed files.  It also allows conversion of each field to keep
;;     to a given data type.
;; 
;; CALLING SEQUENCE:
;; 
;;      STD_2TIMES, Idir, Odir, Itemplate_Sav, Utc_Time_Idx, Server_Time_Idx,
;;                  Keep_Fields
;; 
;; INPUTS:
;; 
;;     Idir:                  Input directory (no trailing separator).
;;     Odir:                  Output directory (no trailing separator).
;;     Itemplate_Sav:         Ascii template to read input files.
;;     Utc_Time_Idx:          Index (in template) where UTC time is.
;;     Server_Time_Idx:       Index (in template) where server time is.
;;     Keep_Fields:           String array or scalar with names of fields
;;                            to keep.
;; 
;; KEYWORD PARAMETERS:
;; 
;;     Keep_Types:            Array of data type codes for each element of
;;                            Keep_Fields.
;;     File_Type:             String scalar with file type (NAV, RMC, RAD, etc.)
;;     OVERWRITE:             Whether to overwrite files in Odir.
;; 
;; SIDE EFFECTS:
;; 
;;     Writes files in Odir.
;; 
;; EXAMPLE:
;; 
;; nav_raw_keep_fields=['prog_version', 'latitude', 'longitude', $
;;                      'sog', 'cog', 'heading', 'pitch', 'roll', $
;;                      'accel_x', 'accel_y', 'accel_z']
;; nav_raw_keep_types=[7, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
;; STD_2TIMES, expand_path('~/tmp/ArcticNet2011/NAV'), $
;;             expand_path('~/tmp/ArcticNet2011/NAV/STD'), $
;;             'nav_raw_template.sav', 6, 1, nav_raw_keep_fields, $
;;             keep_types=nav_raw_keep_types, file_type='NAV', /overwrite
;; 
;;- -----------------------------------------------------------------------
;;; Code:

PRO STD2, IDIR, ODIR, ITEMPLATE_SAV, UTC_TIME_IDX, SERVER_TIME_IDX, $
          KEEP_FIELDS, KEEP_TYPES=KEEP_TYPES, FILE_TYPE=FILE_TYPE, $
          OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 6) THEN $
     message, 'Usage: STD_2TIMES, IDIR, ODIR, ITEMPLATE_SAV, ' + $
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
  n_kf=n_elements(keep_fields)
  IF (n_kf EQ 0) THEN message, 'KEEP_FIELDS is undefined'
  n_kt=n_elements(keep_types)
  IF (n_kt GT 0) AND (n_kt NE n_kf) THEN $
     message, 'KEEP_TYPES and KEEP_FIELDS must have the same number of elements'
  IF ((n_elements(file_type) EQ 0) OR (odir EQ '')) THEN $
     message, 'FILE_TYPE is undefined or is empty string'
  idir_files=file_search(idir + path_sep() + '*', count=nidir_files, $
                         /nosort, /fold_case, /test_regular)
  IF nidir_files LT 1 THEN message, 'No input files found'

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
     idata_times=idata.(where(idata_names EQ $
                              strlowcase(field_names[utc_time_idx])))
     utc_dims=size(idata_times, /dimensions)
     valid_flag=make_array(utc_dims[1], type=2, value=1)
     srv_time_loc=where(idata_names EQ $
                        strlowcase(field_names[server_time_idx]))
     idata_times_srv=idata.(srv_time_loc)
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
              valid_flag[ibad]=0 ; bad UTC times -> invalid record
           ENDIF
        ENDFOREACH
     ENDIF
     IF size(idata_times_srv, /type) EQ 7 THEN BEGIN
        FOREACH fld, indgen((size(idata_times_srv, $
                                  /dimensions))[0]) DO BEGIN
           ok=strsplit(idata_times_srv[fld, *], '" -/:', /extract)
           ok=(temporary(ok)).toArray()
           ok=strjoin(transpose(temporary(ok)))
           idata_times_srv[fld, *]=ok
           is_valid_srv=valid_num(idata_times_srv[fld, *])
           ibad_srv=where(~is_valid_srv, bcount)
           IF bcount GT 0 THEN $
              idata_times_srv[fld, ibad_srv]=''
        ENDFOREACH
     ENDIF
     match2, idata_names, strlowcase(field_names[tags2remove]), is_time
     FOREACH fld, (indgen(n_tags(idata)))[where(is_time LT 0)] DO BEGIN
        IF size(idata.(fld), /type) EQ 7 THEN BEGIN
           ok=strsplit(idata.(fld), '" ', /extract)
           idata.(fld)=ok.toArray()
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
     ;; Re-check validity and subset if needed
     ibad=where(~valid_flag, bcount, complement=ok, ncomplement=nok)
     IF nok LT 1 THEN BEGIN
        message, 'All UTC time stamps are invalid.  Skipping this file', $
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
        idata_times_srv=idata_times_srv[*, ok]
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
     ;; Obtain full UTC time details
     itimes_utc_std=parse_times(idata_times, tnames_last, time_locs)
     ;; Obtain full server time details
     itimes_srv_std=parse_times(idata_times_srv, tnames_last_srv, $
                                time_locs_srv)

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

     ;; Fix things for certain file types
     CASE file_type OF
        'NAV': BEGIN
           ;; Original comment: check for bad compass data showing up as
           ;; missing accel_z
           badcompass=where(~ finite(odata.accel_z), nbadcompass)
           IF nbadcompass GT 0 THEN BEGIN
              badflds=['heading', 'pitch', 'roll', $
                       'accel_x', 'accel_y', 'accel_z']
              FOREACH fld, badflds DO BEGIN
                 match_fld=where(strupcase(fld) EQ tag_names(odata))
                 odata.(match_fld)[badcompass]=!VALUES.F_NAN
              ENDFOREACH
           ENDIF
           badlat=where((odata.latitude LT 0.0) OR $
                        (odata.latitude GT 90.0), nbadlat)
           badlon=where((odata.longitude LT -180.0) OR $
                        (odata.longitude GT 0.0), nbadlon)
           badsog=where((odata.sog LT 0.0) OR $
                        (odata.sog GT 20.0), nbadsog)
           badcog=where((odata.cog LT 0.0) OR $
                        (odata.cog GT 360.0), nbadcog)
           badhead=where((odata.heading LT 0.0) OR $
                         (odata.heading GT 360.0), nbadhead)
           badroll=where((odata.roll LT -3.0) OR $
                         (odata.roll GT 3.0), nbadroll)
           badpitch=where((odata.pitch LT -3.0) OR $
                          (odata.pitch GT 3.0), nbadpitch)
           IF nbadlat GT 0 THEN odata.latitude[badlat]=!VALUES.F_NAN
           IF nbadlon GT 0 THEN odata.longitude[badlon]=!VALUES.F_NAN
           IF nbadsog GT 0 THEN odata.sog[badsog]=!VALUES.F_NAN
           IF nbadcog GT 0 THEN odata.cog[badcog]=!VALUES.F_NAN
           IF nbadhead GT 0 THEN odata.heading[badhead]=!VALUES.F_NAN
           IF nbadroll GT 0 THEN odata.roll[badroll]=!VALUES.F_NAN
           IF nbadpitch GT 0 THEN odata.pitch[badpitch]=!VALUES.F_NAN
        END
        'RMC': BEGIN
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
;;; std2.pro ends here
