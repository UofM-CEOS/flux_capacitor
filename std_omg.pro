;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-09-26T18:38:38+0000
;; Last-Updated: 2013-10-08T18:32:04+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     STD_OMG
;; 
;; PURPOSE:
;; 
;;     Standardize OMG files.  It will likely standardize other files in
;;     the future.
;; 
;; CALLING SEQUENCE:
;; 
;;     STD_OMG, Idir, Odir, Itemplate_Sav, Time_Beg_Idx, Keep_Fields
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
;;     OVERWRITE:             Whether to overwrite files in Odir.
;; 
;; SIDE EFFECTS:
;; 
;;     Writes files in Odir.
;; 
;; EXAMPLE:
;; 
;;     omg_nav_raw_keep_fields=['latitude', 'longitude', 'sog', 'cog']
;;     omg_hdg_raw_keep_fields='heading'
;;     STD_OMG, expand_path('~/tmp/ArcticNet2011/OMG/NAV'), $
;;              expand_path('~/tmp/ArcticNet2011/OMG/NAV/STD'), $
;;              'omg_nav_raw_template.sav', 0, omg_raw_keep_fields
;;     STD_OMG, expand_path('~/tmp/ArcticNet2011/OMG/HDG'), $
;;              expand_path('~/tmp/ArcticNet2011/OMG/HDG/STD'), $
;;              'omg_hdg_raw_template.sav', 0, omg_raw_keep_fields
;; 
;;- -----------------------------------------------------------------------
;;; Code:

PRO STD_OMG, IDIR, ODIR, ITEMPLATE_SAV, TIME_BEG_IDX, KEEP_FIELDS, $
             OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 5) THEN $
     message, 'Usage: STD_OMG, IDIR, ODIR, ITEMPLATE_SAV, ' + $
              'TIME_BEG_IDX, KEEP_FIELDS'
  IF ((n_elements(idir) EQ 0) OR (idir EQ '')) THEN $
     message, 'IDIR is undefined or is empty string'
  IF ((n_elements(odir) EQ 0) OR (odir EQ '')) THEN $
     message, 'ODIR is undefined or is empty string'
  IF ((n_elements(itemplate_sav) EQ 0) OR (itemplate_sav EQ '')) THEN $
     message, 'ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(time_beg_idx) NE 1) OR (time_beg_idx LT 0)) THEN $
     message, 'TIME_BEG_IDX must be a scalar >= zero'
  IF (n_elements(keep_fields) EQ 0) THEN $
     message, 'KEEP_FIELDS is undefined'
  idir_files=file_search(idir + path_sep() + '*', count=nidir_files, $
                         /nosort, /fold_case, /test_regular)
  IF nidir_files LT 1 THEN $
     message, 'No input files found'

  restore, itemplate_sav
  field_names=itemplate.FIELDNAMES
  n_ifields=itemplate.FIELDCOUNT ; N fields in template
  tags2remove=where(field_names EQ field_names[time_beg_idx])
  ;; Times
  tfields=where(itemplate.FIELDGROUPS EQ time_beg_idx, /NULL)
  tnames=strlowcase(field_names[tfields])
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
     time_loc=where(idata_names EQ $
                    strlowcase(field_names[time_beg_idx]))
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
     match2, idata_names, strlowcase(field_names[tags2remove]), is_time
     FOREACH fld, (indgen(n_tags(idata)))[where(is_time LT 0)] DO BEGIN
        IF size(idata.(fld), /type) EQ 7 THEN BEGIN
           ok=strsplit(idata.(fld), '" ', /extract)
           idata.(fld)=ok.toArray()
        ENDIF
     ENDFOREACH

     ;; Obtain full time details
     itimes_std=parse_times(idata_times, tnames_last, time_locs)
     
     odata=remove_structure_tags(idata, field_names[tags2remove])
     ;; Find indices to keep
     match2, strlowcase(tag_names(odata)), keep_fields, toss, keep
     tags2remove_odata=where(toss LT 0, nremove)
     IF nremove GT 0 THEN $
        odata=remove_structure_tags(odata, $
                                    (tag_names(odata))[tags2remove_odata])
     odata=create_struct('year', reform(itimes_std[0, *]), $
                         'month', reform(itimes_std[1, *]), $
                         'day', reform(itimes_std[2, *]), $
                         'hour', reform(itimes_std[3, *]), $
                         'minute', reform(itimes_std[4, *]), $
                         'second', reform(itimes_std[5, *]), odata)
     delvar, idata

     write_csv, ofile_name, odata, header=strlowcase(tag_names(odata))

  ENDFOR

END



;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; std_omg.pro ends here
