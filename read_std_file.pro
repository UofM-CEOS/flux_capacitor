;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-11-12T19:18:12+0000
;; Last-Updated: 2013-11-22T15:56:22+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     READ_STD_FILE
;; 
;; PURPOSE:
;; 
;;     Reads a standardized file having a group of fields with time data,
;;     using an ASCII template structure.  It returns a structure like the
;;     ones returned by READ_ASCII, but does cleaning of quotes,
;;     separators, and spaces, so the structure is simpler to use later.
;; 
;; CALLING SEQUENCE:
;; 
;;     idata=read_std_file(Ifile, Itemplate, Time_Idx)
;; 
;; INPUTS:
;; 
;;     Ifile:         Path of the CSV file to read.
;;     Itemplate:     ASCII template, as returned by ASCII_TEMPLATE.
;;     Time_Idx:      Index (in template) where time matrix is located.
;; 
;; OUTPUTS:
;; 
;;     It returns a structure like the ones returned by READ_ASCII.
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

FUNCTION READ_STD_FILE, IFILE, ITEMPLATE, TIME_IDX

  ;; Parse input template
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
  ;; Obtain full time info
  times_std=parse_times(times, tnames_last, time_locs)
  match2, idata_names, field_names[tags2remove], is_time
  FOREACH fld, (indgen(n_tags(idata)))[where(is_time LT 0)] DO BEGIN
     IF size(idata.(fld), /type) EQ 7 THEN BEGIN
        ok=strsplit(idata.(fld), '" ', /extract)
        idata.(fld)=ok.toArray()
     ENDIF
  ENDFOREACH

  ostruct=create_struct(tnames[0], times_std)
  ;; Add the rest of the data
  FOREACH fld, (indgen(n_tags(idata)))[where(is_time LT 0)] DO BEGIN
     ostruct=create_struct(ostruct, idata_names[fld], idata.(fld))
  ENDFOREACH

  RETURN, OSTRUCT

END


;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     READ_STD2_FILE
;; 
;; PURPOSE:
;; 
;;     This works as READ_STD_FILE, but takes files with matrices of time
;;     data.
;; 
;; CALLING SEQUENCE:
;; 
;;     idata=read_std2_file(Ifile, Itemplate, Time1_Idx, Time2_Idx)
;; 
;; INPUTS:
;; 
;;     Ifile:         Path of the CSV file to read.
;;     Itemplate:     ASCII template, as returned by ASCII_TEMPLATE.
;;     Time1_Idx:     Index (in template) where the first time matrix is
;;                    located.
;;     Time2_Idx:     Index (in template) where the second time matrix is
;;                    located.
;; 
;; OUTPUTS:
;; 
;;     It returns a structure like the ones returned by READ_ASCII.
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

FUNCTION READ_STD2_FILE, IFILE, ITEMPLATE, TIME1_IDX, TIME2_IDX

  ;; Parse input template
  field_names=strlowcase(itemplate.FIELDNAMES)
  field_types=itemplate.FIELDTYPES
  is_time1_field=itemplate.FIELDGROUPS EQ time1_idx
  is_time2_field=itemplate.FIELDGROUPS EQ time2_idx
  non_time_fields=where((~is_time1_field) AND (~is_time2_field))
  non_time_field_names=field_names[non_time_fields]
  tags2remove=where((field_names EQ field_names[time1_idx]) OR $
                    (field_names EQ field_names[time2_idx]))
  ;; Times 1
  tfields1=where(is_time1_field, /NULL)
  tnames1=field_names[tfields1]
  tnamesl1=strsplit(tnames1, '_', /extract)
  tnames1_last=strarr(n_elements(tnamesl1))
  ntbits=n_elements(tnamesl1[0])
  tnames1_id=strjoin((tnamesl1[0])[0:ntbits - 2], '_')
  FOR i=0L, n_elements(tnames1) - 1 DO $
     tnames1_last[i]=tnamesl1[i, n_elements(tnamesl1[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  time1_locs=locate_time_strings(tnames1_last)
  ;; Times 2
  tfields2=where(is_time2_field, /NULL)
  tnames2=field_names[tfields2]
  tnamesl2=strsplit(tnames2, '_', /extract)
  tnames2_last=strarr(n_elements(tnamesl2))
  ntbits=n_elements(tnamesl2[0])
  tnames2_id=strjoin((tnamesl2[0])[0:ntbits - 2], '_')
  FOR i=0L, n_elements(tnames2) - 1 DO $
     tnames2_last[i]=tnamesl2[i, n_elements(tnamesl2[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  time2_locs=locate_time_strings(tnames2_last)
  
  idata=read_ascii(ifile, template=itemplate)
  idata_names=strlowcase(tag_names(idata))
  ;; Obtain times1
  idata_time1_loc=where(idata_names EQ field_names[time1_idx])
  times1=idata.(idata_time1_loc)
  times1_dims=size(times1, /dimensions)
  ;; Remove quotes
  IF size(times1, /type) EQ 7 THEN BEGIN
     FOREACH fld, indgen(times1_dims[0]) DO BEGIN
        ok=strsplit(times1[fld, *], '" -/:', /extract)
        ok=(temporary(ok)).toArray()
        ok=strjoin(transpose(temporary(ok)))
        times1[fld, *]=ok
     ENDFOREACH
  ENDIF
  ;; Obtain full time info
  times1_std=parse_times(times1, tnames1_last, time1_locs)
  ;; Obtain times2
  idata_time2_loc=where(idata_names EQ field_names[time2_idx])
  times2=idata.(idata_time2_loc)
  times2_dims=size(times2, /dimensions)
  ;; Remove quotes
  IF size(times2, /type) EQ 7 THEN BEGIN
     FOREACH fld, indgen(times2_dims[0]) DO BEGIN
        ok=strsplit(times2[fld, *], '" -/:', /extract)
        ok=(temporary(ok)).toArray()
        ok=strjoin(transpose(temporary(ok)))
        times2[fld, *]=ok
     ENDFOREACH
  ENDIF
  times2_std=parse_times(times2, tnames2_last, time2_locs)
  match2, idata_names, field_names[tags2remove], is_time
  FOREACH fld, (indgen(n_tags(idata)))[where(is_time LT 0)] DO BEGIN
     IF size(idata.(fld), /type) EQ 7 THEN BEGIN
        ok=strsplit(idata.(fld), '" ', /extract)
        idata.(fld)=ok.toArray()
     ENDIF
  ENDFOREACH

  ostruct=create_struct(tnames1[0], times1_std)
  ostruct=create_struct(ostruct, tnames2[0], times2_std)
  ;; Add the rest of the data
  FOREACH fld, (indgen(n_tags(idata)))[where(is_time LT 0)] DO BEGIN
     ostruct=create_struct(ostruct, idata_names[fld], idata.(fld))
  ENDFOREACH

  RETURN, OSTRUCT

END



;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; read_std_file.pro ends here
