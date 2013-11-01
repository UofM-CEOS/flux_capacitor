;; Author: Sebastian Luque
;; Created: 2013-10-31T20:49:18+0000
;; Last-Updated: 2013-11-01T18:51:22+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     REMOVE_FIELDS_ASCII
;; 
;; PURPOSE:
;; 
;;     Remove fields from CSV files in input directory, overwriting the
;;     files therein.  This procedure is meant to remove meaningless fields
;;     when averaging files via AVERAGE_SERIES that have already been
;;     averaged with that same procedure.  In these case, the standard
;;     deviation fields lose their meaning so they need to be removed.
;;     This is the case when calculating the 20 min average for processed
;;     MET files, in order to determine the suitability of periods for flux
;;     analyses.
;; 
;; CALLING SEQUENCE:
;; 
;;     REMOVE_FIELDS_ASCII, Idir, Itemplate_Sav, Remove_Field_Names
;; 
;; INPUTS:
;; 
;;     Idir:                  Input directory (no trailing separator).
;;     Itemplate_Sav:         Ascii template to read input files.
;;     Remove_Field_Names:    String array or scalar with names of fields
;;                            to remove.
;; 
;; EXAMPLE:
;; 
;; 
;; 
;;- -----------------------------------------------------------------------
;;; Code:

PRO REMOVE_FIELDS_ASCII, IDIR, ITEMPLATE_SAV, REMOVE_FIELD_NAMES

  ;; Check parameters
  IF (n_params() NE 3) THEN $
     message, 'Usage: REMOVE_FIELDS_ASCII, IDIR, ITEMPLATE_SAV, ' + $
              'REMOVE_FIELD_NAMES'
  IF ((n_elements(idir) EQ 0) OR (idir EQ '')) THEN $
     message, 'IDIR is undefined or is empty string'
  IF ((n_elements(itemplate_sav) EQ 0) OR (itemplate_sav EQ '')) THEN $
     message, 'ITEMPLATE_SAV is undefined or is empty string'
  IF (n_elements(remove_field_names) EQ 0) THEN $
     message, 'REMOVE_FIELD_NAMES is undefined'
  idir_files=file_search(idir + path_sep() + '*', count=nidir_files, $
                         /nosort, /fold_case, /test_regular)
  IF nidir_files LT 1 THEN message, 'No input files found'

  restore, itemplate_sav
  ;; Ignore groups
  itemplate.FIELDGROUPS=indgen(itemplate.FIELDCOUNT)
  field_types=itemplate.FIELDTYPES
  str_flds=where(field_types EQ 7, nstr_flds)

  FOR k=0, nidir_files - 1 DO BEGIN
     ifile=idir_files[k]
     message, 'Processing ' + ifile, /informational
     ;; Read input file
     idata=read_ascii(ifile, template=itemplate)
     ;; Remove quotes and separators from string fields
     IF nstr_flds GT 0 THEN BEGIN
        FOREACH fld, str_flds DO BEGIN
           ok=strsplit(idata.(fld), '" -/:', /extract)
           ok=(temporary(ok)).toArray()
           ok=strjoin(transpose(temporary(ok)))
           idata.(fld)=ok
        ENDFOREACH
     ENDIF
     odata=remove_structure_tags(idata, remove_field_names)
     odata_names=strlowcase(tag_names(odata))
     write_csv, ifile, odata, header=odata_names
  ENDFOR

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; remove_fields_ascii.pro ends here
