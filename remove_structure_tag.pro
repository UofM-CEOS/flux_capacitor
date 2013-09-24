;; $Id: $
;;; remove_structure_tag.pro --- Remove tag from structure
;; Author: Sebastian Luque
;; Created: 2013-09-24T15:38:09+0000
;; Last-Updated: 2013-09-24T15:41:00+0000
;;           By: Sebastian Luque
;; + ----------------------------------------------------------------------
;; NAME:
;; 
;;  remove_structure_tag
;; 
;; PURPOSE:
;; 
;;  Remove tag(s) from a structure
;; 
;; CATEGORY:
;; 
;;  Programming and Control
;; 
;; CALLING SEQUENCE:
;; 
;;  Result=REMOVE_STRUCTURE_TAG(Struct, Tag_Name)
;; 
;; INPUTS:
;; 
;;  Struct: Structure.
;;
;;  Tag_Name: String array or scalar.
;; 
;; - ----------------------------------------------------------------------
;;; Code:
FUNCTION REMOVE_STRUCTURE_TAG, STRUCT, TAG_NAME
  COMPILE_OPT idl2
  IF n_params() NE 2 THEN MESSAGE, 'N_PARAMS=2'
  
  search_tag=strupcase(tag_name)
  tag_names=tag_names(struct)
  a=[-1]
  IF n_elements(tag_name) EQ 1 THEN BEGIN 
     a=[a, where(tag_names NE search_tag)] 
  ENDIF ELSE BEGIN
     FOR i=0, n_elements(tag_names) - 1 DO $
        IF (where(search_tag EQ tag_names[i]))[0] EQ -1 THEN a=[a, i]
  ENDELSE
  IF (n_elements(a) EQ 1) OR (a[1] EQ -1) THEN RETURN, struct
  ostruct=create_struct(tag_names[a[1]], struct.(a[1]))
  IF n_elements(a) GT 2 THEN BEGIN
     FOR i=2, n_elements(a) - 1 DO BEGIN
        ostruct=create_struct(ostruct, tag_names[a[i]], struct.(a[i]))
     ENDFOR
  ENDIF

  RETURN, ostruct
  
END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; remove_structure_tag.pro ends here