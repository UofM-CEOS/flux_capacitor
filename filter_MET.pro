;; $Id: $
;; Author: Sebastian Luque
;; Created: 2013-10-29T14:29:43+0000
;; Last-Updated: 2013-11-01T22:17:39+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     FILTER_MET
;; 
;; PURPOSE:
;; 
;;     This procedure is used to check whether processed MET data within a
;;     particular time period (20 min) meets the conditions required for
;;     eddy covariance calculations.  It generates a file with requested
;;     fields from the processed MET file, consisting of the mean values
;;     within the requested period, starting from a requested start time to
;;     and end time.
;; 
;; CALLING SEQUENCE:
;; 
;; 
;; 
;; INPUTS:
;; 
;; 
;; 
;; OPTIONAL INPUTS:
;; 
;; 
;; 
;; KEYWORD PARAMETERS:
;; 
;; 
;; 
;; OUTPUTS:
;; 
;; 
;; 
;; OPTIONAL OUTPUTS:
;; 
;; 
;; 
;; COMMON BLOCKS:
;; 
;; 
;; 
;; SIDE EFFECTS:
;; 
;; 
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

PRO FILTER_MET, FULL_DIR, AVG_DIR, FULL_ITEMPLATE_SAV, $
                FULL_TIME_IDX, AVG_ITEMPLATE_SAV, AVG_TIME_IDX, $
                AVG_PERIOD, FILTER_IDX, FILTER_THR, OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 9) THEN $
     message, 'Usage: FILTER_MET, FULL_DIR, AVG_DIR, ' + $
              'FULL_ITEMPLATE_SAV, FULL_TIME_IDX, AVG_ITEMPLATE_SAV, ' + $
              'AVG_TIME_IDX, AVG_PERIOD, FILTER_IDX, FILTER_THR'
  IF ((n_elements(full_dir) EQ 0) OR (full_dir EQ '')) THEN $
     message, 'FULL_DIR is undefined or is empty string'
  IF ((n_elements(avg_dir) EQ 0) OR (avg_dir EQ '')) THEN $
     message, 'AVG_DIR is undefined or is empty string'
  IF ((n_elements(full_itemplate_sav) EQ 0) OR $
      (full_itemplate_sav EQ '')) THEN $
         message, 'FULL_ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(full_time_idx) NE 1) OR (full_time_idx LT 0)) THEN $
     message, 'FULL_TIME_IDX must be a scalar >= zero'
  IF ((n_elements(avg_itemplate_sav) EQ 0) OR $
      (avg_itemplate_sav EQ '')) THEN $
         message, 'AVG_ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(avg_time_idx) NE 1) OR (avg_time_idx LT 0)) THEN $
     message, 'AVG_TIME_IDX must be a scalar >= zero'
  IF ((n_elements(avg_period) NE 1) OR (avg_period LT 0)) THEN $
     message, 'AVG_PERIOD must be a scalar >= zero'
  n_fi=n_elements(filter_idx)
  n_ft=n_elements(filter_thr)
  IF n_fi EQ 0 THEN message, 'FILTER_IDX is undefined'
  IF n_fi NE n_ft THEN $
     message, 'FILTER_THR must be supplied along with FILTER_IDX'
  IF n_ft GT 0 THEN BEGIN
     fi_size=size(filter_idx)
     IF (fi_size[0] GT 1) && (fi_size[2] NE 2) THEN $
        message, 'FILTER_IDX must be an integer scalar or vector'
     ft_size=size(filter_thr)
     IF (ft_size[0] GT 1) && (ft_size[2] NE 2) THEN $
        message, 'FILTER_THR must be an integer scalar or vector'
     IF ft_size[1] NE fi_size[1] THEN $
        message, 'FILTER_THR and FILTER_IDX must have the same ' + $
                 'number of elements'
  ENDIF
  idir_files=file_search(idir + path_sep() + '*', count=nidir_files, $
                         /nosort, /fold_case, /test_regular)
  IF nidir_files LT 1 THEN message, 'No input files found'

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; filter_MET.pro ends here
