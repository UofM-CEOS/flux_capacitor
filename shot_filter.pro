;; $Id$
;; Author: Brent Else
;; Created: 2013-10-29T18:55:13+0000
;; Last-Updated: 2013-10-29T20:49:09+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;; 
;; 
;; PURPOSE:
;; 
;; 
;; 
;; CATEGORY:
;; 
;; 
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

FUNCTION SHOT_FILTER, X

    n=n_elements(x)
    x_idx=findgen(n)
    ;; Check for shot noise
    delta=stddev(x, /double, /NAN)
    x_mean=mean(x, /double, /NAN)
    x_ok=where(abs(x - x_mean) LE 3*delta, nok)
    IF nok GT 0 THEN BEGIN
      x_screened=interpol(x[x_ok], x_idx[x_ok], x_idx)
   ENDIF
    ;; Check for 'NaN'
    x_ok=where(finite(x) EQ 1, n_ok)
    IF n_ok GT 0 THEN BEGIN
      x_screened=interpol(x_screened[x_ok], x_idx[x_ok], x_idx)
   ENDIF
    
  RETURN, screened_var
END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; shot_filter.pro ends here
