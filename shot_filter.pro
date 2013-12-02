;; $Id$
;; Author: Brent Else, Sebastian Luque
;; Created: 2013-10-29T18:55:13+0000
;; Last-Updated: 2013-12-02T02:21:36+0000
;;           By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     SHOT_FILTER
;; 
;; PURPOSE:
;; 
;; 
;; 
;; CALLING SEQUENCE:
;; 
;; 
;; 
;; INPUTS:
;; 
;;     X:   Vector
;; 
;; OUTPUTS:
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
    IF nok GT 0 THEN $
       x_screened=interpol(x[x_ok], x_idx[x_ok], x_idx)
    ;; Check for 'NaN'
    x_ok=where(finite(x), n_ok)
    IF n_ok GT 0 THEN $
       x_screened=interpol(x_screened[x_ok], x_idx[x_ok], x_idx)
    
  RETURN, x_screened

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; shot_filter.pro ends here
