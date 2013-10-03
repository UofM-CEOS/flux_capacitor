;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-09-18T20:47:27+0000
;; Last-Updated: 2013-10-03T03:08:26+0000
;;           By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;;  NAME:
;; 
;;   Original comment by BJ (probably) below.
;; 
;;   This program calculates the standard deviation of wind direction using
;;   the Yamartino method (same method used by Campbell, see
;;   http://www.wikidoc.org/index.php/Yamartino_method).
;; 
;; 
;;  PURPOSE:
;; 
;; 
;; 
;;  CATEGORY:
;; 
;; 
;; 
;;  CALLING SEQUENCE:
;; 
;; 
;; 
;;  INPUTS:
;; 
;; 
;; 
;;  OPTIONAL INPUTS:
;; 
;; 
;; 
;;  KEYWORD PARAMETERS:
;; 
;; 
;; 
;;  OUTPUTS:
;; 
;; 
;; 
;;  OPTIONAL OUTPUTS:
;; 
;; 
;; 
;;  COMMON BLOCKS:
;; 
;; 
;; 
;;  SIDE EFFECTS:
;; 
;; 
;; 
;;  RESTRICTIONS:
;; 
;; 
;; 
;;  PROCEDURE:
;; 
;; 
;; 
;;  EXAMPLE:
;; 
;; 
;; 
;;- -----------------------------------------------------------------------
;;; Code:

FUNCTION STDDEV_YAMARTINO, X

  ;;convert to radians
  rads=x * !DTOR
  ;;calculate mean sin and cosine
  goodrad=where(finite(rads) EQ 1, goodcount)
  IF goodcount LT 1 THEN BEGIN
     sd='NaN'
  ENDIF ELSE BEGIN
     mean_sin=mean(sin(rads[goodrad]))
     mean_cos=mean(cos(rads[goodrad]))
     ;;calculate epsilon
     eps=sqrt(1 - (mean_sin ^ 2 + mean_cos ^ 2))
     ;;calculate standard deviation
     sd=(asin(eps) * (1 + 0.1547 * eps ^ 3)) / !DTOR
  ENDELSE

  RETURN, sd

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; stdev_yamartino.pro ends here
