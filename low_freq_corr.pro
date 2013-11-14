;; $Id: $
;; Author: Sebastian Luque
;; Created: 2013-11-14T22:43:46+0000
;; Last-Updated: 2013-11-14T23:43:24+0000
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

FUNCTION LOW_FREQ_CORR, U, V, COG, SOG, HEADING, ZREF

  lines=n_elements(heading)
  corr_wind=fltarr(6,lines)
  ;; [SPL: I think the looping here is just to allow the original TRUEWIND
  ;; routine to process a single record at a time.  This may no longer be
  ;; needed, and we can just process in proper array form.]
  FOR q=0L, lines - 1 DO BEGIN
     ;; The coordinate definitons for a sonic anemometer are a bit strange
     ;; compared to what you would normally consider for converting it into
     ;; a compass direction... So, we do a bit of modification here to get
     ;; it right
     x=v[q]
     y=u[q]
     ;; CALCULATE THE WDIR/VEL AS A COMPASS DIRECTION MEASURED BY THE SONIC
     ;; Code was all borrowed from truewind.pro
     T=(x ^ 2 + y ^ 2) ^ (0.5)          ;vector magnitude
     ;; [SPL: Check all this... I think we should just use the atan2 form
     ;; of arguments here.]
     vec=atan(x / y)
     deg=vec / !DTOR
     IF ((y GE 0) AND (x LE 0)) THEN BEGIN
        compass=180.0 - deg
     ENDIF
     IF ((y LE 0) AND (x LE 0)) THEN BEGIN
        compass=360.0 - deg
     ENDIF
     IF ((y LE 0) AND (x GE 0)) THEN BEGIN
        compass=0.0 - deg
     ENDIF
     IF ((y GE 0) AND (x GE 0)) THEN BEGIN
        compass=180.0 - deg
     ENDIF
     ;; Call TRUEWIND to calculate true wind. Note: I realize this is a bit
     ;; redundant... We could just modify the true wind program... But this
     ;; way I'll be sure it's right
     Utrue=truewind(zref, cog(q), sog(q), heading(q), compass, T)
     ;; calculate proper x/y (u/v) variables
     wspd=Utrue[1]
     Apo=270 - Utrue[0]
     x=wspd * cos(Apo * !PI /180)
     y=wspd * sin(Apo * !PI / 180)
     corr_wind[0, q]=y
     corr_wind[1, q]=x
     corr_wind[2, q]=Utrue[1]
     corr_wind[3, q]=Utrue[0]
     corr_wind[4, q]=T
     corr_wind[5, q]=compass
  ENDFOR

  RETURN, corr_wind

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; low_freq_corr.pro ends here
