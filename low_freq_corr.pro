;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-11-14T22:43:46+0000
;; Last-Updated: 2013-11-27T17:02:28+0000
;;           By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;; 
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
  corr_wind=fltarr(6, lines)
  FOR q=0L, lines - 1 DO BEGIN
     ;; The coordinate definitons for a sonic anemometer are a bit strange
     ;; compared to what you would normally consider for converting it into
     ;; a compass direction... So, we do a bit of modification here to get
     ;; it right
     x=v[q]
     y=u[q]
     ;; CALCULATE THE WDIR/VEL AS A COMPASS DIRECTION MEASURED BY THE SONIC
     ;; Code was all borrowed from truewind.pro
     T=(x ^ 2 + y ^ 2) ^ (0.5)  ; vector magnitude
     ;; [SPL: Check all this... I think we should just use the atan2 form
     ;; of arguments here.]
     deg=atan(x / y) / !DTOR
     IF ((y GE 0) AND (x LE 0)) THEN compass=180.0 - deg
     IF ((y LE 0) AND (x LE 0)) THEN compass=360.0 - deg
     IF ((y LE 0) AND (x GE 0)) THEN compass=0.0 - deg
     IF ((y GE 0) AND (x GE 0)) THEN compass=180.0 - deg
     ;; Call TRUEWIND to calculate truewind. Note: I realize this is a bit
     ;; redundant... We could just modify the true wind program... But this
     ;; way I'll be sure it's right
     Utrue=truewind(zref, cog[q], sog[q], heading[q], compass, T)
     ;; calculate proper x/y (u/v) variables
     wspd=Utrue[1]
     Apo=270 - Utrue[0]
     corr_wind[0, q]=wspd * sin(Apo * !PI / 180)
     corr_wind[1, q]=wspd * cos(Apo * !PI /180)
     corr_wind[2, q]=Utrue[1]
     corr_wind[3, q]=Utrue[0]
     corr_wind[4, q]=T
     corr_wind[5, q]=compass
  ENDFOR

  ;; ;; ;; [SPL: I think the looping here is just to allow the original TRUEWIND
  ;; ;; ;; routine to process a single record at a time.  This may no longer be
  ;; ;; ;; needed, and we can just process in proper array form as below.]
  ;; T=(v ^ 2 + u ^ 2) ^ 0.5
  ;; deg=atan(v / u) / !DTOR
  ;; q1=where((u GE 0) AND (v LE 0), nq1)
  ;; q2=where((u LE 0) AND (v LE 0), nq2)
  ;; q3=where((u LE 0) AND (v GE 0), nq3)
  ;; q4=where((u GE 0) AND (v GE 0), nq4)
  ;; IF nq1 GT 0 THEN deg[q1]=180.0 - deg[q1]
  ;; IF nq2 GT 0 THEN deg[q2]=360.0 - deg[q2]
  ;; IF nq3 GT 0 THEN deg[q3]=0.0 - deg[q3]
  ;; IF nq4 GT 0 THEN deg[q4]=180.0 - deg[q4]
  ;; utrue=truewind(zref, cog, sog, heading, deg, T)
  ;; wspd=utrue[*, 1]
  ;; Apo=270 - utrue[*, 0]
  ;; wind_true=[]
  ;; RETURN, wind_true

  RETURN, corr_wind

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; low_freq_corr.pro ends here
