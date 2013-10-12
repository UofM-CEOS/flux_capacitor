;; Author: Brent Else, Sebastian Luque
;; Created: 2013-10-11T22:38:01+0000
;; Last-Updated: 2013-10-11T22:49:07+0000
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

;; Note that the code that accompanies the smith paper is at:
;; http://www.coaps.fsu.edu/WOCE/truewind/
;;
;; Input variables:
;;   zref: zero reference of anemometer (ususally 0)
;;   cog: course over ground
;;   sog: speed over ground
;;   head: vessel heading
;;   wdir: platform wind direction (meteorological conventions used)
;;   wspd: platform/apparent wind speed
;;
;;   EXAMPLE:
;;     truewind,0,45,5.0,30,250,10

FUNCTION TRUEWIND, ZREF, COG, SOG, HEAD, WDIR, WSPD

  ;; Apparent Wind Direction
  apo=270 - (head + zref + wdir)
  ;; COG in math coordinates
  cpo=90 - cog
  ;; Tu and Tv vectors
  x=wspd * cos(Apo * !PI / 180) + sog * cos(Cpo * !PI / 180)
  y=wspd * sin(Apo * !PI / 180) + sog * sin(Cpo * !PI / 180)
  ;; T magnitude
  T=(x ^ 2 + y ^ 2) ^ 0.5
  calm_flag=1
  ;; The To... Note that the atan thing is backwards in IDL; THIS I STOLE
  ;; STRAIGHT FROM SMITH!
  IF (ABS(x) GT 0.00001) THEN BEGIN
     mtruedir=ATAN(y, x) / !DTOR
  ENDIF ELSE BEGIN
     IF(ABS(y) GT 0.00001) THEN BEGIN
        mtruedir=180.0 - 90.0 * y / ABS(y)
     ENDIF ELSE BEGIN
        mtruedir=270.0
        calm_flag=0.0
     ENDELSE
  ENDELSE
  ;; Convert to meteorological winds
  to=270.0 - mtruedir
  ;; Again, this I stole from Smith... Not to sure what it does, but it
  ;; seems to be kinda important
  WHILE (to LT 0.0) DO to=(to + 360.0) * calm_flag
  WHILE (To GE 360.0) DO to=(to - 360.0) * calm_flag
  ;; Make a meteorological definition of calm winds having direction=0, and
  ;; North winds direction=360
  IF to EQ 0 THEN to=360
  IF T LT 0.01 THEN BEGIN
     to=0
     T=0
  ENDIF
  ;; Put the true wind vector in an array to return out of the program
  true=fltarr(1,2)
  true[0, 0]=T
  true[0, 1]=to
  RETURN, true

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; true_wind.pro ends here
