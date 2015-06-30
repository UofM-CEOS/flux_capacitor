;; Author: Brent Else, Sebastian Luque
;; Created: 2013-10-11T22:38:01+0000
;; Last-Updated: 2015-06-30T19:53:56+0000
;;           By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     TRUEWIND
;; 
;; PURPOSE:
;; 
;;     Calculate true wind.
;;
;;     Note that the code that accompanies the Smith paper is at:
;;     http://www.coaps.fsu.edu/WOCE/truewind/
;; 
;; CALLING SEQUENCE:
;; 
;;     truew=truewind(Zref, Cog, Sog, Head, Wdir, Wspd)
;; 
;; INPUTS:
;; 
;;     Zref:   Zero reference of anemometer (ususally 0).
;;     Cog:    Vector array with course over ground.
;;     Sog:    Vector array with speed over ground.
;;     Head:   Vector array with Vessel heading.
;;     Wdir:   Vector array with platform wind direction (meteorological
;;             conventions used).
;;     Wspd:   Vector array with platform (apparent) wind speed
;; 
;; OUTPUTS:
;; 
;;     A two-row array with true wind direction in the first row and
;;     true wind speed in the second.
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
;; print, truewind(0, 0.0, 0.0, 0.0, 90.0, 5.0)
;; 
;;- -----------------------------------------------------------------------
;;; Code:

FUNCTION TRUEWIND, ZREF, COG, SOG, HEAD, WDIR, WSPD

  ;; Apparent Wind Direction
  apo=270 - (head + zref + wdir)
  ;; COG in math coordinates
  cpo=90 - cog
  ;; Tu and Tv vectors
  x=wspd * cos(apo * !DTOR) + sog * cos(cpo * !DTOR)
  y=wspd * sin(apo * !DTOR) + sog * sin(cpo * !DTOR)
  ;; T magnitude
  T=reform((x ^ 2 + y ^ 2) ^ 0.5)
  calm_flag=make_array(n_elements(T), type=2, value=1)
  ;; The To... Note that the atan thing is backwards in IDL; THIS I STOLE
  ;; STRAIGHT FROM SMITH!
  mtruedir=make_array(n_elements(T), type=4, value=!VALUES.F_NAN)
  is_xok=where(abs(x) GT 1e-5, nxok)
  is_yok=where((abs(x) LE 1e-5) AND (abs(y) GT 1e-5), nyok)
  is_calm=where((abs(x) LE 1e-5) AND (abs(y) LE 1e-5), ncalm)
  IF nxok GT 0 THEN $
     mtruedir[is_xok]=atan(y[is_xok], x[is_xok]) / !DTOR
  IF nyok GT 0 THEN $
     mtruedir[is_yok]=180.0 - 90.0 * y[is_yok] / ABS(y[is_yok])
  IF ncalm GT 0 THEN BEGIN
     mtruedir[is_calm]=270
     calm_flag[is_calm]=0.0
  ENDIF
  ;; Convert to meteorological winds
  to=270.0 - mtruedir
  ;; Again, this I stole from Smith... Not to sure what it does, but it
  ;; seems to be kinda important
  FOR i=0L, n_elements(to) - 1 DO BEGIN
     WHILE (to[i] LT 0.0) DO to[i]=(to[i] + 360.0) * calm_flag[i]
     WHILE (to[i] GE 360.0) DO to[i]=(to[i] - 360.0) * calm_flag[i]
  ENDFOR
  ;; Make a meteorological definition of calm winds having direction=0, and
  ;; North winds direction=360
  is_north=where(to EQ 0, nnorth)
  IF nnorth GT 0 THEN to[is_north]=360
  is_calm1=where(T LT 0.01, nncalm1)
  IF nncalm1 GT 0 THEN BEGIN
     to[is_calm1]=0
     T[is_calm1]=0
  ENDIF
  ;; Put the true wind vector in an array to return out of the program
  RETURN, reform([to, T], n_elements(to), 2)

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
