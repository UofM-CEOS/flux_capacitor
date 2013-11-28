;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-11-14T22:43:46+0000
;; Last-Updated: 2013-11-28T15:59:46+0000
;;           By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     TRUEWIND_SONIC
;; 
;; PURPOSE:
;;
;;     Original comments: The coordinate definitons for a sonic anemometer
;;     are a bit strange compared to what you would normally consider for
;;     converting it into a compass direction... So, we do a bit of
;;     modification here to get it right.
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

FUNCTION TRUEWIND_SONIC, U, V, COG, SOG, HEADING, ZREF

  ;; Calculate the wdir/vel as a compass direction measured by the sonic
  ;; anemometer.  Code was all borrowed from truewind.pro
  T=reform((v ^ 2 + u ^ 2) ^ 0.5)
  deg=reform(atan(v / u) / !DTOR)
  q1=where((u GE 0) AND (v LE 0), nq1)
  q2=where((u LE 0) AND (v LE 0), nq2)
  q3=where((u LE 0) AND (v GE 0), nq3)
  q4=where((u GE 0) AND (v GE 0), nq4)
  IF nq1 GT 0 THEN deg[q1]=180.0 - deg[q1]
  IF nq2 GT 0 THEN deg[q2]=360.0 - deg[q2]
  IF nq3 GT 0 THEN deg[q3]=0.0 - deg[q3]
  IF nq4 GT 0 THEN deg[q4]=180.0 - deg[q4]
  ;; Call TRUEWIND to calculate truewind. Note: I realize this is a bit
  ;; redundant... We could just modify the true wind program... But this
  ;; way I'll be sure it's right
  utrue=truewind(zref, cog, sog, heading, deg, T)
  ;; calculate proper x/y (u/v) variables
  wspd=utrue[*, 1]
  Apo=270 - utrue[*, 0]
  wind_true=[[wspd * sin(Apo * !DTOR)], [wspd * cos(Apo * !DTOR)], $
             [utrue[*, 1]], [utrue[*, 0]], [T], [deg]]

  RETURN, transpose(wind_true)

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; truewind_sonic.pro ends here
