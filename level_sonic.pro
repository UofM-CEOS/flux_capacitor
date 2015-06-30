;; Author: Brent Else, Sebastian Luque
;; Created: 2013-11-07T19:11:00+0000
;; Last-Updated: 2015-06-30T19:51:26+0000
;;           By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     LEVEL_SONIC
;; 
;; PURPOSE:
;; 
;;     This function digitally levels the sonic anemometer wind speed
;;     components, given mean roll pitch angles.  This could be done using
;;     the mean roll/pitch output from the "level_motionpak" function, but
;;     because the anemometer is on a different tower section than the
;;     Motion Pak, an independent correction is required.
;; 
;; CALLING SEQUENCE:
;; 
;;     level=level_sonic(Wind, Roll, Pitch)
;; 
;; INPUTS:
;; 
;;     Wind:   3-column array with raw u, v, and w wind speeds, in that
;;             order.
;;     Pitch:  Mean pitch angle of the sonic anemometer platform, in
;;             radians and in a LHS.
;;     Roll:   Mean roll angle of the sonic anemometer platform, in radians
;;             and in a LHS .
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

FUNCTION LEVEL_SONIC, WIND, ROLL, PITCH

  wind_dims=size(wind, /dimensions)

  ;;  Decompose angles
  cr=rebin([cos(roll)], wind_dims[1]) ; cosine of roll (phi)
  sr=rebin([sin(roll)], wind_dims[1]) ; sine of roll (phi)
  sp=rebin([sin(pitch)], wind_dims[1]) ; sine of pitch (theta)
  cp=rebin([cos(pitch)], wind_dims[1]) ; cosine of pitch (theta)
  cy=1.0                               ; yaw angle is 0, cos(0)=1.
  sy=0.0                               ; yaw angle is 0, sin(0)=0.
  ;; Now do the simple rotation
  u_corr=wind[0, *] * cp * cy + wind[1, *] * (sp * sr * cy - cr * sy) + $
         wind[2, *] * (sy * sr + sp * cr * cy)                                                 ;
  v_corr=wind[0, *] * cp * sy + wind[1, *] * (cr * cy + sr * sp * sy) + $
         wind[2, *] * (-sr * cy + cr * sy * sp)
  w_corr=-wind[0, *] * sp + wind[1, *] * cp * sr + wind[2, *] * cp * cr
  ;; Reconstruct array
  u_corr=reform(u_corr)
  v_corr=reform(v_corr)
  w_corr=reform(w_corr)

  RETURN, transpose([[u_corr], [v_corr], [w_corr]])

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
