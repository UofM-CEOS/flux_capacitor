;; $Id$
;; Author: Brent Else, Sebastian Luque
;; Created: 2013-11-07T20:41:26+0000
;; Last-Updated: 2013-11-20T00:24:53+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     LEVEL_MOTIONPAK
;; 
;; PURPOSE:
;; 
;;     Given raw Motion Pak data in a proper RH coordinate system
;;     (i.e. following Anctil et al. 1994), an iterative approach is used
;;     to determine the roll and pitch angles that should be applied to
;;     minimize mean x and y acceleration.  Yaw angle is assumed to be
;;     correct.
;; 
;; CALLING SEQUENCE:
;; 
;;     level=level_motionpak(Acc, Rate, R_Range, P_Range)
;; 
;; INPUTS:
;; 
;;      Acc:      3-column by N array with acceleration in the x, y, and z
;;                axes, in that order and in a right-hand coordinate system
;;                (as per Edson or Drennan).
;;      Rate:     3-column by N array with angular rates phi, theta, and
;;                shi, in that order and in a true right-hand coordinate
;;                system (as per Drennan).  Therefore yaw rate *-1,
;;                compared to Edson.
;;      R_Range:  Range (degrees) that you wish to attempt to roll the
;;                Motion Pak by (i.e. +/- 5 deg).
;;      P_Range:  Range (degrees) that you wish to attempt to pitch the
;;                Motion Pak by (i.e. +/- 5 deg).
;; 
;; OUTPUTS:
;; 
;;      Nx6 array with the following rows: levelled x acceleration,
;;      levelled y acceleration, levelled z acceleration, levelled phi
;;      angular rate, levelled theta angular rate, levelled shi angular
;;      rate.
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

FUNCTION LEVEL_MOTIONPAK, ACC, RATE, R_RANGE, P_RANGE, STATUS=STATUS

  acc_dims=size(acc, /dimensions)
  rate_dims=size(rate, /dimensions)

  ;; Create an array which will test the different roll/pitch angles
  roll=(dindgen(r_range * 100) / 100.0) * !DTOR
  pitch=(dindgen(p_range * 100) / 100.0) * !DTOR
  IF mean(acc[0, *]) GT 0 THEN pitch=-pitch
  IF mean(acc[1, *]) LT 0 THEN roll=-roll

 ;; Minimize roll
  min_roll=dindgen(r_range * 100)
  cr=cos(roll)
  sr=sin(roll)
  sp=sin(0)
  cp=cos(0)
  cy=cos(0)
  sy=sin(0)
  FOR r=0, r_range * 100 - 1 DO BEGIN
     ACCy=ACC[0, *] * cp * sy + ACC[1, *] * $
          (cr[r] * cy + sr[r] * sp * sy) + ACC[2, *] * $
          (-sr[r] * cy + cr[r] * sy * sp)
     min_roll[r]=mean(ACCy)
  ENDFOR
  min_mean_r=min(abs(min_roll))
  min_mean_r_index=where(abs(min_roll) EQ min_mean_r)
  min_roll_angle=roll[min_mean_r_index]
  ;; Set the status error check
  status=1                      ; failed
  IF min_mean_r_index EQ (r_range * 100 - 1) THEN BEGIN
     corrpak=make_array(6, type=4, value=!VALUES.F_NAN)
     RETURN, corrpak
  ENDIF

  ;; Minimize pitch
  min_pitch=dindgen(p_range * 100)
  cr=cos(0)
  sr=sin(0)
  sp=sin(pitch)
  cp=cos(pitch)
  cy=cos(0)
  sy=sin(0)
  FOR p=0, p_range * 100 - 1 DO BEGIN
     ACCx=ACC[0, *] * cp[p] * cy + ACC[1, *] * $
          (sp[p] * sr * cy - cr * sy) + ACC[2, *] * $
          (sy * sr + sp[p] * cr * cy)
     min_pitch[p]=mean(ACCx)
  ENDFOR
  min_mean_p=min(abs(min_pitch))
  min_mean_p_index=where(abs(min_pitch) EQ min_mean_p)
  min_pitch_angle=pitch[min_mean_p_index]
  ;; Reset status error check
  status=1                      ; failed
  IF min_mean_p_index EQ (p_range * 100 - 1) THEN BEGIN
     corrpak=make_array(6, type=4, value=!VALUES.F_NAN)
     RETURN, corrpak
  ENDIF

  ;; TEMPORARY: set min_roll_angle/min_pitch_angle manually.
  ;; min_roll_angle=1.9*!DTOR
  ;; min_pitch_angle=(-1.16)*!DTOR

  ;; Perform the rotations to get ACC/RATE into the proper ship frame
  cr=rebin([cos(min_roll_angle)], acc_dims[1])
  sr=rebin([sin(min_roll_angle)], acc_dims[1])
  sp=rebin([sin(min_pitch_angle)], acc_dims[1])
  cp=rebin([cos(min_pitch_angle)], acc_dims[1])
  cy=rebin([cos(0.0)], acc_dims[1])
  sy=rebin([sin(0.0)], acc_dims[1])

  rate_phi_ship=reform(rate[0, *] * cp * cy + rate[1, *] * $
                       (sp * sr * cy - cr * sy) + rate[2, *] * $
                       (sy*sr+sp*cr*cy))
  rate_theta_ship=reform(rate[0, *] * cp * sy + rate[1, *] * $
                         (cr * cy + sr * sp * sy) + rate[2, *] * $
                         (-sr * cy + cr * sy * sp))
  rate_shi_ship=reform(-rate[0, *] * sp + rate[1, *] * cp * sr + $
                       rate[2, *] * cp * cr)

  acc_x_ship=reform(acc[0, *] * cp * cy + acc[1, *] * $
                    (sp * sr * cy - cr * sy) + acc[2, *] * $
                    (sy * sr + sp * cr * cy))
  acc_y_ship=reform(acc[0, *] * cp * sy + acc[1, *] * $
                    (cr * cy + sr * sp * sy) + acc[2, *] * $
                    (-sr * cy + cr * sy * sp))
  acc_z_ship=reform(-acc[0, *] * sp + acc[1, *] * cp * sr + $
                    acc[2, *] * cp * cr)

  message, 'Motion Pak level - Roll: ' + $
           strcompress(min_roll_angle / !DTOR) + 'degrees, ' + $
           'Pitch: ' + strcompress(min_pitch_angle / !DTOR) + ' degrees', $
           /informational
  corr_pak=[[ACC_x_ship], [ACC_y_ship], [ACC_z_ship], $
            [RATE_phi_ship], [RATE_theta_ship], [RATE_shi_ship]]

  status=0                      ; success
  RETURN, transpose(corr_pak)

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; level_motionpak.pro ends here
