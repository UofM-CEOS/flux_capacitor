;; $Id$
;; Author: Will Drennan, Brent Else, Sebastian Luque
;; Created: 2013-11-14T21:10:00+0000
;; Last-Updated: 2013-11-14T22:32:24+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     MOTCORR
;; 
;; PURPOSE:
;; 
;;     Will Drennan's Matlab code, coverted into IDL by B. Else
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
;;     Wind:   3-column array of u,v,w wind speeds in a right-hand coordinate
;;             system.  This is the same for Edson or Drennan.
;;     Acc:    3-column array of x,y,z accelerations (m/s2) in a right-hand
;;             coordinate system (as per Edson or Drennan)
;;     Angle:  3-column array of pitch,roll,yaw angles in radians in a true
;;             true right-hand coordinate system (as per Drennan, therefore
;;             yaw *-1 compared to Edson).
;;     L:      3-element array with x,y,z displacement between Motion Pak
;;             and anemometer.
;;     Fs:     Scalar indicating sampling frequency
;;     Lf:     Low pass cutoff period to be used when we integrate (for
;;             10Hz data use 0.04).
;;     Hf:     High frequency cutoff period to be used when we
;;             differentiate (for 10Hz data use hf=4.5).
;;     G:      Scalar: gravity.
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

FUNCTION MOTCORR, WIND, ACC, ANGLE, L, FS, LF, HF, G

  ;; BE: filter the computed ANGLES before going any further.  This removes
  ;; some of the hf noise.
  ANGLE[0, *]=lowpass_filter(ANGLE[0, *], fs, hf)
  ANGLE[1, *]=lowpass_filter(ANGLE[1, *], fs, hf)
  ANGLE[2, *]=lowpass_filter(ANGLE[2, *], fs, hf)

  ;; Calculate global variables
  cr=cos(ANGLE[1, *])           ; cosine of roll (phi)
  sr=sin(ANGLE[1, *])           ; sine of roll (phi)
  sp=sin(ANGLE[0, *])           ; sine of pitch (theta)
  cp=cos(ANGLE[0, *])           ; cosine of pitch (theta)

  ;; Mean yaw from the angles (shi) and subtract the mean value from the
  ;; angle derived by filtering.  THIS IS PROBABLY NO LONGER NECESSARY
  myaw=atan(mean(sin(ANGLE[2, *])), mean(cos(ANGLE[2, *])))
  cy=cos(ANGLE[2, *] - myaw)
  sy=sin(ANGLE[2, *] - myaw)

  ;; BE: in my version of the motion correction routine, I got the
  ;; accelerations in to the earth. Frame first... I can't remember why I
  ;; thought that was important....  But, I was really convinced that doing
  ;; it in this order was a mistake.

  ;; Remove the gravitational component from the linear accelerations, then
  ;; perform the integration in frequency space while high pass filtering.
  surge_b=int1byf((ACC[0, *] + g * sp), $
                  fs, lf)       ; boat reference frame (hence surge_b)
  sway_b=int1byf((ACC[1, *] - g * sr * cp), $
                 fs, lf)        ; boat reference frame (hence sway_b)
  heave_b=int1byf((ACC[2, *] - g *cp * cr), $
                  fs, lf)       ; boat reference frame (hence heave_b)

 ;; BE: low pass filter these to get rid of some of the high frequency
 ;; noise.
  surge_b=lowpass_filter(surge_b, fs, hf)
  sway_b=lowpass_filter(sway_b, fs, hf)
  heave_b=lowpass_filter(heave_b, fs, hf)

  ;; Multiply the translational velocities by the transformation matrix (T)
  ;; to derive velocities in the earth reference frame
  surge_e=surge_b * cp * cy + sway_b * (sp * sr * cy - cr * sy) + $
          heave_b * (sy * sr + sp * cr * cy) ; earth reference frame
  sway_e=surge_b * cp * sy + sway_b * (cr * cy + sr * sp * sy) + $
         heave_b * (-sr * cy + cr * sy * sp) ; earth reference frame
  heave_e=-surge_b * sp + sway_b * cp * sr + $
          heave_b * cp * cr     ;  earth reference frame

 ;; Multiply the wind velocities by the transformation matrix (T) to derive
 ;; velocities in the earth reference frame
  u_e=WIND[0, *] * cp * cy + WIND[1, *] * (sp * sr * cy - cr * sy) + $
      WIND[2, *] * (sy * sr + sp * cr * cy)
  v_e=WIND[0, *] * cp * sy + WIND[1, *] * (cr * cy + sr * sp * sy) + $
      WIND[2, *] * (-sr * cy + cr * sy * sp)
  w_e=-WIND[0, *] * sp + WIND[1, *] * cp * sr + WIND[2, *] * cp *cr

  ;; Calculate rotational motion... First differentiate the angles to
  ;; derive rates
  dpdt=diff1byf(ANGLE[0, *], fs, hf)
  drdt=diff1byf(ANGLE[1, *], fs, hf)
  yaw=ANGLE[2, *] - myaw
  replaceyaw1=where((yaw LT - !PI), yawcount1)
  IF yawcount1 GT 0 THEN yaw[replaceyaw1]=yaw[replaceyaw1] + (2 * !PI)
  replaceyaw2=where((yaw GT !PI), yawcount2)
  IF yawcount2 GT 0 THEN yaw[replaceyaw2]=yaw[replaceyaw2] - (2 * !PI)
  dydt=diff1byf(yaw, fs, hf)

  Omx=- dpdt * sy + drdt * cp * cy
  Omy=dpdt * cy + drdt * cp * sy
  Omz=dydt - drdt * sp

  ;; Multiply the positional vector by T
  Lex=L[0] * cp * cy + L[1] * (sp * sr * cy - cr * sy) + $
      L[2] * (sy * sr + sp * cr * cy)
  Ley=L[0] * cp * sy + L[1] * (cr * cy + sr * sp * sy) + $
      L[2] * (-sr * cy + cr * sy * sp)
  Lez=-L[0] * sp + L[1] * cp * sr + L[2] * cp * cr

  ;; Multiply the angular rates by TL
  x_rot=(Omy * Lez - Omz * Ley)
  y_rot=-(Omx * Lez - Omz * Lex)
  z_rot=(Omx * Ley - Omy * Lex)

  ;; Add up the motion correction terms to derive the true wind velocities
  u=reform(u_e + surge_e + x_rot)
  v=reform(v_e + sway_e + y_rot)
  w=reform(w_e + heave_e + z_rot)

  RETURN, transpose([[u], [v], [w]])

END



;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; motcorr.pro ends here
