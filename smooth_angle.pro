;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-11-14T00:15:07+0000
;; Last-Updated: 2013-11-26T22:09:39+0000
;;           By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     SMOOTH_ANGLE
;; 
;; PURPOSE:
;; 
;;     Smooth angles by decomposing them, applying a boxcar average with a
;;     given window, and recomposing back into angles and vectors.
;; 
;; CALLING SEQUENCE:
;; 
;;     vector=Smooth_Angle(Iangle, VMag, Width)
;; 
;; INPUTS:
;; 
;;     Iangle:     Array of angles (degrees).
;;     VMag:       Array of corresponding magnitude.
;;     Width:      Window for the boxcar average smoothing.
;; 
;; OUTPUTS:
;; 
;;     A 2-column array with the smoothed angle and magnitude,
;;     respectively.
;; 
;; EXAMPLE:
;; 
;;     
;; 
;;- -----------------------------------------------------------------------
;;; Code:

FUNCTION SMOOTH_ANGLE, IANGLE, VMAG, WIDTH

  xy=decompose(iangle, vmag)
  x=smooth(xy[0, *], width, /nan, /edge_truncate)
  y=smooth(xy[1, *], width, /nan, /edge_truncate)

  RETURN, recompose(x, y)

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; smooth_angle.pro ends here
