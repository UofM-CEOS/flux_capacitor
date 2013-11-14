;; Author: Brent Else, Sebastian Luque
;; Created: 2013-11-14T00:31:53+0000
;; Last-Updated: 2013-11-14T17:08:09+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     DECOMPOSE
;; 
;; PURPOSE:
;; 
;;     Simple function to decompose angles and associated magnitudes into x
;;     and y vectors.
;; 
;; CALLING SEQUENCE:
;; 
;;     coords=decompose(Angle, Vmag)
;; 
;; INPUTS:
;; 
;;     Angle:    Vector of angles in degrees.
;;     VMag:     Vector of associated magnitudes (units irrelevant).
;; 
;; OUTPUTS:
;; 
;;     Outputs an array with x (sine) and y (cosine) vector components.
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

FUNCTION DECOMPOSE, ANGLE, VMAG

  ;; Reform just to remove degenerate dimensions in input
  x=reform(vmag * sin(angle * !DTOR))
  y=reform(vmag * cos(angle * !DTOR))

  RETURN, transpose([[x], [y]])

END


;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     RECOMPOSE
;; 
;; PURPOSE:
;; 
;;     Simple function to recompose angles and associated magnitudes from x
;;     and y vector components.
;; 
;; CALLING SEQUENCE:
;; 
;;     angmag=recompose(X, Y)
;; 
;; INPUTS:
;; 
;;     X:    Vector of x coordinates (sine).
;;     Y:    Vector of y coordinates (cosine).
;; 
;; OUTPUTS:
;; 
;;     Outputs an array with angles and associated magnitudes.
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

FUNCTION RECOMPOSE, X, Y

  ;; Reform just to remove degenerate dimensions in input
  vmag=reform(sqrt(x ^ 2 + y ^ 2))
  ang=reform(atan(x, y) / !DTOR)
  negs=where(ang LE 0, nneg)
  IF nneg GT 0 THEN ang[negs]=ang[negs] + 360
  ;; When magnitude is 0, then angle should also be 0.  IDL is funny with
  ;; numerical representation issues and comparisons...  It says that an
  ;; angle that printed as 0.0000000000 is less than 0 in the comparison
  ;; above, so we have to fix it here (after having added 360).
  zeros=where(vmag EQ 0, nzeros)
  IF nzeros GT 0 THEN ang[zeros]=0

  RETURN, transpose([[ang], [vmag]])

END

  
;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; decompose.pro ends here
