;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-11-14T00:15:07+0000
;; Last-Updated: 2013-11-26T21:54:49+0000
;;           By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     SMOOTH_ANGLE
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
