;; $Id$
;; Author: Brent Else, Sebastian Luque
;; Created: 2013-11-14T20:33:42+0000
;; Last-Updated: 2013-11-20T15:43:29+0000
;;	     By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;;
;;
;;
;; PURPOSE:
;;
;;     A simple linear detrend function.
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

FUNCTION DETREND, Y

 ;; Create an x-axis
  x=findgen(n_elements(y))
  y_fit=linfit(x, y, /DOUBLE)

  RETURN, y - (y_fit[0] + y_fit[1] * x)

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; detrend.pro ends here
