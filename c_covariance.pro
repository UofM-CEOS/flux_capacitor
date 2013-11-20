;; $Id$
;; Author: Brent Else, Sebastian Luque
;; Created: 2013-11-15T20:04:23+0000
;; Last-Updated: 2013-11-18T16:30:31+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     C_COVARIANCE
;; 
;; PURPOSE:
;; 
;;     This program calculates the cross covariance of two variables (X, Y)
;;     lagged by a value L The program is described in the file
;;     C_Correlate.pdf, and is designed to fix mistakes in the IDL program
;;     C_Correlate.
;;
;; CALLING SEQUENCE:
;; 
;; 
;; 
;; INPUTS:
;; 
;;     X: an n-element vector (n must be >2)
;;     Y: a second n-element vector
;;     L: a vector of size Ln which provides the lag values for calculating
;;        the cross covariance
;; 
;; OUTPUTS:
;; 
;;     an Ln-element vector of covariance, solving for the different lags.
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

FUNCTION C_COVARIANCE, X, Y, L

  ;; Check a couple of things
  Xn=n_elements(X) & Yn=n_elements(Y) & nrecs=Xn
  IF (Xn NE Yn) THEN $
     message, 'X and Y arrays must have the same number of elements.'
  IF (Xn LT 2) THEN $
     message, 'X and Y arrays must contain 2 or more elements.'

  ;; first, find out how many elements are contained by L, then start a
  ;; loop to go through that
  Ln=n_elements(L)
  out=findgen(Ln)               ;this is the vector that will be output
  FOR z=0L, Ln-1 DO BEGIN
     Lag=L[z]          ;derive the lag value that we will use for this loop
     ;; If L is less than 0, shift and truncate the signals appropriately
     IF Lag LT 0 THEN BEGIN
        X_shift=X[0 + abs(lag):nrecs - 1]
        Y_shift=Y[0:nrecs - 1 - abs(Lag)]
        ;; if L is greater than or equal to 0, shift and truncate signals
        ;; appropriately
     ENDIF ELSE BEGIN
        X_shift=X[0:nrecs - 1 - Lag]
        Y_shift=Y[0 + lag:nrecs - 1]
     ENDELSE
     ;; Calculate covariance 
     out[z]=mean((Y_shift - mean(Y_shift, /nan)) * $
                 (X_shift - mean(X_shift, /nan)), /nan)
  ENDFOR

  RETURN, out

END
    


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; c_covariance.pro ends here
