;; $Id: $
;; Author: Will Drennan, Brent Else, Sebastian Luque
;; Created: 2013-11-14T21:45:39+0000
;; Last-Updated: 2013-11-14T22:35:18+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     DIFF1BYF
;; 
;; PURPOSE:
;; 
;;     This source code is designed by Will Drenan of Uni of Miami and is
;;     by the Author copyright
;; 
;;     This function differentiates a signal in the time domain by an
;;     equivalent product in the frequency program, and cut off energy of
;;     periods less than cut-off.
;; 
;; CALLING SEQUENCE:
;; 
;; 
;; 
;; INPUTS:
;; 
;;       Si:  Input signal
;;       Fs:  Sampling frequency [Hz]
;;       Lf:  Low frequency cutoff for high pass filter operation
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

FUNCTION DIFF1BYF, si, fs,lf

  dt=double(1) / fs   ; convert to sampling rate in s (as per Drennan's code)
  si=transpose(si)

  ;; Frequency range estimation.
  lsi=float(n_elements(si))   ; size of the column vector
  IF (lsi MOD 2) EQ 1 THEN si=[si, si[lsi - 1]]

  clsi=ceil(lsi / double(2), /L64)
  f=findgen(1, clsi * 2)
  iter=n_elements(f)
  ;; Estimating frequency range to perform FFT operation
  FOR i=0, (iter - 1) DO f[i]=f[i] - clsi
  f=transpose(f) / lsi / dt     ; estimate the frequency range
  j=complex(0, 1)               ; Define a complex variable J=0 + j

  ;; % Signal transformation to the frequency domain.
  SI=FFT(detrend(si), /INVERSE)
  SI=complex(real_part(SI), -imaginary(SI))
  ;; Performs a FFTShift so that the zero frequency component moves to the
  ;; middle of the spectrum
  SIshift=shift(SI, clsi)

  ;; % Differentiation: From the FT properties, differation in time domain
  ;; is equivalent to multiplication in frequency domain
  SI[0:clsi - 1]=SIshift[0:clsi - 1] * (j * 2 * !pi * f[0:clsi - 1])
  SI[(clsi)]=0 ; sets first element to zero for the double integrations
  SI[clsi + 1:lsi - 1]=SIshift[(clsi + 1):lsi - 1] * $
                       (j * 2 * !pi * f(clsi + 1:lsi - 1))

  ;; % Energy elimination at frequencies smaller than 1 / cutoff, i.e. it
  ;; is a
  iter=n_elements(SI)

  FOR i=0, (iter - 1) DO $
     IF (f(i) LT -lf OR f[i] GT lf) THEN SI[i]=0

  ;; Signal transformation back to the temporal domain.
  si=shift(si, clsi)
  si=complex(real_part(si), -imaginary(si))
  si=FFT(si)
  sireal=real_part(si)

  IF (lsi MOD 2) EQ 1 THEN BEGIN
     temp=findgen(lsi)
     temp[0:lsi - 1]=sireal[0:lsi - 1]
     sireal=temp
  ENDIF

  RETURN, sireal

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; diff1byf.pro ends here
