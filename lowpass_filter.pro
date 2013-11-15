;; $Id: $
;; Author: Brent Else, Sebastian Luque
;; Created: 2013-11-07T22:51:27+0000
;; Last-Updated: 2013-11-15T16:50:33+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;; LOWPASS_FILTER
;; 
;; PURPOSE:
;; 
;;     Given a time series of data, a low pass filter operation is
;;     performed that cuts off frequencies above the user defined cut off
;;     frequency.  Code largely follows the low pass filter operation
;;     performed in Will Drennan's dif1byf.m program.
;; 
;; CALLING SEQUENCE:
;; 
;;     result=high_pass_filter(Signal, Fs, Hf)
;; 
;; INPUTS:
;; 
;;     Signal:    Array with input time series signal.
;;     Fs:        Scalar with sampling interval of the time series in Hz.
;;     Hf:        Scalar with frequency above which frequencies will be cutoff.
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

FUNCTION LOWPASS_FILTER, SIGNAL, FS, HF

  dt=double(1) / fs
  signal=transpose(signal)

  ;; Frequency range estimation.
  lsi=float(n_elements(signal))
  IF (lsi MOD 2) EQ 1 THEN signal=[signal, signal(lsi - 1)]
  corrsi=signal - detrend(signal)
  clsi=ceil(lsi / double(2), /L64)
  f=findgen(1, clsi * 2)
  iter=n_elements(signal)

  ;; Estimating frequency range to perform FFT operation
  FOR i=0, (iter - 1) DO f[i]=f(i) - clsi
  f=transpose(f) / lsi / dt

  ;; Signal transformation to the frequency domain.
  ;; Do I wish to detrend the signal?? I think so...
  SI=FFT(detrend(signal), /INVERSE)
  SI=complex(real_part(SI), -imaginary(SI))
  ;; Performs a FFTShift so that the zero frequency component moves to the
  ;; middle of the spectrum
  SIshift=shift(SI, clsi)

  SI[0:clsi-1]=SIshift[0:clsi - 1]
  SI[clsi]=0 ; Sets first element to zero to perform double integrations
  SI[clsi + 1:lsi - 1]=SIshift[(clsi + 1):lsi - 1]

  ;; Perform the filter
  FOR i=0, (iter - 1) DO $
     IF (f[i] LT -hf OR f[i] GT hf) THEN SI[i]=0

  ;; Return to the signal domain
  SI=shift(SI, clsi)
  SI=complex(real_part(SI), -imaginary(SI))
  si=FFT(SI)
  sireal=real_part(si) + corrsi

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
;;; lowpass_filter.pro ends here
