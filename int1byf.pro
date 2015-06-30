;; Author: Will Drennan, Brent Else, Sebastian Luque
;; Created: 2013-11-07T22:22:17+0000
;; Last-Updated: 2015-06-30T19:51:11+0000
;;	     By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;;
;;     INT1BYF
;;
;; PURPOSE:
;;
;;     This function performs a high pass filtering to the acceleration
;;     data.  It double-integrates a signal in the time domain by an
;;     equivalent product in the frequency program, and cut off energy of
;;     periods greater than cut-off.
;;
;;     Source provided by Will Drennan.  Matlab code from Mohammad Mahfuz
;;     Al-Mamun, and edited by Brent Else.
;;
;; CALLING SEQUENCE:
;;
;;     res=int1byf(Si, Fs, Lf)
;;
;; INPUTS:
;;
;;     Si:     Input signal.
;;     Fs:     Sampling frequency [Hz].
;;     Lf:     Low frequency cutoff for high pass filter operation [Hz].
;;
;; OUTPUTS:
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

FUNCTION INT1BYF, si, fs,lf

  dt=1.0 / fs	      ; convert sampling rate to s (as per Drennan's code)
  si=transpose(si)

  ;; Frequency range estimation
  lsi=float(N_ELEMENTS(si))	; size of column vector
  IF (lsi MOD 2) EQ 1 THEN si=[si, si(lsi - 1)]

  clsi=ceil(lsi / 2.0, /L64)
  f=FINDGEN(1, clsi * 2)
  iter=N_ELEMENTS(f)
  ;; Estimating frequency range to perform FFT operation
  FOR i=0, (iter - 1) DO f[i]=f(i)-clsi
  f=transpose(f) / lsi / dt	; estimate the frequency range
  j=COMPLEX(0, 1)		; define a complex variable J=0 + j

  ;; % Signal transformation to the frequency domain
  SI=FFT(detrend(si), /INVERSE)
  SI=complex(real_part(SI), -imaginary(SI))
  ;; Performs a FFTShift so that the zero frequency component moves to the
  ;; middle of the spectrum
  SIshift=shift(SI, clsi)

  ;; % Double integration: From the FT properties, integration in time
  ;; domain is equivalent to devision in frequency domain
  SI[0:clsi - 1]=SIshift(0:clsi - 1) / (j * 2 * !pi * f(0:clsi - 1))
  ;; Sets the first element to zero to perform the double integrations
  SI[(clsi)]=0
  SI[clsi + 1:lsi - 1]=SIshift((clsi + 1):lsi - 1) / $
		       (j * 2 * !pi * f(clsi + 1:lsi - 1))

  ;; % Energy elimination at frequencies smaller than 1/cut-off
  iter=N_ELEMENTS(SI)
  FOR i=0, (iter - 1) DO BEGIN
     IF ( f(i) GT -lf AND f(i) LT lf ) THEN SI[i]=0;
  ENDFOR

  ;; Signal transformation back to the temporal domain
  si=SHIFT(si, clsi)
  si=complex(real_part(si), -imaginary(si))
  si=FFT(si)
  sireal=real_part(si)

  IF (lsi MOD 2) EQ 1 THEN BEGIN
     temp=findgen(lsi)
     temp(0:lsi - 1)=sireal(0:lsi - 1)
     sireal=temp
  ENDIF

  RETURN, sireal

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
