;; $Id$
;; Author: Brent Else, Sebastian P. Luque
;; Created: 2013-11-29T14:57:17+0000
;; Last-Updated: 2013-11-29T16:03:26+0000
;;           By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;; 
;; 
;; PURPOSE:
;; 
;;     This program is designed to make ogive plots of input flux data,
;;     with data points at intervals of +10 second integration, starting
;;     from an averaging time of 1 minute.
;; 
;; CALLING SEQUENCE:
;; 
;; 
;; 
;; INPUTS:
;; 
;;     X:              Vector of the scalar of interest.
;;     XName:          String which can have a value of: 'wco2_op',
;;                     'wh2o_op', 'wco2_cl', 'wh2o_cl', 'wu', or 'wTair'.
;;     Wrot:           Vector of rotated vertical wind velocity.
;;     Period:         Scalar float indicating the averaging period of
;;                     input data (in minutes, eg. 30 minute flux run).
;;     Isample_Rate:   Scalar float which describes the sampling frequency of
;;                     the input data (in Hz, e.g. 10 Hz).
;;     Maxc:           Scalar integer variable indicating the maximum
;;                     number of records than can be lagged for digital
;;                     time shifting OR a 2-element vector.
;;     Output:         2-element string vector containing the flux run
;;                     stamp, and the ogive output directory.
;; 
;; KEYWORD PARAMETERS:
;; 
;; 
;; 
;; OUTPUTS:
;; 
;;     A nice ogive plot!
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

PRO OGIVE, STAMP, WROT, SCALAR, PERIOD, FREQ, MAXC, OUTPUT

  ;; extract output data
  fluxrun=output[0]
  odir=output[1]

  ;; Calculate the number of iterations we will do at 10 sec intervals,
  ;; starting from 1 minute
  ini=60. / (1. / freq)             ;initial sample chunk
  iters=(period - 1) * 60. / freq
  datarr=fltarr(2, iters + 1)

  ;; Here is how we do it for ones which do not need to be digitally
  ;; shifted (mom, heat flux)
  IF stamp EQ 'wu' OR 'wTair' THEN BEGIN
     FOR x=0, iters DO BEGIN
        ;; Extract the chunk of w data we will use for this
        use_w=wrot[0:ini + (x * (10. * freq)) - 3]
        ;; Extract the chunk of scalar data we will use for this
        use_scalar=scalar[0:ini + (x * (10. * freq)) - 3]
        ;; now, we calculate the covariance for this period using the usual
        ;; approach, and save it into our final array.
        cov_w_scalar=correlate(use_w, use_scalar, /COVARIANCE, /double)
        ;; Frequency value for this sample
        datarr[0, x]=1. / ((ini + (x * (10. * freq))) / 10.)
        ;; Covariance value for this samples
        datarr[1, x]=cov_w_scalar
     ENDFOR
  ENDIF

  ;; Here is how we do it for the open path sensors
  IF stamp EQ 'wco2_op' OR 'wh2o_op' THEN BEGIN
     FOR x=0, iters DO BEGIN
        ;; Extract the chunk of w data we will use for this
        use_w=wrot[0:ini + (x * (10. * freq)) - 3]
        ;; Extract the chunk of scalar data we will use for this
        use_scalar=scalar[0:ini + (x * (10. * freq)) - 3]
        ;; Now, we calculate the covariance for this period using the usual
        ;; approach, and save it into our final array.
        n_lag_op=(maxc * 2) + 1
        lag_op=findgen(n_lag_op) - maxc
        w_scalar=c_covariance(scalar, wrot, lag_op)
        ii=where(abs(w_scalar) EQ max(abs(w_scalar)))
        cov_w_scalar=w_scalar[ii[0]]
        ;; Frequency value for this sample
        datarr[0, x]=1. / ((ini + (x * (10. * freq))) / 10.)
        ;; Covariance for this sample
        datarr[1, x]=cov_w_scalar
     ENDFOR
  ENDIF 

  ;; Here is how we do it for the closed path sensors
  IF stamp EQ 'wco2_cl' OR 'wh2o_cl' THEN BEGIN
     FOR x=0, iters DO BEGIN
        ;; Extract the chunk of w data we will use for this
        use_w=wrot[0:ini + (x * (10. * freq)) - 3]
        ;; Extract the chunk of scalar data we will use for this
        use_scalar=scalar[0:ini + (x * (10. * freq)) - 3]
        ;; Now, we calculate the covariance for this period using the usual
        ;; approach, and save it into our final array.
        n_lag_cl=(maxc_c[1] - maxc_c[0]) + 1
        lag_cl=findgen(n_lag_cl)
        FOR lagpop=(maxc_c[0]), (maxc_c[1]) DO $
           lag_cl[lagpop-maxc_c[0]]=-lagpop
        w_scalar=c_covariance(scalar, wrot, lag_cl)
        ii=where(abs(w_scalar) EQ max(abs(w_scalar)))
        cov_w_scalar=w_scalar[ii[0]]
        ;; Frequency value for this sample
        datarr[0, x]=1. / ((ini + (x * (10. * freq))) / 10.)
        ;; Covariance for this sample
        datarr[1, x]=cov_w_scalar
     ENDFOR
  ENDIF 

  ;; Now we do the plotting
  !P.MULTI=[0, 0, 0, 0]         ; set up a page for a single plot
  set_plot, 'ps'                ; set up the postscript file 
  plot_file=strcompress(output[1] + 'Ogive_' + stamp + $
                        strmid(output[0], 7, 26) + '.ps', /REMOVE_ALL)
  device, file=plot_file
  device, /inches, xsize=11.0, ysize=8.
  device, /times
  device, /landscape
  plot, datarr[0, *], datarr[1, *], xtitle='f(H!Dz!N)', $
        ytitle='Covariance', title='Ogive Plot', /XLOG
  device, /close

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; ogive.pro ends here
