;; $Id: $
;;; bearing_avg.pro --- Average bearing and magnitude
;; Author: Bruce Johnson, Sebastian Luque
;; Created: 2013-09-18T20:11:56+0000
;; Last-Updated: 2013-09-18T22:22:21+0000
;;           By: Sebastian Luque
;; ------------------------------------------------------------------------
;;; Commentary: 
;; 
;; This needs some work and checking of arguments.
;;
;; Original comment by BJ (likely) below.
;;
;; DESCRIPTION:
;; This program calculates an average bearing from an input array of
;; bearings (e.g. heading or COG) It operates in two modes: one in which
;; only an average compass heading is calculated (e.g. for heading) and a
;; second in which a compass value AND a magnitude is calculated (e.g. COG
;; and SOG) It operates by decomposing the input array (iarr) into vectors.
;; When only calculating a compass heading, the magnitude of the vectors is
;; assigned the value of unity.
;; 
;; USAGE
;; the argument "iarr" is the input array of bearings to be averaged
;; (i.e. iarr=heading) For calculating bearing average only, set the
;; argument "vecarr" to 1 When calculating vector magintude (e.g. SOG/COG),
;; vecarr is the magnitude ;(e.g. iarr=COG, vecarr=SOG)
;;  
;; OUTPUT
;; 1 x 2 array.  Row 1 = avg heading, Row 2 = vector magnitude (ignore if
;; not using vector mag.)
;; 
;; NOTE:
;; 
;; inspiration was from this website: http://www.ndbc.noaa.gov/wndav.shtml
;; has been checked against all the examples proposed in this website.
;; ------------------------------------------------------------------------
;;; Code:

FUNCTION BEARING_AVG, X, VMAG, dimension=dimension

  ;; Convert to radians
  radarr=x * !DTOR
  ;; calculate x and y coordinates
  yarr=cos(radarr) * vmag
  xarr=sin(radarr) * vmag
  x_mean=mean(xarr, dimension=dimension, /NAN)    ; x_mean is in the v_wind direction
  y_mean=mean(yarr, dimension=dimension, /NAN)    ; y_mean is in the u_wind direction
  mag=sqrt(x_mean ^ 2 + y_mean ^ 2)
  vec_mean=atan(x_mean / y_mean)
  vec_mean1=atan(x_mean, y_mean)
  deg_mean=vec_mean / !DTOR
  deg_mean1=vec_mean1 / !DTOR
  IF ((x_mean GE 0) AND (y_mean GE 0)) THEN $
     compass_mean=0 + deg_mean
  IF ((x_mean GE 0) AND (y_mean LE 0)) THEN $
     compass_mean=180 + deg_mean
  IF ((x_mean LE 0) AND (y_mean LE 0)) THEN $
     compass_mean=180 + deg_mean
  IF ((x_mean LE 0) AND (y_mean GE 0)) THEN $
     compass_mean=360 + deg_mean
  RETURN, compass_mean
  ;; RETURN, compass_mean
  ;; output=fltarr(n_elements(,2)
  ;; output[0,0]=compass_mean
  ;; output[0,1]=mag
  ;; RETURN, output

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; bearing_avg.pro ends here
