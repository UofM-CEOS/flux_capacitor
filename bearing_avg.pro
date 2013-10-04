;; $Id$
;; Author: Bruce Johnson, Sebastian Luque
;; Created: 2013-09-18T20:11:56+0000
;; Last-Updated: 2013-10-04T15:45:18+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;; 
;; 
;; PURPOSE:
;; 
;;  Original comment by BJ (likely) below.
;;
;;  This program calculates an average bearing from an input array of
;;  bearings (e.g. heading or COG) It operates in two modes: one in which
;;  only an average compass heading is calculated (e.g. for heading) and a
;;  second in which a compass value AND a magnitude is calculated (e.g. COG
;;  and SOG) It operates by decomposing the input array (iarr) into
;;  vectors.  When only calculating a compass heading, the magnitude of the
;;  vectors is assigned the value of unity.
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

FUNCTION BEARING_AVG, X, VMAG, DIMENSION=dimension

  ;; Convert to radians
  rads=x * !DTOR
  ;; calculate x and y coordinates
  ys=cos(rads) * vmag
  xs=sin(rads) * vmag
  x_mean=mean(xs, dimension=dimension, /NAN)    ; x_mean is in the v_wind direction
  y_mean=mean(ys, dimension=dimension, /NAN)    ; y_mean is in the u_wind direction
  mag=sqrt(x_mean ^ 2 + y_mean ^ 2)
  vec_mean=atan(x_mean, y_mean)
  deg_mean=vec_mean / !DTOR
  negs=where(deg_mean LE 0, nneg)
  IF nneg GT 0 THEN $
     deg_mean[negs]=temporary(deg_mean[negs]) + 360
  RETURN, reform([deg_mean, mag], n_elements(deg_mean), 2)

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; bearing_avg.pro ends here
