;; $Id$
;; Author: Brent Else, Sebastian Luque
;; Created: 2013-11-15T19:12:08+0000
;; Last-Updated: 2013-12-02T02:48:15+0000
;;	     By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;;
;;     YAWPITCH
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

FUNCTION YAWPITCH, U1, V1, W1, N

  ;; first rotation for yaw correction (cross wind to zero)
  U=mean(u1, /nan)
  V=mean(v1, /nan)
  W=mean(w1, /nan)

  ce=U / sqrt(U ^ 2 + V ^ 2)
  se=V / sqrt(U ^ 2 + V ^ 2)

  yaw=[[ce, se, 0], [-1 * se, ce, 0],[0, 0, 1]]

  dum=fltarr(3,n)
  dum[0, 0:n - 1]=u1[0:n - 1]
  dum[1, 0:n - 1]=v1[0:n - 1]
  dum[2, 0:n - 1]=w1[0:n - 1]
  dum1=transpose(dum)

  xyz=yaw ## dum1
  txyz=transpose(xyz)

  u2=findgen(n)
  w2=findgen(n)

  u2[0:n - 1]=txyz[0, 0:n - 1]
  w2[0:n - 1]=txyz[2, 0:n - 1]

  U=mean(u2, /nan)
  W=mean(w2, /nan)

  ;; Second rotation for pitch corection (average w to zero)
  ct=U / sqrt(U ^ 2 + W ^ 2)
  st=W / sqrt(U ^ 2 + W ^ 2)

  pitch=[[ct, 0, st], [0, 1, 0], [-1 * st, 0, ct]]

  xyz1=transpose(fltarr(3, n))
  xyz1=pitch ## xyz
  txyz1=transpose(xyz1)

  ;; Third rotation for roll correction (cov_vw to zero)
  V=mean(txyz1[1, 0:n - 1], /double, /nan)
  W=mean(txyz1[2, 0:n - 1], /double, /nan)
  vprime=findgen(n)
  wprime=findgen(n)
  vprime=(txyz[1, 0:n - 1] - V)
  wprime=(txyz[2, 0:n - 1] - W)
  cov_vw=correlate(vprime[where(finite(vprime))], $
		   wprime[where(finite(wprime))], /covariance, /double)
  var_v=variance(vprime, /double, /nan)
  var_w=variance(wprime, /double, /nan)
  bb=(2 * cov_vw) / (var_v - var_w)
  bbeta=0.5 * atan(bb)
  roll=[[1, 0, 0], $
	[0, cos(bbeta), sin(bbeta)], $
	[0, -1 * sin(bbeta), cos(bbeta)]]
  xyz2=roll ## xyz1

  RETURN, transpose(xyz2)

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; yawpitch.pro ends here
