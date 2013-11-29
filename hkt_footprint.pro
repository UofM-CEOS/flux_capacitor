;; $Id$
;; Author: Brent Else, Sebastian Luque
;; Created: 2013-11-19T16:31:01+0000
;; Last-Updated: 2013-11-20T18:12:07+0000
;;	     By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;;
;;     HKT_FOOTPRINT
;;
;; PURPOSE:
;;
;;     This program does footprint modeling according to Hsieh et al. 2000,
;;     "An Approximate analytical model for footprint estimation of scalar
;;     fluxes in thermally stratified atmospheric flows" Advances in Water
;;     Resources 23 (2000) 765-777.
;;
;; CALLING SEQUENCE:
;;
;;
;;
;; INPUTS:
;;
;;     Mol:	  Monin-Obukhov length.
;;     Zo:	  Roughness length.
;;     Zm:	  Height from surface of flux package.
;;     Interval:  Interval (in m) to calculate footprint contribution.
;;     Maxx:	  Maximum upwind distance for which the footprint will be
;;		  calculated (in km)... set VERY high if you don't want it
;;		  to ever crash here.
;;     MaxRat:	  Maximum CNF ratio that we will calculate up to... To
;;		  avoid ending due to this keyword, set to 100.
;;
;; KEYWORD PARAMETERS:
;;
;;     PLOT_FILE: Set this keyword equal to a .ps file where the graphs
;;		  should be printed... if not set, graphs are printed to
;;		  screen.
;;
;; OUTPUTS:
;;
;;     One graph page, containing CNF and footprint plots, printed to
;;     screen (default) or to file (if PLOT_FILE keyword is set).
;;
;;     A 2-element array: result(0)=distance to maximum contributin source
;;     area (km) result(1)=distance to point where 90% of the total flux is
;;     contributed to the measurement (km).
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
;;     result=hkt_footprint(-50.0, 0.04, 4.0, 10.0, 2.0, 0.95)
;;
;;- -----------------------------------------------------------------------
;;; Code:

FUNCTION HKT_FOOTPRINT, MOL, ZO, ZM, INTERVAL, MAXX, MAX_RATIO, $
			PLOT_FILE=PLOT_FILE

  ;; Calculate the zu term
  zu=zm *(alog(zm / zo) - 1 + zo / zm)
  ;; Figure out the stability, and set the appropriate parameters from my
  ;; understanding of the paper, unstable: zu/L < -0.04, neutral: -0.04 <
  ;; zu/L < 0.04, stable: zu/L > 0.04
  thresh=0.04
  IF zu / MOL LT -thresh THEN BEGIN
     Dval=0.28 & Pval=0.59
  ENDIF
  IF zu / MOL GE -thresh AND zu / MOL LE thresh THEN BEGIN
     Dval=0.97 & Pval=1.0
  ENDIF
  IF zu / MOL GE thresh THEN BEGIN
     Dval=2.44 & Pval=1.33
  ENDIF
  ;; Going to run this doing a while loop... if I can...  first, we'll make
  ;; two that we'll slowly add to:
  CNFcurr=0.0 & CNF=[CNFcurr]	;fraction
  XLcurr=0.0 & XL=[XLcurr]	;distance (km)
  ;; Footprint (1/m)... not exactly sure what this means, but it identifies
  ;; the point at which maximum flux contribution occurs.
  FTcurr=0.0 & FT=[FTcurr]

  ;; Now we'll make a counter, which starts at 1, and we'll initiate the
  ;; while loop
  cnt=1.0
  WHILE (CNFcurr LT max_ratio) AND (XLcurr LT maxx * 1000.0) DO BEGIN
     ;; Distance that we're using for this point
     XLcurr=(cnt*interval)
     ;; Evaluate eq'n 16
     CNFcurr=exp((-1 / (0.4 ^ 2 * XLcurr)) * $
		 Dval * zu ^ Pval * abs(MOL) ^(1.0 - Pval) )
     ;; Evaluate eq'n 17
     FTcurr=((1 / (0.4 ^ 2 * XLcurr ^ 2)) * $
	     Dval * zu ^ Pval * abs(MOL) ^ (1.0 - Pval) ) * CNFcurr

     CNF=[CNF, CNFcurr]
     XL=[XL, XLcurr / 1000.0]	; in km
     FT=[FT, FTcurr]
     cnt=cnt++
  ENDWHILE

  ;; Calculate the peak location of the footprint
  peak=((Dval * zu ^ Pval * abs(MOL) ^ (1 - Pval)) / $
	(2.0 * 0.4 ^ 2.0)) / 1000.0
  ;; Calculate the location of the 90% flux contribution
  sub90=abs(CNF - 90.0) & min90=min(sub90)
  dist90=XL(where(sub90 EQ min90))

  ;; Plotting
  !P.MULTI=[0, 1, 2, 0] 	   ;set up a page with 2 graphs
  IF keyword_set(plot_file) THEN BEGIN
     set_plot, 'ps'
     device, file=plot_file, /inches, xsize=11.0, ysize=8.0, /times, $
	     /landscape
  ENDIF ELSE BEGIN
     set_plot, 'X'
     window, 0
  ENDELSE
  ;; Plot the cumulative flux
  CNFtitle='Cumulative Flux: zo=' + strmid(strcompress(zo), 0, 7) + $
	   ', 90% flux at' + strmid(strcompress(dist90), 0, 5) + ' km'
  plot, XL, CNF, TITLE=CNFtitle, XTITLE='Upwind Distance (km)', $
	YTITLE='Cumulative Flux Fraction'
  ;; Plot the footprint... here we'll truncate the distance a bit, because
  ;; it's much shorter distance to get to the peak than the 90% val or
  ;; whatever
  endpt=(where(FT EQ max(FT), endcount)) * 10.0
  IF (endpt[0] GT n_elements(FT)) OR (endcount GT 1) THEN $
     endpt=n_elements(FT) - 1.0
  FTtitle='Footprint: Peak Located at' + strcompress(peak) + 'km'
  plot, XL[0:endpt], FT[0:endpt], TITLE=FTtitle, $
	XTitle='Upwind Distance (km)', YTITLE='Footprint (1/m)'

  RETURN, [peak, dist90]

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; hkt_footprint.pro ends here
