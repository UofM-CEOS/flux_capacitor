;; $Id$
;; Author: Brent Else, Sebastian Luque
;; Created: 2013-11-18T22:20:50+0000
;; Last-Updated: 2013-12-05T19:28:31+0000
;;	     By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;;
;;     EC_CLOSED
;;
;; PURPOSE:
;;
;;     This program calculates CO2 and H2O flux from a closed path system.
;;
;; CALLING SEQUENCE:
;;
;;
;;
;; INPUTS:
;;
;;	Wind:	       3xn element array of unrotated wind [u,v,w].
;;	XCO2_m:        1xn array of mole fraction of CO2 (umol co2 / mol
;;		       moist air).
;;	XH2O_m:        1xn array of mole fraction of H2O (mmol h2o / mol
;;		       moist air).
;;	IRGA_P:        1xn array of pressure (kPa) in the IRGA cell.
;;      IRGA_T:        If /PSEUDO_WPL is set: a 1xn array of temperature
;;                     (C) in the IRGA cell. If /PSEUDO_WPL is not set: a
;;                     single floating point value giving the average IRGA
;;                     cell temperature for the period.
;;	MET_T:	       Scalar float variable of the mean air temperature
;;		       measured by the MET logger (C).
;;      MET_P:         Scalar float variable of the mean pressure measured
;;                     by the MET logger (KpA).
;;	MET_RH:        Scalar float variable of the mean relative humidity
;;		       measured by the MET logger (percentage).
;;	MAXC:	       Scalar integer variable indicating the maximum number
;;		       of records that can be lagged for digital time
;;		       shifting.
;;	Flow:	       A variable indicating the gas flow through the tube in
;;		       LPM.
;;      EC_Period:     Scalar float variable of the averaging period under
;;                     consideration (seconds) (e.g. 30min=1200 s flux
;;                     runs).
;;      Isample_Freq:  Scalar float variable of the data sampling frequency
;;                     in Hz (e.g. 10 Hz).
;;
;; KEYWORD PARAMETERS:
;;
;;     PSEUDO_WPL:     When this keyword is NOT set, point-by-point
;;		       conversions to mixing ratio (for h2o/co2) are made
;;		       from mole fraction.  This only is available if high
;;		       res'n IRGA_T & IRGA_P are available with Xco2_m and
;;		       Xh2o_m If this keyword is set, the flux of co2 will
;;		       be calculated with a pseudo WPL correction (as per
;;		       eqn3b of Ibrom et al. 2007).
;;     CORR_MASSMAN:   Set this keyword to run a spectral correction
;;		       following Massman this keyword must be used to set a
;;		       number of parameters for the massman correction.
;;		       ALL fields must be filled.  If they are not
;;		       applicable to the specific application, set them to
;;		       the string value 'NAN' (for not applicable): [0:
;;		       sample rate (sec), 1: sonic path length (m), 2:
;;		       scalar sensor separation in x coord (m), 3: scalar
;;		       sensor separation in y coord (m), 4: line length for
;;		       scalar sensor with line averaging (m), 5: diameter
;;		       of scalar sensor with volume averaging (m), 6:
;;		       length of scalar sensor with volume averaging (m),
;;		       7: tube flow (LPM), 8: tube diameter (m), 9: tube
;;		       length (m)].
;;     OGIVE_OFILES:   Set this key word to produce an ogive plot for all
;;                     fluxes which are calculated using this routine.
;;                     This keyword must be set to a 2-element string
;;                     array, indicating the directory where the ogive plot
;;                     files for CO2 and H2O (in that order) will be
;;                     placed.
;;
;; OUTPUTS:
;;
;;     A 10-element array:
;;
;;     [0]: w/Tair covariance, [1]: correction factor for heat flux,
;;     [2]: heat flux (W/m2), [3]: w/co2 covariance
;;     [4]: correction factor for co2 flux, [5]: CO2 flux (mmol/m2/day)
;;     [6]: w/h2o covariance, [7]: correction factor for h2o flux
;;     [8]: h2o flux (mol/m2/day), [9]: Latent heat flux (W/m2)
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

FUNCTION EC_CLOSED, WIND, XCO2_M, XH2O_M, IRGA_P, IRGA_T, MET_T, MET_RH, $
		    MET_P, MAXC_C, EC_PERIOD, ISAMPLE_FREQ, $
		    PSEUDO_WPL=PSEUDO_WPL, CORR_MASSMAN=CORR_MASSMAN, $
		    OGIVE_OFILES=OGIVE_OFILES

  ;; CONSTANTS:
  r=8.31451		 ; j/mol/k universal gas constant
  rv=461.5		 ; J/kg/K gas constant for water vapour Stull(1995)
  mv=18.02		 ; g/mol molecular weight for water, Stull(1995)
  ma=28.96		 ; g/mol molecular weight for dry air, Stull(1995)
  mc=44.0098		 ; g/mol molecular weight co2 Stull(1995)
  ;; J/g/K specific heat for dry air at constant pressure Stull(1995)
  cpd=1004.67 / 1000
  ;; J/g/K specific heat for water vapour at constant pressure Stull(1995)
  cpv=1875 / 1000
  mu=0.622			; ratio, Rd/Rv, Stull(1988)
  g=9.807			; m2/s
  ;; ratio of molar mass c to molar mass co2.  Multiplier to convert from
  ;; mass flux co2 to mass flux c
  conv=0.2729
  nrecs=n_elements(c_co2)

  ;; Convert moist air mixing ratios into mol/mol
  Xco2_m=Xco2_m / 1000000.0
  Xh2o_m=Xh2o_m / 1000.0

  ;; Calculate the mean horizontal wind speeds (unrotated)... this is
  ;; required in the spectral correction
  horwind=[mean(WIND[0, *], /nan), mean(WIND[1, *], /nan)]
  ;; do the wind rotations
  WINDrot=yawpitch(WIND[0, *], WIND[1, *], WIND[2, *], $
		   n_elements(WIND[0, *]))
  urot=WINDrot[0, *] & vrot=WINDrot[1, *] & wrot=WINDrot[2, *]

  IF keyword_set(pseudo_wpl) THEN BEGIN
     ;; CALCULATE FLUXES FROM MOIST AIR MIXING RATIOS (REQUIRES A
     ;; PSEUDO-WPL CORRECTION FOR H2O EFFECTS)
     ;; Calculate mean components
     mean_Xco2_m=mean(Xco2_m, /NAN, /DOUBLE)
     mean_Xh2o_m=mean(Xh2o_m, /NAN, /DOUBLE)
     fluc_Xco2_m=Xco2_m - mean_Xco2_m
     fluc_Xh2o_m=Xh2o_m - mean_Xh2o_m
     ;; Calculate covariances using max covariance
     n_lag_cl=(maxc_c[1] - maxc_c[0]) + 1
     lag_cl=findgen(n_lag_cl)

     FOR lagpop=maxc_c[0], maxc_c[1] DO $
	lag_cl[lagpop-maxc_c[0]]=-lagpop

     w_Xco2_m=c_covariance(Xco2_m, wrot, lag_cl)
     w_Xh2o_m=c_covariance(Xh2o_m, wrot, lag_cl)
     ii_co2=where(abs(w_Xco2_m) EQ max(abs(w_Xco2_m)))
     ii_h2o=where(abs(w_Xh2o_m) EQ max(abs(w_Xh2o_m)))
     lag_co2_cl=lag_cl[ii_co2(0)]
     lag_h2o_cl=lag_cl[ii_h2o(0)]
     cov_w_Xco2_m=w_Xco2_m[ii_co2[0]]
     ;; Covariance of water vapour, lagged by CO2 time constant
     cov_w_Xh2o_m_forWPL=w_Xh2o_m[ii_co2[0]]
     cov_w_Xh2o_m=w_Xco2_m[ii_h2o[0]] ;covariance of water vapour

     IF keyword_set(ogive_ofiles) THEN BEGIN
	ogive, 'wco2_cl', wrot, Xco2_m, ec_period, isample_freq, maxc_c, $
	       ogive_ofiles[0]
	ogive, 'wh2o_cl', wrot, Xh2o_m, ec_period, isample_freq, maxc_c, $
	       ogive_ofiles[1]
     ENDIF
     cf_wXco2m=!VALUES.D_NAN
     cf_wXh2om=!VALUES.D_NAN
     ;; Now do spectral correction for the CO2 and H2O
     IF keyword_set(corr_massman) THEN BEGIN
	;; First, do CO2 and H2O, then do H2O lagged by CO2 time constant.
	;; First, shift the CO2 and H2O signals appropriately
	co2_shift=Xco2_m[0 + abs(lag_co2_cl):nrecs - 1]
	w_shift_co2=wrot[0:nrecs - 1 - abs(lag_co2_cl)]
	h2o_shift=Xh2o_m[0 + abs(lag_h2o_cl):nrecs - 1]
	w_shift_h2o=wrot[0:nrecs - 1 - abs(lag_h2o_cl)]
	h2o_shift_WPL=Xh2o_m[0 + abs(lag_co2_cl):nrecs - 1]
	w_shift_h2o_WPL=wrot[0:nrecs - 1 - abs(lag_co2_cl)]
	;; Then, run the spectral correction
	cf_wXco2m=spec_massman(w_shift_co2[where(finite(w_shift_co2))], $
			       co2_shift[where(finite(w_shift_co2))], $
			       horwind, corr_massman[0], corr_massman[1], $
			       GEOM=[corr_massman[2], corr_massman[3]], $
                               TUBEATT=[corr_massman[7], corr_massman[8], $
                                        corr_massman[9], $
					mean(IRGA_P, /NAN), met_T], $
			       GAS='CO2')
	cf_wXh2om=spec_massman(w_shift_h2o[where(finite(w_shift_h2o))], $
			       h2o_shift[where(finite(w_shift_h2o))], $
			       horwind, corr_massman[0], corr_massman[1], $
			       GEOM=[corr_massman[2], corr_massman[3]], $
                               TUBEATT=[corr_massman[7], corr_massman[8], $
                                        corr_massman[9], $
					mean(IRGA_P, /NAN), met_T], $
			       GAS='H2O')
	ok=where(finite(w_shift_h2o_WPL))
	cf_wXco2m_forWPL=spec_massman(w_shift_h2o_WPL[ok], $
				      h2o_shift_WPL[ok], $
				      horwind, 0.1, 0.145, $
                                      GEOM=[corr_massman[2], $
                                            corr_massman[3]], $
                                      TUBEATT=[corr_massman[7], $
                                               corr_massman[8], $
					       corr_massman[9], $
					       mean(IRGA_P, /NAN), met_T], $
				      GAS='H2O')
	cov_w_Xco2_m=cov_w_Xco2_m * cf_wXco2m
	cov_w_Xh2o_m=cov_w_Xh2o_m * cf_wXh2om
	cov_w_Xh2o_m_forWPL=cov_w_Xh2o_m_forWPL * cf_wXco2m_forWPL
     ENDIF

     ;; Calculate the flux of CO2, applying a pseudo-WPL correction as per
     ;; eqn 3b of Ibrom et al. 2007.  Mean molar air concentration (moist
     ;; air).
     mean_c_m=(mean(IRGA_P, /nan) * 1000.0) / $
	      (R * (mean(IRGA_T, /nan) + 273.15))
     ;; Mean molar concentration of water vapour
     mean_c_v=mean(Xh2o_m, /nan) * mean_c_m
     mean_c_c=mean(Xco2_m, /nan) * mean_c_m ; mean molar concentration of co2
     ;; Mean dry air mixing ratio of co2
     mean_Xco2_d=mean_c_c / (mean_c_m - mean_c_v)

     c_m_atm=(met_P * double(1000)) / (R * (met_T + 273.15))
     ;; Flux of CO2, units mol/m2/s
     Fco2_cl=c_m_atm*[cov_w_Xco2_m + (mean_Xco2_d * cov_w_Xh2o_m_forWPL)]
     ;; Calculate the flux of H2O, no WPL correction necessary (assumes, as
     ;; in the previous step that temperature fluctuations are completely
     ;; attenuated)
     E_cl=c_m_atm * cov_w_Xh2o_m ; Flux of H2O, units mol/m2/s

     ;; Vaporization
     IF ((met_T) GE 0) THEN  Lv=(2.50057 - 0.00245715 * (met_T)) * 1000.0
     ;; Sublimation
     IF ((met_T) LT 0)	THEN Lv=(2.83539 - 0.000135713 * (met_T)) * 1000.0
     Qe_cl=double(E_cl * mv * Lv)

     RETURN, [cov_w_Xco2_m, cf_wXco2m, lag_co2_cl, Fco2_cl * 86400000.0, $
              cov_w_Xh2o_m, cf_wXh2om, lag_h2o_cl, E_cl * 86400.0, Qe_cl, $
              mean_Xco2_m * 1000000.0, mean(c_v, /nan) / mean(c_m, /nan)]

  ENDIF ELSE BEGIN

     ;; CALCULATE FLUXES BASED ON POINT-BY-POINT CONVERSIONS TO DRY AIR
     ;; MIXING RATIOS (NO WPL REQUIRED)

     ;; Do point-by-point conversions to dry air mixing ratio
     c_m=(IRGA_P * double(1000)) / $
	 (R * (IRGA_T + 273.15)) ; molar concentration of moist air
     c_v=Xh2o_m * c_m		 ; molar concentration of water vapour
     c_d=c_m - c_v		   ; molar concentration of dry air
     c_c=Xco2_m * c_m		   ; molar concentration of co2
     Xco2_d=c_c / (c_m - c_v)	   ; dry air mixing ratio of co2
     Xh2o_d=c_v / (c_m - c_v)	       ; dry air mixing ratio of h2o
     ;; Calculate means
     mean_Xco2_d=mean(Xco2_d, /NAN, /DOUBLE)
     mean_Xh2o_d=mean(Xh2o_d, /NAN, /DOUBLE)
     ;; Calculate covariances using max covariance (NOTE: h2o lag is
     ;; determined independently of co2)
     n_lag_cl=(maxc_c[1] - maxc_c[0]) + 1
     lag_cl=findgen(n_lag_cl)

     FOR lagpop=maxc_c[0], maxc_c[1] DO $
	lag_cl[lagpop - maxc_c[0]]=-lagpop

     w_Xco2_d=c_covariance(Xco2_d, wrot, lag_cl)
     w_Xh2o_d=c_covariance(Xh2o_d, wrot, lag_cl)
     ii_co2=where(abs(w_Xco2_d) EQ max(abs(w_Xco2_d)))
     ii_h2o=where(abs(w_Xh2o_d) EQ max(abs(w_Xh2o_d)))
     lag_co2_cl=lag_cl[ii_co2(0)]
     lag_h2o_cl=lag_cl[ii_h2o(0)]
     cov_w_Xco2_d=w_Xco2_d[ii_co2[0]]
     cov_w_Xh2o_d=w_Xh2o_d[ii_h2o[0]]

     IF keyword_set(ogive_ofile) THEN BEGIN
	ogive, 'wco2_cl', wrot, Xco2_d, ec_period, isample_freq, maxc_c, $
	       ogive_ofile
	ogive, 'wh2o_cl', wrot, Xh2o_d, ec_period, isample_freq, maxc_c, $
	       ogive_ofile
     ENDIF

     cf_wXco2d=!VALUES.D_NAN
     cf_wXh2od=!VALUES.D_NAN
     ;; Now do spectral correction for the CO2 and H2O
     IF keyword_set(corr_massman) THEN BEGIN
	;; First, do CO2 and H2O, then do H2O lagged by CO2 time constant.
	;; First, shift the CO2 and H2O signals appropriately.
	co2_shift=Xco2_m[0 + abs(lag_co2_cl):nrecs - 1]
	w_shift_co2=wrot[0:nrecs - 1 - abs(lag_co2_cl)]
	h2o_shift=Xh2o_m[0 + abs(lag_h2o_cl):nrecs - 1]
	w_shift_h2o=wrot[0:nrecs - 1 - abs(lag_h2o_cl)]
	;; Then, run the spectral correction
	cf_wXco2d=spec_massman(w_shift_co2[where(finite(w_shift_co2))], $
			       co2_shift[where(finite(w_shift_co2))], $
			       horwind, corr_massman[0], corr_massman[1], $
			       GEOM=[corr_massman[2], corr_massman[3]], $
                               TUBEATT=[corr_massman[7], corr_massman[8], $
                                        corr_massman[9], $
					mean(IRGA_P, /NAN), met_T], $
			       GAS='CO2')
	cf_wXh2od=spec_massman(w_shift_h2o[where(finite(w_shift_h2o))], $
			       h2o_shift[where(finite(w_shift_h2o))], $
			       horwind, corr_massman[0], corr_massman[1], $
			       GEOM=[corr_massman[2], corr_massman[3]], $
                               TUBEATT=[corr_massman[7], corr_massman[8], $
                                        corr_massman[9], $
					mean(IRGA_P, /NAN), met_T], $
			       GAS='H2O')
	cov_w_Xco2_d=cov_w_Xco2_d * cf_wXco2d
	cov_w_Xh2o_d=cov_w_Xh2o_d * cf_wXh2od
     ENDIF
     ;; Calculate fluxes as per equation 6.22 of Leuning 2004.	Calculate
     ;; the vapour pressure in the atmosphere (not in the IRGA) in Pa
     ev=mean(Xh2o_m, /NAN) * met_P * double(1000)
     ;; Calculate the dry air molar concentration in the atmosphere, mol/m3.
     c_d_atm=(met_P * 1000.0 - ev) / (R * (met_T + 273.15))
     Fco2_cl=c_d_atm * cov_w_Xco2_d ;mol/m2s
     E_cl=c_d_atm * cov_w_Xh2o_d    ;mol/m2s

     ;; Calculate latent heat flux
     IF keyword_set(open) THEN $
	checkT=mean_Tair - 273.15 $
     ELSE checkT=met_T
     IF ((checkT) GE 0) THEN $
	Lv=(2.50057 - 0.00245715 * (checkT)) * 1000.0 ; vaporization
     IF ((checkT) LT 0) THEN $
	Lv=(2.83539 - 0.000135713 * (checkT)) * 1000.0 ; sublimation
     Qe_cl=double(E_cl * mv * Lv)

     RETURN, [cov_w_Xco2_d, cf_wXco2d, lag_co2_cl, Fco2_cl * 86400000.0, $
              cov_w_Xh2o_d, cf_wXh2od, lag_h2o_cl, E_cl * 86400.0, Qe_cl, $
              mean_Xco2_d * 1000000.0, mean(c_v, /nan) / mean(c_m, /nan)]

  ENDELSE

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; ec_closed.pro ends here
