;; Author: Brent Else, Sebastian Luque
;; Created: 2013-11-15T21:32:25+0000
;; Last-Updated: 2015-06-30T19:49:40+0000
;;	     By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;;
;;     EC_OPEN
;;
;; PURPOSE:
;;
;;     This program calculates heat flux, latent heat flux and CO2 flux
;;     from an open path EC system.
;;
;; CALLING SEQUENCE:
;;
;;
;;
;; INPUTS:
;;
;;	 Wind:	        3xn element array of unrotated wind [u,v,w].
;;	 Ts:	        1xn array of sonic temperature (C).
;;	 CO2:	        1xn array of co2 values in units of mass
;;		        concentration (mmol/m3).
;;	 H2O:	        1xn array of h2o values in units of mass
;;		        concentration (mmol/m3).
;;	 P:	        1xn array of atmospheric pressure (kPa).
;;	 MaxC:	        Scalar integer indicating maximum number of records
;;                      that can be lagged for digital time shifting.
;;       EC_Period:     Scalar float with the averaging period under
;;                      consideration (seconds) (e.g. 30min=1200 s flux
;;                      runs).
;;       Isample_Freq:  Scalar float variable of the data frequency in Hz
;;                      (e.g. 10Hz).
;;
;; KEYWORD PARAMETERS:
;;
;;     CORR_MASSMAN:	 Set this keyword to run a spectral correction
;;			 following Massman this keyword must be used to set
;;			 a number of parameters for the massman correction.
;;			 ALL fields must be filled.  If they are not
;;			 applicable to the specific application, set them
;;			 to the string value 'NAn' (for not
;;			 applicable). [0-sample rate (in seconds), 1-sonic
;;			 path length (in m), 2-scalar sensor seaparation in
;;			 x coord (in m), 3-scalar sensor separation in y
;;			 coord (in m), 4-line length for scalar sensor with
;;			 line averaging (in m), 5-diameter of scalar sensor
;;			 with volume averaging (in m), 6-length of scalar
;;			 sensor with volume averaging (in m)].
;;     OGIVE_OFILES:	 Set this key word to produce an ogive plot for all
;;                       fluxes which are calculated using this routine.
;;                       This keyword must be set to a 3-element string
;;                       array indicating the paths where the ogive plot
;;                       files for air temperature, CO2, and H2O (in thar
;;                       order) will be saved.
;;     BURBA:		 Set this keyword to perform the Burba correction
;;			 (Burba et al. 2008, Glob. Ch. Biol.)  this keyword
;;			 must be used to set a number of parameters for the
;;			 Burb correction.  ALL fields must be filled.  If
;;			 they are not applicable to the specific
;;			 application, set them to the string value 'NAn'
;;			 (for not applicable).	[0 - Incoming shortwave
;;			 radiation (w/m2), 1 - Incoming LW radiation
;;			 (w/m2), 2 - Mean raw sonic wind speed (NO
;;			 MOTION CORRECTION) (m/s)].
;;     PKT:		 Set this keyword to run a correction as per
;;			 Prytherch et al. 2010 (doi: 10.1029/2009GL041482)
;;			 keyword MUST BE SET EQUAL TO u* (which can be
;;			 calculated in a previous call to EC_momentum.
;;
;; OUTPUTS:
;;
;;     20-element array:
;;
;;     [0]=w/Tair covariance
;;     [1]=correction factor for heat flux
;;     [2]=heat flux (W/m2)
;;     [3]=w/co2 covariance
;;     [4]=correction factor for co2 flux
;;     [5]=CO2 flux (mmol/m2/day), with WPL correction
;;     [6]=w/h2o covariance
;;     [7]=correction factor for h2o flux
;;     [8]=h2o flux (mol/m2/day)
;;     [9]Latent heat flux (W/m2)
;;     [10]=max covariance lag
;;     [11]=mean dry air co2 mixing ratio (ppm)
;;     [12]=mean specific humidity (g/g)
;;     [13]=WPL contribution to CO2 flux (mmol/m2/d)
;;     [14]=Burba correction value for CO2 using mutlivariate approach (IS
;;	    NOT ADDED TO CO2 FLUX... MUST BE ADDED AFTER)
;;     [15]=Burba correction value for CO2 using linear approach (IS NOT
;;	    ADDED TO CO2 FLUX... MUST BE ADDED AFTER)
;;     [16]=Burba correction value for H2O using multivariate approach (IS
;;	    NOT ADDED TO CO2 FLUX... MUST BE ADDED AFTER)
;;     [17]=Burba correction value H2O using linear approach (IS NOT ADDED
;;	    TO CO2 FLUX... MUST BE ADDED AFTER)
;;     [18]=CO2 flux (mmol/m2/day) calculated with PKT correction (includes
;;	    WPL correction, but no BURBA correction)
;;     [19]=Number of iterations required for PKT correction to converge
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

FUNCTION EC_OPEN, WIND, TS, CO2, H2O, P_ATM, MAXC, EC_PERIOD, ISAMPLE_FREQ, $
                  CORR_MASSMAN=CORR_MASSMAN, OGIVE_OFILE=OGIVE_OFILE, $
                  BURBA=BURBA, PKT=PKT

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
  nrecs=n_elements(co2)

  ;; Calculate the mean horizontal wind speed (unrotated)... this is
  ;; required in the spectral correction
  horwind=[mean(WIND[0, *], /nan), mean(WIND[1, *], /nan)]
  ;; do the wind rotations
  WINDrot=yawpitch(WIND[0, *], WIND[1, *], WIND[2, *], $
		   n_elements(WIND[0, *]))
  urot=WINDrot[0, *] & vrot=WINDrot[1, *] & wrot=WINDrot[2, *]

  ;;========CALCULATE NECESSARY TERMS=================

  ;; Get units of co2 and h2o in mol/m3
  c_co2=co2 / double(1000)
  c_h2o=h2o / double(1000)
  mean_c_h2o=mean(c_h2o, /DOUBLE, /NAN)
  mean_c_co2=mean(c_co2, /DOUBLE, /NAN)

  ;; Get Ts into Tair
  P=P_atm * double(1000)		  ; Atmospheric pressure in Pa
  Tair=sonicT2airT(Ts, c_h2o, P)
  Tair=Tair + 273.15		; Air temperature into K
  mean_Tair=mean(Tair, /nan)

  ev=c_h2o * R * Tair		 ; water vapour pressure (Pa)
  rho_v=c_h2o * mv		 ; water vapour density (g/m3)
  mean_rho_v=mean(rho_v, /DOUBLE, /NAN) ; mean water vapour density (g/m3)
  ;; mean dry air density (g/m3)
  mean_rho_d=mean((P - ev) / (R * Tair), /DOUBLE, /NAN) * ma
  mean_rho=mean_rho_d + mean_rho_v ; mean moist air density (g/m3)
  ;; mean specific humidity (g_h2o/g_moist air)
  mean_qh2o=mean_rho_v / mean_rho
  ;; mean specific heat capacity of moist air (J/g/K)
  mean_cp=cpd * (double(1) + 0.84 * mean_qh2o)
  ;; mean air molar concentration (mol/m3)
  mean_c=(mean_rho_d) / (ma + mean_rho_v / mv)
  mean_c_dry=mean_rho_d / ma	;mean dry air molar concentration
  mean_Xv=mean_c_h2o / mean_c_dry ; mean dry air mixing ratio of water
  ;; mean dry air mixing mixing ration of co2 (in ppm)
  co2ppm=mean_c_co2 / mean_c_dry * double(1000000)

  ;;========DO H CALCULATIONS UTILIZING THE OPEN PATH=========

  ;;calculate w/Tair covariance
  cov_w_Tair=correlate(wrot, Tair, /COVARIANCE, /DOUBLE)

  IF keyword_set(ogive_ofiles) THEN $
     ogive, 'wTair', wrot, Tair, ec_period, isample_freq, maxc, $
            ogive_ofiles[0]
  ;; Do spectral corrrection
  cf_wTair=!VALUES.D_NAN
  IF keyword_set(corr_massman) THEN BEGIN
     cf_wTair=spec_Massman(wrot[where(finite(wrot))], $
			   Tair[where(finite(wrot))], horwind, $
			   corr_massman[0], corr_massman[1], $
                           GEOM=[corr_massman[2], corr_massman[3]], $
                           SCALARLINE=corr_massman[1])
     cov_w_Tair=cf_wTair * cov_w_Tair
  ENDIF

  ;;calculate heat flux
  H=mean_rho * mean_cp * cov_w_Tair ;  Units: W/m2

  ;;===============CALCULATE CO2/H2O COVARIANCES==================
  n_lag_op=(maxc * 2) + 1
  lag_op=findgen(n_lag_op) - maxc
  w_c_co2=c_covariance(c_co2, wrot, lag_op)
  w_c_h2o=c_covariance(c_h2o, wrot, lag_op)
  ii_co2=where(abs(w_c_co2) EQ max(abs(w_c_co2)))
  lag_co2_op_use=ii_co2[0]
  ii_h2o=where(abs(w_c_h2o) EQ max(abs(w_c_h2o)))
  lag_h2o_op_use=ii_h2o[0]

  ;; For some additional analysis purposes, we're going to calculate the
  ;; lagged w_Ts and w_u covariances
  w_Tair=c_covariance(Tair, wrot, lag_op)
  w_u=c_covariance(urot, wrot, lag_op)

  ;; Now, decide what we're going to use for the lags of each of
  ;; these... it's a bit complicated.  [SPL: removed all these keywords,
  ;; since no longer used
  lag_h2o_op_use=lag_co2_op_use

  ;; Calculate the covariances, using the appropriate lags
  cov_w_c_co2=w_c_co2[lag_co2_op_use]
  cov_w_c_h2o=w_c_h2o[lag_h2o_op_use]
  lag_co2_op=lag_op[lag_co2_op_use]
  lag_h2o_op=lag_op[lag_h2o_op_use]

  ;;==========CALL OGIVE==========================================
  IF keyword_set(ogive_ofiles) THEN BEGIN
     ogive, 'wco2_op', wrot, c_co2, ec_period, isample_freq, maxc, $
            ogive_ofiles[1]
     ogive, 'wh2o_op', wrot, c_h2o, ec_period, isample_freq, maxc, $
            ogive_ofiles[2]
  ENDIF

  ;;==============SPECTRAL CORRECTION=============================
  IF keyword_set(corr_massman) THEN BEGIN
     ;; first, the signals need to be digitally shifted by the lags
     ;; calculated above.
     IF lag_co2_op LT 0 THEN BEGIN
	co2_shift=c_co2[0 + abs(lag_co2_op):nrecs - 1]
	w_shift_co2=wrot[0:nrecs - 1 - abs(lag_co2_op)]
     ENDIF
     IF lag_co2_op GE 0 THEN BEGIN
	co2_shift=c_co2[0:nrecs - 1 - lag_co2_op]
	w_shift_co2=wrot[0 + lag_co2_op:nrecs - 1]
     ENDIF
     IF lag_co2_op LT 0 THEN BEGIN
	h2o_shift=c_h2o[0 + abs(lag_co2_op):nrecs - 1]
	w_shift_h2o=wrot[0:nrecs - 1 - abs(lag_co2_op)]
     ENDIF
     IF lag_co2_op GE 0 THEN BEGIN
	h2o_shift=c_h2o[0:nrecs - 1 - lag_co2_op]
	w_shift_h2o=wrot[0 + lag_co2_op:nrecs - 1]
     ENDIF
     ;; now do the spectral correction on the shifted signals
     cf_wco2=spec_massman(w_shift_co2[where(finite(w_shift_co2))], $
			  co2_shift[where(finite(co2_shift))], horwind, $
			  corr_massman[0], corr_massman[1], $
			  GEOM=[corr_massman[2], corr_massman[3]], $
			  SCALARVOL=[corr_massman[5], corr_massman[6]])
     cf_wh2o=spec_massman(w_shift_h2o[where(finite(w_shift_h2o))], $
			  h2o_shift[where(finite(w_shift_h2o))], $
			  horwind, corr_massman[0], corr_massman[1], $
			  GEOM=[corr_massman[2], corr_massman[3]], $
			  SCALARVOL=[corr_massman[5], corr_massman[6]])
     cov_w_c_co2=cov_w_c_co2 * cf_wco2
     cov_w_c_h2o=cov_w_c_h2o * cf_wh2o
  ENDIF ELSE BEGIN
     cf_wco2=!VALUES.D_NAN
     cf_wh2o=!VALUES.D_NAN
  ENDELSE

  ;; ;;=====================================================================
  ;; ;;Write out the lags and covariances to a temporary file (this is for
  ;; ;;visualizing the maximum covariance signal strength)
  ;; lag_arr=fltarr(5, n_lag_op)
  ;; lag_arr[0, *]=lag_op & lag_arr[1, *]=w_c_co2 & lag_arr[2, *]=w_c_h2o
  ;; lag_arr[3, *]=w_Tair & lag_arr[4, *]=w_u
  ;; openw, 7, 'lag_op.txt'
  ;; printf, 7, lag_arr, format='(5A)'
  ;; close, 7

  ;;===========BURBA CORRECTION=========================================
  IF keyword_set(burba) THEN BEGIN
     SW=burba[0] & LW=burba[1] & rawU=burba[2]
     meanT_C=mean_Tair - 273.15
     ;; first, calculate temperatures of all of the different LI7500
     ;; components using the multivariate and linear approaches.  These are
     ;; the calculations for "daytime" conditions, which we'll define as SW
     ;; > 5 w/m2
     IF SW GE 5 THEN BEGIN
	Tbot_mt=meanT_C + 2.8 - 0.0681 * meanT_C + $
		0.0021 * SW - 0.334 * rawU ; multiple regression approach
	Ttop_mt=meanT_C - 0.1 - 0.0044 * meanT_C + $
		0.0011 * SW - 0.022 * rawU ; multiple regression approach
	Tspar_mt=meanT_C + 0.3 - 0.00077 * meanT_C + 0.0006 * $
		 SW - 0.044 * rawU ; multiple regression approach
	Tbot_ln=0.944 * meanT_C + 2.57 ;linear regression approach
	Ttop_ln=1.005 * meanT_C + 0.24 ;linear regression approach
	Tspar_ln=1.01 * meanT_C + 0.36 ;linear regression approach
     ENDIF ELSE BEGIN
	Tbot_mt= meanT_C + 0.5 - 0.116 * meanT_C + $
		 0.0087 * LW - 0.206 * rawU ; multiple regression approach
	Ttop_mt=meanT_C - 1.7 - 0.016 * meanT_C + $
		0.0051 * LW - 0.029 * rawU ; multiple regression approach
	Tspar_mt=meanT_C - 2.1 - 0.02 * meanT_C + $
		 0.007 * LW + 0.026 * rawU ; multiple regression approach
	Tbot_ln=0.883 * meanT_C + 2.17	   ; linear regression approach
	Ttop_ln=1.008 * meanT_C - 0.41	   ;linear regression approach
	Tspar_ln=1.01 * meanT_C - 0.17	   ;linear regression approach
     ENDELSE
     ;; This empirical approach is far from perfect.... If Tbot or Tspar
     ;; are lt air T, set those to air T (it should never get cooler than
     ;; air T)
     IF (Tbot_mt LT meanT_C) AND (Ttop_mt LT meanT_C) THEN Tspar_mt=meanT_C
     IF Tbot_mt LT meanT_C THEN Tbot_mt=meanT_C
     ;;   if Ttop_mt lt meanT_C then begin
     ;;     Ttop_mt=meanT_C
     ;;   endif
     ;;
     ;;   if (Tbot_ln lt meanT_C) and (Ttop_ln lt meanT_C) then begin
     ;;     Tspar_ln=meanT_C
     ;;   endif
     ;;   if Tbot_ln lt meanT_C then begin
     ;;     Tbot_ln=meanT_C
     ;;   endif
     ;;   if Ttop_ln lt meanT_C then begin
     ;;     Ttop_ln=meanT_C
     ;;   endif

     ;; now, calculate the heat fluxes for these various
     ;; components. Dimensions of the LI7500 (m)
     lspar=0.005 & rtop=0.0225 & rspar=0.0025 & lbot=0.065 & ltop=0.045
     ;; Boundary layer thickness
     del_bot=0.004 + 0.004 * sqrt(lbot / rawU)
     del_top=0.0045 + 0.0028 * sqrt(ltop / rawU) + 0.00025 / rawU
     del_spar=0.0058 * sqrt(lspar / rawU)
     ;; Individual heat fluxes
     kair=0.0243		; Thermal conductivity of air
     Sbot_mt=kair * ((Tbot_mt - meanT_C) / del_bot)
     Stop_mt=kair * ((rtop + del_top) * (Ttop_mt - meanT_C) / $
		     (rtop * del_top))
     Sspar_mt=kair *  (Tspar_mt - meanT_C) / $
	      (rspar * alog((rspar + del_spar) / rspar))
     Ssensor_mt=Sbot_mt + Stop_mt + (0.15 * Sspar_mt) ; (W/m2)

     Sbot_ln=kair * ((Tbot_ln - meanT_C) / del_bot)
     Stop_ln=kair * ((rtop + del_top) * (Ttop_ln - meanT_C) / $
		     (rtop * del_top))
     Sspar_ln=kair * (Tspar_ln - meanT_C) / $
	      (rspar * alog((rspar + del_spar) / rspar))
     Ssensor_ln=Sbot_ln + Stop_ln + (0.15 * Sspar_ln) ; (W/m2)
  ENDIF

  ;;==================CALCULATE FLUXES===================================
  ;;Calculate the water vapour flux, applying the WPL correction in terms
  ;;of H [EQ'N 6.26, Handbook of Micrometeorology] THIS SHOULD HAVE THE
  ;;BURBA CORRECTION TOO....
  E_op=(double(1) + mean_Xv) * $
       [cov_w_c_h2o + (mean_c_h2o / mean_Tair) * $
	(H / (mean_rho * mean_cp))] ;  Units: mol/m2s
  ;; Calculate the CO2 flux, pre-WPL correction
  noWPL_Fco2_op=cov_w_c_co2	; Units: mol/m2
  ;; Calculate CO2 flux, applying the WPL correction in terms of H and E
  WPL_Fco2_op=cov_w_c_co2 + mean_c_co2 * $
	      [(E_op / mean_c) + $
	       (H / (mean_rho * mean_cp * mean_Tair))] ;  Units: mol/m2s
  ;; Calculate contribution of WPL term to the overall flux (for interest's
  ;; sake, really)
  WPLcont=WPL_Fco2_op * double(86400000) - noWPL_Fco2_op * double(86400000)

  IF keyword_set(burba) THEN BEGIN
     BURBA_E_op_mt=(double(1) + mean_Xv) * $
		   [cov_w_c_h2o + (mean_c_h2o / mean_Tair) * $
		    ([H + Ssensor_mt] / (mean_rho * mean_cp))] ; mol/m2s
     BURBA_E_op_ln=(double(1) + mean_Xv) * $
		   [cov_w_c_h2o + (mean_c_h2o / mean_Tair) * $
		    ([H + Ssensor_ln] / (mean_rho * mean_cp))] ; mol/m2s
     H2OBurba_mt=BURBA_E_op_mt * double(86400) - E_op
     H2OBurba_ln=BURBA_E_op_ln * double(86400) - E_op

     BURBA_Fco2_op_mt=cov_w_c_co2 + mean_c_co2 * $
		      [(BURBA_E_op_mt / mean_c) + $
		       ([H + Ssensor_mt] / (mean_rho*mean_cp * mean_Tair))]
     CO2Burba_mt=BURBA_Fco2_op_mt * double(86400000) - $
		 WPL_Fco2_op * double(86400000)

     BURBA_Fco2_op_ln=cov_w_c_co2 + mean_c_co2 * $
		      [(BURBA_E_op_ln / mean_c) + $
		       ([H + Ssensor_ln] / $
			(mean_rho * mean_cp * mean_Tair))]
     CO2Burba_ln=BURBA_Fco2_op_ln * double(86400000) - $
		 WPL_Fco2_op * double(86400000)
  ENDIF ELSE BEGIN
     FINALco2_op=WPL_Fco2_op
     CO2Burba_mt=!VALUES.D_NAN & CO2Burba_ln=!VALUES.D_NAN
     H2OBurba_mt=!VALUES.D_NAN & H2OBurba_ln=!VALUES.D_NAN
     FINAL_E_op=E_op
  ENDELSE

  ;; Calculate latent heat flux
  IF ((mean_Tair - 273.15) GE 0) THEN $ ; vaporization
     Lv=(2.50057 - 0.00245715 * (mean_Tair - 273.15)) * double(1000)
  IF ((mean_Tair - 273.15) LT 0) THEN $ ; sublimation
     Lv=(2.83539 - 0.000135713 * (mean_Tair - 273.15)) * double(1000)
  Qe_op=double(E_op * mv * Lv)

  ;;========DO PKT CORRECTION==============
  IF keyword_set(pkt) THEN BEGIN
     ustar=pkt
     ;; First, calculate mixing ratio of co2 on point-by-point basis, as
     ;; per Fairall et al. (2000) (doi: 10.1023/A:1002662826020) eq'n 62:
     rhoc=c_co2 * mc		; co2 density in g/m3
     rhov=c_h2o * mv		; h2o density in g/m3
     ;; Shift the rhoc/rhov signals by the lag to make sure they line up
     ;; with T properly.
     IF lag_co2_op LT 0 THEN BEGIN
	rhoc_shift=rhoc[0 + abs(lag_co2_op):nrecs - 1]
	rhov_shift=rhov[0 + abs(lag_co2_op):nrecs - 1]
	;; calculated water vapour pressure... need for calc'n rh
	ev_shift=ev[0 + abs(lag_co2_op):nrecs - 1]
	w_shift=wrot[0:nrecs - 1 - abs(lag_co2_op)]
	T_shift=Tair[0:nrecs - 1 - abs(lag_co2_op)]
     ENDIF
     IF lag_co2_op GE 0 THEN BEGIN
	rhoc_shift=rhoc[0:nrecs - 1 - lag_co2_op]
	rhov_shift=rhov[0:nrecs - 1 - lag_co2_op]
	ev_shift=ev[0:nrecs - 1 - lag_co2_op]
	w_shift=wrot[0 + lag_co2_op:nrecs - 1]
	T_shift=Tair[0 + lag_co2_op:nrecs - 1]
     ENDIF
     ;; calculate the point-by-point mixing ratio of CO2, RH for H2O and
     ;; Mixing Ratio (kg/kg) of H2O.
     ;; Fluctuations of rhoc shift
     rhoc_shift_pr=rhoc_shift - (mean(rhoc_shift, /nan))
     ;; Fluctuations of rhov shift
     rhov_shift_pr=rhov_shift - (mean(rhov_shift, /nan))
     T_shift_pr=T_shift - (mean(T_shift, /nan)) ; fluctuations of T shift

     ;; High frequency dry air density (from LICOR/SONIC)
     rhod_shift=((P - ev_shift) / (R * T_shift)) * ma
     ;; calculate fluctuating term of the co2 dry air mixing ratio:
     Xc_shift_pr=(rhoc_shift_pr + $
		  [(ma * rhov_shift_pr) / $
		   (mv * rhod_shift) + (double(1) + $
					(ma * rhov_shift) / $
					(mv * rhod_shift)) * $
		   (T_shift_pr / T_shift)] * (mean_c_co2 * mc)) / $
		 mean_rho_d
     ;; Not sure if we should be using dry air density or "moist" air
     ;; density... shouldn't make a big difference, but Prytech and Fairall
     ;; are not very specific
     Xc_shift=Xc_shift_pr + mean(rhoc_shift, /nan) / mean(rhod_shift, /nan)
     ;; Reconstruct the full term. NOTE: the units here are g of C/g of Dry
     ;; Air... this is not really a mixing ratio the way I usually
     ;; calculate it (i.e. mols C/mol dry air).
     ;; Calculate the various water vapor terms they need.
     ;; Sat'n vapour pressure (Pa)
     es_shift=[6.112 * exp((17.67 * (T_shift - 273.15)) / $
			   ((T_shift - 273.15) + 243.5))] * double(100)
     ;; High frequency relative humidity (from LICOR)
     RH_shift=(ev_shift / es_shift) * 100
     ;; High frequency H2O mixing ratio... This does not include the
     ;; dilution effect on this term due to temperature fluctuations... May
     ;; want to add that!
     Xv_shift=rhov_shift / rhod_shift
     ;; Calculate q*
     E_op_kg=E_op * mv / double(1000) ; convert mol/m2s to kg/m2s
     qstar=E_op_kg[0] / ustar

     ;; Calculate RH-q relationship... we need to do things slightly
     ;; different to recreate the exact results of Matlab polyfit.
     Xv_Mu1=mean(Xv_shift, /NAN)
     Xv_Mu2=stddev(Xv_shift, /NAN)
     Xv_hat=(Xv_shift - Xv_Mu1) / Xv_Mu2
     fitXv=poly_fit(Xv_hat, RH_shift, 1)
     drh_by_dq=fitXv[1] / Xv_Mu2
     ;; Detrend CO2 mixing ratio with respect to RH
     RH_Mu1=mean(RH_shift, /NAN)
     RH_Mu2=stddev(RH_shift, /NAN)
     RH_hat=(RH_shift - RH_Mu1) / RH_Mu2
     fitRH=poly_fit(RH_hat, Xc_shift, 3)
     slope=fitRH[3] / RH_Mu2
     ;; + fitRH(0); for some reason, the Prytech code leaves off the last
     ;; coefficient... not sure why
     response=fitRH[3] * RH_hat ^ double(3) + $
	      fitRH[2] * RH_hat ^ double(2) + $
	      fitRH[1] * RH_hat
     Xc_det=Xc_shift-response
     ;; Calculate CO2 flux using the detrended Xc signal
     cflux=c_covariance(w_shift, Xc_det, 0) * mean_rho_d ; flux in g/m2s
     cflux=cflux[0] / double(1000)			 ; kg/m3s
     cstar=cflux / ustar				 ; kg/m3s

     dc_by_dq=cstar / qstar
     ;; Determine c*/RH* ratio (Eq 2, Prytherch et al, 2009)
     dc_by_drh=dc_by_dq / drh_by_dq

     cflux=cflux * double(1000) / mc ; mols/(m2.s)
     cflux=cflux * double(60) * double(60) * $
	   double(24) * double(365)	; mols/(m2 yr)

     cflux_det=cflux * 2.737909263
     dc_by_drh_initial=dc_by_drh
     dc_by_dq_initial=dc_by_dq

     ;;  Iterate CO2 until we get some sort of convergence.  I think this
     ;; is done so that the while loop doesn't immediately die due to the
     ;; cfluxold-cflux < 1 criteria
     cfluxold=cflux + 2
     loop=double(0)
     WHILE (abs(cflux - cfluxold) GT 1) AND $
	(loop LT 100) AND (abs(cflux) LT 1000) DO BEGIN
	loop=loop + 1
	cfluxold=cflux
	Xc_nu=Xc_det + $
	      ((RH_shift - mean(RH_shift, /nan)) * dc_by_drh) / double(2)
	;; Now redo the flux calculatioon
	cflux=c_covariance(w_shift, Xc_nu, 0) * mean_rho_d ; flux in g/m2s
	cflux=cflux[0] / double(1000)			   ; kg/m3s
	cstar=cflux / ustar				   ; kg/m3s

	dc_by_dq=cstar / qstar
	dc_by_drh=dc_by_dq / drh_by_dq ; new c*/rh* ratio

	cflux=cflux * 1000 / mc  ; % mols/(m2.s)
	cflux=cflux * double(60) * double(60) * $
	      double(24) * double(365) ; % mols/(m2 yr)
     ENDWHILE

     dc_by_drh_final=dc_by_drh
     dc_by_dq_final=dc_by_dq

     pkt_FCO2_op=cflux * 2.737909263
     pkt_loop=loop

  ENDIF ELSE BEGIN
     pkt_FCO2_op=!VALUES.D_NAN & pkt_loop=!VALUES.D_NAN
     cflux_det=!VALUES.D_NAN & dRH_by_dq=!VALUES.D_NAN
     dc_by_drh_initial=!VALUES.D_NAN & dc_by_drh_final=!VALUES.D_NAN
     dc_by_dq_initial=!VALUES.D_NAN & dc_by_dq_final=!VALUES.D_NAN
  ENDELSE
  ;; Mean dry air density (g/m3)
  mean_rho_d=mean((P - ev) / (R * Tair), /DOUBLE, /NAN) * ma

  RETURN, [cov_w_Tair, cf_wTair, H, cov_w_c_co2, cf_wco2, $
           WPL_Fco2_op * double(86400000), cov_w_c_h2o, cf_wh2o, $
           E_op * double(86400), Qe_op, lag_co2_op, co2ppm, mean_qh2o, $
           WPLcont, CO2Burba_mt, CO2Burba_ln, H2OBurba_mt, H2OBurba_ln, $
           pkt_FCO2_op, pkt_loop, cflux_det, dRH_by_dq, dc_by_drh_initial, $
           dc_by_drh_final, dc_by_dq_initial, dc_by_dq_final, lag_h2o_op]
END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
