;; $Id$
;; Author: Brent Else, Sebastian Luque
;; Created: 2013-11-15T18:01:54+0000
;; Last-Updated: 2013-11-27T22:09:28+0000
;;           By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     SPEC_MASSMAN
;; 
;; PURPOSE:
;; 
;;     This program performs spectral correction by using the Massman 2000
;;     time constants, and applying the first-order transfer functions to
;;     the cospectra and doing the integration.  NOTE: digital shifting for
;;     maximum covariance must be done BEFORE calling this script.
;; 
;; CALLING SEQUENCE:
;; 
;; 
;; 
;; INPUTS:
;; 
;;     VertWind:   Vector of rotated vertical wind velocity.
;; 
;;     Scalar:     Vector with scalar to perform eddy covariance with (CAN
;;                 ALSO BE THE HORIZONTAL WIND FOR MOM. FLUX).
;;     HorWind:    2-element vector [mean_u,mean_v] with average of u and v
;;                 from sonic anemometer over the EC period.
;;     Rate:       Scalar float with sampling rate in s (e.g. for 10Hz,
;;                 rate=0.1).
;;     SonicLine:  Scalar float describing the sonic anemometer path length
;;                 in m.
;; 
;; KEYWORD PARAMETERS:
;; 
;;     GEOM:       Set keyword equal to a 2-element vector [x,y] describing
;;                 the scalar sensor separation from the sonic anemometer,
;;                 using a RH coordinate system with its origin equal to
;;                 the axis of the sonic [REQUIRED FOR SCALAR CORRECTION,
;;                 but can be set to [0,0]].
;;     MOMENTUM:   Set keyword if calculating spectral correction for
;;                 momentum flux.
;;     SCALARLINE  Set keyword as a scalar to define the path length of the
;;                 scalar sensor being corrected.
;;     SCALARVOL:  Set keyword as a 2-element vector [diameter,length] to
;;                 describe the diameter and length of a right circular
;;                 cylinder of the scalar sensor being used to calculate
;;                 volume averaging effects.
;;     TUBEATT:    Set keyword to do theoretical tube attenuation.  Set a
;;                 5-element vector [tube flow,tube diam,tube L,irga P,air
;;                 T] that describes the average tube flow (in LPM), the
;;                 tube diameter (in m), the tube length (in m), irga
;;                 pressure (kPa) and air temperature (deg C) during the
;;                 eddy covariance period.
;;     GAS:        Set keyword as a string, either 'H2O' or 'CO2' to define
;;                 if spectral correction is for H2O or CO2 [REQUIRED IF
;;                 TUBE ATT SET].
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
;; SOURCES:
;; 
;;     Massman (2000) A simple method for estimating frequency response
;;     corrections for eddy covariance systems, Agricultural and forest
;;     Meteorology 104, p 185-198.
;;                    
;;- -----------------------------------------------------------------------
;;; Code:

FUNCTION SPEC_MASSMAN, VERTWIND, SCALAR, HORWIND, RATE, SONICLINE, $
                       GEOM=GEOM, MOMENTUM=MOMENTUM, SCALARLINE=SCALARLINE, $
                       SCALARVOL=SCALARVOL, TUBEATT=TUBEATT, GAS=GAS

  ;; Mean wind velocity from sonic anemometer
  meanu=sqrt(horwind[0] ^ 2 + horwind[1] ^ 2)
  ;; Frequency range we're dealing with. Check if number of samples is even
  ;; or odd (if n is even, ev_odd=0, otherwise ev_odd=1)
  N=float(n_elements(vertwind))
  ev_odd=0
  IF (N AND 1) THEN ev_odd=1

  ;; Calculate the nyquist frequency (not to sure how to designate nyquist
  ;; frequency when N is odd... I think this is ok).
  nf=N / double(2)
  IF ev_odd GT 0 THEN nf=(N + double(1)) / double(2)

  ;; Calculate frequency range
  P=(N + double(1)) * rate      ; units=s
  xaxis_n=findgen(1, N + 1) + 1  
  xaxis_f=xaxis_n[1:nf] / P

  ;; Calculate the first order time constants.  Define an overall time
  ;; constant, and set it to zero.  We will build on this time constant
  ;; throughout the script as per equation 9 of Massman 2000
  t_overall=0.0

  ;; First, check if we're doing momentum flux
  IF keyword_set(momentum) THEN BEGIN
     ;; Calculate time constants for sonic anemometer line averaging
     t_w=sonicline / (5.7 * meanu) ; for the vertical wind
     ;; For horizontal wind... we will consider this the scalar
     t_u=sonicline / (2.8 * meanu)
     t_overall=t_overall + t_u ^ 2 + t_w ^ 2
  ENDIF ELSE BEGIN
     ;; Not doing momentum flux, proceed with scalar transfer function
     ;; calculations.  Calculate lateral and longitudinal separation
     ;; distances
     scalar_x=geom[0] & scalar_y=geom[1]
     IF (scalar_x EQ 0.0) AND (scalar_y EQ 0.0) THEN BEGIN
        lat=0.0 & lon=0.0
     ENDIF ELSE BEGIN
        angle=atan(scalar_y / scalar_x)
        ;; calculate hypotenuse and angle if wind direction=0
        hyp=SQRT(scalar_x ^ 2 + scalar_y ^ 2)
        windangle=atan(horwind[1] / horwind[0]) ; wind angle
        lat=abs(sin(angle + windangle) * hyp)   ; longitudinal separation
        lon=abs(cos(angle + windangle) * hyp)   ; lateral separation
     ENDELSE
     ;; Calculate time constants for longitudinal and lateral separation
     t_lat=lat / (1.1 * meanu)
     t_lon=lon / (1.05 * meanu)
     t_overall=t_overall + t_lat ^ 2 + t_lon ^ 2
     ;; Calculate sonic anemometer line averaging transfer function
     t_wline=sonicline / (8.4 * meanu)
     t_overall=t_overall + t_wline ^ 2
     ;; Calculate scalar sensor line averaging transfer function
     IF keyword_set(scalarline) THEN BEGIN
        t_scline=scalarline / (4.0 * meanu)
        t_overall=t_overall + t_scline ^ 2
     ENDIF
     ;; Calculate scalar sensor volume averaging transfer function
     IF keyword_set(scalarvol) THEN BEGIN
        t_scvol=(0.2 + 0.4 * (scalarvol[0] / scalarvol[1])) * $
                (scalarvol[1] / meanu)
        t_overall=t_overall + t_scvol ^ 2
     ENDIF
     ;; Calculate tube attenuation transfer function
     IF keyword_set(tubeatt) THEN BEGIN
        ;; Reynolds number calculation.  Problems: is IRGAP representative
        ;; of P through the line?  What is the temperature in the line?
        ;; Dynamic viscosity, units kg/ms SOURCE: Sutherland's formula,
        ;; http://en.wikipedia.org/wiki/Viscosity
        dynvisc=0.00001827 * [(291.15 + double(120)) / $
                              ((tubeatt[4] + 273.15) + double(120))] * $
                ((tubeatt[4] + 273.15) / 291.15) ^ (double(3) / double(2))
        ;; Estimate density of air travelling through the tube (it's an
        ;; estimate because we don't know what the air temperature is in
        ;; the line).  Note R=287 is gas constant in terms of J/K/kg
        rho_a=(tubeatt[3] * double(1000)) / (double(287) * $
                                             (tubeatt[4] + 273.15))
        ;; kinematic viscosity, units m2/s SOURCE:
        ;; http://en.wikipedia.org/wiki/Viscosity
        kinvisc=dynvisc / rho_a
        Qtube=tubeatt[0] * 0.0000167 ; convert from LPM to m3/s
        Re=(Qtube * tubeatt[1]) / $
           (kinvisc * !PI * (tubeatt[1] / double(2)) ^ 2) ; Reynolds number
        ;; Differentiate between a couple of parameters for CO2 and H2O
        IF keyword_set(gas) THEN BEGIN
           IF gas EQ 'CO2' THEN BEGIN
              ;; Coefficients for calculating Lambda, as per Shimizu 2007
              ;; (eq'n 8), and CO2 diffusivity (m2/s) as per
              ;; http://www.cco.caltech.edu/~brokawc/Bi145/Diffusion.html
              c1=0.75 & c2=0.04 & diff=0.000016
           ENDIF ELSE BEGIN
              ;; Coefficients for calculating Lambda, as per Shimizu 2007
              ;; (eq'n 8), and H2O diffusivity (m2/s) as per
              ;; http://www.cco.caltech.edu/~brokawc/Bi145/Diffusion.html
              c1=0.76 & c2=0.039 & diff=0.000025
           ENDELSE
        ENDIF ELSE BEGIN
           message, 'When performing tube attenuation spectral correction, ' + $
                    'keyword GAS must be set'
        ENDELSE
        ;; Calculate time constant, and transfer function for tube
        ;; attenuation. Shimizu 2007, eq'n 8
        lambda=0.639 / [8 * !PI ^ 2 * (alog(c1 * Re ^ c2)) ^ 2]
        t_tube=sqrt(lambda * tubeatt[2] * (tubeatt[1] / double(2))) / $
               (0.83 * (Qtube / (!PI * (tubeatt[1] / double(2)) ^ 2)))
        t_overall=t_overall + t_tube ^ 2
     ENDIF
     ;; Calculate the low frequency attenuation
  ENDELSE

  ;; Calculate the overall time constant, and associated first order
  ;; transfer function for the cospectra.
  t_overall=sqrt(t_overall)
  Tfcospec=double(1) / (double(1) + $
                        (double(2) * !PI * xaxis_f * t_overall) ^ 2)

  ;; Calculate cospectra while applying transfer functions

  ;; Calculate the FFTS for both variables
  fft_vert=fft([vertwind]) & fft_scal=fft([scalar])
  real_vert=real_part(fft_vert) & imag_vert=imaginary(fft_vert)
  real_scal=real_part(fft_scal) & imag_scal=imaginary(fft_scal)
  ;; Calculate the uncorrected cospectra, and sum it up to get the
  ;; uncorrected covariance
  uncorrCo=real_vert * real_scal + imag_vert * imag_scal
  uncorr_cov=total(uncorrCo(1:n_elements(uncorrCo)-1))

  ;; I'm not 100% sure how I want to apply the transfer functions to the
  ;; "two-sided" spectra and cospectra.  METHOD 1: apply a "mirror" image
  ;; of the transfer function so that the frequencies above the nyquist
  ;; frequency get the same transfer function as they should get when
  ;; they're "folded" over.  I believe this is the best way to do it.
  fcount=n_elements(xaxis_f)
  index=findgen(fcount)
  invindex=abs(index(1:fcount - 1) - index(fcount - 1))
  Tfcospec_full=findgen(fcount * 2 - 1)
  Tfcospec_full[index]=Tfcospec[index]
  Tfcospec_full[fcount:fcount * 2 - 2]=Tfcospec_full[invindex]

  ;; Now, calculate the cospectrum and apply the transfer function
  corrCo=uncorrCo[1:n_elements(uncorrCo) - 1] / Tfcospec_full

  ;; Now sum up the corrected cospectrum, and calculate the correction
  ;; factor
  corr_cov=total(corrCo)
  cf=corr_cov / uncorr_cov

  RETURN, cf

END


;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     EC_MOMENTUM
;; 
;; PURPOSE:
;; 
;;     This program calculates momentum flux and stability.
;; 
;; CALLING SEQUENCE:
;; 
;; 
;; 
;; INPUTS:
;; 
;;     Wind:         3-column array of unrotated wind [u,v,w].
;;     Ts:           Vector of sonic temperature (in deg C).
;;     Met_T:        Scalar float with mean air temperature measured by the
;;                   MET logger (C).
;;     Met_P:        Scalar float with mean pressure measured by the MET
;;                   logger (KpA).
;;     Met_RH:       Scalar float with mean relative humidity measured by
;;                   the MET logger (%).
;;     Avg_Period:   Scalar float with the averaging period under
;;                   consideration (minutes) (e.g. 30min flux runs, 20min
;;                   flux runs, 10min flux runs, etc.)
;;     Data_Freq:    Scalar float with the sampling frequency in Hz
;;                   (e.g. 10Hz).
;; 
;; KEYWORD PARAMETERS:
;; 
;;     CORR_MASSMAN: Set this keyword to run a spectral correction
;;                   following Massman this keyword must be used to set a
;;                   number of parameters for the Massman correction.  ALL
;;                   fields must be filled.  If they are not applicable to
;;                   the specific application, set them to the string value
;;                   'NAN' (for not applicable).  CORR_MASSMAN=[sample rate
;;                   (in seconds), sonic path length (in m)].
;;     OGIVE:        Set this keyword to produce an ogive plot for all
;;                   fluxes which are calculated using this routine.  this
;;                   keyword must be set to a string which indicates the
;;                   directory where the ogive should be placed.
;; 
;; OUTPUTS:
;; 
;;     A 5 element array, [0]=w/u covariance (corrected or uncorrected)
;;     [1]=correction factor for momentum flux [2]=friction velocity, u*
;;     (m/s) [3]=momentum flux (N/m2) [4]=Obukhov length (m).
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

FUNCTION EC_MOMENTUM, WIND, TS, MET_T, MET_RH, MET_P, AVG_PERIOD, $
                      DATA_FREQ, CORR_MASSMAN=CMASS, OGIVE=O_OUTPUT

  ;; Constants
  r=8.31451                ; j/mol/k universal gas constant
  mv=18.02                 ; g/mol molecular weight for water, Stull(1995)
  ma=28.96                 ; g/mol molecular weight for dry air, Stull(1995)
  ;; J/g/K specific heat for dry air at constant pressure Stull(1995)
  cpd=1004.67 / 1000.0

  ;; Calculate the mean horizontal wind velocities (unrotated)... this is
  ;; required in the spectral correction
  horwind=[mean(WIND[0, *], /nan), mean(WIND[1, *], /nan)]
  ;; do the wind rotations
  WINDrot=yawpitch(WIND[0, *], WIND[1, *], $
                   WIND[2, *], n_elements(WIND[0, *]))
  urot=WINDrot[0, *] & vrot=WINDrot[1, *] & wrot=WINDrot[2, *]

  ;; calculate w/u covariance
  cov_w_u=c_covariance(wrot, urot, 0)
  cf_wu=!VALUES.D_NAN          ; default in case no spectral correction run

  IF keyword_set(o_output) THEN BEGIN
     maxc=0.
     ogive,'wu', wrot, urot, avg_period, data_freq, maxc, o_output
  ENDIF

  ;; Run spectral correction to calculate correction factor --> sensor
  ;; geometry, specifications are currently hardcoded for the 2007/08
  ;; Amundsen cruises
  IF keyword_set(cmass) THEN BEGIN
     cf_wu=spec_massman(wrot[where(finite(wrot))], $
                        urot[where(finite(urot))], horwind, $
                        cmass[0], cmass[1], /MOMENTUM)
     cov_w_u=cov_w_u * cf_wu    ; calculate the corrected covariance
  ENDIF

  ;; calculate u*
  Ustar=sqrt(sqrt(cov_w_u ^ double(2))) ; as per eq'n 2.10b of Stull

  ;; calculate Tau.  Need to calculate rho_d and rho_v... we'll use the met
  ;; data for convenience... should be ok.  Sat'n vapour pressure (Pa)
  es_met=[6.112 * exp((17.67 * met_T) / (met_T + 243.5))] * double(100)
  ev_met=(met_RH / double(100)) * es_met ; vapour pressure (Pa)
  c_h2o_met=ev_met / (R * (met_T + 273.15)) ;mol 
  rho_v_met=c_h2o_met * mv                               ;g/m3
  rho_d_met=((met_P * double(1000) - ev_met) / $
             (R * (met_T + 273.15))) * ma 
  rho_met=rho_v_met + rho_d_met ;g/m3
  ;; mean specific humidity (g_h2o/g_moist air)
  mean_qh2o=rho_v_met / rho_met
  ;; mean specific heat capacity of moist air (J/g/K)
  mean_cp=cpd * (double(1) + 0.84 * mean_qh2o)

  Tau=(rho_met / double(1000)) * (cov_w_u) ; as per eq'n 2.10a of Stull

  ;; To get the stability parameter out, need to calculate sonic potential
  ;; temperature. This is as per eq'n 1.5.1c of Stull, which is a bit
  ;; strange... it uses the met Patm
  Ts_K=Ts + 273.15
  Pot_multiplier=findgen(n_elements(Ts_k))
  Pot_multiplier[*]=[(double(100) / met_P) ^ (0.286)]
  ;; Should probably use open path P instead of met P, but I doubt it's a
  ;; big deal
  Ts_Pot=Ts_K * Pot_multiplier
  cov_w_Ts_Pot=c_covariance(wrot, Ts_Pot, 0)

  ;; Run the spectral correction to calculate the correction factor for
  ;; this covariance
  IF keyword_set(cmass) THEN BEGIN
     cf_wT=spec_massman(wrot[where(finite(wrot))], $
                        Ts_Pot[where(finite(wrot))], horwind, $
                        cmass[0], cmass[1], $
                        GEOM=[double(0), double(0)], SCALARLINE=cmass[1])
     cov_w_Ts_Pot=cov_w_Ts_Pot * cf_wT
  ENDIF

  ;; calculate L
  grav=9.81 & vonk=0.4
  mean_Ts_Pot=mean(Ts_Pot, /nan)
  ;;  As per eq'n 5.7c of Stull
  L=(-mean_Ts_Pot * Ustar ^ 3) / (vonk * grav * cov_w_Ts_Pot)

  ;; We will also calculate the sonic temperature flux, which we will
  ;; output here just as an additional parameter
  cov_w_Ts=c_covariance(wrot, Ts_K, 0)
  IF keyword_set(cmass) THEN BEGIN
     cf_wTs=spec_massman(wrot[where(finite(wrot))], $
                         Ts_K[where(finite(wrot))], horwind, $
                         cmass[0], cmass[1], $
                         GEOM=[double(0), double(0)], SCALARLINE=cmass[1])
     cov_w_Ts=cov_w_Ts * cf_wTs
  ENDIF

  ;; and we'll calculate a modified true heat flux, using the met RH/T as
  ;; input for the sonic2airT program
  c_h2o_arr=findgen(n_elements(wrot)) & P_arr=c_h2o_arr
  c_h2o_arr[*]=c_h2o_met
  P_arr[*]=met_P * double(1000)
  ;; Calculate a 'pseudo' Tair using the low F RH data
  psTair=sonicT2airT(Ts, c_h2o_arr, P_arr)
  psTair=psTair + 273.15        ; Air temperature into K 
  cov_w_psTair=c_covariance(wrot, psTair, 0)
  IF keyword_set(cmass) THEN BEGIN
     cf_psTair=spec_massman(wrot[where(finite(wrot))], $
                            psTair[where(finite(wrot))], horwind, $
                            cmass[0], cmass[1], $
                            GEOM=[double(0), double(0)], $
                            SCALARLINE=cmass[1])
     cov_w_psTair=cov_w_psTair * cf_psTair
  ENDIF
  psH=rho_met * mean_cp * cov_w_psTair ;  Units: W/m2


  returnvec=[cov_w_u, cf_wu, Ustar, Tau, L, cov_w_Ts, cov_w_psTair, psH]
  ;;         0        1      2      3    4  5         6             7

  RETURN, returnvec

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; ec_momentum.pro ends here
