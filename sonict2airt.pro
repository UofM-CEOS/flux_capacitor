;; $Id$
;; Author: Brent Else, Sebastian Luque
;; Created: 2013-11-15T20:12:42+0000
;; Last-Updated: 2013-11-26T23:11:33+0000
;;	     By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;;
;;     SONICT2AIRT
;;
;; PURPOSE:
;;
;;     This little function tries to bootstrap the calculation procedure
;;     for conversion of sonic temperature (Tair_c in degC) from virtual to
;;     absolute (Tair_a).  It's bootstrapping because we use an open path
;;     sensor to get the chi and chi already needs the absolute temperature.
;;
;; CALLING SEQUENCE:
;;
;;     airt=sonict2airt(Tair_C, C_H2O, P_barometric)
;;
;; INPUTS:
;;
;;     Tair_C:	     Sonic virtual temperature; degC
;;     C_H2O:	     Water vapour concentration; mol/m3
;;     P_barometric: Barometric pressure; Pa
;;
;; OUTPUTS:
;;
;;     Tair_a = actual air temperature; degC.
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
;; HISTORY:
;;
;;     May 26, 2003 - Fixed chi calculation which was mistakenly assuming
;;		      H_density in mmol/m3 was equiv. to Ma Note that this
;;		      meant that chi_calc = chi_actual./0.622 at all times.
;;		      And fixed Tair calculation to use T in Kelvins rather
;;		      than deg C.
;;
;;     Sept 23, 2003 - Fixed Tair calculation which had chi (g h2o /g air)
;;		       erroneously divided by 1000.  I've now moved to
;;		       Kai's formulas with mole fraction rather than
;;		       specific humidity.
;;
;;     Sept 2008 (BE) - Modified input water vapour term to be mol/m3 (to
;;			coincide with the AMD_EC program) - Made the
;;			bootstrapping procedure a loop - easy to add more
;;			iterations of necessary.
;;
;;     References - Stull, R.B. (1995) Meteorology Today for Scientists and
;;		    Engineers, West Publ. Co.  Stull, R.B. (1988) Boundary
;;		    Layer Meteorology, Kluwer Academic Publ.  Wallace and
;;		    Hobbs(1977) Atmospheric Science - An Introductory
;;		    Survey, Academic Press Inc.
;;
;;- -----------------------------------------------------------------------
;;; Code:

FUNCTION SONICT2AIRT, TAIR_C, C_H2O, P_BAROMETRIC

  ZeroK=273.15	   ; zero deg. C in deg K.
  R=8.31451	   ; J/mol/K universal gas constant
  Tair_K=Tair_C + ZeroK 	  ; convert to Kelvins
  FOR x=0, 5 DO BEGIN
     ;; (mol h2o/m3 air * J/mol/K * K ) / Pa=mol h2o/mol wet air
     chi=c_h2o * [(R * Tair_K) / P_barometric]
     Tair_a=Tair_K / (1 + 0.32 * chi)
     Tair_K=Tair_a
  ENDFOR
  ;; convert back to degree C
  Tair_a=Tair_a - ZeroK

  RETURN, Tair_a

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; sonict2airt.pro ends here
