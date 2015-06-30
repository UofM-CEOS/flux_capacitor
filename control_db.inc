;; Author: Sebastian Luque
;; Created: 2013-08-27T17:52:26+0000
;; Last-Updated: 2015-06-30T19:48:09+0000
;;           By: Sebastian P. Luque
;; ------------------------------------------------------------------------
;;; Commentary: 
;;
;; This is a control file where variables for flux_capacitor_2013.pro are
;; defined.
;;; -----------------------------------------------------------------------
;;; Code:

;; Temporarily setting these here; they will be set from main GUI program.
root_dir=expand_path('~/Data/ArcticNet/2013/Tower')
templates_dir=expand_path('~/Scripts/Projects/Flux_Capacitor/Templates')

;;;_ : Flux (eddy-covariance)

ec_daily_rate=0.1               ; sampling rate (s)
ec_stamp='EC'
ec_period=1200                  ; duration of flux period (s)

;; Raw files
ec_idir=root_dir + path_sep() + 'EC'
;; Period files
ec_period_dir=ec_idir + path_sep() + 'FromDB'
;; Output directory for motion corrected files
ec_motcorr_dir=ec_period_dir + path_sep() + 'Motion_Corrected'
;; Output directory for footprint plots
ec_footprint_dir=ec_period_dir + path_sep() + 'Footprint'
;; Output file for fluxes
ec_fluxes_file=ec_idir + path_sep() + 'fluxes.csv'

;; See DB_FLUX docs for meaning of variables below
motpak_offset=[1.700, 0, 2.725]
sog_thr=0.25
lfreq_thr=0.02
hfreq_thr=1.0
xover_freq_thr=0.02

;;;_  . Templates

ec_std_template=templates_dir + path_sep() + 'ec_from_db.sav'



;;;_ * Emacs local variables
;; Local Variables:
;; mode: idlwave
;; allout-layout: (-2 + : 0)
;; End:
