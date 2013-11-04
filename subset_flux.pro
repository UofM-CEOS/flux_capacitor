;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-11-03T18:49:19+0000
;; Last-Updated: 2013-11-04T04:07:15+0000
;;           By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     SUBSET_FLUX
;; 
;; PURPOSE:
;; 
;;     This procedure subsets the (standardized) daily EC flux files,
;;     extracting the periods with suitable data, based on the diagnostic
;;     flag from MET data.  It also subsets the corresponding RMC and GYRO
;;     data for these periods.
;; 
;; CALLING SEQUENCE:
;; 
;;     SUBSET_FLUX, Idir, Odir, Itemplate_Sav, Time_Beg_Idx, Isample_Rate, $
;;                  DIAG_DIR, Diag_Itemplate_Save, Diag_Time_Idx, Diag_Idx, $
;;                  Flux_Period
;; 
;; INPUTS:
;; 
;;     Idir:                  Input directory (no trailing separator) for
;;                            (standardized) daily EC flux files.
;;     Odir:                  Output directory (no trailing separator).
;;     Itemplate_Sav:         Ascii template to read input files.
;;     Time_Beg_Idx:          Index (in template) where time is.
;;     Isample_Rate:          Scalar indicating the frequency (s) with
;;                            which input data were sampled.
;;     Diag_Dir:              Directory (no trailing separator) containing
;;                            processed MET files, averaged across flux
;;                            periods, and with a diagnostic field
;;                            indicating whether the period is suitable for
;;                            flux analyses (flag=0).
;;     Diag_Itemplate_Sav:    Ascii template to read Diag_Dir files.
;;     Diag_Time_Idx:         Index (in template) where time is in files in
;;                            Diag_Dir.
;;     Diag_Idx:              Index (in template) where the diagnostic
;;                            field is in files in Diag_Dir.
;;     Flux_Period:           Scalar indicating the duration (s) of flux
;;                            study periods.
;;     Rmc_Dir:               Directory (no trailing separator) containing
;;                            standardized daily RMC files.
;;     Rmc_Itemplate_Sav:     Ascii template to read Rmc_Dir files.
;;     Rmc_Time_Idx:          Index (in template) where time is in files in
;;                            Rmc_Dir.
;;     Gyro_Dir:              Directory (no trailing separator) containing
;;                            standardized daily GYRO files.
;;     Gyro_Itemplate_Sav:    Ascii template to read Gyro_Dir files.
;;     Gyro_Time_Idx:         Index (in template) where time is in files in
;;                            Gyro_Dir.
;; 
;; KEYWORD PARAMETERS:
;; 
;; 
;; 
;; SIDE EFFECTS:
;; 
;;     Files are written in Odir.
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

PRO SUBSET_FLUX, IDIR, ODIR, ITEMPLATE_SAV, TIME_BEG_IDX, ISAMPLE_RATE, $
                 DIAG_DIR, DIAG_ITEMPLATE_SAV, DIAG_TIME_IDX, DIAG_IDX, $
                 FLUX_PERIOD, RMC_DIR, RMC_ITEMPLATE_SAV, RMC_TIME_IDX, $
                 GYRO_DIR, GYRO_ITEMPLATE_SAV, GYRO_TIME_IDX, /OVERWRITE

  ;; Check parameters
  IF (n_params() NE 16) THEN $
     message, 'Usage: SUBSET_FLUX, IDIR, ODIR, ITEMPLATE_SAV, ' + $
              'TIME_BEG_IDX, ISAMPLE_RATE, DIAG_DIR, DIAG_ITEMPLATE_SAV, ' + $
              'DIAG_TIME_IDX, DIAG_IDX, FLUX_PERIOD, RMC_DIR, ' + $
              'RMC_ITEMPLATE_SAV, RMC_TIME_IDX, GYRO_DIR, ' + $
              'GYRO_ITEMPLATE_SAV, GYRO_TIME_IDX'
  IF ((n_elements(idir) EQ 0) OR (idir EQ '')) THEN $
     message, 'IDIR is undefined or is empty string'
  IF ((n_elements(odir) EQ 0) OR (odir EQ '')) THEN $
     message, 'ODIR is undefined or is empty string'
  IF ((n_elements(itemplate_sav) EQ 0) OR (itemplate_sav EQ '')) THEN $
         message, 'ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(time_beg_idx) NE 1) OR (time_beg_idx LT 0)) THEN $
     message, 'TIME_BEG_IDX must be a scalar >= zero'
  IF ((n_elements(isample_rate) NE 1) OR (isample_rate LT 0)) THEN $
     message, 'ISAMPLE_RATE must be a scalar >= zero'
  IF ((n_elements(diag_dir) EQ 0) OR (diag_dir EQ '')) THEN $
     message, 'DIAG_DIR is undefined or is empty string'
  IF ((n_elements(diag_itemplate_sav) EQ 0) OR $
      (diag_itemplate_sav EQ '')) THEN $
         message, 'DIAG_ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(diag_time_idx) NE 1) OR (diag_time_idx LT 0)) THEN $
     message, 'DIAG_TIME_IDX must be a scalar >= zero'
  IF ((n_elements(diag_idx) NE 1) OR (diag_idx LT 0)) THEN $
     message, 'DIAG_IDX must be a scalar >= zero'
  IF ((n_elements(flux_period) NE 1) OR (flux_period LT 0)) THEN $
     message, 'FLUX_PERIOD must be a scalar >= zero'
  IF ((n_elements(rmc_dir) EQ 0) OR (rmc_dir EQ '')) THEN $
     message, 'RMC_DIR is undefined or is empty string'
  IF ((n_elements(rmc_itemplate_sav) EQ 0) OR $
      (rmc_itemplate_sav EQ '')) THEN $
         message, 'RMC_ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(rmc_time_idx) NE 1) OR (rmc_time_idx LT 0)) THEN $
     message, 'RMC_TIME_IDX must be a scalar >= zero'
  IF ((n_elements(gyro_dir) EQ 0) OR (gyro_dir EQ '')) THEN $
     message, 'GYRO_DIR is undefined or is empty string'
  IF ((n_elements(gyro_itemplate_sav) EQ 0) OR $
      (gyro_itemplate_sav EQ '')) THEN $
         message, 'GYRO_ITEMPLATE_SAV is undefined or is empty string'
  IF ((n_elements(gyro_time_idx) NE 1) OR (gyro_time_idx LT 0)) THEN $
     message, 'GYRO_TIME_IDX must be a scalar >= zero'
  idir_files=file_search(idir + path_sep() + '*', count=nidir_files, $
                         /nosort, /fold_case, /test_regular)
  diag_files=file_search(diag_dir + path_sep() + '*', count=ndiag_files, $
                         /nosort, /fold_case, /test_regular)
  rmc_files=file_search(rmc_dir + path_sep() + '*', count=nrmc_files, $
                        /nosort, /fold_case, /test_regular)
  gyro_files=file_search(gyro_dir + path_sep() + '*', count=ngyro_files, $
                         /nosort, /fold_case, /test_regular)
  IF (nidir_files EQ 0) OR (ndiag_files EQ 0) OR $
     (nrmc_files EQ 0) OR (ngyro_files EQ 0) THEN $
        message, 'No files in at least one of IDIR, DIAG_DIR, ' + $
                 'RMC_DIR, or GYRO_DIR'


END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; subset_flux.pro ends here
