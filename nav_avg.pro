;; $Id$
;;; nav_avg.pro --- calculate averages in NAV files
;; Author: Bruce Johnson, Sebastian Luque
;; Created: 2013-09-17T14:59:07+0000
;; Last-Updated: 2013-09-24T20:34:42+0000
;;           By: Sebastian Luque
;; ------------------------------------------------------------------------
;;; Commentary: 
;; 
;; Example call:
;;
;; nav_avg, expand_path('~/tmp/NAV/Daily'), expand_path('~/tmp/NAV/1min'), $
;;          1, 60, 15, 16, 17, 'NAV', expand_path('nav_template.sav'), $
;;          'headera, headerb, ...'
;;
;; Original comments from BJ below.
;; 
;; This program takes daily 1 second files and makes 1 minute averages
;; NOTE: will need a new one for higher freq data (i.e. 10Hz), but should
;; be easy to modify.
;; 
;; THIS PROGRAM HAS BEEN MODIFIED FOR EXCLUSIVE USE WITH NAV DATA FROM
;; CR1000 LOGGER -- in addition to calculating the averages of SOG and COG,
;; it also calculates the StdDev.  For SOG, this is done in a straight
;; stdev calc, for COG, it is done using the Yamartino method
;; 
;; VARIABLE
;; idir       - input directory
;; odir       - output directory
;; temp_sav   - ASCII template save file
;; n_fields   - number of fields in the data being processed (i.e. input
;;              file)
;; bear_field - set the number of a field which has bearing (compass) data
;;              (set to 0 for no data)
;; vec_field  - set the number of a field which has vector magnitude data
;;              associated with the bearing field (set to 0 for no data)
;; stamp      - stamp you wish to add
;; header     - header for the file
;; 
;; for NAV data
;; -------------

;; idir       = daily_nav_dir   --> file format "NAV_YYYYDOY.dat"
;; odir       = min_nav_dir     --> file format "NAV_YYYYDOY_min.dat
;; n_fields   = 23
;; vec_field  = 15
;; bear_field = 16
;; head_field = 17
;; stamp      = 'NAV'
;; header     = 'Year, Month, Day, Hour, Minute, Second, ProgVers,
;;               Latitude, Longitude, SOG(kts), COG(deg), Heading, Pitch,
;;               Roll, Accelx, Accely, Accelz, SOG_stdev, COG_stdev,
;;               Heading_stdev'
;;
;; ------------------------------------------------------------------------
;;; Code:

PRO nav_avg, IDIR, ODIR, ISAMPLE_RATE, OSAMPLE_RATE, BMAG_FIELD, $
             BEAR_FIELD, HEAD_FIELD, STAMP, ITEMPLATE_SAV, HEADER

  ;; Check parameters
  IF (n_params() NE 10) THEN $
     message, 'Usage: ONE_MIN_NAV, IDIR, ODIR, ISAMPLE_RATE, ' + $
              'OSAMPLE_RATE, BMAG_FIELD, BEAR_FIELD, HEAD_FIELD, ' + $
              'STAMP, ITEMPLATE_SAV, HEADER'
  IF ((n_elements(idir) EQ 0) OR (idir EQ '')) THEN $
     message, 'IDIR is undefined or is empty string'
  IF ((n_elements(odir) EQ 0) OR (odir EQ '')) THEN $
     message, 'ODIR is undefined or is empty string'
  IF ((n_elements(bmag_field) EQ 0) OR (bmag_field EQ '')) THEN $
     message, 'BMAG_FIELD is undefined or is empty string'
  IF ((n_elements(isample_rate) EQ 0) OR (isample_rate EQ '')) THEN $
     message, 'ISAMPLE_RATE is undefined or is empty string'
  IF ((n_elements(osample_rate) EQ 0) OR (osample_rate EQ '')) THEN $
     message, 'OSAMPLE_RATE is undefined or is empty string'
  IF ((osample_rate MOD isample_rate) NE 0) THEN $
     message, 'ISAMPLE_RATE must be an integer divisor of OSAMPLE_RATE'
  IF ((n_elements(bear_field) EQ 0) OR (bear_field EQ '')) THEN $
     message, 'BEAR_FIELD is undefined or is empty string'
  IF ((n_elements(head_field) EQ 0) OR (head_field EQ '')) THEN $
     message, 'HEAD_FIELD is undefined or is empty string'
  IF ((n_elements(stamp) EQ 0) OR (stamp EQ '')) THEN $
     message, 'STAMP is undefined or is empty string'
  IF ((n_elements(itemplate_sav) EQ 0) OR (itemplate_sav EQ '')) THEN $
     message, 'ITEMPLATE_SAV is undefined or is empty string'
  IF (n_elements(header) EQ 0) THEN $
     message, 'HEADER is undefined'

  idir_files=file_search(idir + path_sep() + '*.dat', count=nidir_files, $
                         /nosort)
  IF nidir_files LT 1 THEN BEGIN
     print, 'No input files found'
     RETURN
  ENDIF

  restore, itemplate_sav
  n_fields=nav_template.FIELDCOUNT
  ang_cols=[bmag_field, bear_field, head_field]
  avg_cols=(cgSetDifference(indgen(n_fields - 13) + 13, ang_cols))
  n_avg_cols=size(avg_cols, /n_elements)
  ncols_osample=osample_rate / isample_rate
  nrows_osample=86400 / osample_rate
  header=strsplit(header, ', ', /extract)
  header_out=strarr(n_fields - 3)
  header_out[0:6]=header[0:6]
  header_out[7:(n_fields - 3 - 4)]=header[13:(n_fields - 1)]
  stdevstr=['SOG_stdev', 'COG_stdev', 'Heading_stdev']
  header_out[(n_fields - 6):(n_fields - 4)]=stdevstr
  header_out=strjoin(header_out, ', ')
  FOR k=0, nidir_files - 1 DO BEGIN
     ifile=idir_files[k]
     print, 'Producing 1-min average for file: ' + ifile
     data=read_ascii(ifile, count=n_inputfile, template=nav_template)
     beg_jd=julday(data.FIELD02[0], data.FIELD03[0], data.FIELD01[0], 0)
     mins=timegen(start=beg_jd, final=beg_jd + (86399.0 / 86400), $
                  step_size=osample_rate, units='seconds')
     mins=jul2timestamp(temporary(mins))
     ;; Make a working array, populate with 'NaN'
     all_arr=fltarr(n_fields - 3, size(mins, /n_elements))
     all_arr[*, *]='NaN'
     all_arr[0, *]=strmid(mins, 0, 4)
     all_arr[1, *]=strmid(mins, 5, 2)
     all_arr[2, *]=strmid(mins, 8, 2)
     all_arr[3, *]=strmid(mins, 11, 2)
     all_arr[4, *]=strmid(mins, 14, 2)
     all_arr[5, *]=strmid(mins, 17, 2)
     all_arr[6,*]=data.FIELD07[0] ; place Program Version
     FOREACH col, avg_cols DO BEGIN
        col_data=reform(data.(col), ncols_osample, nrows_osample)
        ;; This will throw "floating point illegal operand" arithmetic
        ;; error warning when all elements in each row are NaN, which
        ;; should should just be ignored.  In these cases, the result is
        ;; -NaN.
        avgs=mean(col_data, dimension=1, /nan)
        all_arr[col - 6, *]=avgs
     ENDFOREACH
     ;; Means for SOG, COG, and heading
     cog=data.(bear_field)
     brg2d=reform(cog, ncols_osample, nrows_osample)
     sog=data.(bmag_field)
     bmag2d=reform(sog, ncols_osample, nrows_osample)
     head=data.(head_field)
     head2d=reform(head, ncols_osample, nrows_osample)
     brgmag_mean=bearing_avg(brg2d, bmag2d, dimension=1)
     all_arr[bmag_field - 6, *]=brgmag_mean[*, 1]
     all_arr[bear_field - 6, *]=brgmag_mean[*, 0]
     head_mean=bearing_avg(head2d, 1, dimension=1)
     all_arr[head_field - 6, *]=head_mean[*, 0]
     ;; Standard deviations
     FOR i=0, nrows_osample - 1 DO BEGIN
        sog_ok=where(finite(bmag2d[*, i]) GT 0, n_sog_ok, complement=sog_bad)
        cog_ok=where(finite(brg2d[*, i]) GT 0, n_cog_ok, complement=cog_bad)
        head_ok=where(finite(head2d[*, i]) GT 0, n_head_ok, $
                      complement=head_bad)
        IF n_sog_ok GT 10 THEN $
           all_arr[17, i]=stddev(bmag2d[sog_ok, i])
        IF n_cog_ok GT 10 THEN $
           all_arr[18, i]=stddev_yamartino(brg2d[cog_ok, i])
        IF n_head_ok GT 10 THEN $
           all_arr[19, i]=stdev_yamartino(head2d[head_ok, i])
     ENDFOR
     
     file_stamp=file_basename(ifile, '.dat') + '_' + $
                strtrim(osample_rate, 2) + 's.dat'
     fmt_nfields=strtrim(n_fields - 4)
     fmt_str='(' + fmt_nfields + '(a,","),a)'
     all_arr=strcompress(temporary(all_arr), /remove_all)

     print, 'Writing file: ' + file_stamp
     openw, lun, odir + path_sep() + file_stamp, /get_lun
     printf, lun, header_out
     printf, lun, all_arr, format=fmt_str
     free_lun, lun
  ENDFOR
        
END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; nav_avg.pro ends here
