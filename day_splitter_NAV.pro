;; $Id$
;;; day_splitter_NAV.pro --- split NAV data into daily files
;; Author: Bruce Johnson, Sebastian Luque
;; Created: 2013-08-29T03:29:07+0000
;; Last-Updated: 2013-09-20T16:00:43+0000
;;           By: Sebastian Luque
;; ------------------------------------------------------------------------
;;; Commentary: 
;;
;; This procedure needs to be improve by removing those loops all over the
;; place.  We also rely on the ITEMPLATE_SAV input argument with the path
;; to a save file that contains a template built via ascii_template().
;;
;; Example call:
;;
;; day_splitter_nav, '20111015', '20111025', $
;;                   expand_path('~/tmp/NAV/STD'), $
;;                   expand_path('~/tmp/NAV/Daily'), $
;;                   expand_path('nav_template.sav'), $
;;                   'headera, headerb, ...', 60, 'NAV', 8
;;
;; Comments from BJ below.
;; 
;; This program takes all of the .dat files in the std_nav directory, and
;; breaks them up into individual day files.  It will create a day file for
;; each day between the startdate and enddate, and fill any missing data
;; with 'NaN'
;;
;; NOTE: in it's present configuration, you CANNOT process more than 1 year
;; of data However, you can process across a calendar year (e.g. from
;; 20071001 to 20080701) It may slog a bit across a calendar year... try to
;; avoid this.
;; 
;; NOTE: for this to work, the input files should have the dates/times in
;; their first columns: for 1min or 1sec data: year,month,day,hour,min,sec
;; 
;; Input variables:
;; 
;;   startdate     - first day for generating daily files (YYYYMMDD) -->
;;                   User prompted to enter date
;;   
;;   enddate       - last day for generating daily files (YYYYMMDD) -->
;;                   User prompted to enter date
;;   
;;   idir          - input directory of files to be processed
;;   
;;   odir          - output directory for new daily files
;;   
;;   itemplate_sav - directory of .sav file of template being used NOTE:
;;                   the template must have been saved with the variable
;;                   name "mytemplate"
;;   
;;   header        - text string that defines the header
;;   
;;   samplerate    - The sample rate which data was collected in seconds
;;                   (e.g. 1, 60)
;;   
;;   stamp         - gives a prefix to the output files (e.g. if set to
;;                   pCO2, files will be pCO2_YYYY_JD_MMDD.dat)
;;   
;;   n_prefix      - gives the number of characters BEFORE the date of the
;;                   input file name (e.g. for pAMD2011_pCO2_071028.dat,
;;                   the number is 14)
;;  
;; For NAV data files:
;; 
;;  idir          = nav_std_dir -> '/ArcticNet2010/TowerData/NAV/std-NAV/'
;;                  file format "std_AMD2011_NAV_mmdd_hhmm.dat"
;;  
;;  odir          = daily_nav_dir ->
;;                  '/ArcticNet2010/TowerData/NAV/daily-NAV/' file format
;;                  "NAV_yyyy_JDxxx_mmdd.dat"
;;  
;;  itemplate_sav = nav_template ->
;;                  '/ArcticNet2010/Templates/nav_template.sav'
;;  
;;  header        = nav_header -> 'Year, Month, Day, Hour, Minute, Second,
;;                                 ProgVers, GPSDate, GPSTime, Latitude,
;;                                 Longitude, SOG(kts), COG(deg), Heading,
;;                                 Pitch, Roll, Accelx, Accely, Accelz'
;;                                      
;;  samplerate    = met_timing     --> 1
;;  stamp         = nav_stamp      --> 'NAV'
;;  n_prefix      = nav_prefix     --> 16
;;
;; ------------------------------------------------------------------------
;;; Code:

PRO DAY_SPLITTER_NAV, STARTDATE, ENDDATE, IDIR, ODIR, ITEMPLATE_SAV, $
                      HEADER, SAMPLERATE, STAMP, OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 8) THEN $
     message, 'Usage: DAY_SPLITTER_NAV, STARTDATE, ENDDATE, IDIR, ' + $
              'ODIR, ITEMPLATE_SAV, HEADER, SAMPLERATE, STAMP, [/OVERWRITE]'
  IF ((n_elements(startdate) EQ 0) OR (idir EQ '')) THEN $
     message, 'STARTDATE is undefined or is empty string'
  IF ((n_elements(enddate) EQ 0) OR (idir EQ '')) THEN $
     message, 'ENDDATE is undefined or is empty string'
  IF (long(enddate) LT long(startdate)) THEN $
     message, 'ENDDATE must be greater than or equal to STARTDATE'
  IF ((n_elements(idir) EQ 0) OR (idir EQ '')) THEN $
     message, 'IDIR is undefined or is empty string'
  IF ((n_elements(odir) EQ 0) OR (odir EQ '')) THEN $
     message, 'ODIR is undefined or is empty string'
  IF ((n_elements(itemplate_sav) EQ 0) OR (itemplate_sav EQ '')) THEN $
     message, 'ITEMPLATE_SAV is undefined or is empty string'
  IF (n_elements(header) EQ 0) THEN $
     message, 'HEADER is undefined'
  IF (n_elements(samplerate) EQ 0) THEN $
     message, 'SAMPLERATE is undefined'
  IF ((n_elements(stamp) EQ 0) OR (stamp EQ '')) THEN $
     message, 'STAMP is undefined or is empty string'

  idir_files=file_search(idir + path_sep() + '*.dat', $
                         count=nidir_files, /nosort)
  IF nidir_files LT 1 THEN BEGIN
     print, 'No input files found'
     RETURN
  ENDIF

  ;; Extract the year, month, and day from the startdate.
  beg_year=long(strmid(startdate, 0, 4))
  beg_mon=long(strmid(startdate, 4, 2))
  beg_day=long(strmid(startdate, 6, 2))
  ;; Extract the year, month, and day from the enddate
  end_year=long(strmid(enddate, 0, 4))
  end_mon=long(strmid(enddate, 4, 2))
  end_day=long(strmid(enddate, 6, 2))
  ;; Find the start and end day-of-year
  startDOY=calendar2doy(beg_year, beg_mon, beg_day)
  endDOY=calendar2doy(end_year, end_mon, end_day)
  ;; Find how many days we have, starting from 0
  beg_jd=julday(beg_mon, beg_day, beg_year, 0) ; start at 00:00:00
  end_jd=julday(end_mon, end_day, end_year, 0) ; start at 00:00:00
  n_days=long(end_jd - beg_jd)
  n_recs=86400 / samplerate     ; records per day
  restore, itemplate_sav
  n_fields=nav_template.FIELDCOUNT ; hard-coding this for now (MUST FIX)

  ;; Add 1 to the last day to set the breaks for value_locate()
  days=timegen(start=beg_jd, final=end_jd + 1, step_size=1, units='days')
  days=jul2timestamp(temporary(days))
  ;; Add up to 23:59:59 of the last requested day
  times=timegen(start=beg_jd, final=end_jd + (86399.0 / 86400), $
                step_size=samplerate, units='seconds')
  times=jul2timestamp(temporary(times))
  ;; ;; Index times in days vector.  Note that value_locate() does work with
  ;; ;; strings, which in this case is appropriate given jul2timestamp()
  ;; pos=value_locate(days, times)

  all_arr=fltarr(n_fields, size(times, /n_elements))
  all_arr[*, *]='NaN'
  all_arr[0, *]=strmid(times, 0, 4)
  all_arr[1, *]=strmid(times, 5, 2)
  all_arr[2, *]=strmid(times, 8, 2)
  all_arr[3, *]=strmid(times, 11, 2)
  all_arr[4, *]=strmid(times, 14, 2)
  all_arr[5, *]=strmid(times, 17, 2)

  ;; Build file stamps and check output files already there.  We could use
  ;; the vectorized strsplit(), instead of strmid() in IDL >= 8.0
  file_yyyydoy=strmid(days, 0, 4) + $
               calendar2doy(strmid(days, 0, 4), $
                            strmid(days, 5, 2), $
                            strmid(days, 8, 2))
  file_stamps=strcompress(stamp + '_' + file_yyyydoy + '.dat', $
                          /remove_all)
  odir_files=file_basename(file_search(odir + path_sep() + '*.dat', $
                                       count=nodir_files, /nosort))

  skipdays=intarr(size(days, /n_elements))
  FOR i=0L, n_days, 1L DO BEGIN
     is_ofile=where(strmatch(odir_files, file_stamps[i]), n_ofiles)
     IF n_ofiles GT 0 THEN BEGIN
        IF keyword_set(overwrite) THEN BEGIN
           print, 'Daily NAV file ' + file_stamps[i] + ' already exists.  ' + $
                  'Overwriting'
        ENDIF ELSE BEGIN
           print, 'Daily NAV file' + file_stamps[i] + ' already exists.  ' + $
                  'Not overwriting'
           skipdays[i]=1
           CONTINUE
        ENDELSE
     ENDIF
  ENDFOR
  ;; We find the rows for the days that had to be skipped, if requested.
  skip=where(skipdays GT 0, nskip) * n_recs
  IF nskip GT 0 THEN BEGIN
     file_stamps=file_stamps[where(~skipdays)]
     ;; We use rebin() to get an array where the rows are the skip values
     ;; repeated for as many columns as there are records in a day
     ;; (n_recs). We use matrix multiplication to obtain an array with the
     ;; same dimensions, containing the sequence 1-n_recs in each column,
     ;; and then we add these two arrays.  We use reform to turn the result
     ;; into a simple vector.
     skips=reform(transpose(rebin(skip, nskip, n_recs, /sample) + $
                            (replicate(1, nskip) # indgen(n_recs))), $
                  nskip * n_recs)
     idx=histogram(skips, min=0, max=(size(all_arr, /dimensions))[1] - 1)
     ;; Subset the array and other objects
     all_arr=all_arr[*, where(~idx)]
     times=times[where(~idx)]
  ENDIF

  ;; Read each file
  FOR k=0L, nidir_files - 1 DO BEGIN
     ifile=idir_files[k]
     print, 'Processing File: ' + ifile
     data=read_ascii(ifile, count=n_inputfile, template=nav_template)
     n_krecs=n_elements(data.(0))
     ;; Check each line and match against array
     FOR i=0L, n_krecs - 1, 1L DO BEGIN
        IF (data.(0)[i] LE 0) THEN BEGIN
           print, i, $
                  format='("Skipping record with unintelligible ' + $
                  'time stamp at: ", i)'
           CONTINUE
        ENDIF
        i_ts=string(data.(0)[i], format='(i4)') + '-' + $
             string(data.(1)[i], format='(i02)') + '-' + $
             string(data.(2)[i], format='(i02)') + ' ' + $
             string(data.(3)[i], format='(i02)') + ':' + $
             string(data.(4)[i], format='(i02)') + ':' + $
             string(data.(5)[i], format='(i02)')
        match=where(times EQ i_ts, mcount)
        IF mcount EQ 1 THEN BEGIN
           FOR col=6, n_fields - 1, 1L DO BEGIN $
              all_arr[col, match]=data.(col)[i]
           ENDFOR
        ENDIF
     ENDFOR
  ENDFOR

  ;; Write each daily array
  fmt_nfields=strtrim(n_fields - 1, 2)
  fmt_str='(' + fmt_nfields + '(a,","),a)'
  all_arr=strcompress(temporary(all_arr), /remove_all)
  FOR begi=0L, n_elements(times) - 1, n_recs DO BEGIN
     file_stamp=file_stamps[begi / n_recs]
     openw, lun, odir + path_sep() + file_stamp, /get_lun
     printf, lun, header
     printf, lun, all_arr[*, begi:(begi + n_recs - 1)], format=fmt_str
     free_lun, lun
  ENDFOR

END


;; TESTS




;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; day_splitter_NAV.pro ends here
