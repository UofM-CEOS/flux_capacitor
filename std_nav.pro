;; $Id$
;;; std_nav_AN11.pro --- Standardize NAV files
;; Author: Bruce Johnson, Sebastian Luque
;; Created: 2013-08-28T17:48:36+0000
;; Last-Updated: 2013-09-17T00:00:34+0000
;;           By: Sebastian Luque
;; ------------------------------------------------------------------------
;;; Commentary: 
;;
;; Check ChangeLog for history.
;;
;; Below is from Bruce Johnson's comments.
;; 
;; raw_nav_dir = input directory for raw files
;; out_columns = number of columns in output file
;; out_dir     = output directory for standardized files
;; header      = header information for output file
;; 
;; parameters for standardizing Nav data
;; -------------------------------------
;;
;; Input directory for raw files
;; raw_nav_dir = '/ArcticNet2011/TowerData/NAV/raw-NAV'
;; number of columns in the output file
;; out_columns = 23
;; output directory for plots
;; out_dir     = '/ArcticNet2011/TowerData/NAV/std-NAV'
;; header      = 'Year, Month, Day, Hour, Minute, Second, ProgVers, GPSDate,
;;                GPSTime, Latitude, Longitude, SOG, COG, 
;;                Heading, Roll, Pitch, Accelx, Accely, Accelz'
;;
;; ------------------------------------------------------------------------
;;; Code:

PRO STD_NAV, RAW_NAV_DIR, OUT_DIR, HEADER, OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 3) THEN $
     message, 'Usage: STD_NAV, RAW_NAV_DIR, HEADER'
  IF ((n_elements(raw_nav_dir) EQ 0) OR (raw_nav_dir EQ '')) THEN $
     message, 'RAW_NAV_DIR is undefined or is empty string'
  IF ((n_elements(out_dir) EQ 0) OR (out_dir EQ '')) THEN $
     message, 'OUT_DIR is undefined or is empty string'
  IF (n_elements(header) EQ 0) THEN $
     message, 'HEADER is undefined'

  cd, raw_nav_dir
  navlist=file_search('*.dat', count=nav_recs, /nosort)
  ;; Code currently works under these assumptions (silly state of affairs):
  in_cols=20
  out_cols=23
     
  ;; Loop through files in the raw-NAV directory
  FOR k=0, nav_recs - 1 DO BEGIN
     iname=strsplit(navlist[k], '.', /extract)
     ;; Get a name for the file, check if it already exists
     outname=strcompress(out_dir + path_sep() + iname[0] + '_std.' + $
                         iname[1], /remove_all)
     file_stamp=strcompress(iname[0] + '_std.' + iname[1], /remove_all)
     
     cd, out_dir
     out_list=file_search('*.dat', /nosort)
     matchfiles=where(file_stamp EQ out_list, matchfilecount)
     IF matchfilecount GT 0 THEN BEGIN
        IF keyword_set(overwrite) THEN BEGIN
           print, 'Standardized NAV file ' + file_stamp + $
                  ' already exists.  Overwriting'
        ENDIF ELSE BEGIN
        print, 'Standardized NAV file ' + file_stamp + $
               ' already exists.  Not overwriting'
        CONTINUE
     ENDELSE
     ENDIF

     ifile_nav=strcompress(raw_nav_dir + path_sep() + navlist[k], $
                           /remove_all)
     print, 'Processing file: ' + navlist[k]

     ;; Designate number of lines and columns for output file
     lines=file_lines(ifile_nav)
     work_arr=fltarr(out_cols, lines)
     work_arr[*, *]='NAN' 

     ;; Open input file as a string array
     openr, 1, ifile_nav
     ifile_str=strarr(1, lines) 
     readf, 1, ifile_str
     close, 1

     ;; Need to modify this time array so it works, should be (0-year,
     ;; 1-month, 2-day, 3-hour, 4-min, 5-sec)
     time_arr=fltarr(6, lines)
     time_arr[*, *]='NaN'
     
     ;; Loop through records (lines)
     FOR exd=0L, lines[0] - 1  DO BEGIN
        ;; Extract date/time info and fix NAN problem.
        splitstr=strsplit(ifile_str[exd], ",", /extract) ; split input
                                ; columns
        ;; Skip if wrong number of records
        IF (n_elements(splitstr) NE in_cols) THEN CONTINUE
        findnan=where(splitstr EQ '"NAN"', foundnan) 
        IF foundnan GT 0 THEN BEGIN
           splitstr[findnan]='NaN' ;look for "NAN" and replace with "NaN"
        ENDIF
        
        ;; Get the DOY and then the date and time
        calendar=doy2calendar(long(splitstr[1]), long(splitstr[2]))
        time_arr[0, exd]=long(strmid(calendar, 0, 4))
        time_arr[1, exd]=long(strmid(calendar, 4, 2))
        time_arr[2, exd]=long(strmid(calendar, 6, 2))

        hrmin=string(splitstr[3], format='(I04)')
        time_arr[3, exd]=long(strmid(hrmin, 0, 2))
        time_arr[4, exd]=long(strmid(hrmin, 2, 2))
        time_arr[5, exd]=float(splitstr[4])
        
        ;; Reconstruct the string, with the fixed NaN's and new timestamp
        ;; Place year, month, day, hour, minutes
        work_arr[0:4, exd]=strcompress(long(time_arr[0:4, exd]), $
                                       /remove_all)
        ;; Place seconds
        work_arr[5, exd]=strcompress(time_arr[5, exd], /remove_all)
        work_arr[6, exd]=(splitstr[5]) ;place program version
        
        ;; Check length of GPSDate, if not 8 characters, substitute NaN
        IF strlen(splitstr[6]) EQ 8 THEN BEGIN
           work_arr[7, exd]=long(strmid(splitstr[6], 5, 2)) + 2000 ;place GPSYear
           work_arr[8, exd]=long(strmid(splitstr[6], 3, 2))  ;place GPSMonth
           work_arr[9, exd]=long(strmid(splitstr[6], 1, 2))  ;place GPSDay
        ENDIF ELSE BEGIN
           work_arr[7:9, exd]='NaN'
        ENDELSE
        
        ;; Check length of GPSTime, if not 8 characters, substitute NaN
        IF strlen(splitstr[7]) EQ 8 THEN BEGIN
           work_arr[10, exd]=long(strmid(splitstr[7], 1, 2)) ;place GPSHour
           work_arr[11, exd]=long(strmid(splitstr[7], 3, 2)) ;place GPSMinute
           work_arr[12, exd]=long(strmid(splitstr[7], 5, 2)) ;place GPSSecond
        ENDIF ELSE BEGIN
           work_arr[10:12, exd]='NaN'
        ENDELSE
        
        work_arr[13:14, exd]=splitstr[8:9] ;place Latitude and Longitude
        
        ;; Check length of SOG, if not 7 characters keep value as 'NaN'
        IF strlen(splitstr[10]) EQ 7 THEN BEGIN
           work_arr[15, exd]=float(strmid(splitstr[10], 1, 5)) ;place SOG
        ENDIF ELSE BEGIN
           work_arr[15, exd]='NaN'
        ENDELSE
        
        ;; Check length of COG, if not 7 characters keep value as 'NaN'
        IF strlen(splitstr[11]) EQ 7 THEN BEGIN
           work_arr[16, exd]=float(strmid(splitstr[11], 1, 5)) ;place COG
        ENDIF ELSE BEGIN
           work_arr[16, exd]='NaN'
        ENDELSE
        
        ;; Check for bad data line from digital compass, shows up as 'NaN'
        ;; in last field
        IF splitstr[19] EQ 'NaN' THEN BEGIN
           work_arr[17:22, exd]='NaN'
        ENDIF ELSE BEGIN
           work_arr[17:22, exd]= splitstr[14:19] ;place Heading, Pitch, Roll, Accel(3)
        ENDELSE

     ENDFOR

     ;; Filtering for out of range values
     ;; ---------------------------------
     ;; latitude  -- < 0 or > 90    - northern hemisphere values
     ;; longitude -- < -180 or > 0  - western hemisphere values
     ;; SOG       -- < 0 or > 20    - assuming ship can't be > 20 kts
     ;; COG       -- < 0 or > 360  
     ;; heading   -- < 0 or > 360  
     ;; pitch     -- < -3 or > 3    - check to ensure these are proper bounds
     ;; roll      -- < -3 or > 3

     badlat=where((work_arr[13, *] LT 0.0) OR $
                  (work_arr[13, *] GT 90.0), badlatcount)
     badlong=where((work_arr[14, *] LT -180.0) OR $
                   (work_arr[14, *] GT 0.0), badlongcount)
     badsog=where((work_arr[15, *] LT 0.0) OR $
                  (work_arr[15, *] GT 20.0), badsogcount)
     badcog=where((work_arr[16, *] LT 0.0) OR $
                  (work_arr[16, *] GT 360.0), badcogcount)
     badhead=where((work_arr[17, *] LT 0.0) OR $
                   (work_arr[17, *] GT 360.0), badheadcount)
     badroll=where((work_arr[18, *] LT -3.0) OR $
                   (work_arr[18, *] GT 3.0), badrollcount)
     badpitch=where((work_arr[19, *] LT -3.0) OR $
                    (work_arr[19, *] GT 3.0), badpitchcount)
     
     IF badlatcount GT 0 THEN work_arr[13, badlat]='NaN'
     IF badlongcount GT 0 THEN work_arr[14, badlong]='NaN'
     IF badsogcount GT 0 THEN work_arr[15, badsog]='NaN'
     IF badcogcount GT 0 THEN work_arr[16, badcog]='NaN'
     IF badheadcount GT 0 THEN work_arr[17, badhead]='NaN'
     IF badrollcount GT 0 THEN work_arr[18, badroll]='NaN'
     IF badpitchcount GT 0 THEN work_arr[19, badpitch]='NaN'

     ;;----------------Now print it out-----------------------------------
     ;;make the array a string for easy printing
     formatnum=string(out_cols - 1, format='(i2)')
     formatstring='(' + formatnum + '(a,", "),a)'
     openw, unit1, outname, /get_lun
     printf, unit1, header
     printf, unit1, strcompress(work_arr, /remove_all), format=formatstring
     close, /all

  ENDFOR

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; std_nav_AN11.pro ends here
