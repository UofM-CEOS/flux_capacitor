;; $Id: $
;;; day_splitter_MET.pro --- split MET data into daily files
;; Author: Sebastian P. Luque
;; Created: 2013-09-20T03:54:03+0000
;; Last-Updated: 2013-09-20T12:30:47+0000
;;           By: Sebastian P. Luque
;; ------------------------------------------------------------------------
;;; Commentary: 
;; 
;; Original comment from BJ (likely) below.
;; 
;; This program takes all of the .dat files in an input directory, and
;; breaks them up into individual day files.  It will create a day file for
;; each day input between startdate-enddate and fill any missing data with
;; 'NaN'
;; 
;; NOTE: in it's present configuration, you CANNOT process more than 1 year
;; of data However, you can process across a calendar year (e.g. from
;; 20071001 to 20080701) It may slog a bit across a calendar year... try to
;; avoid this.
;; 
;; NOTE: this program might be a bit inefficient when it comes to
;; processing large flux files
;; 
;; NOTE: for this to work, the input files should have the dates/times in
;; their first columns: For 1min or 1sec data: year,month,day,hour,min,sec
;; For high frequency data: year,month,day,hour,min,sec,tenth of sec
;; 
;; NOTE: this program has NOT been thoroughly tested for anything except
;; 60s data.
;; 
;; NOTE: you need to double check that the time stamp on the files is
;; correct (i.e. the proper day (UTC) it was downloaded)
;; 
;; Input variables:
;; 
;;   startdate     - first day for generating daily files (YYYYMMDD)
;;   enddate       - last day for generating daily files (YYYYMMDD)
;;   idir          - input directory of files to be processed
;;   odir          - output directory for new daily files
;;   itemplate_sav - directory of .sav file of template being used NOTE:
;;                   the template must have been saved with the variable
;;                   name "mytemplate"
;;   header        - text string that defines the header
;;   samplerate    - The sample rate which data was collected in seconds
;;                   (e.g. 1, 60)
;;   nfields       - The number of fields that are in the data file to be
;;                   created
;;   stamp         - gives a prefix to the output files (e.g. if set to
;;                   pCO2, files will be pCO2_YYYY_JD_MMDD.dat)
;;   n_prefix      - gives the number of characters BEFORE the date of the
;;                   input file name (e.g. for pAMD_pCO2_071028.dat, the
;;                   number is 10)
;;  
;; For MET data files:
;; 
;;  idir           = met_idir --> 'raw-MET/' -- file format
;;                                              "AMD_MET_mmdd_hhmm.dat"
;;  odir           = met_dailydir --> 'daily-MET/' -- file format
;;                                                    "MET_yyyy_JDxxx_mmdd.dat"
;;  itemplate_sav  = met_template --> 'met_template.sav'
;;  header         = met_header --> 'Year, Month, Day, Hour, Minute,
;;                                      Sec, Prog Version, Battery(V),
;;                                      Panel_T(C), Pressure(kPa),
;;                                      Air_T(C),' +$ 'RH(%), Surface_T(C),
;;                                      Raw Wind Speed(m/s), Raw Wind
;;                                      Direction(deg), Wind StdDev,
;;                                      PAR(umol/m2/s), Pitch(deg),' +$
;;                                      'Roll (deg), Battery StdDev,
;;                                      Panel_T StdDev, Pressure StDev,
;;                                      Air_T StDev, RH StDev, Surface_T
;;                                      StDev, PAR StDev'
;;  samplerate     = met_timing --> 60
;;  nfields        = met_outcolumns --> 26
;;  stamp          = met_stamp --> 'MET'
;;  n_prefix       = met_prefix     --> 8
;; 
;; For RAD data files:
;; 
;;  idir           = raw_rad_dir --> '/ArcticNet2011/RadData/raw-RAD'
;;  odir           = daily_rad_dir --> '/ArcticNet2011/RadData/daily-RAD/'
;;  itemplate_sav  = Rad_template --> '/ArcticNet2011/Templates/RAD_template.sav'
;;  header         = RAD_head --> 'Year, Month, Day, Hour, Minute, Sec,
;;                                      Battery(V), Panel_T(C), Batt_stdev,
;;                                      PanelT_stdev, Kdown(W/m2),
;;                                      Thermopile,' +$
;;                                      'Tcase(K),Tdome_avg, LWin(W/m2),
;;                                      PAR(umol/m2/s), Tuv(C), UVb(W/m2),
;;                                      UVa(W/m2), UVbroad(W/m2),
;;                                      Kdown_stdev,' +$ 'Thermopile_stdev,
;;                                      Tcase_stdev, Tdome_stdev,
;;                                      LWin_stdev, PAR_stdev, Tuv_stdev,
;;                                      Uvb_stdev, Uva_stdev,
;;                                      UVbroad_stdev'
;;  samplerate     = rad_timing     --> 60
;;  nfields        = rad_outcolumns --> 30
;;  stamp          = rad_stamp      --> 'RAD'
;;  n_prefix       = rad_prefix     --> 8
;; 
;; ------------------------------------------------------------------------
;;; Code:

PRO DAY_SPLITTER_MET, STARTDATE, ENDDATE, IDIR, ODIR, ITEMPLATE_SAV, $
                      HEADER, SAMPLERATE, STAMP, N_PREFIX, $
                      OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 9) THEN $
     message, 'Usage: DAY_SPLITTER_NAV, STARTDATE, ENDDATE, IDIR, ' + $
              'ODIR, ITEMPLATE_SAV, HEADER, SAMPLERATE, STAMP, N_PREFIX'
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
  IF (n_elements(n_prefix) EQ 0) THEN $
     message, 'N_PREFIX is undefined'

  idir_files=file_search(idir + path_sep() + '*.dat', $
                         count=nidir_files, /nosort)
  IF nidir_files LT 1 THEN BEGIN
     print, 'No input files found'
     RETURN
  ENDIF

;Extract the year, month, and day from the input dates
s_y = LONG(STRMID(startdate,0,4))
s_m = LONG(STRMID(startdate,4,2))
s_d = LONG(STRMID(startdate,6,2))
e_y = LONG(STRMID(enddate,0,4))
e_m = LONG(STRMID(enddate,4,2))
e_d = LONG(STRMID(enddate,6,2))

s_y=s_y[0]
s_m=s_m[0]
s_d=s_d[0]
e_y=e_y[0]
e_m=e_m[0]
e_d=e_d[0]

;Calculate the Julian Days under consideration
startJD = calendar_to_jd(s_y,s_m,s_d)
endJD   = calendar_to_jd(e_y,e_m,e_d)

;check if start/end years are leap years
s_y_leap = 0
IF ((s_y MOD 4) EQ 0) THEN BEGIN
  s_y_leap = 1
ENDIF

e_y_leap = 0
IF ((e_y MOD 4) EQ 0) THEN BEGIN
  e_y_leap = 1
ENDIF

;Now calculate how big of a loop you need
IF (s_y EQ e_y) THEN BEGIN
  n_days = endJD - startJD
ENDIF ELSE BEGIN
IF s_y_leap EQ 0 THEN BEGIN
 s_days = 365 - startJD
ENDIF ELSE BEGIN
  s_days = 366 - startJD
ENDELSE
  n_days = s_days + endJD
ENDELSE

n_days=n_days[0]
FOR loopdays = 0, n_days DO BEGIN

  ;Calculate what Julian Day and year you're working with
  current_JD   = startJD + loopdays
  current_year = s_y
  IF s_y_leap EQ 0 THEN BEGIN
    IF (current_JD GT 365) THEN BEGIN
      current_JD   = current_JD - 365
      current_year = s_y + 1
   ENDIF
 ENDIF ELSE BEGIN
    IF (current_JD GT 366) THEN BEGIN
      current_JD   = current_JD - 366
      current_year = s_y + 1
   ENDIF
 ENDELSE

  ;now calculate what calendar day we're working with AN09 --> don't need, re: JD is in the file
  current_date = jd_to_calendar(current_year,current_JD)
  current_m    = long(strmid(current_date,4,2))
  current_day  = long(strmid(current_date,6,2))

  print, 'Processing Day:  ' + current_date

  ;check to see if this file has already been created
  CD, odir
  out_list=FILE_SEARCH('*.dat',COUNT=nouts,/NOSORT)
  IF current_JD LT 100 THEN BEGIN
    JD_stamp = strcompress(current_JD)
    JD_stamp = '0'+JD_stamp
    IF current_JD LT 10 THEN BEGIN
      JD_stamp = strcompress(current_JD)
      JD_stamp = '00'+JD_stamp
   ENDIF
 ENDIF ELSE BEGIN
    JD_stamp = strcompress(current_JD, /REMOVE_ALL)
 ENDELSE

  date_stamp = strmid(current_date,4,4)
  year_stamp = STRCOMPRESS(current_year)
  file_stamp = STRCOMPRESS(stamp+'_'+year_stamp+'_JD'+JD_stamp+'_'+date_stamp+'.dat', /REMOVE_ALL)
  matchfiles = where(file_stamp EQ out_list,matchfilecount)
  IF matchfilecount GT 0 THEN BEGIN
    print, 'Daily file '+file_stamp+' already exists'
    print, 'To overwrite this file, manually remove it from the directory and re-run script'
    GOTO, skipday
 ENDIF

  ;Now we create an array which will be initially populated with dates, times, and 'NaN'
  n_recs        = 86400/samplerate               ;determines how many records go in each array
  work_arr      = fltarr(n_field,n_recs)
  work_arr(*,*) = 'NaN'
  work_arr(0,*) = current_year
  work_arr(1,*) = current_m
  work_arr(2,*) = current_day

  ;place the hours
  ;---------------
  FOR r=0,23 DO BEGIN
    num   = n_recs/24
    start = r*num
    en    = start+num-1
    work_arr(3,start:en) = r 
 ENDFOR

  IF samplerate EQ 60 THEN GOTO, min_sample
  IF samplerate EQ 1 THEN GOTO, sec_sample
  IF samplerate LT 1 THEN GOTO, highf_sample

  ;place the minutes when the sampling time is 1 minute
  ;----------------------------------------------------
  min_sample:
  FOR r=0,59 DO BEGIN
    x = r
    REPEAT BEGIN
      work_arr(4,x) = r
      x = x + 60      
   ENDREP UNTIL (x EQ (n_recs+r))
 ENDFOR
 
  ;place the seconds when the sampling time is 1 minute
  ;----------------------------------------------------
  work_arr(5,*)=0
  GOTO, time_finished

  ;place the minutes when the sampling time is 1 second
  ;----------------------------------------------------
  sec_sample:
  FOR r=0,59 DO BEGIN
    x = long(r*60)
    REPEAT BEGIN
      start = x
      en    = long(x+60-1)
      work_arr(4,start:en) = r
      x = x + (60*60)
   ENDREP UNTIL (x EQ n_recs+(r*60))
 ENDFOR

  ;place the seconds when the sampling time is 1 second
  ;----------------------------------------------------
  FOR r=0,59 DO BEGIN
    x = LONG(r)
    REPEAT BEGIN
      work_arr(5,x) = r
      x = x + 60      
   ENDREP UNTIL (x EQ (n_recs+r))
 ENDFOR
  GOTO, time_finished

  ;place the minutes when the sampling time is high frequency
  ;----------------------------------------------------------
  highf_sample:
  FOR r=0,59 DO BEGIN
    x = long(r*60*(1/samplerate))
    REPEAT BEGIN
      start = x
      en    = long(x+60*(1/samplerate)-1)
      work_arr(4,start:en) = r
      x = x + (60*60*(1/samplerate))
   ENDREP UNTIL (x EQ (n_recs+(r*60*(1/samplerate))))
 ENDFOR

  ;place the seconds when the sampling time is high frequency
  ;----------------------------------------------------------
  FOR r=0,59 DO BEGIN
    x = float(r*(1/samplerate))
    REPEAT BEGIN
      start = x
      en    = long(x+(1/samplerate)-1)
      work_arr(5,start:en) = r
      x = x + (60*(1/samplerate))
   ENDREP UNTIL (x EQ (n_recs+(r*(1/samplerate))))
 ENDFOR

  ;place the tenths of seconds when the sampling time is high frequency
  ;--------------------------------------------------------------------
  interv=1/samplerate
  FOR r=0,interv-1 DO BEGIN
    x = LONG(r)
    REPEAT BEGIN
      work_arr(6,x) = r
      x = x + interv      
   ENDREP UNTIL (x EQ (n_recs+r))
 ENDFOR
  GOTO, time_finished

  time_finished:

  work_str = strcompress(work_arr, /REMOVE_ALL)

  ;Now start a loop that will go through each of the input files, and then it should be NO PROBLEM
  ;to match the data with the array using the WHERE function.

  CD, idir
  list = FILE_SEARCH('*.dat', COUNT=recs, /NOSORT)
  FOR k=0,recs-1 DO BEGIN
    ifile=list(k)
    ;decode the year and JD of the file being examined
    ex_year = 2011
    ex_mon  = LONG(STRMID(ifile,n_prefix,2))
    ex_day  = LONG(STRMID(ifile,n_prefix+2,2))
    ex_JD   = calendar_to_jd(ex_year,ex_mon,ex_day)
    
    IF ex_year LT current_year THEN GOTO, skipdoom
    
    ;if the file's JD is less than the current JD, we can go to the next file
    IF ((ex_JD LT current_JD) AND (ex_year EQ current_year)) THEN GOTO, skipdoom
  
    restore, file = itemplate_sav
    print, 'Processing Day:  ' + current_date + '     Opening File:   ' + ifile
    print, systime()
    print, '------------------------------------------'
    data = READ_ASCII(ifile, COUNT=n_inputfile, TEMPLATE=mytemplate);

    ;some extra processing for the MET data
    ;--------------------------------------
    ;for 2011 removed RH data when sensor was not working
    ;from 0931 UTC on JD204 until 1435 on JD207

    IF stamp EQ 'MET' THEN BEGIN 
    
    ;for RH data on julian days 204, replace with NaN if time greater than 0930
      rh_nc = where(data.(1)[*] EQ 2011 AND data.(2)[*] EQ 204 AND long(data.(3)[*] GT 930), rh_nc_count)
      IF rh_nc_count GT 0 THEN BEGIN
        data.(10)[rh_nc] = 'NaN'
     ENDIF
      
    ;for RH data on julian days 205 and 206, replace with NaN
      rh_nc = where(data.(1)[*] EQ 2011 AND data.(2)[*] GT 204 AND data.(2)[*] LT 207, rh_nc_count)
      IF rh_nc_count GT 0 THEN BEGIN
        data.(10)[rh_nc] = 'NaN'
     ENDIF
      
    ;for RH data on julian day 207, replace with NaN if time less than 1435
      rh_nc = where(data.(1)[*] EQ 2011 AND data.(2)[*] EQ 207 AND long(data.(3)[*] LT 1436), rh_nc_count)
      IF rh_nc_count GT 0 THEN BEGIN
        data.(10)[rh_nc] = 'NaN'
     ENDIF
   ENDIF
    
    stop
    ;Some extra processing for the RAD data
    ;--------------------------------------
    ;converting time of 2400 to time of 0000, and adding 1 to Julian Day
    IF stamp EQ 'RAD' THEN BEGIN
      tfclock = where(long(data.(3)[*]) EQ 2400, tfclock_count)
      IF tfclock_count GT 0 THEN BEGIN
        data.(3)[tfclock] = '0'
        data.(2)[tfclock] = data.(2)[tfclock]+1
     ENDIF
      
      ;for 2011 substituting 'NaN' when the UVS-AB-T sensor was not installed
      ;no data for all days before Julian Day 213 (Aug 1)
      no_uv_sensor = where((data.(1)[*] EQ 2011) AND (data.(2)[*] LT 213), no_uv_count)
      IF no_uv_count GT 0 THEN BEGIN
        data.(14)[no_uv_sensor] = 'NaN'
        data.(15)[no_uv_sensor] = 'NaN'
        data.(16)[no_uv_sensor] = 'NaN'
        data.(24)[no_uv_sensor] = 'NaN'
        data.(25)[no_uv_sensor] = 'NaN'
        data.(26)[no_uv_sensor] = 'NaN'
     ENDIF 
      
      ;no data for all times before 1515 UTC on Julian Day 213 (Aug 1)
      no_uv_sensor = where((data.(1)[*] EQ 2011) AND (data.(2)[*] EQ 213) AND (long(data.(3)[*] LT 1515)), no_uv_count)
      IF no_uv_count GT 0 THEN BEGIN
        data.(14)[no_uv_sensor] = 'NaN'
        data.(15)[no_uv_sensor] = 'NaN'
        data.(16)[no_uv_sensor] = 'NaN'
        data.(24)[no_uv_sensor] = 'NaN'
        data.(25)[no_uv_sensor] = 'NaN'
        data.(26)[no_uv_sensor] = 'NaN'
     ENDIF
   ENDIF 

    ;decode our date/time information
    input_hour = lonarr(1,n_elements(data.(3)[*]))
    input_min  = input_hour
    input_0    = where(long(data.(3)[*]) LT 60, match0)
    input_1    = where((long(data.(3)[*]) GT 59) AND (long(data.(3)[*] LT 1000)),match1)
    input_2    = where(long(data.(3)[*]) GT 999,match2)

    IF match0 GT 0 THEN BEGIN
      input_min(input_0) = long(data.(3)[input_0])
   ENDIF

    IF match1 GT 0 THEN BEGIN
      input_hour(input_1) = long(strmid(data.(3)[input_1],0,1))
      input_min(input_1)  = long(strmid(data.(3)[input_1],1,2))
   ENDIF

    IF match2 GT 0 THEN BEGIN
      input_hour(input_2) = long(strmid(data.(3)[input_2],0,2))
      input_min(input_2)  = long(strmid(data.(3)[input_2],2,2))
   ENDIF

    JD_arr=fltarr(1,n_recs)
    JD_arr(0,*)=current_JD


    ;START ANALYSIS FOR LOW FREQUENCY DAT
    ;------------------------------------
    IF samplerate LT 1 THEN GOTO, highfreqrun

    ;Loop through each file to see if they match
    FOR y=0L,n_inputfile-1L DO BEGIN

      IF stamp EQ 'RAD' THEN BEGIN
        match = where((data.(1)[y] EQ work_arr(0,*)) AND (data.(2)[y] EQ JD_arr(0,*)) AND (input_hour(y) EQ work_arr(3,*)) AND $
                (input_min(y) EQ work_arr(4,*)),matchcount)
     ENDIF ELSE BEGIN
        match = where((data.(1)[y] EQ work_arr(0,*)) AND (data.(2)[y] EQ JD_arr(0,*)) AND (input_hour(y) EQ work_arr(3,*)) AND $
                (input_min(y) EQ work_arr(4,*)) AND (data.(4)[y] EQ work_arr(5,*)),matchcount)
     ENDELSE

      IF match LT 0 THEN GOTO, skip_pop
      IF stamp EQ 'RAD' THEN BEGIN
        FOR poploop=6,n_field-1 DO BEGIN
          work_str(poploop,match) = strcompress(data.(poploop-2)[y], /REMOVE_ALL)
       ENDFOR
     ENDIF ELSE BEGIN
        FOR poploop=6,n_field-1 DO BEGIN
          work_str(poploop,match)=strcompress(data.(poploop-1)[y],/REMOVE_ALL)
       ENDFOR
     ENDELSE

      skip_pop:

  ENDFOR
    GOTO, donerun


    ;START ANALYSIS FOR HIGH FREQUENCY DATA
    ;--------------------------------------
    highfreqrun:
  
    FOR y=0L,n_inputfile-1L DO BEGIN
      match = where((data.(0)[y] EQ work_arr(0,*)) AND (data.(1)[y] EQ work_arr(1,*)) AND (data.(2)[y] EQ work_arr(2,*)) AND $
              (data.(3)[y] EQ work_arr(3,*)) AND (data.(4)[y] EQ work_arr(4,*)) AND (data.(5)[y] EQ work_arr(5,*)) AND $
              (data.(6)[y] EQ work_arr(6,*)),matchcount)
      IF match LT 0 THEN GOTO, skip_pop2
      FOR poploop=7,n_field-1 DO BEGIN
        work_str(poploop,match)=data.(poploop)[y]
     ENDFOR
      skip_pop2:
   ENDFOR
    GOTO, donerun

    donerun:
    close, /all; not sure if this is necessary....

    ;if the current filename JD was greater than the current JD, then we can end the loop.
    IF ((ex_JD GT current_JD) AND (ex_year EQ current_year)) THEN GOTO, skipdoom2
    skipdoom:
    close, /all; not sure if this is necessary....
 ENDFOR
  skipdoom2:

  ;Print it out
  print_file   = STRCOMPRESS(odir + stamp + '_' + year_stamp + '_JD' + JD_stamp + '_' + date_stamp + '.dat', /REMOVE_ALL)
  formatnum    = strcompress(n_field-1)
  formatstring = strcompress('('+formatnum+'(a,","),a)', /REMOVE_ALL)

  ;make the array a string for easy printing, and DO IT
  ;work_str=strcompress(work_arr, /REMOVE_ALL)
  openw, unit1, print_file, /get_lun
  printf, unit1, header
  printf, unit1, work_str, FORMAT = formatstring

  skipday:
  close, /all
ENDFOR


END

;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; day_splitter_MET.pro ends here
