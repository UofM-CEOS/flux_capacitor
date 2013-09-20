;This program takes all of the .dat files in an input directory, and breaks them up into 
;individual day files.  It will create a day file for each day input between startdate-enddate
;and fill any missing data with 'NaN'

;NOTE: in it's present configuration, you CANNOT process more than 1 year of data
;  However, you can process across a calendar year (e.g. from 20071001 to 20080701)
;  It may slog a bit across a calendar year... try to avoid this.

;NOTE: this program might be a bit inefficient when it comes to processing large flux files

;NOTE: for this to work, the input files should have the dates/times in their first columns:
;  For 1min or 1sec data: year,month,day,hour,min,sec
;  For high frequency data: year,month,day,hour,min,sec,tenth of sec

;NOTE: this program has NOT been thoroughly tested for anything except 60s data.

;NOTE: you need to double check that the time stamp on the files is correct (i.e. the proper day (UTC) it was downloaded)

;Input variables:
;  startdate     - first day for generating daily files (YYYYMMDD) --> User prompted to enter date
;  enddate       - last day for generating daily files (YYYYMMDD)  --> User prompted to enter date
;  idir          - input directory of files to be processed
;  odir          - output directory for new daily files
;  itemplate_sav - directory of .sav file of template being used NOTE: the template must have been saved with the variable name "mytemplate"
;  header        - text string that defines the header
;  samplerate    - The sample rate which data was collected in seconds (e.g. 1, 60)
;  nfields       - The number of fields that are in the data file to be created 
;  stamp         - gives a prefix to the output files (e.g. if set to pCO2, files will be pCO2_YYYY_JD_MMDD.dat)
;  n_prefix      - gives the number of characters BEFORE the date of the input file name (e.g. for pAMD_pCO2_071028.dat, the number is 10)
; 
;For MET data files:
; idir           = met_idir       --> '/ArcticNet2011/TowerData/MET/raw-MET/'     -- file format "AMD_MET_mmdd_hhmm.dat"
; odir           = met_dailydir   --> '/ArcticNet2011/TowerData/MET/daily-MET/'   -- file format "MET_yyyy_JDxxx_mmdd.dat"
; itemplate_sav  = met_template   --> '/ArcticNet2011/Templates/met_template.sav'
; header         = met_header     --> 'Year, Month, Day, Hour, Minute, Sec, Prog Version, Battery(V), Panel_T(C), Pressure(kPa), Air_T(C),' +$
;                                     'RH(%), Surface_T(C), Raw Wind Speed(m/s), Raw Wind Direction(deg), Wind StdDev, PAR(umol/m2/s), Pitch(deg),' +$
;                                     'Roll (deg), Battery StdDev, Panel_T StdDev, Pressure StDev, Air_T StDev, RH StDev, Surface_T StDev, PAR StDev'
; samplerate     = met_timing     --> 60
; nfields        = met_outcolumns --> 26
; stamp          = met_stamp      --> 'MET'
; n_prefix       = met_prefix     --> 8
;
;For RAD data files:
; idir           = raw_rad_dir    --> '/ArcticNet2011/RadData/raw-RAD'
; odir           = daily_rad_dir  --> '/ArcticNet2011/RadData/daily-RAD/'
; itemplate_sav  = Rad_template   --> '/ArcticNet2011/Templates/RAD_template.sav'
; header         = RAD_head       --> 'Year, Month, Day, Hour, Minute, Sec, Battery(V), Panel_T(C), Batt_stdev, PanelT_stdev, Kdown(W/m2), Thermopile,' +$
;                                     'Tcase(K),Tdome_avg, LWin(W/m2), PAR(umol/m2/s), Tuv(C), UVb(W/m2), UVa(W/m2), UVbroad(W/m2), Kdown_stdev,' +$
;                                     'Thermopile_stdev, Tcase_stdev, Tdome_stdev, LWin_stdev, PAR_stdev, Tuv_stdev, Uvb_stdev, Uva_stdev, UVbroad_stdev' 
; samplerate     = rad_timing     --> 60
; nfields        = rad_outcolumns --> 30
; stamp          = rad_stamp      --> 'RAD'
; n_prefix       = rad_prefix     --> 8


PRO day_splitter_MET, startdate, enddate, idir, odir, itemplate_sav, header, samplerate, n_field, stamp, n_prefix

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
IF ((s_y mod 4) eq 0) then begin
  s_y_leap = 1
ENDIF

e_y_leap = 0
IF ((e_y mod 4) eq 0) then begin
  e_y_leap = 1
ENDIF

;Now calculate how big of a loop you need
IF (s_y eq e_y) THEN BEGIN
  n_days = endJD - startJD
ENDIF ELSE BEGIN
IF s_y_leap eq 0 then begin
 s_days = 365 - startJD
ENDIF ELSE BEGIN
  s_days = 366 - startJD
ENDELSE
  n_days = s_days + endJD
ENDELSE

n_days=n_days[0]
FOR loopdays = 0, n_days do begin

  ;Calculate what Julian Day and year you're working with
  current_JD   = startJD + loopdays
  current_year = s_y
  IF s_y_leap eq 0 then begin
    IF (current_JD gt 365) then begin
      current_JD   = current_JD - 365
      current_year = s_y + 1
    ENDIF
  ENDIF ELSE BEGIN
    IF (current_JD gt 366) then begin
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
  if current_JD lt 100 then begin
    JD_stamp = strcompress(current_JD)
    JD_stamp = '0'+JD_stamp
    if current_JD lt 10 then begin
      JD_stamp = strcompress(current_JD)
      JD_stamp = '00'+JD_stamp
    endif
  endif else begin
    JD_stamp = strcompress(current_JD, /REMOVE_ALL)
  endelse

  date_stamp = strmid(current_date,4,4)
  year_stamp = STRCOMPRESS(current_year)
  file_stamp = STRCOMPRESS(stamp+'_'+year_stamp+'_JD'+JD_stamp+'_'+date_stamp+'.dat', /REMOVE_ALL)
  matchfiles = where(file_stamp eq out_list,matchfilecount)
  if matchfilecount gt 0 then begin
    print, 'Daily file '+file_stamp+' already exists'
    print, 'To overwrite this file, manually remove it from the directory and re-run script'
    goto, skipday
  endif

  ;Now we create an array which will be initially populated with dates, times, and 'NaN'
  n_recs        = 86400/samplerate               ;determines how many records go in each array
  work_arr      = fltarr(n_field,n_recs)
  work_arr(*,*) = 'NaN'
  work_arr(0,*) = current_year
  work_arr(1,*) = current_m
  work_arr(2,*) = current_day

  ;place the hours
  ;---------------
  for r=0,23 do begin
    num   = n_recs/24
    start = r*num
    en    = start+num-1
    work_arr(3,start:en) = r 
  endfor

  if samplerate eq 60 then goto, min_sample
  if samplerate eq 1 then goto, sec_sample
  if samplerate lt 1 then goto, highf_sample

  ;place the minutes when the sampling time is 1 minute
  ;----------------------------------------------------
  min_sample:
  for r=0,59 do begin
    x = r
    Repeat begin
      work_arr(4,x) = r
      x = x + 60      
    Endrep until (x eq (n_recs+r))
  endfor
 
  ;place the seconds when the sampling time is 1 minute
  ;----------------------------------------------------
  work_arr(5,*)=0
  goto, time_finished

  ;place the minutes when the sampling time is 1 second
  ;----------------------------------------------------
  sec_sample:
  for r=0,59 do begin
    x = long(r*60)
    repeat begin
      start = x
      en    = long(x+60-1)
      work_arr(4,start:en) = r
      x = x + (60*60)
    endrep until (x eq n_recs+(r*60))
  endfor

  ;place the seconds when the sampling time is 1 second
  ;----------------------------------------------------
  for r=0,59 do begin
    x = LONG(r)
    Repeat begin
      work_arr(5,x) = r
      x = x + 60      
    Endrep until (x eq (n_recs+r))
  endfor
  goto, time_finished

  ;place the minutes when the sampling time is high frequency
  ;----------------------------------------------------------
  highf_sample:
  for r=0,59 do begin
    x = long(r*60*(1/samplerate))
    repeat begin
      start = x
      en    = long(x+60*(1/samplerate)-1)
      work_arr(4,start:en) = r
      x = x + (60*60*(1/samplerate))
    endrep until (x eq (n_recs+(r*60*(1/samplerate))))
  endfor

  ;place the seconds when the sampling time is high frequency
  ;----------------------------------------------------------
  for r=0,59 do begin
    x = float(r*(1/samplerate))
    repeat begin
      start = x
      en    = long(x+(1/samplerate)-1)
      work_arr(5,start:en) = r
      x = x + (60*(1/samplerate))
    endrep until (x eq (n_recs+(r*(1/samplerate))))
  endfor

  ;place the tenths of seconds when the sampling time is high frequency
  ;--------------------------------------------------------------------
  interv=1/samplerate
  for r=0,interv-1 do begin
    x = LONG(r)
    Repeat begin
      work_arr(6,x) = r
      x = x + interv      
    Endrep until (x eq (n_recs+r))
  endfor
  goto, time_finished

  time_finished:

  work_str = strcompress(work_arr, /REMOVE_ALL)

  ;Now start a loop that will go through each of the input files, and then it should be NO PROBLEM
  ;to match the data with the array using the WHERE function.

  CD, idir
  list = FILE_SEARCH('*.dat', COUNT=recs, /NOSORT)
  for k=0,recs-1 do begin
    ifile=list(k)
    ;decode the year and JD of the file being examined
    ex_year = 2011
    ex_mon  = LONG(STRMID(ifile,n_prefix,2))
    ex_day  = LONG(STRMID(ifile,n_prefix+2,2))
    ex_JD   = calendar_to_jd(ex_year,ex_mon,ex_day)
    
    if ex_year lt current_year then goto, skipdoom
    
    ;if the file's JD is less than the current JD, we can go to the next file
    if ((ex_JD lt current_JD) AND (ex_year eq current_year)) then goto, skipdoom
  
    restore, file = itemplate_sav
    print, 'Processing Day:  ' + current_date + '     Opening File:   ' + ifile
    print, systime()
    print, '------------------------------------------'
    data = READ_ASCII(ifile, COUNT=n_inputfile, TEMPLATE=mytemplate);

    ;some extra processing for the MET data
    ;--------------------------------------
    ;for 2011 removed RH data when sensor was not working
    ;from 0931 UTC on JD204 until 1435 on JD207

    if stamp eq 'MET' then begin 
    
    ;for RH data on julian days 204, replace with NaN if time greater than 0930
      rh_nc = where(data.(1)[*] eq 2011 AND data.(2)[*] eq 204 AND long(data.(3)[*] gt 930), rh_nc_count)
      if rh_nc_count gt 0 then begin
        data.(10)[rh_nc] = 'NaN'
      endif
      
    ;for RH data on julian days 205 and 206, replace with NaN
      rh_nc = where(data.(1)[*] eq 2011 AND data.(2)[*] gt 204 AND data.(2)[*] lt 207, rh_nc_count)
      if rh_nc_count gt 0 then begin
        data.(10)[rh_nc] = 'NaN'
      endif
      
    ;for RH data on julian day 207, replace with NaN if time less than 1435
      rh_nc = where(data.(1)[*] eq 2011 AND data.(2)[*] eq 207 AND long(data.(3)[*] lt 1436), rh_nc_count)
      if rh_nc_count gt 0 then begin
        data.(10)[rh_nc] = 'NaN'
      endif
    endif
    
    stop
    ;Some extra processing for the RAD data
    ;--------------------------------------
    ;converting time of 2400 to time of 0000, and adding 1 to Julian Day
    if stamp eq 'RAD' then begin
      tfclock = where(long(data.(3)[*]) eq 2400, tfclock_count)
      if tfclock_count gt 0 then begin
        data.(3)[tfclock] = '0'
        data.(2)[tfclock] = data.(2)[tfclock]+1
      endif
      
      ;for 2011 substituting 'NaN' when the UVS-AB-T sensor was not installed
      ;no data for all days before Julian Day 213 (Aug 1)
      no_uv_sensor = where((data.(1)[*] eq 2011) AND (data.(2)[*] lt 213), no_uv_count)
      if no_uv_count gt 0 then begin
        data.(14)[no_uv_sensor] = 'NaN'
        data.(15)[no_uv_sensor] = 'NaN'
        data.(16)[no_uv_sensor] = 'NaN'
        data.(24)[no_uv_sensor] = 'NaN'
        data.(25)[no_uv_sensor] = 'NaN'
        data.(26)[no_uv_sensor] = 'NaN'
      endif 
      
      ;no data for all times before 1515 UTC on Julian Day 213 (Aug 1)
      no_uv_sensor = where((data.(1)[*] eq 2011) AND (data.(2)[*] eq 213) AND (long(data.(3)[*] lt 1515)), no_uv_count)
      if no_uv_count gt 0 then begin
        data.(14)[no_uv_sensor] = 'NaN'
        data.(15)[no_uv_sensor] = 'NaN'
        data.(16)[no_uv_sensor] = 'NaN'
        data.(24)[no_uv_sensor] = 'NaN'
        data.(25)[no_uv_sensor] = 'NaN'
        data.(26)[no_uv_sensor] = 'NaN'
      endif
    endif 

    ;decode our date/time information
    input_hour = lonarr(1,n_elements(data.(3)[*]))
    input_min  = input_hour
    input_0    = where(long(data.(3)[*]) lt 60, match0)
    input_1    = where((long(data.(3)[*]) gt 59) AND (long(data.(3)[*] lt 1000)),match1)
    input_2    = where(long(data.(3)[*]) gt 999,match2)

    if match0 gt 0 then begin
      input_min(input_0) = long(data.(3)[input_0])
    endif

    if match1 gt 0 then begin
      input_hour(input_1) = long(strmid(data.(3)[input_1],0,1))
      input_min(input_1)  = long(strmid(data.(3)[input_1],1,2))
    endif

    if match2 gt 0 then begin
      input_hour(input_2) = long(strmid(data.(3)[input_2],0,2))
      input_min(input_2)  = long(strmid(data.(3)[input_2],2,2))
    endif

    JD_arr=fltarr(1,n_recs)
    JD_arr(0,*)=current_JD


    ;START ANALYSIS FOR LOW FREQUENCY DAT
    ;------------------------------------
    if samplerate lt 1 then goto, highfreqrun

    ;Loop through each file to see if they match
    for y=0L,n_inputfile-1L do begin

      if stamp eq 'RAD' then begin
        match = where((data.(1)[y] eq work_arr(0,*)) AND (data.(2)[y] eq JD_arr(0,*)) AND (input_hour(y) eq work_arr(3,*)) AND $
                (input_min(y) eq work_arr(4,*)),matchcount)
      endif else begin
        match = where((data.(1)[y] eq work_arr(0,*)) AND (data.(2)[y] eq JD_arr(0,*)) AND (input_hour(y) eq work_arr(3,*)) AND $
                (input_min(y) eq work_arr(4,*)) AND (data.(4)[y] eq work_arr(5,*)),matchcount)
      endelse

      if match lt 0 then goto, skip_pop
      if stamp eq 'RAD' then begin
        for poploop=6,n_field-1 do begin
          work_str(poploop,match) = strcompress(data.(poploop-2)[y], /REMOVE_ALL)
        endfor
      endif else begin
        for poploop=6,n_field-1 do begin
          work_str(poploop,match)=strcompress(data.(poploop-1)[y],/REMOVE_ALL)
        endfor
      endelse

      skip_pop:

    endfor
    goto, donerun


    ;START ANALYSIS FOR HIGH FREQUENCY DATA
    ;--------------------------------------
    highfreqrun:
  
    for y=0L,n_inputfile-1L do begin
      match = where((data.(0)[y] eq work_arr(0,*)) AND (data.(1)[y] eq work_arr(1,*)) AND (data.(2)[y] eq work_arr(2,*)) AND $
              (data.(3)[y] eq work_arr(3,*)) AND (data.(4)[y] eq work_arr(4,*)) AND (data.(5)[y] eq work_arr(5,*)) AND $
              (data.(6)[y] eq work_arr(6,*)),matchcount)
      if match lt 0 then goto, skip_pop2
      for poploop=7,n_field-1 do begin
        work_str(poploop,match)=data.(poploop)[y]
      endfor
      skip_pop2:
    endfor
    goto, donerun

    donerun:
    close, /all; not sure if this is necessary....

    ;if the current filename JD was greater than the current JD, then we can end the loop.
    if ((ex_JD gt current_JD) AND (ex_year eq current_year)) then goto, skipdoom2
    skipdoom:
    close, /all; not sure if this is necessary....
  endfor
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
endfor


end
