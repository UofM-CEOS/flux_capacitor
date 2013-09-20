;; $Id: $
;;; day_splitter_MET.pro --- split MET data into daily files
;; Author: Sebastian P. Luque
;; Created: 2013-09-20T03:54:03+0000
;; Last-Updated: 2013-09-20T22:21:54+0000
;;           By: Sebastian Luque
;; ------------------------------------------------------------------------
;;; Commentary: 
;;
;; Lots to be done here.  Why are MET files not standardized, as NAV are?
;;
;; Example call:
;;
;; day_splitter_met, '20111015', '20111025', $
;;                   expand_path('~/tmp/MET'), $
;;                   expand_path('~/tmp/MET/Daily'), $
;;                   expand_path('met_template.sav'), $
;;                   'headera, headerb, ...', 60, 'MET'
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

  idir_files=file_search(idir + path_sep() + '*.dat', count=nidir_files, $
                         /nosort)
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
  n_ifields=met_raw_template.FIELDCOUNT ; hard-coding this for now (MUST FIX)
  header=strsplit(temporary(header), ', ', /extract)
  n_ofields=n_elements(header)

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

  all_arr=fltarr(n_ofields, size(times, /n_elements))
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
           print, 'Daily MET file ' + file_stamps[i] + ' already exists.  ' + $
                  'Overwriting'
        ENDIF ELSE BEGIN
           print, 'Daily MET file' + file_stamps[i] + ' already exists.  ' + $
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
     data=read_ascii(ifile, count=n_inputfile, template=met_raw_template)
     n_krecs=n_elements(data.(0))
     ;; Check each line and match against array
     FOR i=0L, n_krecs - 1, 1L DO BEGIN
        IF (data.(0)[i] LE 0) THEN BEGIN
           print, i, $
                  format='("Skipping record with unintelligible ' + $
                  'time stamp at: ", i)'
           CONTINUE
        ENDIF
        yyyymmdd=doy2calendar(data.(1)[i], data.(2)[i])
        hhmm=string(data.(3)[i], format='(i04)')
        ss=string(data.(4)[i], format='(i02)')
        i_ts=strmid(yyyymmdd, 0, 4) + '-' + strmid(yyyymmdd, 4, 2) + $
             '-' + strmid(yyyymmdd, 6, 2) + ' ' + strmid(hhmm, 0, 2) + $
             ':' + strmid(hhmm, 2) + ':' + ss
        match=where(times EQ i_ts, mcount)
        IF mcount EQ 1 THEN BEGIN
           FOR col=6, n_ifields - 1, 1L DO BEGIN $
              all_arr[col, match]=data.(col)[i]
           ENDFOR
        ENDIF
     ENDFOR
  ENDFOR

  ;; Write each daily array
  fmt_nfields=strtrim(n_ofields - 1, 2)
  fmt_str='(' + fmt_nfields + '(a,","),a)'
  all_arr=strcompress(temporary(all_arr), /remove_all)
  FOR begi=0L, n_elements(times) - 1, n_recs DO BEGIN
     file_stamp=file_stamps[begi / n_recs]
     openw, lun, odir + path_sep() + file_stamp, /get_lun
     printf, lun, strjoin(header, ', ')
     printf, lun, all_arr[*, begi:(begi + n_recs - 1)], format=fmt_str
     free_lun, lun
  ENDFOR


;;   ;Now start a loop that will go through each of the input files, and then it should be NO PROBLEM
;;   ;to match the data with the array using the WHERE function.

;;     ;some extra processing for the MET data
;;     ;--------------------------------------
;;     ;for 2011 removed RH data when sensor was not working
;;     ;from 0931 UTC on JD204 until 1435 on JD207

;;     IF stamp EQ 'MET' THEN BEGIN 
    
;;     ;for RH data on julian days 204, replace with NaN if time greater than 0930
;;       rh_nc = where(data.(1)[*] EQ 2011 AND data.(2)[*] EQ 204 AND long(data.(3)[*] GT 930), rh_nc_count)
;;       IF rh_nc_count GT 0 THEN BEGIN
;;         data.(10)[rh_nc] = 'NaN'
;;      ENDIF
      
;;     ;for RH data on julian days 205 and 206, replace with NaN
;;       rh_nc = where(data.(1)[*] EQ 2011 AND data.(2)[*] GT 204 AND data.(2)[*] LT 207, rh_nc_count)
;;       IF rh_nc_count GT 0 THEN BEGIN
;;         data.(10)[rh_nc] = 'NaN'
;;      ENDIF
      
;;     ;for RH data on julian day 207, replace with NaN if time less than 1435
;;       rh_nc = where(data.(1)[*] EQ 2011 AND data.(2)[*] EQ 207 AND long(data.(3)[*] LT 1436), rh_nc_count)
;;       IF rh_nc_count GT 0 THEN BEGIN
;;         data.(10)[rh_nc] = 'NaN'
;;      ENDIF
;;    ENDIF
    
;;     stop
;;     ;Some extra processing for the RAD data
;;     ;--------------------------------------
;;     ;converting time of 2400 to time of 0000, and adding 1 to Julian Day
;;     IF stamp EQ 'RAD' THEN BEGIN
;;       tfclock = where(long(data.(3)[*]) EQ 2400, tfclock_count)
;;       IF tfclock_count GT 0 THEN BEGIN
;;         data.(3)[tfclock] = '0'
;;         data.(2)[tfclock] = data.(2)[tfclock]+1
;;      ENDIF
      
;;       ;for 2011 substituting 'NaN' when the UVS-AB-T sensor was not installed
;;       ;no data for all days before Julian Day 213 (Aug 1)
;;       no_uv_sensor = where((data.(1)[*] EQ 2011) AND (data.(2)[*] LT 213), no_uv_count)
;;       IF no_uv_count GT 0 THEN BEGIN
;;         data.(14)[no_uv_sensor] = 'NaN'
;;         data.(15)[no_uv_sensor] = 'NaN'
;;         data.(16)[no_uv_sensor] = 'NaN'
;;         data.(24)[no_uv_sensor] = 'NaN'
;;         data.(25)[no_uv_sensor] = 'NaN'
;;         data.(26)[no_uv_sensor] = 'NaN'
;;      ENDIF 
      
;;       ;no data for all times before 1515 UTC on Julian Day 213 (Aug 1)
;;       no_uv_sensor = where((data.(1)[*] EQ 2011) AND (data.(2)[*] EQ 213) AND (long(data.(3)[*] LT 1515)), no_uv_count)
;;       IF no_uv_count GT 0 THEN BEGIN
;;         data.(14)[no_uv_sensor] = 'NaN'
;;         data.(15)[no_uv_sensor] = 'NaN'
;;         data.(16)[no_uv_sensor] = 'NaN'
;;         data.(24)[no_uv_sensor] = 'NaN'
;;         data.(25)[no_uv_sensor] = 'NaN'
;;         data.(26)[no_uv_sensor] = 'NaN'
;;      ENDIF
;;    ENDIF 

;;     ;decode our date/time information
;;     input_hour = lonarr(1,n_elements(data.(3)[*]))
;;     input_min  = input_hour
;;     input_0    = where(long(data.(3)[*]) LT 60, match0)
;;     input_1    = where((long(data.(3)[*]) GT 59) AND (long(data.(3)[*] LT 1000)),match1)
;;     input_2    = where(long(data.(3)[*]) GT 999,match2)

;;     IF match0 GT 0 THEN BEGIN
;;       input_min(input_0) = long(data.(3)[input_0])
;;    ENDIF

;;     IF match1 GT 0 THEN BEGIN
;;       input_hour(input_1) = long(strmid(data.(3)[input_1],0,1))
;;       input_min(input_1)  = long(strmid(data.(3)[input_1],1,2))
;;    ENDIF

;;     IF match2 GT 0 THEN BEGIN
;;       input_hour(input_2) = long(strmid(data.(3)[input_2],0,2))
;;       input_min(input_2)  = long(strmid(data.(3)[input_2],2,2))
;;    ENDIF

;;     JD_arr=fltarr(1,n_recs)
;;     JD_arr(0,*)=current_JD


;;     ;START ANALYSIS FOR LOW FREQUENCY DAT
;;     ;------------------------------------
;;     IF samplerate LT 1 THEN GOTO, highfreqrun

;;     ;Loop through each file to see if they match
;;     FOR y=0L,n_inputfile-1L DO BEGIN

;;       IF stamp EQ 'RAD' THEN BEGIN
;;         match = where((data.(1)[y] EQ work_arr(0,*)) AND (data.(2)[y] EQ JD_arr(0,*)) AND (input_hour(y) EQ work_arr(3,*)) AND $
;;                 (input_min(y) EQ work_arr(4,*)),matchcount)
;;      ENDIF ELSE BEGIN
;;         match = where((data.(1)[y] EQ work_arr(0,*)) AND (data.(2)[y] EQ JD_arr(0,*)) AND (input_hour(y) EQ work_arr(3,*)) AND $
;;                 (input_min(y) EQ work_arr(4,*)) AND (data.(4)[y] EQ work_arr(5,*)),matchcount)
;;      ENDELSE

;;       IF match LT 0 THEN GOTO, skip_pop
;;       IF stamp EQ 'RAD' THEN BEGIN
;;         FOR poploop=6,n_field-1 DO BEGIN
;;           work_str(poploop,match) = strcompress(data.(poploop-2)[y], /REMOVE_ALL)
;;        ENDFOR
;;      ENDIF ELSE BEGIN
;;         FOR poploop=6,n_field-1 DO BEGIN
;;           work_str(poploop,match)=strcompress(data.(poploop-1)[y],/REMOVE_ALL)
;;        ENDFOR
;;      ENDELSE

;;       skip_pop:

;;   ENDFOR
;;     GOTO, donerun


;;     ;START ANALYSIS FOR HIGH FREQUENCY DATA
;;     ;--------------------------------------
;;     highfreqrun:
  
;;     FOR y=0L,n_inputfile-1L DO BEGIN
;;       match = where((data.(0)[y] EQ work_arr(0,*)) AND (data.(1)[y] EQ work_arr(1,*)) AND (data.(2)[y] EQ work_arr(2,*)) AND $
;;               (data.(3)[y] EQ work_arr(3,*)) AND (data.(4)[y] EQ work_arr(4,*)) AND (data.(5)[y] EQ work_arr(5,*)) AND $
;;               (data.(6)[y] EQ work_arr(6,*)),matchcount)
;;       IF match LT 0 THEN GOTO, skip_pop2
;;       FOR poploop=7,n_field-1 DO BEGIN
;;         work_str(poploop,match)=data.(poploop)[y]
;;      ENDFOR
;;       skip_pop2:
;;    ENDFOR
;;     GOTO, donerun

;;     donerun:
;;     close, /all; not sure if this is necessary....

;;     ;if the current filename JD was greater than the current JD, then we can end the loop.
;;     IF ((ex_JD GT current_JD) AND (ex_year EQ current_year)) THEN GOTO, skipdoom2
;;     skipdoom:
;;     close, /all; not sure if this is necessary....
;;  ENDFOR
;;   skipdoom2:

;;   ;Print it out
;;   print_file   = STRCOMPRESS(odir + stamp + '_' + year_stamp + '_JD' + JD_stamp + '_' + date_stamp + '.dat', /REMOVE_ALL)
;;   formatnum    = strcompress(n_field-1)
;;   formatstring = strcompress('('+formatnum+'(a,","),a)', /REMOVE_ALL)

;;   ;make the array a string for easy printing, and DO IT
;;   ;work_str=strcompress(work_arr, /REMOVE_ALL)
;;   openw, unit1, print_file, /get_lun
;;   printf, unit1, header
;;   printf, unit1, work_str, FORMAT = formatstring

;;   skipday:
;;   close, /all
;; ENDFOR


END

;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; day_splitter_MET.pro ends here
