;; $Id$
;;; day_splitter.pro --- split input data into daily files
;; Author: Sebastian P. Luque
;; Created: 2013-09-20T03:54:03+0000
;; Last-Updated: 2013-09-28T19:37:07+0000
;;           By: Sebastian Luque
;; ------------------------------------------------------------------------
;;; Commentary: 
;;
;; This is a complete mess; lots to be done here.  Why are MET files not
;; standardized, as NAV are?
;;
;; Example call:
;;
;; day_splitter, '20110723', '20111022', $
;;               expand_path('~/tmp/ArcticNet2011/NAV/STD'), $
;;               expand_path('~/tmp/Arcticnet2011/NAV/Daily'), $
;;               expand_path('nav_std_template.sav'), 0, 1, 'NAV', /overwrite
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

PRO DAY_SPLITTER, STARTDATE, ENDDATE, IDIR, ODIR, ITEMPLATE_SAV, $
                  TIME_BEG_IDX, ISAMPLERATE, STAMP, OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 8) THEN $
     message, 'Usage: DAY_SPLITTER, STARTDATE, ENDDATE, IDIR, ' + $
              'ODIR, ITEMPLATE_SAV, TIME_BEG_IDX, ISAMPLERATE, ' + $
              'STAMP, [/OVERWRITE]'
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
  IF ((n_elements(time_beg_idx) NE 1) OR $
      ((size(time_beg_idx, /type) NE 2) || time_beg_idx LT 0)) THEN $
         message, 'TIME_BEG_IDX must be an integer scalar >= zero'
  IF (n_elements(isamplerate) EQ 0) THEN $
     message, 'ISAMPLERATE is undefined'
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
  n_recs=86400 / isamplerate     ; records per day
  restore, itemplate_sav
  n_fields=itemplate.FIELDCOUNT
  header=itemplate.FIELDNAMES
  is_time_field=itemplate.FIELDGROUPS EQ time_beg_idx
  time_fields=where(is_time_field, /NULL)
  time_field_names=strlowcase(header[time_fields])
  year_subfield=where(time_field_names EQ 'year') ; locate year
  month_subfield=where(time_field_names EQ 'month') ; locate month
  day_subfield=where(time_field_names EQ 'day') ; locate day
  hour_subfield=where(time_field_names EQ 'hour') ; locate hour
  minute_subfield=where(time_field_names EQ 'minute') ; locate minute
  second_subfield=where(time_field_names EQ 'second') ; locate second
  ;; Check if we have a name with "second" as substring somewhere after the
  ;; first character; matches e.g.: "decisecond", "millisecond"
  subsecond_exists=strpos(time_field_names, 'second', 1)
  subsecond_subfield=where(subsecond_exists GE 0)

  ;; Add 1 to the last day to set the breaks for value_locate()
  days=timegen(start=beg_jd, final=end_jd + 1, step_size=1, units='days')
  days=jul2timestamp(temporary(days))
  ;; Add up to 23:59:59 of the last requested day
  times=timegen(start=beg_jd, final=end_jd + (86399.0 / 86400), $
                step_size=isamplerate, units='seconds')
  times=jul2timestamp(temporary(times))
  ;; ;; Index times in days vector.  Note that value_locate() does work with
  ;; ;; strings, which in this case is appropriate given jul2timestamp()
  ;; pos=value_locate(days, times)
  
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
           print, 'Daily ' + stamp + ' file ' + file_stamps[i] + $
                  ' already exists.  Overwriting'
        ENDIF ELSE BEGIN
           print, 'Daily ' + stamp + ' file ' + file_stamps[i] + $
                  ' already exists.  Not overwriting'
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
     idx=histogram(skips, min=0, max=n_elements(times) - 1)
     ;; Subset times for matching
     times=times[where(~idx)]
  ENDIF

  ;; We create a hash to hold full time series
  ts_all=hash(header)           ; just the keys here
  ;; Iterate through time fields, and fill the hash with time data and
  ;; holders for the rest. We do not care about matching fractions of
  ;; seconds, for the high frequency data.  It is treated as a non-time
  ;; field.
  FOREACH value, ts_all, fld DO BEGIN
     CASE fld OF
        'year': ts_all[fld]=strmid(times, 0, 4)
        'month': ts_all[fld]=strmid(times, 5, 2)
        'day': ts_all[fld]=strmid(times, 8, 2)
        'hour': ts_all[fld]=strmid(times, 11, 2)
        'minute': ts_all[fld]=strmid(times, 14, 2)
        'second': ts_all[fld]=strmid(times, 17, 2)
        ELSE: BEGIN
           fld_idx=where(itemplate.FIELDNAMES EQ fld)
           ifield_type=itemplate.FIELDTYPES[fld_idx]
           val=(ifield_type EQ 4) ? !VALUES.F_NAN : ''
           ts_all[fld]=make_array(n_elements(times), type=ifield_type, $
                                  value=val)
        END
     ENDCASE
  ENDFOREACH
  ;; Separate times from full hash
  ts_times=ts_all.remove(time_field_names)
  non_time_fields=where(~is_time_field)
  non_time_field_names=strlowcase(header[non_time_fields])

  ;; Read each file
  FOR k=0L, nidir_files - 1 DO BEGIN
     ifile=idir_files[k]
     print, 'Processing File: ' + ifile
     idata=read_ascii(ifile, count=n_inputfile, template=itemplate)
     ;; Extract times and remove quotes or spaces from strings
     idata_times=idata.(time_beg_idx)
     idata=remove_structure_tag(idata, (tag_names(idata))[time_beg_idx])
     idata_names=strlowcase(tag_names(idata))
     IF size(idata_times, /type) EQ 7 THEN BEGIN
        FOREACH fld, indgen((size(idata_times, /dimensions))[0]) DO BEGIN
           ok=strsplit(idata_times[fld, *], '" ', /extract)
           idata_times[fld, *]=ok.toArray()
        ENDFOREACH
     ENDIF
     FOREACH fld, indgen(n_tags(idata)) DO BEGIN
        IF size(idata.(fld), /type) EQ 7 THEN BEGIN
           ok=strsplit(idata.(fld), '" ', /extract)
           idata.(fld)=ok.toArray()
        ENDIF
     ENDFOREACH
     n_krecs=n_elements(idata_times[0, *])
     ;; Check each line and match against array
     is_valid=valid_num(idata_times[year_subfield, *])
     iok=where(is_valid, vcount) ; work on valid data only
     IF vcount LT 1 THEN CONTINUE
     idata_times=idata_times[*, iok]
     yyyy=string(idata_times[year_subfield, *], format='(i4)')
     mo=string(idata_times[month_subfield, *], format='(i02)')
     dd=string(idata_times[day_subfield, *], format='(i02)')
     hh=string(idata_times[hour_subfield, *], format='(i02)')
     mm=string(idata_times[minute_subfield, *], format='(i02)')
     ss=string(idata_times[second_subfield, *], format='(i02)')
     i_ts=yyyy + '-' + mo + '-' + dd + ' ' + hh + ':' + mm + ':' + ss
     match2, times, i_ts, times_in_its, i_ts_in_times
     t_matches=where(times_in_its GE 0, mcount, /null)
     IF mcount LT 1 THEN CONTINUE
     its_matches=where(i_ts_in_times GE 0)
     FOREACH value, ts_times, fld DO BEGIN ; loop time hash
        match_fld=where(time_field_names EQ strlowcase(fld))
        ts_times[fld, t_matches]=(idata_times[match_fld, *])[its_matches]
     ENDFOREACH
     FOREACH value, ts_all, fld DO BEGIN ; loop non-time hash
        match_fld=where(non_time_field_names EQ strlowcase(fld))
        ;; Note we have to subset the field in the structure, as structure
        ;; elements cannot have their size/type changed.
        ts_all[fld, t_matches]=(idata.(match_fld)[iok])[its_matches]
     ENDFOREACH
  ENDFOR

  ;; Prepare output
  ts_full=ts_times + ts_all
  delvar, ts_times, ts_all
  ;; Write each daily array
  FOR begi=0L, n_elements(times) - 1, n_recs DO BEGIN
     endi=(begi + n_recs - 1)
     file_stamp=file_stamps[begi / n_recs]
     ;; ts=ts_full.toStruct(/no_copy)
     ;; There must be a better way to re-order tags
     ts=create_struct(header[0], ts_full[header[0], begi:endi])
     FOREACH fld, header[1:*] DO BEGIN
        ts=create_struct(ts, header[where(header EQ fld)], $
                         ts_full[fld, begi:endi])
     ENDFOREACH
     write_csv, odir + path_sep() + file_stamp, ts, header=header
  ENDFOR

END

;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; day_splitter.pro ends here
