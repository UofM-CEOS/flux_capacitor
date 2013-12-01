;; $Id$
;; Author: Sebastian P. Luque
;; Created: 2013-09-20T03:54:03+0000
;; Last-Updated: 2013-12-01T02:30:39+0000
;;           By: Sebastian P. Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;;
;;     DAY_SPLITTER
;;
;; PURPOSE:
;;
;;     Generates a regular time series from a given start date through a
;;     given end date, with data for each time step from input files.  It
;;     produces one file per day, consisting of the entire day's time
;;     series.
;;
;; CALLING SEQUENCE:
;;
;;     DAY_SPLITTER, StartDate, EndDate, Idir, Odir, Itemplate_Sav, $
;;     Time_Beg_Idx, Step_Time, Stamp
;;
;; INPUTS:
;;
;;     StartDate:       String scalar for starting date.
;;     EndDate:         String scalar for ending date.
;;     Idir:            Input directory (no trailing separator).
;;     Odir:            Output directory (no trailing separator).
;;     Itemplate_Sav:   Ascii template to read input files.
;;     Time_Idx:        Index (in template) where time is.
;;     Step_Time:       Scalar for step time in seconds.
;;     Stamp:           Scalar for string to preprend to output file name.
;;
;; KEYWORD PARAMETERS:
;;
;;     OVERWRITE:             Whether to overwrite files in Odir.
;;
;; EXAMPLE:
;;
;;     DAY_SPLITTER, '20110719', '20111022', $
;;                   expand_path('~/tmp/ArcticNet2011/NAV/STD'), $
;;                   expand_path('~/tmp/ArcticNet2011/NAV/Daily'), $
;;                   'nav_std_template.sav', 0, nav_daily_rate, $
;;                   nav_stamp, /overwrite
;;
;;- -----------------------------------------------------------------------
;;; Code:

PRO DAY_SPLITTER, STARTDATE, ENDDATE, IDIR, ODIR, ITEMPLATE_SAV, $
                  TIME_IDX, STEP_TIME, STAMP, OVERWRITE=OVERWRITE

  ;; Check parameters
  IF (n_params() NE 8) THEN $
     message, 'Usage: DAY_SPLITTER, STARTDATE, ENDDATE, IDIR, ' + $
              'ODIR, ITEMPLATE_SAV, TIME_IDX, STEP_TIME, ' + $
              'STAMP, [/OVERWRITE]'
  idir_info=file_info(idir)
  itpl_info=file_info(itemplate_sav)
  IF ((n_elements(startdate) EQ 0) OR (idir EQ '')) THEN $
     message, 'STARTDATE is undefined or is empty string'
  IF ((n_elements(enddate) EQ 0) OR (idir EQ '')) THEN $
     message, 'ENDDATE is undefined or is empty string'
  IF (long(enddate) LT long(startdate)) THEN $
     message, 'ENDDATE must be greater than or equal to STARTDATE'
  IF (~idir_info.directory) THEN $
     message, 'IDIR must be a string pointing to an existing directory'
  IF (~itpl_info.read) THEN $
     message, 'ITEMPLATE_SAV must be a string pointing to a readable file'
  IF ((n_elements(time_idx) NE 1) OR $
      ((size(time_idx, /type) NE 2) || time_idx LT 0)) THEN $
         message, 'TIME_IDX must be an integer scalar >= zero'
  IF ((n_elements(odir) NE 1) OR (size(odir, /type) NE 7)) THEN $
         message, 'ODIR must be a string scalar'
  IF ((n_elements(step_time) EQ 0) OR (step_time LT 0)) THEN $
     message, 'STEP_TIME must be a scalar >= zero'
  IF ((n_elements(stamp) EQ 0) OR (stamp EQ '')) THEN $
     message, 'STAMP is undefined or is empty string'

  idir_files=file_search(idir + path_sep() + '*std*', count=nidir_files, $
                         /nosort, /fold_case, /test_regular)
  IF nidir_files LT 1 THEN message, 'No input files found'

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
  n_recs=86400 / step_time     ; records per day
  restore, itemplate_sav
  n_fields=itemplate.FIELDCOUNT
  field_names=strlowcase(itemplate.FIELDNAMES)
  field_groups=itemplate.FIELDGROUPS
  ;; Times
  is_time_field=field_groups EQ time_idx
  ;; Ignore other groups when reading the data
  itemplate.FIELDGROUPS=indgen(itemplate.FIELDCOUNT)
  itemplate.FIELDGROUPS[where(is_time_field)]=time_idx
  non_time_fields=where(~is_time_field)
  non_time_field_names=field_names[non_time_fields]
  tfields=where(is_time_field, /NULL)
  tnames=field_names[tfields]
  tnamesl=strsplit(tnames, '_', /extract)
  tnames_last=strarr(n_elements(tnamesl))
  tnames_id=strjoin((tnamesl[0])[0:n_elements(tnamesl[0]) - 2], '_')
  FOR i=0L, n_elements(tnames) - 1 DO $
     tnames_last[i]=tnamesl[i, n_elements(tnamesl[i]) - 1]
  ;; Determine where in these names we're supposed to get each time field
  ;; (year, month, day, hour, minute, second, subsecond)
  time_locs=locate_time_strings(tnames_last)

  ;; Add 1 to the last day to set the breaks for value_locate()
  days=timegen(start=beg_jd, final=end_jd + 1, step_size=1, units='days')
  days=jul2timestamp(temporary(days))
  ;; Add up to 23:59:59 of the last requested day
  step_d=float(step_time) / 86400
  times=timegen(start=beg_jd, $
                final=end_jd + 1 - (step_d / 2), $
                step_size=step_time, units='seconds')
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
                                       count=nodir_files, /nosort, $
                                       /fold_case, /test_regular))

  skipdays=intarr(size(days, /n_elements))
  FOR i=0L, n_days, 1L DO BEGIN
     is_ofile=where(strmatch(odir_files, file_stamps[i]), n_ofiles)
     IF n_ofiles GT 0 THEN BEGIN
        IF keyword_set(overwrite) THEN BEGIN
           message, 'Daily ' + stamp + ' file ' + file_stamps[i] + $
                    ' already exists.  Overwriting', /informational
        ENDIF ELSE BEGIN
           message, 'Daily ' + stamp + ' file ' + file_stamps[i] + $
                    ' already exists.  Not overwriting', /informational
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
  ts_all=hash(field_names)           ; just the keys here
  ;; Note we use our time_locs variable to locate each time info
  ts_all[tnames[time_locs[0]]]=strmid(times, 0, 4)
  ts_all[tnames[time_locs[1]]]=strmid(times, 5, 2)
  ts_all[tnames[time_locs[2]]]=strmid(times, 8, 2)
  ts_all[tnames[time_locs[3]]]=strmid(times, 11, 2)
  ts_all[tnames[time_locs[4]]]=strmid(times, 14, 2)
  ;; Becareful here with the decimal places formatting
  ts_all[tnames[time_locs[5]]]=strmid(times, 17)
  ;; Iterate through time fields, and fill the hash with time data and
  ;; holders for the rest. We do not care about matching fractions of
  ;; seconds, for the high frequency data.  It is treated as a non-time
  ;; field.
  FOREACH fld, field_names[where(~is_time_field)] DO BEGIN
     fld_idx=where(field_names EQ fld)
     ifield_type=itemplate.FIELDTYPES[fld_idx]
     val=(ifield_type EQ 4) ? !VALUES.F_NAN : ''
     ts_all[fld]=make_array(n_elements(times), type=ifield_type, $
                            value=val)
  ENDFOREACH
  ;; Separate times from full hash
  ts_times=ts_all.remove(field_names[tfields])

  ;; Read each file
  is_match=intarr(n_elements(times)) ; to check how many matches
  FOR k=0L, nidir_files - 1 DO BEGIN
     ifile=idir_files[k]
     message, 'Processing ' + ifile, /informational
     idata=read_ascii(ifile, count=n_inputfile, template=itemplate)
     ;; Extract times and remove quotes or spaces from strings
     idata_times=idata.(time_idx)
     idata=remove_structure_tags(idata, (tag_names(idata))[time_idx])
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
     ;; Obtain full time details array
     itimes_std=parse_times(idata_times, tnames_last, time_locs)
     n_krecs=n_elements(itimes_std[0, *])
     ;; Check each line and match against array
     is_valid=valid_num(itimes_std[0, *])
     iok=where(is_valid, vcount) ; work on valid data only
     IF vcount LT 1 THEN CONTINUE
     itimes_std=itimes_std[*, iok]
     yyyy=string(itimes_std[0, *], format='(i4)')
     mo=string(itimes_std[1, *], format='(i02)')
     dd=string(itimes_std[2, *], format='(i02)')
     hh=string(itimes_std[3, *], format='(i02)')
     mm=string(itimes_std[4, *], format='(i02)')
     ss=string(itimes_std[5, *], format='(f06.3)')
     i_ts=yyyy + '-' + mo + '-' + dd + ' ' + hh + ':' + mm + ':' + ss
     match2, times, i_ts, times_in_its
     t_matches=where(times_in_its GE 0, mcount, /null)
     IF mcount LT 1 THEN CONTINUE
     is_match[t_matches]=1      ; mark matches
     its_matches=times_in_its[t_matches]
     FOREACH value, ts_all, fld DO BEGIN ; loop non-time hash
        match_fld=where(non_time_field_names EQ fld)
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
     day_idcs=lindgen(n_recs, start=begi)
     day_matches=where(is_match[day_idcs] GT 0, n) ; output only matches
     IF n LT 1 THEN CONTINUE
     file_stamp=file_stamps[begi / n_recs]
     ;; ts=ts_full.toStruct(/no_copy)
     ;; There must be a better way to re-order tags.
     ;; We could use day_idcs[day_matches] instead of begi:endi to output
     ;; only the matches for the day, but the end the average_series.pro
     ;; needs to change.
     ts=create_struct(field_names[0], $
                      ts_full[field_names[0], begi:endi])
     FOREACH fld, field_names[1:*] DO BEGIN
        ts=create_struct(ts, field_names[where(field_names EQ fld)], $
                         ts_full[fld, begi:endi])
     ENDFOREACH
     write_csv, odir + path_sep() + file_stamp, ts, $
                header=strlowcase(tag_names(ts))
  ENDFOR

END

;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; day_splitter.pro ends here
