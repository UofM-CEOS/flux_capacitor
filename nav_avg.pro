;; $Id$
;;; nav_avg.pro --- calculate averages in NAV files
;; Author: Bruce Johnson, Sebastian Luque
;; Created: 2013-09-17T14:59:07+0000
;; Last-Updated: 2013-10-02T17:42:20+0000
;;           By: Sebastian Luque
;; ------------------------------------------------------------------------
;;; Commentary: 
;; 
;; Example call:
;;
;; nav_avg, expand_path('~/tmp/ArcticNet2011/NAV/Daily'), $
;;          expand_path('~/tmp/ArcticNet2011/NAV/1min'), $
;;          1, 60, 15, 16, 17, 'NAV', expand_path('nav_std_template.sav'), 0
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

PRO NAV_AVG, IDIR, ODIR, ISAMPLE_RATE, OSAMPLE_RATE, BMAG_FIELD, $
             BEAR_FIELD, HEAD_FIELD, STAMP, ITEMPLATE_SAV, TIME_BEG_IDX

  ;; Check parameters
  IF (n_params() NE 10) THEN $
     message, 'Usage: ONE_MIN_NAV, IDIR, ODIR, ISAMPLE_RATE, ' + $
              'OSAMPLE_RATE, BMAG_FIELD, BEAR_FIELD, HEAD_FIELD, ' + $
              'STAMP, ITEMPLATE_SAV, TIME_BEG_IDX'
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
  IF ((n_elements(time_beg_idx) NE 1) OR $
      ((size(time_beg_idx, /type) NE 2) || time_beg_idx LT 0)) THEN $
         message, 'TIME_BEG_IDX must be an integer scalar >= zero'

  idir_files=file_search(idir + path_sep() + '*.dat', count=nidir_files, $
                         /nosort)
  IF nidir_files LT 1 THEN BEGIN
     message, 'No input files found', /informational
     RETURN
  ENDIF

  restore, itemplate_sav
  header=itemplate.FIELDNAMES
  is_time_field=itemplate.FIELDGROUPS EQ time_beg_idx
  time_fields=where(is_time_field, /NULL)
  time_field_names=strlowcase(header[time_fields])
  non_time_fields=where(~is_time_field)
  non_time_field_names=strlowcase(header[non_time_fields])
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

  ang_cols=[bmag_field, bear_field, head_field]
  stdev_str=['sog_stdev', 'cog_stdev', 'heading_stdev']
  avg_cols_maybe=(cgSetDifference(non_time_fields, ang_cols))
  ;; Only average those fields that are not strings (BECAREFUL HERE)
  avg_cols=avg_cols_maybe[where(itemplate.FIELDTYPES[avg_cols_maybe] NE 7)]
  n_avg_cols=size(avg_cols, /n_elements)
  ncols_osample=osample_rate / isample_rate
  nrows_osample=86400 / osample_rate
  oheader=[[time_field_names, header[avg_cols], stdev_str]]
  FOR k=0, nidir_files - 1 DO BEGIN
     ifile=idir_files[k]
     message, 'Producing 1-min average for ' + ifile, /informational
     idata=read_ascii(ifile, count=n_inputfile, template=itemplate)
     ;; Extract times and remove quotes or spaces from strings
     idata_times=idata.(time_beg_idx)
     idata=remove_structure_tags(idata, (tag_names(idata))[time_beg_idx])
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
     beg_jd=julday(idata_times[month_subfield, 0], $
                   idata_times[day_subfield, 0], $
                   idata_times[year_subfield, 0], 0)
     otimes=timegen(start=beg_jd, final=beg_jd + (86399.0 / 86400), $
                    step_size=osample_rate, units='seconds')
     otimes=jul2timestamp(otimes)
     ;; Set up hash that will contain output data
     ohash=hash(oheader)
     FOREACH value, ohash, fld DO BEGIN
        CASE fld OF
           'year': ohash[fld]=strmid(otimes, 0, 4)
           'month': ohash[fld]=strmid(otimes, 5, 2)
           'day': ohash[fld]=strmid(otimes, 8, 2)
           'hour': ohash[fld]=strmid(otimes, 11, 2)
           'minute': ohash[fld]=strmid(otimes, 14, 2)
           'second': ohash[fld]=strmid(otimes, 17, 2)
           ELSE: BEGIN
              fld_idx=where(itemplate.FIELDNAMES EQ fld)
              ifield_type=itemplate.FIELDTYPES[fld_idx]
              val=(ifield_type EQ 4) ? !VALUES.F_NAN : ''
              ohash[fld]=make_array(n_elements(otimes), type=ifield_type, $
                                    value=val)
           END
        ENDCASE
     ENDFOREACH

     FOREACH col, avg_cols DO BEGIN
        col_idata=reform(idata.(where(idata_names EQ header[col])), $
                         ncols_osample, nrows_osample)
        ;; This will throw "floating point illegal operand" arithmetic
        ;; error warning when all elements in each row are NaN, which
        ;; should should just be ignored.  In these cases, the result is
        ;; -NaN.
        avgs=mean(col_idata, dimension=1, /nan)
        ohash[header[col]]=avgs
     ENDFOREACH
     ;; Means for SOG, COG, and heading
     cog=idata.(where(idata_names EQ header[bear_field]))
     brg2d=reform(cog, ncols_osample, nrows_osample)
     sog=idata.(where(idata_names EQ header[bmag_field]))
     bmag2d=reform(sog, ncols_osample, nrows_osample)
     head=idata.(where(idata_names EQ header[head_field]))
     head2d=reform(head, ncols_osample, nrows_osample)
     brgmag_mean=bearing_avg(brg2d, bmag2d, dimension=1)
     ohash[header[bmag_field]]=brgmag_mean[*, 1]
     ohash[header[bear_field]]=brgmag_mean[*, 0]
     head_mean=bearing_avg(head2d, 1, dimension=1)
     ohash[header[head_field]]=head_mean[*, 0]
     ;; Standard deviations
     FOR i=0, nrows_osample - 1 DO BEGIN
        sog_ok=where(finite(bmag2d[*, i]) GT 0, n_sog_ok, complement=sog_bad)
        cog_ok=where(finite(brg2d[*, i]) GT 0, n_cog_ok, complement=cog_bad)
        head_ok=where(finite(head2d[*, i]) GT 0, n_head_ok, $
                      complement=head_bad)
        IF n_sog_ok GT 10 THEN $
           ohash[stdev_str[0], i]=stddev(bmag2d[sog_ok, i])
        IF n_cog_ok GT 10 THEN $
           ohash[stdev_str[1], i]=stddev_yamartino(brg2d[cog_ok, i])
        IF n_head_ok GT 10 THEN $
           ohash[stdev_str[2], i]=stddev_yamartino(head2d[head_ok, i])
     ENDFOR
     
     file_stamp=file_basename(ifile, '.dat') + '_' + $
                strtrim(osample_rate, 2) + 's.dat'
     ts=create_struct(oheader[0], ohash[oheader[0]])
     FOREACH fld, oheader[1:*] DO BEGIN
        ts=create_struct(ts, oheader[where(oheader EQ fld)], ohash[fld])
     ENDFOREACH
     write_csv, odir + path_sep() + file_stamp, ts, header=oheader
  ENDFOR
        
END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; nav_avg.pro ends here
