;; $Id$
;;; nav_avg.pro --- calculate averages in NAV files
;; Author: Bruce Johnson, Sebastian Luque
;; Created: 2013-09-17T14:59:07+0000
;; Last-Updated: 2013-09-18T19:22:24+0000
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

PRO nav_avg, IDIR, ODIR, ISAMPLE_RATE, OSAMPLE_RATE, VEC_FIELD, $
             BEAR_FIELD, HEAD_FIELD, STAMP, ITEMPLATE_SAV, HEADER

  ;; Check parameters
  IF (n_params() NE 10) THEN $
     message, 'Usage: ONE_MIN_NAV, IDIR, ODIR, VEC_FIELD, BEAR_FIELD' + $
              'HEAD_FIELD, STAMP, ITEMPLATE_SAV, HEADER'
  IF ((n_elements(idir) EQ 0) OR (idir EQ '')) THEN $
     message, 'IDIR is undefined or is empty string'
  IF ((n_elements(odir) EQ 0) OR (odir EQ '')) THEN $
     message, 'ODIR is undefined or is empty string'
  IF ((n_elements(vec_field) EQ 0) OR (vec_field EQ '')) THEN $
     message, 'VEC_FIELD is undefined or is empty string'
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

  idir_files=file_search(idir + path_sep() + '*.dat', $
                         count=nidir_files, /nosort)
  IF nidir_files LT 1 THEN BEGIN
     print, 'No input files found'
     RETURN
  ENDIF

  restore, itemplate_sav
  n_fields=nav_template.FIELDCOUNT
  ang_cols=[bear_field, head_field]
  avg_cols=(cgSetDifference(indgen(n_fields - 13) + 13, ang_cols))
  n_avg_cols=size(avg_cols, /n_elements)
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
        col_data=reform(data.(col), osample_rate / isample_rate, $
                        86400 / osample_rate)
        ;; This will throw "floating point illegal operand" arithmetic
        ;; error warning when all elements in each row are NaN, which
        ;; should should just be ignored.  In these cases, the result is
        ;; -NaN.
        avgs=mean(col_data, dimension=1, /nan)
        all_arr[col - 6, *]=avgs
     ENDFOREACH

     ;; FOR fld=13L, n_fields - 1 DO BEGIN
     ;;    print, data.(fld)
     ;;    ;; IF fld EQ bear_field THEN BEGIN
     ;;    ;;    bearings=data.(fld)
     ;;    ;; ENDIF
     ;; ENDFOR
  ENDFOR
;; ;create loops to calculate 60s averages
;; ;--------------------------------------
;;      FOR y=0L, 1440-1L DO BEGIN
;;         st=y * 60
;;         en=st + 59
;;         FOR f=13, n_fields-1 DO BEGIN
           
;;                                 ;calculate averages for COG and SOG
;;            IF f EQ bear_field THEN BEGIN
;;               bear_data=data.FIELD01(f,st:en)
;;               vec_data =1
              
;;               IF vec_field GT 0 THEN BEGIN
;;                  vec_data=data.FIELD01(vec_field,st:en)
;;               ENDIF
              
;;               values=bearing_avg(bear_data, vec_data)
;;               work_arr[bear_field-6,y)=values(0,0)
              
;;               IF vec_field GT 0 THEN BEGIN
;;                  work_arr[vec_field-6,y)=values(0,1)
;;               ENDIF
;;               GOTO, skip
;;            ENDIF
           
;;                                 ;calculate average for Heading
;;            IF f EQ head_field THEN BEGIN
;;               head_data=data.FIELD01(f,st:en)
;;               vec_data =1
;;               values=bearing_avg(head_data, vec_data)
;;               work_arr[head_field-6,y)=values(0,0)
;;               GOTO, skip
;;            ENDIF

;;            work_arr[f-6,y)=mean(data.FIELD01(f,st:en), /NAN)
;;            skip:
;;         ENDFOR
        
;;         sog=data.FIELD01(15,st:en)
;;         cog=data.FIELD01(16,st:en)
;;         head=data.FIELD01(17,st:en)
;;         good_sog=where(finite(sog) GT 0, n_good_sog)
;;         good_cog=where(finite(cog) GT 0, n_good_cog)
;;         good_head=where(finite(head) GT 0, n_good_head)

;;         IF n_good_sog GT 10 THEN BEGIN
;;            work_arr[17,y)=stdev(sog[good_sog])
;;         ENDIF

;;         IF n_good_cog GT 10 THEN BEGIN
;;            work_arr[18,y)=yamartino(cog[good_cog])
;;         ENDIF

;;         IF n_good_head GT 10 THEN BEGIN
;;            work_arr[19,y)=yamartino(head[good_head])
;;         ENDIF

;;      ENDFOR

;; ;create a name for the output file
;; ;---------------------------------
;;      current_JD  =LONG(calendar_to_jd(work_arr[0,0),work_arr[1,0),work_arr[2,0)))
;;      current_date=jd_to_calendar(long(work_arr[0,0)),long(current_JD))

;;      IF current_JD LT 100 THEN BEGIN
;;         JD_stamp=strcompress(current_JD)
;;         JD_stamp='0'+JD_stamp
;;         IF current_JD LT 10 THEN BEGIN
;;            JD_stamp=strcompress(current_JD)
;;            JD_stamp='00'+JD_stamp
;;         ENDIF
;;      ENDIF ELSE BEGIN
;;         JD_stamp  =strcompress(current_JD, /REMOVE_ALL)
;;      ENDELSE

;;      date_stamp=strmid(current_date,4,4)
;;      year_stamp=STRCOMPRESS(LONG(work_arr[0,0)))

;; ;create a file name for the output file with JD stamp and date stamp
;; ;format for NAV output file name --> 1min_NAV_yyyy_JDxxx_hhmm.dat
;;      outname=STRCOMPRESS(odir + '1min_' + stamp + '_' + year_stamp + '_JD' + JD_stamp + '_' + date_stamp + '.dat', /REMOVE_ALL)

;;      formatnum   =strcompress(n_fields-4)
;;      formatstring=strcompress('(' + formatnum + '(a,","),a)', /REMOVE_ALL)

;; ;make the array a string for easy printing, and DO IT
;;      work_str=strcompress(work_arr, /REMOVE_ALL)
;;      header  ='Year, Month, Day, Hour, Minute, Second, ProgVers, Latitude, Longitude, SOG(kts), COG(deg), Heading, Pitch, ' +$
;;               'Roll, Accelx, Accely, Accelz, SOG_stdev, COG_stdev, Heading_stdev'
;;      print, 'Creating file:  ' + outname
;;      print, '---------------------------'

;;      openw, unit1, outname, /get_lun
;;      printf, unit1, header
;;      printf, unit1, work_arr, FORMAT=formatstring

;; skiploop:
;;      close, /all
;;   ENDFOR

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; nav_avg.pro ends here
