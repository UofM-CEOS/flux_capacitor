;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-10-03T19:23:21+0000
;; Last-Updated: 2013-10-24T18:19:13+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     PARSE_TIMES
;; 
;; PURPOSE:
;; 
;;     This function takes a numeric array with columns containing time
;;     data such as year, month, day, etc., along with a string array
;;     corresponding to the names of these columns (possible names are the
;;     same as those used in locate_time_strings()), and an integer array
;;     with the indices corresponding to the location of the fields from
;;     which to extract year, month, day, hour, minute, second, and
;;     possibly sub-second.  It returns a six-column array with these time
;;     data, but where sub-seconds are included as the fractional part of
;;     seconds.
;; 
;; CALLING SEQUENCE:
;; 
;;     times_std=parse_times(times, time_names, locations)
;; 
;; INPUTS:
;; 
;;     Times:         Numeric array where columns represent time data.
;;     Time_Names:    String array with the names of the columns in Times
;;                    (possible names are those used in the function
;;                    locate_time_strings).
;;     Locations:     Integer array with indices corresponding to the
;;                    location of the fields from which to extract year,
;;                    month, day, hour, minute, second, and possibly
;;                    sub-second.
;; 
;; RESTRICTIONS:
;; 
;;     Data must be valid numbers, otherwise this may fail or provide
;;     erroneous output.
;; 
;; PROCEDURE:
;; 
;; 
;; 
;; EXAMPLE:
;; 
;;     times_std=parse_times()
;; 
;;- -----------------------------------------------------------------------
;;; Code:

FUNCTION PARSE_TIMES, TIMES, TIME_NAMES, LOCATIONS

  times_dims=size(times, /dimensions)

  CASE time_names[locations[0]] OF
     'year': $
        yyyy=string(times[locations[0], *], format='(i04)')
     'yyyymmdd': BEGIN
        ystr=string(times[locations[0], *], format='(i08)')
        yyyy=strmid(ystr, 0, 4)
     END
     'mmddyyyy': BEGIN
        ystr=string(times[locations[0], *], format='(i08)')
        yyyy=strmid(ystr, 4, 4)
     END
     'ddmmyyyy': BEGIN
        ystr=string(times[locations[0], *], format='(i08)')
        yyyy=strmid(ystr, 4, 4)
     END
     'ddmmyy': BEGIN
        message, 'Assuming current century', /informational
        tstamp=jul2timestamp(systime(/julian))
        ystr=string(times[locations[0], *], format='(i06)')
        yyyy=strmid(tstamp, 0, 2) + strmid(ystr, 4, 2)
     END
     ELSE: message, 'Do not know how to extract year from this field'
  ENDCASE
  CASE time_names[locations[1]] OF
     'month': mo=string(times[locations[1], *], format='(i02)')
     'yyyymmdd': BEGIN
        mostr=string(times[locations[1], *], format='(i08)')
        mo=strmid(mostr, 4, 2)
     END
     'mmddyyyy': BEGIN
        mostr=string(times[locations[1], *], format='(i08)')
        mo=strmid(mostr, 0, 2)
     END
     'ddmmyyyy': BEGIN
        mostr=string(times[locations[1], *], format='(i08)')
        mo=strmid(mostr, 2, 2)
     END
     'ddmmyy': BEGIN
        mostr=string(times[locations[1], *], format='(i06)')
        mo=strmid(mostr, 2, 2)
     END
     'doy': BEGIN
        calendar=doy2calendar(yyyy, times[locations[1], *])
        mo=strmid(calendar, 4, 2)
     END
     ELSE: message, 'Do not know how to extract month from this field'
  ENDCASE
  CASE time_names[locations[2]] OF
     'day': dd=string(times[locations[2], *], format='(i02)')
     'yyyymmdd': BEGIN
        dstr=string(times[locations[2], *], format='(i08)')
        dd=strmid(dstr, 6, 2)
     END
     'mmddyyyy': BEGIN
        dstr=string(times[locations[2], *], format='(i08)')
        dd=strmid(dstr, 2, 2)
     END
     'ddmmyyyy': BEGIN
        dstr=string(times[locations[2], *], format='(i08)')
        dd=strmid(dstr, 0, 2)
     END
     'ddmmyy': BEGIN
        dstr=string(times[locations[2], *], format='(i06)')
        dd=strmid(dstr, 0, 2)
     END
     'doy': dd=strmid(calendar, 6, 2) ; we already have calendar
     ELSE: message, 'Do not know how to extract day from this field'
  ENDCASE
  CASE time_names[locations[3]] OF
     'hour': hh=string(times[locations[3], *], format='(i02)')
     'hhmmss': BEGIN
        hstr=string(times[locations[3], *], format='(i06)')
        hh=strmid(hstr, 0, 2)
     END
     'hhmm': BEGIN
        hstr=string(times[locations[3], *], format='(i04)')
        hh=strmid(hstr, 0, 2)
     END
     ELSE: message, 'Do not know how to extract hour from this field'
  ENDCASE
  CASE time_names[locations[4]] OF
     'minute': mm=string(times[locations[4], *], format='(i02)')
     'hhmmss': BEGIN
        mstr=string(times[locations[4], *], format='(i06)')
        mm=strmid(mstr, 2, 2)
     END
     'hhmm': BEGIN
        mstr=string(times[locations[4], *], format='(i04)')
        mm=strmid(mstr, 2, 2)
     END
     ELSE: message, 'Do not know how to extract minute from this field'
  ENDCASE
  IF locations[5] LT 0 THEN BEGIN ; assume 0 seconds if no seconds found
     ss=string(rebin([0], times_dims[1]), format='(f06.3)')
  ENDIF ELSE BEGIN
     CASE time_names[locations[5]] OF
        'second': ss=string(times[locations[5], *], format='(f06.3)')
        ;; Becareful with the format here - we need to know how much this
        ;; varies in input files
        'hhmmss': BEGIN
           sstr=string(times[locations[5], *], format='(f010.3)')
           ss=strmid(sstr, 4)
        END
        ELSE: message, 'Do not know how to extract second from this field'
     ENDCASE
  ENDELSE
  IF locations[6] GE 0 THEN $   ; concatenate if we have fractional ss
     ss=string(temporary(ss), format='(i02)') + '.' + $
        string(times[locations[6], *], format='(i03)')

  otimes=[reform(yyyy, 1, times_dims[1]), $
          reform(mo, 1, times_dims[1]), $
          reform(dd, 1, times_dims[1]), $
          reform(hh, 1, times_dims[1]), $
          reform(mm, 1, times_dims[1]), $
          reform(ss, 1, times_dims[1])]

  RETURN, otimes

END



;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; parse_times.pro ends here
