;; $Id$
;; Author: Sebastian P. Luque
;; Created: 2013-10-01T02:19:03+0000
;; Last-Updated: 2013-10-24T21:18:51+0000
;;           By: Sebastian Luque
;;+ -----------------------------------------------------------------------
;; NAME:
;; 
;;     LOCATE_TIME_STRINGS
;; 
;; PURPOSE:
;; 
;;     Given a string array corresponding to field names such as 'year',
;;     'month', 'day', etc., this function returns an integer array with
;;     seven indices indicating where the year, month, day, hour, minute,
;;     seconds, and subseconds are available in the input array.
;; 
;; CALLING SEQUENCE:
;; 
;;     time_locs=locate_time_strings(names)
;; 
;; INPUTS:
;; 
;;     Names:    A string array with the following possible values (case is
;;               irrelevant): 'year', 'month', 'day', 'hour', 'minute',
;;               'second', 'doy', 'yyyymmdd', 'mmddyyyy', 'ddmmyy',
;;               'hhmmss', 'hhmm'.  For subseconds, any string containing
;;               the substring 'second' anywhere after the first character
;;               is elligible.
;; 
;; RESTRICTIONS:
;; 
;; 
;; 
;; PROCEDURE:
;; 
;;     If there are no sub-seconds, then the last element of output array
;;     is -1.
;; 
;; EXAMPLE:
;; 
;;     time_locs=locate_time_strings(['year', 'DOY', 'hhmmss'])
;; 
;;- -----------------------------------------------------------------------
;;; Code:

FUNCTION LOCATE_TIME_STRINGS, NAMES

  ;; Define a vector with possible names
  time_names=['year', $  ; string that can be formatted to 4-digit year
              'month', $ ; string that can be formatted to 2-digit month
              'day', $   ; string that can be formatted to 2-digit day
              'hour', $  ; string that can be formatted to 2-digit hour
              'minute', $       ; string that can be formatted to 2-digit minute
              'second', $       ; string that can be formatted to 2-digit second
              'doy', $          ; string that can be formatted to 3-digit DOY
              ;; String with 4-digit year, 2-digit month, and 2-digit
              ;; day, in that order, with any one-character separator
              ;; between them
              'yyyymmdd', $
              'mmddyyyy', $     ; same with different ordering
              'ddmmyyyy', $     ; same with different ordering
              'ddmmyy', $       ; we have to assume a year
              ;; String with 2-digit hour, 2-digit minute, and 2-digit
              ;; second, in that order, with any one-character separator
              ;; between them
              'hhmmss', $
              'hhmm' $          ; same with different ordering
  ]

  ;; Consider only the last string after underscore (allow for some
  ;; identifiers preceding the string)
  inames=strlowcase(names)
  match2, inames, time_names, a, b
  year_maybe=b[[0, 7, 8, 9, 10]]
  month_maybe=b[[1, 7, 8, 9, 10, 6]]
  day_maybe=b[[2, 7, 8, 9, 10, 6]]
  hour_maybe=b[[3, 11, 12]]
  minute_maybe=b[[4, 11, 12]]
  second_maybe=b[[5, 11]]

  year_candidates=where(year_maybe GE 0, year_count)
  month_candidates=where(month_maybe GE 0, month_count)
  day_candidates=where(day_maybe GE 0, day_count)
  hour_candidates=where(hour_maybe GE 0, hour_count)
  minute_candidates=where(minute_maybe GE 0, minute_count)
  second_candidates=where(second_maybe GE 0, second_count)

  IF year_count LT 1 THEN $
     message, 'Cannot determine year with these strings'
  IF month_count LT 1 THEN $
     message, 'Cannot determine month with these strings'
  IF day_count LT 1 THEN $
     message, 'Cannot determine day with these strings'
  IF hour_count LT 1 THEN $
     message, 'Cannot determine hour with these strings'
  IF minute_count LT 1 THEN $
     message, 'Cannot determine minute with these strings'

  year=year_maybe[year_candidates[0]]
  month=month_maybe[month_candidates[0]]
  day=day_maybe[day_candidates[0]]
  hour=hour_maybe[hour_candidates[0]]
  minute=minute_maybe[minute_candidates[0]]
  second=second_maybe[second_candidates[0]]
  ;; Look for subsecond name
  subsecond_exists=strpos(inames, 'second', 1)
  subsecond=where(subsecond_exists GE 0)
  RETURN, [year, month, day, hour, minute, second, subsecond]

END


;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; locate_time_strings.pro ends here
