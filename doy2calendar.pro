;; $Id$
;; Author: Sebastian Luque
;; Created: 2013-08-28T17:54:17+0000
;; Last-Updated: 2013-10-03T03:03:24+0000
;;           By: Sebastian P. Luque
;; ------------------------------------------------------------------------
;;; Commentary: 
;;
;; This is meant to replace the function 'JD_to_calendar'.
;; ------------------------------------------------------------------------
;;; Code:

FUNCTION DOY2CALENDAR, YEAR, DOY
  
  year_int=fix(year)
  doy_int=fix(doy)
  caldat, julday(1, doy_int, year_int), mm, dd, yyyy
  calendar=strcompress(string(yyyy, format='(i04)') + $
                       string(mm, format='(i02)') + $
                       string(dd, format='(i02)'), /REMOVE_ALL)
  RETURN, calendar

END

FUNCTION CALENDAR2DOY, YEAR, MONTH, DAY
  
  jd=julday(month, day, year)
  caldat, jd, Null0, Null1, year
  doy=string(jd - julday(12, 31, year - 1), format='(i03)')
  RETURN, doy

END

FUNCTION JUL2TIMESTAMP, JUL
  caldat, jul, mo, dd, yyyy, hh, mm, ss
  calendar=strcompress(string(yyyy, format='(i04, "-")') + $
                       string(mo, format='(i02, "-")') + $
                       string(dd, format='(i02, " ")') + $
                       string(hh, format='(i02, ":")') + $
                       string(mm, format='(i02, ":")') + $
                       string(ss, format='(i02)'))
  RETURN, calendar
END



;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; doy2calendar.pro ends here
