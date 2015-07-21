### uv_spectral.R --- calibrate UV based on independent data
## Author: Sebastian Luque
## Created: 2014-01-28T16:45:03+0000
## Last-Updated: 2015-04-06T20:55:44+0000
##           By: Sebastian P. Luque
## copyright (c) 2014 Sebastian P. Luque
## 
## This program is Free Software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3, or (at your option)
## any later version.
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
## You should have received a copy of the GNU General Public License
## along with GNU Emacs; see the file COPYING. If not, write to the
## Free Software Foundation, Inc., 59 Temple Place - Suite 330,
## Boston, MA 02111-1307, USA.
## ------------------------------------------------------------------------
### Commentary: 
## 
## ------------------------------------------------------------------------
### Code:

rad <- read.csv("~/Dropbox/CEOS/RAD_geo/Resolute_UV_2010_2011.csv")
rad <- within(rad, {
    uv <- cut(wavelength, breaks=c(0, 315, 400), labels=c("B", "A"))
    timer <- as.POSIXct(as.character(time), format="%H:%M:%S", tz="GMT")
    date.time <- as.POSIXct(strptime(paste(as.character(date),
                                           as.character(time)),
                                     format="%Y-%m-%d %H:%M:%S"),
                            tz="GMT")
})

junk <- subset(rad, date == "2010-04-01" & time == "03:52:31")
junk[] <- lapply(junk, function(x) x[, drop=TRUE])

by(junk, junk$uv, function(x) {
    ## simple linear interpolating function
    irr.linfun <- approxfun(x$wavelength, x$s_irradiance)
    ## Integrate between the corresponding end points
    irr.linint <- integrate(irr.linfun,
                            min(x$wavelength), max(x$wavelength))
    irr.linint
})

"wave.integrate" <- function(x) {
    ## nrow(x)
    ## simple linear interpolating function
    irr.linfun <- approxfun(x$wavelength, x$s_irradiance)
    ## Integrate between the corresponding end points
    irr.linint <- integrate(irr.linfun,
                            min(x$wavelength), max(x$wavelength))
    data.frame(x[1, ], s_irradiance_integral=irr.linint$value)
}


rad.int <- by(rad, list(rad$date, rad$time, rad$uv), wave.integrate)
rad.int <- Filter(Negate(is.null), rad.int)
rad.int <- do.call(rbind, rad.int)
rad.int <- rad.int[order(rad.int$date.time), ]
tt <- by(rad.int, list(rad.int$uv, rad.int$timer), function(x) {
    data.frame(x[1, ], s_irradiance_avg=mean(x$s_irradiance_integral))
})
tt <- do.call(rbind, tt)
xyplot(s_irradiance_avg ~ timer |
       cut(date.time, breaks="year", labels=c("2010", "2011")),
       data=tt, groups=uv, type="p", pch=19, cex=0.2, auto.key=TRUE)

xyplot(s_irradiance_integral ~ timer |
       cut(date.time, breaks="year", labels=c("2010", "2011")),
       data=rad.int, groups=uv, type="p", pch=19, cex=0.2, auto.key=TRUE,
       panel=panel.superpose,
       panel.groups=function(x, y, ...) {
           panel.xyplot(x, y, ...)
           panel.smooth(x, y, ...)
       })


### uv_spectral.R ends here
