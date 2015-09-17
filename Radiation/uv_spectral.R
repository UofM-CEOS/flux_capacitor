### uv_spectral.R --- calibrate UV based on independent data
## Author: Sebastian Luque
## Created: 2014-01-28T16:45:03+0000
## Last-Updated: 2015-09-17T20:06:53+0000
##           By: Sebastian P. Luque
## copyright (c) 2014-2015 Sebastian P. Luque
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
## Analyses for Virgine Galindo on UVa and UVb radiation.
## ------------------------------------------------------------------------
### Code:

rad <- read.csv("~/Dropbox/CEOS/RAD_geo/Resolute_UV_2010_2011.csv")
rad <- within(rad, {
    uv <- cut(wavelength, breaks=c(0, 315, 400), labels=c("B", "A"))
    ## We know they're all negative, so ignore the "-", convert to seconds,
    ## and then add offset
    UTC_offset <- substring(as.character(UTC_offset), 2)
    UTC_offset.s <- (as.numeric(substr(UTC_offset, 1, 2)) * 3600) +
        (as.numeric(substr(UTC_offset, 4, 5)) * 60) +
        as.numeric(substr(UTC_offset, 7, 8))
    ## The time provided, as POSIXct
    solar_time <- as.POSIXct(paste(date, time),
                             format="%Y-%m-%d %H:%M:%S", tz="GMT")
    ## Add the offset seconds to express as UTC
    UTC <- solar_time + UTC_offset.s
    ## Only the time portion (date arbitrarily set to today)
    UTC_time <- as.POSIXct(format(UTC, format="%H:%M:%S", tz="GMT"),
                           format="%H:%M:%S", tz="GMT")
})
rad$SO2 <- rad$Err_SO2 <- NULL          # nothing here
xyplot(log(s_irradiance) ~ wavelength, data=rad,
       subset=date == "2011-06-21" & time == "11:59:56",
       type=c("p", "smooth"), degree=2)
plot(s_irradiance ~ UTC, data=rad, subset=date == "2011-06-21" & uv == "A")

junk <- subset(rad, date == "2011-06-21" & time == "11:59:56")
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
    ## Build a simple linear interpolating function for lambda <= 315
    irr.linfun <- approxfun(x$wavelength, x$s_irradiance)
    ## Integrate between the ends of the data we do have
    irr.linint <- integrate(irr.linfun, min(x$wavelength),
                            max(x$wavelength))
    ## Fit a smoothing spline function.  We'll use this for lambda > 315
    irr.spline <- smooth.spline(x$wavelength, x$s_irradiance, spar=0.6,
                                all.knots=TRUE)
    irr.splinefun <- function(y) predict(irr.spline, y)$y
    ## Integrate between the corresponding end points.  We extrapolate for
    ## UVa, because we don't have anything above 315
    lambda.max <- ifelse(unique(x$uv) == "B", max(x$wavelength), 400)
    ## ## Below just for debugging/checking
    ## irr.hat <- predict(irr.spline,
    ##                    seq(min(x$wavelength), lambda.max, by=0.5))
    ## plot(x$wavelength, x$s_irradiance)
    ## lines(irr.hat, col="gray")
    ## plot(x$wavelength, x$s_irradiance,
    ##      xlim=c(min(x$wavelength), lambda.max),
    ##      ylim=c(min(x$s_irradiance), max(irr.hat$y)))
    ## lines(irr.hat, col="gray")
    ## Integrate from upper end to 400 if on UVa; this is zero if on UVb
    irr.splint <- integrate(irr.splinefun, max(x$wavelength), lambda.max)
    ## Sum the integrals from availale data and estimated extrapolation.
    irr.totint <- irr.linint$value + irr.splint$value
    data.frame(x[1, ], s_irradiance_integral=irr.totint)
}

rad.int <- by(rad, list(rad$uv, rad$solar_time), wave.integrate)
## rad.int <- Filter(Negate(is.null), rad.int)
rad.int <- do.call(rbind, rad.int)
rad.int <- rad.int[order(rad.int$UTC), ]
rad.int[] <- lapply(rad.int, function(x) x[, drop=TRUE])
## xyplot(s_irradiance ~ wavelength, data=rad,
##        subset=date == "2010-04-01" & time == "03:52:31",
##        type=c("p", "smooth"), degree=2)

trellis.device(pdf, file="daily_irradiance.pdf", width=11, height=5,
               color=FALSE)
par.sets=list(superpose.symbol=list(col=c("black", "gray")))
for (d in as.character(levels(rad.int$date))) {
    if (nrow(subset(rad.int, date == d)) > 2) {
        pp <- xyplot(s_irradiance_integral ~
                     as.POSIXct(format(solar_time, format="%H:%M:%S"),
                                format="%H:%M:%S"), data=rad.int,
                     subset=date == d, groups=uv, type="b", pch=19,
                     par.settings=par.sets,
                     xlim=range(as.POSIXct(format(rad.int$solar_time,
                       format="%H:%M:%S"), format="%H:%M:%S")),
                     scales=list(alternating=1, rot=0, tck=c(0.5, 0),
                       x=list(format="%H:%M:%S")),
                     xlab=paste("Solar time (", d, ")", sep=""),
                     ylab=expression(paste("S-Irradiance integral ",
                         (W / m ^ 2))),
                     panel=panel.superpose,
                     panel.groups=function(x, y, ...) {
                         panel.xyplot(x, y, ...)
                     },
                     key=list(text=list(levels(rad.int$uv)),
                       points=list(pch=19, col=c("black", "gray")),
                       lines=Rows(trellis.par.get("superpose.line"),
                         c(1:2)), columns=2))
        print(pp)
    }
}
dev.off()

trellis.device(pdf, file="irradiance.pdf", width=11, height=5, color=FALSE)
xyplot(s_irradiance_integral ~
       as.POSIXct(format(solar_time, format="%H:%M:%S"),
                  format="%H:%M:%S") |
       cut(solar_time, breaks="year", labels=c("2010", "2011")),
       data=rad.int, groups=uv, pch=".", lty=2, col.line=c("red", "blue"),
       par.settings=list(superpose.symbol=list(col=c("black", "gray"))),
       scales=list(alternating=1, rot=0, tck=c(0.5, 0),
         x=list(format="%H:%M:%S")),
       xlab="Solar time",
       ylab=expression(paste("S-Irradiance integral ", (W / m ^ 2))),
       panel=panel.superpose,
       panel.groups=function(x, y, lwd, ...) {
           panel.xyplot(x, y, ...)
           panel.loess(x, y, span=0.2, degree=1, lwd=2, ...)
       },
       key=list(text=list(levels(rad.int$uv)),
         points=list(pch="."),
         lines=list(lty=2, col=c("red", "blue")), columns=2))
dev.off()

## Data we measured
rad.ours <- read.csv("~/Data/UV_Spectral/Rad_Res2010_work.csv",
                     na.strings=c("-99999"))
rad.ours <- subset(rad.ours, select=-c(TC_1:TC_25))
## Strangely, and as was done in 2011, these data have "2400" as the last
## time stamp of the day, which doesn't work for conversion to POSIXct
## times, so we fix it
rad.ours$hhmm <- sprintf("%04d", rad.ours$hhmm)
bad.doy <- which(rad.ours$hhmm == "2400")
rad.ours$JulDay[bad.doy] <- rad.ours$JulDay[bad.doy] + 1
rad.ours$hhmm[bad.doy] <- "0000"
rad.ours <- within(rad.ours, {
    UV.A <- UV.A * 0.70266              # re-scaling as per Tim's email
    UV.B <- UV.B * 0.7345
    date <- format(strptime(paste(year, JulDay), format="%Y %j"))
    UTC <- as.POSIXct(paste(date, substr(hhmm, 1, 2), substr(hhmm, 3, 4)),
                      format="%Y-%m-%d %H %M", tz="GMT")
    UTC_time <- as.POSIXct(format(UTC, format="%H:%M:%S", tz="GMT"),
                           format="%H:%M:%S", tz="GMT")
    local.time <- UTC - (60 * 60 * 5)   # UTC - 5
    days <- cut(local.time, "day")
})
rad.ours <- rad.ours[order(rad.ours$UTC), ]

save(rad, rad.int, rad.ours, file="rad.RData") # to save time now that
                                               # we've got the process
                                               # streamlined...

plot(UV.B ~ local.time, data=rad.ours,
     subset=days == levels(days)[2],
     cex=0.2, type="l")

trellis.device(pdf, file="daily_rad_ours.pdf", width=11, height=5,
               color=FALSE)
par.sets=list(superpose.symbol=list(col=c("black", "gray")))
for (d in levels(rad.ours$days)) {
    if (nrow(subset(rad.ours, days == d)) > 2) {
        pp <- xyplot(UV.A + UV.B ~
                     as.POSIXct(format(local.time, format="%H:%M:%S"),
                                format="%H:%M:%S"), data=rad.ours,
                     subset=days == d, type="b", pch=".",
                     par.settings=par.sets,
                     xlim=range(as.POSIXct(format(rad.ours$local.time,
                       format="%H:%M:%S"), format="%H:%M:%S")),
                     scales=list(alternating=1, rot=0, tck=c(0.5, 0),
                       x=list(format="%H:%M:%S")),
                     xlab=paste("Local time (", d, ")", sep=""),
                     ylab=expression(paste("Radiation ",
                         (W / m ^ 2))),
                     panel=panel.superpose,
                     panel.groups=function(x, y, ...) {
                         panel.xyplot(x, y, ...)
                     },
                     key=list(text=list(c("B", "A")),
                       points=list(pch=19, col=c("black", "gray")),
                       lines=Rows(trellis.par.get("superpose.line"),
                         c(1:2)), columns=2))
        print(pp)
    }
}
dev.off()

## Plot both time series to compare

## UV B
pdf("rad_series.pdf", width=17, height=5)
## plot(UV.B ~ local.time, data=rad.ours, type="l", ylim=c(0, 2),
##      xlab="Local (CEOS) or solar (government) time",
##      ylab=expression(paste("UV B ", (W / m ^ 2))))
plot(UV.B ~ UTC, data=rad.ours, type="l", ylim=c(0, 2),
     xlab="UTC time",
     ylab=expression(paste("UV B ", (W / m ^ 2))))
## with(subset(rad.int, uv == "B" & cut(solar_time, "year") == "2010-01-01"),
##      lines(UTC, s_irradiance_integral, type="l",
##            col="red"))
with(subset(rad.int, uv == "B" & cut(UTC, "year") == "2010-01-01"),
     lines(UTC, s_irradiance_integral, type="l",
           col="red"))
legend("topright", c("CEOS", "Government"), lty=1, col=c("black", "red"),
       bty="n")
## UV A
## plot(UV.A ~ local.time, data=rad.ours, type="l", ylim=c(0, 80),
##      xlab="Local (CEOS) or solar (government) time",
##      ylab=expression(paste("UV A ", (W / m ^ 2))))
plot(UV.A ~ UTC, data=rad.ours, type="l", ylim=c(0, 80),
     xlab="Local (CEOS) or solar (government) time",
     ylab=expression(paste("UV A ", (W / m ^ 2))))
## with(subset(rad.int, uv == "A" & cut(solar_time, "year") == "2010-01-01"),
##      lines(solar_time, s_irradiance_integral, type="l",
##            col="red"))
with(subset(rad.int, uv == "A" & cut(UTC, "year") == "2010-01-01"),
     lines(UTC, s_irradiance_integral, type="l",
           col="red"))
legend("topright", c("CEOS", "Government"), lty=1, col=c("black", "red"),
       bty="n")
dev.off()

## Now calculate the differences when PAR < 1 (night time).  We might have
## to merge the two data sets

load("rad.RData")

par.smooth <- with(subset(rad.ours, !is.na(PAR)),
                   smooth.spline(UTC, PAR, spar=0.3))
pdf("par_swdown_series.pdf", width=17, height=5)
plot(PAR ~ UTC, data=rad.ours, type="l", ylim=c(0, 300))
lines(par.smooth, col="red")
plot(predict(par.smooth, as.numeric(rad.ours$UTC), deriv=1), type="l")
plot(predict(par.smooth, as.numeric(rad.ours$UTC), deriv=2), type="l",
      col="red")
lines(predict(par.smooth, as.numeric(rad.ours$UTC), deriv=1), type="l")
with(subset(rad.ours, !is.na(PAR)),
     lines(smooth.spline(UTC, PAR, spar=0.3), col="red"))
plot(SWdown ~ UTC, data=rad.ours, type="l", ylim=c(0, 200))
dev.off()

rad.nopar <- subset(rad.ours, PAR < 1)  # only 1 record below 1
summary(rad.ours$PAR)
summary(rad.ours$SWdown)
rad.nokd <- subset(rad.ours, SWdown < 1)

## Let's just "eyeball" it...

## 2011
rad.ours2011 <- read.csv("~/Data/UV_Spectral/RES2011_RAD_april_may20.csv")
rad.ours2011 <- subset(rad.ours2011, select=-c(Ice.T220:Snow.T21))
rad.ours2011$hhmm <- sprintf("%04d", rad.ours2011$hhmm)
bad.doy <- which(rad.ours2011$hhmm == "2400")
rad.ours2011$JulDay[bad.doy] <- rad.ours2011$JulDay[bad.doy] + 1
rad.ours2011$hhmm[bad.doy] <- "0000"
rad.ours2011 <- within(rad.ours2011, {
    ## UV.A <- UV.A * 0.70266              # re-scaling as per Tim's email
    ## UV.B <- UV.B * 0.7345
    date <- format(strptime(paste(Year, JulDay), format="%Y %j"))
    UTC <- as.POSIXct(paste(date, substr(hhmm, 1, 2), substr(hhmm, 3, 4)),
                      format="%Y-%m-%d %H %M", tz="GMT")
    UTC_time <- as.POSIXct(format(UTC, format="%H:%M:%S", tz="GMT"),
                           format="%H:%M:%S", tz="GMT")
    local.time <- UTC - (60 * 60 * 5)   # UTC - 5
    days <- cut(local.time, "day")
})
rad.ours2011 <- rad.ours2011[order(rad.ours2011$UTC), ]

pdf("uvb_offsets.pdf", width=17, height=5)
op <- par(no.readonly = TRUE) # the whole list of settable par's.
plot(I(UV.B - 0.215) ~ UTC, data=rad.ours, type="l", ylim=c(0, 0.5),
     main="Offset: 0.215", xlab="UTC time (2010)",
     ylab=expression(paste("UV B ", (W / m ^ 2))), las=1)
with(subset(rad.int, uv == "B" & cut(UTC, "year") == "2010-01-01"),
     lines(UTC, s_irradiance_integral, type="l",
           col="red"))
par(xpd=TRUE)
inset <- c(0, -0.18)
legend("topright", c("CEOS", "Government"), lty=1, col=c("black", "red"),
       bty="n", inset=inset)
par(op)
plot(I(UV.B - 0.3) ~ UTC, data=rad.ours, type="l", ylim=c(-0.1, 0.5),
     main="Offset: 0.3", xlab="UTC time (2010)",
     ylab=expression(paste("UV B ", (W / m ^ 2))), las=1)
with(subset(rad.int, uv == "B" & cut(UTC, "year") == "2010-01-01"),
     lines(UTC, s_irradiance_integral, type="l",
           col="red"))
par(xpd=TRUE)
legend("topright", c("CEOS", "Government"), lty=1, col=c("black", "red"),
       bty="n", inset=inset)
par(op)
plot(I(UV.B) ~ UTC, data=rad.ours2011, type="l", ylim=c(0, 0.5),
     main="Offset: 0", xlab="UTC time (2011)",
     ylab=expression(paste("UV B ", (W / m ^ 2))), las=1)
with(subset(rad.int, uv == "B" & cut(UTC, "year") == "2011-01-01"),
     lines(UTC, s_irradiance_integral, type="l",
           col="red"))
par(xpd=TRUE)
inset <- c(0, -0.18)
legend("topright", c("CEOS", "Government"), lty=1, col=c("black", "red"),
       bty="n", inset=inset)
par(op)
dev.off()

## Daily average from government data
rad.int.daily <- by(rad.int,
                    list(cut(rad.int$solar_time, "days"), rad.int$uv),
                    function(x) {
                        irr.mu <- mean(x$s_irradiance_integral, na.rm=TRUE)
                        irr.sd <- sd(x$s_irradiance_integral, na.rm=TRUE)
                        data.frame(solar_date=unique(trunc(x$solar_time,
                                     "days")),
                                   uv=unique(x$uv), s_irradiance_mu=irr.mu,
                                   s_irradiance_sd=irr.sd)
                    })
rad.int.daily <- do.call(rbind, rad.int.daily)

write.csv(rad.int, file="uv_integrals.csv", row.names=FALSE)
write.csv(rad.int.daily, file="uv_daily_means.csv", row.names=FALSE)


### uv_spectral.R ends here
