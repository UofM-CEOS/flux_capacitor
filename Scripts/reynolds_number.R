"reynolds.tube" <- function(diam, rate, kvisc)
{
    ## Value: numeric; Reynolds number (unit-less)
    ## --------------------------------------------------------------------
    ## Arguments: diam=diameter of tube (m); rate=volumetric flow rate
    ## (m3/s) of fluid in tube; kvisc=kinematic viscosity of fluid (m2/s)
    ## --------------------------------------------------------------------
    ## Purpose: Calculate Reynolds number for eddy covariance studies,
    ## where choice of pipe diameter and flow rates are key variables
    ## because we want to achieve turbulent flow (Re > 4000).  We use the
    ## simplified Reynolds number calculation using kinematic viscosity,
    ## rather than the more complex formula using dynamic viscosity and
    ## density of fluid.
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    vel <- (4 / pi) * (rate / (diam ^ 2)) # m/s
    (vel * diam) / kvisc
}

## Some unit conversions
## 1 m3/s = 60000 L/min
## 1 in = 0.0254 m
## Kinematic viscosity of air at sea level, 15C = 14.55e-6 m2/s
kvisc <- 14.55e-6
## Diameter and rate vectors of same length for square grid
rates.lpm <- seq(6, 30)                 # L/min
rates <- rates.lpm / 60000              # m3/s
diameters.in <- seq(1/8, 1/2, 1/64)     # inches
diameters <- diameters.in * 0.0254      # m
Re <- outer(diameters, rates, reynolds.tube, kvisc)
## Plot
pdf("reynolds_dry_air_15C.pdf")
image(diameters.in, rates.lpm, Re, col=heat.colors(100),
      xlab="Tube diameter (in)", ylab="Flow rate (LPM)",
      main="Reynolds # for dry air at sea level, 15C", las=1)
contour(diameters.in, rates.lpm, Re,
        levels=c(seq(1000, 3000, 1000), seq(5000, 13000, 1000)),
        add=TRUE)
contour(diameters.in, rates.lpm, Re, levels=4000, lwd=3, add=TRUE)
dev.off()
