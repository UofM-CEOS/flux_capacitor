"""Processing of underway pCO2 system data."""

import os.path as osp
import glob
import numpy as np
import pandas as pd
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
plt.style.use('ggplot')

idir="/home/sluque/Data/ArcticNet/2010/FromDB/LowFreq_1w30min"
iheader = ["time_30min", "time_study", "longitude", "longitude_30min",
           "latitude", "latitude_30min", " atmospheric_pressure",
           "atmospheric_pressure_30min", "air_temperature",
           "air_temperature_30min", "relative_humidity",
           "relative_humidity_30min", "surface_temperature",
           "surface_temperature_30min", "true_wind_speed",
           "true_wind_direction", "PAR", "PAR_30min", "K_down",
           "K_down_30min", "LW_down", "LW_down_30min",
           "UV_sensor_temperature", "UV_sensor_temperature_30min",
           "UV_b", "UV_b_30min", "UV_a", "UV_a_30min", "UV_broad",
           "UV_broad_30min", "equ_temperature", "uw_CO2_cellA",
           "uw_CO2_cellB", "uw_CO2_fraction", "uw_H2O_cellA",
           "uw_H2O_cellB", "uw_H2O_fraction", "uw_temperature_analyzer",
           "uw_pressure_analyzer", "equ_pressure", "H2O_flow",
           "air_flow_analyzer", "equ_speed_pump", "ventilation_flow",
           "condensation_atm", "condensation_equ", "drip_1", "drip_2",
           "condenser_temperature", "temperature_dry_box", "ctd_pressure",
           "ctd_temperature", "ctd_conductivity", "ctd_O2_saturation",
           "ctd_O2_concentration", "uw_pH", "uw_redox_potential",
           "temperature_external", "temperature_in"]
odir = osp.join(idir, "Processed")
ifiles = glob.glob(osp.join(idir, "BULK*.csv"))
gasR = 8.31451                    # J/mol/K universal gas constant
salinity_default = 30.0
use_default_salinity = True
anemometer_height = 16          # in m

for ifile in ifiles[0:10]:
        # Get a file name prefix to be shared by the output files from this
        # period.  Note iname is THE SAME AS THE INDEX IN OSUMMARY
        iname = osp.basename(ifile)
        iname_prefix = osp.splitext(iname)[0]
        # Read into data frame
        bulk = pd.read_csv(ifile, header=None, names=iheader)
        if bulk.shape[0] < 30:
                continue

        # Convert from molar concentration to dry air mixing ratio - for
        # both underway and MET data
        
        # Underway - seawater

        # Bulk air molar concentration, CO2 molar concentration, water
        # vapor molar concentration
        air_sw = ((bulk.uw_pressure_analyzer * 100.0) /
                  (gasR * (bulk.uw_temperature_analyzer + 273.15)))
        CO2_sw = (bulk.uw_CO2_fraction / 1.0e6) * air_sw
        H2O_sw = (bulk.uw_H2O_fraction / 1000.0) * air_sw
        # Dry air mixing ratio
        pCO2_sw = (CO2_sw / (air_sw - H2O_sw)) * 1.0e6

        # MET - atmosphere

        # Bulk air molar concentration, CO2 molar concentration, water
        # vapor molar concentration
        air_atm = ((bulk.cp_pressure * 100.0) /
                   (gasR * (bulk.cp_temperature + 273.15)))
        CO2_atm = (bulk.cp_CO2_fraction / 1.0e6) * air_atm
        H2O_sw = (bulk.cp_H2O_fraction / 1000.0) * air_atm
        # Dry air mixing ratio
        pCO2_atm = (CO2_atm / (air_atm - H2O_atm)) * 1.0e6
        # Extrapolate wind to 10 m height; using eq. 4.14b from Stull -
        # Meteorology for Scientists and Engineers
        wind_speed10 = (bulk.true_wind_speed * log(10 / 0.0002) /
                        log(anemometer_height / 0.0002))

        # Calculate salinity - from: UNESCO Technical Papers in Marine
        # Science #44: Algorithms for computation of fundamental properties
        # of seawater
        ctd_pressure_dbar = bulk.ctd_pressure / 10.0 # dbar
        # This is the conductivity of seawater at S=35,T=15,p=0 (hard
        # number to dig up)
        bigR = bulk.ctd_conductivity / 42.914

        # Define coefficients
        a0, a1, a2, a3, a4, 15 = (0.008, -0.1692, 25.3851,
                                  14.0941, -7.0261, 2.7081)
        b0, b1, b2, b3, b4, b5 = (0.0005, -0.0056, -0.0066,
                                  -0.0375, 0.0636, -0.0144)
        c0, c1, c2, c3, c4 = (0.6766097, 2.00564e-2, 1.104259e-4,
                              -6.9698e-7, 1.0031e-9)
        d1, d2, d3, d4 = (3.426e-2, 4.464e-4, 4.215e-1, -3.107e-3)
        e1, e2, d3 = (2.070e-5, -6.370e-10, 3.989e-15)
        K = 0.0162

        # Salinity calculation
        rt = (c0 + (c1 * bulk.ctd_temperature) +
              (c2 * bulk.ctd_temperature ** 2.0) +
              (c3 * bulk.ctd_temperature ** 3.0) +
              (c4 * bulk.ctd_temperature ** 4.0)) # [eqn 3]
        Rp = 1.0 + ((ctd_pressure_dbar *
                     (e1 + e2 * ctd_pressure_dbar +
                      e3 * ctd_pressure_dbar ** 2.0)) /
                    (1.0 + d1 * bulk.ctd_temperature +
                     d2 * bulk.ctd_temperature ** 2.0 +
                     (d3 + d4 * bulk.ctd_temperature) * bigR)) #  [eqn 4]
        bigRt = bigR / (Rp * rt) # eqn after eq'n 4
        delS = (((bulk.ctd_temperature - 15.0) /
                 (1.0 + K * (bulk.ctd_temperature - 15.0))) *
                (b0 + b1 * bigR ** (1 / 2.0) + b2 * bigR + b3 *
                 bigR ** (3 / 2.0) + b4 * bigR ** 2.0 + b5 *
                 bigR ** (5 / 2.0))) # [eqn 2]
        sal = (a0 + a1 * bigRt ** (1 / 2.0) + a2 * bigRt + a3 *
               bigRt ** (3 / 2.0) + a4 * bigRt ** 2.0 + a5 *
               bigRt ** (5 / 2.0) + delS)

        # Saturated pCO2eq calculation
        # (standardized to 1 atm, not pressure from MET tower)

        # For now, we'll use H2O T (not tank T)... will resolve that
        # problem later.
        T_eq = bulk.equ_temperature + 273.15
        T_air = bulk.air_temperature + 273.15
        PW_sw = (np.exp(24.4543 - 67.4509 * (100.0 / T_eq) - 4.8489 *
                        log(T_eq / 100.0) - 0.000544 * sal))
        PW_atm = (np.exp(24.4543 - 67.4509 * (100.0 / T_air) - 4.8489 *
                         alog(T_air / 100.0) - 0.000544 * sal))
        # Replace dry air mixing ration with pCO2eq
        pCO2_sw = pCO2_sw * (1 - PW_sw) # units: uatm, here 1 -> 1 atm
        pCO2_atm = pCO2_atm * (1 - PW_atm) # units: uatm, here 1 -> 1 atm

        # Calculate pCO2_sw (apply temperature corrections)

        # 2014 we had an inline temperature sensor to measure seawater
        # temperature continuously
        Tsw = bulk.temperature_external + 273.15
        # If above absent, we need regression coefficients from data (below
        # is from 2013 - as per Tonya):
        Tsw = (0.8866 * (T_eq - 273.15) - 0.7466) + 273.15
        # Now apply the temperature correction (Takahashi et al. 1993)
        pCO2_sw = pCO2_sw * np.exp(0.0423 * (Tsw - T_eq))

        # Calculate solubility (mmol/m3/atm) (Weiss, 1974)
        A1, A2, A3 = -58.0931, 90.5069, =22.2940
        B1, B2, B3 = 0.027766, -0.025888, 0.0050578
        sol = (np.exp(A1 + A2 * (100.0 / Tsw) + A3 * np.log(Tsw / 100.0) +
                      sal * (B1 + B2 * (Tsw / 100.0) + B3 *
                             (Tsw / 100.0) ** 2.0)) * 1e6)
        # and Schmidt number (Jahne et al. 1987)
        Tsw = Tsw - 273.15
        sc = (2073.1 - 125.62 * Tsw + 3.6276 * Tsw ** 2 - 0.043219 *
              Tsw ** 3.0)
