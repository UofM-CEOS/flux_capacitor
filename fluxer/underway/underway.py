#! /usr/bin/env python

"""Processing of underway pCO2 system data.

We take a single file output from the database, containing the 1-min
resolution meteorology, radiation, and underway system data.  This input
file, and other processing options, are specified in a configuration file,
which is the single input parameter for this module.

"""

__all__ = ["main", "underway_pCO2"]

import os.path as osp
import glob
import numpy as np
import pandas as pd
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
plt.style.use("ggplot")
from fluxer import parse_config

def underway_pCO2(period_file, config):
    """"""
    # Extract all the config pieces
    colnames = config["UW Inputs"]["colnames"]
    dtypes = config["UW Inputs"]["dtypes"]
    # Read, specifying the options matching what we get in our database
    # output files
    uw = pd.read_csv(period_file, dtype=dtypes,
                     parse_dates=[0, 1], index_col=1, names=colnames,
                     na_values=["NAN"], true_values=["t"],
                     false_values=["f"])
    uw_nrows = len(ec.index)

    # Convert from molar concentration to dry air mixing ratio - for both
    # underway and MET data

    # Underway - seawater

    # Bulk air molar concentration, CO2 molar concentration, water vapor
    # molar concentration
    air_sw = ((bulk.uw_pressure_analyzer * 100.0) /
              (gasR * (bulk.uw_temperature_analyzer + 273.15)))
    CO2_sw = (bulk.uw_CO2_fraction / 1.0e6) * air_sw
    H2O_sw = (bulk.uw_H2O_fraction / 1000.0) * air_sw
    # Dry air mixing ratio
    pCO2_sw = (CO2_sw / (air_sw - H2O_sw)) * 1.0e6

    # MET - atmosphere

    # Bulk air molar concentration, CO2 molar concentration, water vapor
    # molar concentration
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

    # Calculate salinity - from: UNESCO Technical Papers in Marine Science
    # #44: Algorithms for computation of fundamental properties of seawater
    ctd_pressure_dbar = bulk.ctd_pressure / 10.0 # dbar
    # This is the conductivity of seawater at S=35,T=15,p=0 (hard
    # number to dig up)
    bigR = bulk.ctd_conductivity / 42.914

    # Define coefficients
    a0, a1, a2, a3, a4, a5 = (0.008, -0.1692, 25.3851,
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

    # For now, we'll use H2O T (not tank T)... will resolve that problem
    # later.
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
    # If above absent, we need regression coefficients from data (below is
    # from 2013 - as per Tonya):
    Tsw = (0.8866 * (T_eq - 273.15) - 0.7466) + 273.15
    # Now apply the temperature correction (Takahashi et al. 1993)
    pCO2_sw = pCO2_sw * np.exp(0.0423 * (Tsw - T_eq))

    # Calculate solubility (mmol/m3/atm) (Weiss, 1974)
    A1, A2, A3 = -58.0931, 90.5069, 22.2940
    B1, B2, B3 = 0.027766, -0.025888, 0.0050578
    sol = (np.exp(A1 + A2 * (100.0 / Tsw) + A3 * np.log(Tsw / 100.0) +
                  sal * (B1 + B2 * (Tsw / 100.0) + B3 *
                         (Tsw / 100.0) ** 2.0)) * 1e6)
    # and Schmidt number (Jahne et al. 1987)
    Tsw = Tsw - 273.15
    sc = (2073.1 - 125.62 * Tsw + 3.6276 * Tsw ** 2 - 0.043219 *
          Tsw ** 3.0)


def main(config_file):
    """"""
    config = parse_config(config_file)
    uw_idir = config["UW Inputs"]["input_directory"]
    uw_files = config["UW Inputs"]["input_files"]
    colnames = config["UW Inputs"]["colnames"]
    summary_file = config["UW Outputs"]["summary_file"]
    # Stop if we don't have any files
    if (len(uw_files) < 1):
        raise Exception("There are no input files")

    for uw_file in uw_files:
        print(uw_file)             # REMOVE FOR PRODUCTION
        # Get a file name prefix to be shared by the output files from this
        # period.  Note iname is THE SAME AS THE INDEX IN OSUMMARY
        iname = osp.basename(uw_file)
        iname_prefix = osp.splitext(iname)[0]
        # ec_wind_corr, ec_flags = flux_period(ec_file, config)
        # # Save to file with suffix '_mc.csv'
        # ec_wind_corr.to_csv(osp.join(ec_idir, iname_prefix + '_mc.csv'),
        #                     index_label=colnames[1])
        # for k, v in ec_flags.iteritems():
        #     osummary.loc[iname, k] = v

    # Now we have the summary DataFrame filled up and can work with it.
    osummary.to_csv(summary_file, index_label="input_file")
    print("Summary of fluxes written to " + summary_file)


if __name__ == "__main__":
    import argparse
    description = ("Perform underway pCO2 calculations, " +
                   "given a configuration file.")
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("config_file", type=str,
                        help="Path to configuration file")
    args = parser.parse_args()
    main(args.config_file)
