#! /usr/bin/env python

"""Processing of underway pCO2 system data.

We take a single file output from the database, containing the 1-min
resolution meteorology, radiation, and underway system data.  This input
file, and other processing options, are specified in a configuration file,
which is the single input parameter for this module.

The following variable names are expected in input with units:

uw_pressure_analyzer (mbar)
uw_temperature_analyzer (C)
uw_CO2_fraction (umol/mol)
uw_H2O_fraction (mmol/mol)
ctd_pressure (dbar)
ctd_temperature (C)
ctd_conductivity (S/m; Siemens/m)
equ_temperature (C)
temperature_external (C)
bad_CO2_flag (boolean)
bad_H2O_flag (boolean)
bad_temperature_external (boolean)
bad_temperature_analyzer (boolean)
bad_pressure_analyzer_flag (boolean)
bad_equ_temperature_flag (boolean)
"""

import os.path as osp
import glob
import numpy as np
import pandas as pd
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
plt.style.use("ggplot")
from fluxer import parse_config

__all__ = ["main", "underway_pCO2"]

def underway_pCO2(period_file, config):
    """Perform pCO2 computations on period."""
    R_u = 8.31451              # j/mol/k universal gas constant
    # Extract all the config pieces
    colnames = config["UW Inputs"]["colnames"]
    dtypes = config["UW Inputs"]["dtypes"]
    uw_depth = config["UW Inputs"]["uw_intake_depth"]
    calc_ext_temp = config["UW Inputs"]["uw_regress_temperature_external"]
    ext_temp_coefs = config["UW Inputs"]["uw_temperature_external_coefs"]
    # Read, specifying the options matching what we get in our database
    # output files
    uw = pd.read_csv(period_file, dtype=dtypes, header=1,
                     parse_dates=[0, 1], index_col=1, names=colnames,
                     na_values=["NAN"], true_values=["t"],
                     false_values=["f"])
    uw_nrows = len(uw.index)

    # Convert from molar concentration to dry air mixing ratio

    # In most cases, this should be very minor... provided that the
    # condensors are running properly.

    # Underway - seawater
    # Bulk air molar concentration, CO2 molar concentration, water vapor
    # molar concentration
    air_sw = ((uw.uw_pressure_analyzer * 100.0) /
              (R_u * (uw.uw_temperature_analyzer + 273.15)))
    CO2_sw = (uw.uw_CO2_fraction / 1.0e6) * air_sw
    H2O_sw = (uw.uw_H2O_fraction / 1000.0) * air_sw
    # Dry air mixing ratio
    xCO2 = (CO2_sw / (air_sw - H2O_sw)) * 1.0e6

    # # MET - atmosphere
    # # Bulk air molar concentration, CO2 molar concentration, water vapor
    # # molar concentration
    # air_atm = ((uw.cp_pressure * 100.0) /
    #            (R_u * (uw.cp_temperature + 273.15)))
    # CO2_atm = (uw.cp_CO2_fraction / 1.0e6) * air_atm
    # H2O_atm = (uw.cp_H2O_fraction / 1000.0) * air_atm
    # # Dry air mixing ratio
    # pCO2_atm = (CO2_atm / (air_atm - H2O_atm)) * 1.0e6
    # # Extrapolate wind to height of 2D anemometer; using eq. 4.14b from
    # # Stull - Meteorology for Scientists and Engineers
    # anemometer_height = config["UW Inputs"]["anemometer2d_height"]
    # wind_speed_adj = (uw.true_wind_speed * np.log(10 / 0.0002) /
    #                   np.log(anemometer_height / 0.0002))

    # Calculate salinity - from: UNESCO Technical Papers in Marine Science
    # #44: Algorithms for computation of fundamental properties of seawater
    # This is the conductivity of seawater at S=35, T=15, p=0 (hard number
    # to dig up)
    bigR = uw.ctd_conductivity / 42.914
    # Tim points out that this computation should use hydrostatic pressure
    # at the underway's intake depth, so we need to calculate it. The
    # formula that gives the P pressure (N/m2, Pa) on an object submerged
    # in a fluid is: P = r * g * h, where r (rho) is the density of the
    # fluid, g is the acceleration of gravity, h is the height of the fluid
    # above the object (NASA).  Formula are developed for salinity at
    # depth, where 1 decibar is equall to approximately 1 m depth, assuming
    # pressure at the surface to be zero.  The density of sea water is
    # 1.03e3 kg/m3. Underway intake depth is in m. Convert result from Pa
    # to dbar.
    pressure_sw = (1.03e3 * 9.8 * uw_depth) * 1.0e-4

    # Define coefficients
    a0, a1, a2, a3, a4, a5 = (0.008, -0.1692, 25.3851,
                              14.0941, -7.0261, 2.7081)
    b0, b1, b2, b3, b4, b5 = (0.0005, -0.0056, -0.0066,
                              -0.0375, 0.0636, -0.0144)
    c0, c1, c2, c3, c4 = (0.6766097, 2.00564e-2, 1.104259e-4,
                          -6.9698e-7, 1.0031e-9)
    d1, d2, d3, d4 = (3.426e-2, 4.464e-4, 4.215e-1, -3.107e-3)
    e1, e2, e3 = (2.070e-5, -6.370e-10, 3.989e-15)
    K = 0.0162

    # Salinity calculation
    rt = (c0 + (c1 * uw.ctd_temperature) +
          (c2 * uw.ctd_temperature ** 2.0) +
          (c3 * uw.ctd_temperature ** 3.0) +
          (c4 * uw.ctd_temperature ** 4.0)) # [eqn 3]
    Rp = 1.0 + ((pressure_sw *
                 (e1 + e2 * pressure_sw +
                  e3 * pressure_sw ** 2.0)) /
                (1.0 + d1 * uw.ctd_temperature +
                 d2 * uw.ctd_temperature ** 2.0 +
                 (d3 + d4 * uw.ctd_temperature) * bigR)) #  [eqn 4]
    bigRt = bigR / (Rp * rt) # eqn after eq'n 4
    delS = (((uw.ctd_temperature - 15.0) /
             (1.0 + K * (uw.ctd_temperature - 15.0))) *
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
    T_eq = uw.equ_temperature + 273.15
    # T_air = uw.air_temperature + 273.15
    # PW_sw = (np.exp(24.4543 - 67.4509 * (100.0 / T_eq) - 4.8489 *
    #                 np.log(T_eq / 100.0) - 0.000544 * sal))
    # PW_atm = (np.exp(24.4543 - 67.4509 * (100.0 / T_air) - 4.8489 *
    #                  np.log(T_air / 100.0) - 0.000544 * sal))
    Pw = (np.exp(24.4543 - 67.4509 * (100.0 / T_eq) - 4.8489 *
                 np.log(T_eq / 100.0) - 0.000544 * sal))
    # Calculate pCO2 in the equilibrator
    pCO2_eq = xCO2 * (1 - Pw) # units: uatm, here 1 -> 1 atm

    # Calculate pCO2_sw (apply temperature corrections)

    # Since 2014 we have an inline temperature sensor to measure seawater
    # temperature continuously.  Otherwise, we calculate a surrogate from
    # equilibrator temperature
    if calc_ext_temp:
        Tsw = (ext_temp_coefs[0] + ext_temp_coefs[1] *
               uw.equ_temperature) + 273.15
    else:
        Tsw = uw.temperature_external + 273.15
    # Now apply the temperature correction (Takahashi et al. 1993)
    pCO2_sw = pCO2_eq * np.exp(0.0423 * (Tsw - T_eq))

    # Calculate fugacity at base of MDL as per Weiss 1974
    fA = 12.0408 * Tsw
    fB = 3.27957e-2 * Tsw ** 2
    fC = 3.16528e-5 * Tsw ** 3
    fBBII = -1636.75 + fA - fB + fC
    dCO2_air = 57.7 - (0.118 * Tsw)
    fCO2 = ((pCO2_sw * 1.0e-6 * 101325) *
            np.exp((((fBBII * 1.0e-6) + 2 * (dCO2_air * 1.0e-6)) * 101325) /
                   (R_u * Tsw)))
    fCO2 = (fCO2 / 101325) * 1.0e6

    # Calculate solubility (mmol/m3/atm) (Weiss, 1974)
    A1, A2, A3 = -58.0931, 90.5069, 22.2940
    B1, B2, B3 = 0.027766, -0.025888, 0.0050578
    sol = (np.exp(A1 + A2 * (100.0 / Tsw) + A3 * np.log(Tsw / 100.0) +
                  sal * (B1 + B2 * (Tsw / 100.0) + B3 *
                         (Tsw / 100.0) ** 2.0)) * 1.0e6)
    # and Schmidt number (Jahne et al. 1987)
    Tsw = Tsw - 273.15          # back to C
    sc = (2073.1 - 125.62 * Tsw + 3.6276 * Tsw ** 2 - 0.043219 *
          Tsw ** 3.0)

    # NULL data according to flags (these were setup in database)
    sal[uw.bad_ctd_flag.fillna(False)] = np.NaN
    sol[uw.bad_ctd_flag.fillna(False)] = np.NaN
    bad_pCO2 = (uw.bad_CO2_flag.fillna(False) |
                uw.bad_H2O_flag.fillna(False) |
                uw.bad_temperature_external_flag.fillna(False) |
                uw.bad_temperature_analyzer_flag.fillna(False) |
                uw.bad_pressure_analyzer_flag.fillna(False) |
                uw.bad_equ_temperature_flag.fillna(False))
    xCO2[bad_pCO2] = np.NaN
    pCO2_eq[bad_pCO2] = np.NaN
    pCO2_sw[bad_pCO2] = np.NaN
    fCO2[bad_pCO2] = np.NaN
    sc[uw.bad_temperature_external_flag.fillna(False)] = np.NaN

    # Return object
    uw_new = pd.DataFrame({'dry_air_mixing_ratio': xCO2,
                           'pCO2_equilibrator': pCO2_eq,
                           'pCO2_seawater': pCO2_sw,
                           'fCO2': fCO2,
                           'solubility': sol,
                           'salinity': sal,
                           'schmidt_number': sc},
                          index=uw.index)
    return pd.concat([uw, uw_new], axis=1)


def main(config_file):
    config = parse_config(config_file)
    uw_idir = config["UW Inputs"]["input_directory"]
    uw_files = config["UW Inputs"]["input_files"]
    colnames = config["UW Inputs"]["colnames"]
    pCO2_dir = config["UW Outputs"]["pco2_directory"]
    # Stop if we don't have any files
    if len(uw_files) < 1:
        raise Exception("There are no input files")

    for uw_file in uw_files:
        # print(uw_file)             # REMOVE FOR PRODUCTION
        # Get a file name prefix to be shared by the output files from this
        # period.  Note iname is THE SAME AS THE INDEX IN OSUMMARY
        iname = osp.basename(uw_file)
        iname_prefix = osp.splitext(iname)[0]
        uw_pCO2 = underway_pCO2(uw_file, config)
        # Save to file with suffix '_mc.csv'
        ofile = osp.join(pCO2_dir, iname_prefix + '_pCO2.csv')
        uw_pCO2.to_csv(ofile, index_label=colnames[1])
        print("pCO2 file written to " + ofile)


if __name__ == "__main__":
    import argparse
    description = ("Perform underway pCO2 calculations, " +
                   "given a configuration file.")
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("config_file", type=str,
                        help="Path to configuration file")
    args = parser.parse_args()
    main(args.config_file)
