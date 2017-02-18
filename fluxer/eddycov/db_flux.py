# pylint: disable=too-many-locals,invalid-name,no-member

"""Steps towards CO2 flux analyses, using data files output from PostgreSQL
database

Main function takes a single argument, which should be the path to a
configuration file containing necessary set up information such as location
of input files and variables.

"""

import os.path as osp
import logging
import numpy as np
import pandas as pd
from .parse_ecfile import (prepare_period, NavigationError,
                           SonicError, FluxError)
from ..flux_config import parse_config
from .flux import (wind3D_correct)
from .tilt_windows import TiltWindows

__all__ = ["main", "wind3D_correct_period"]

logger = logging.getLogger(__name__)
# Add the null handler if importing as library; whatever using this library
# should set up logging.basicConfig() as needed
logger.addHandler(logging.NullHandler())


def wind3D_correct_period(ec_prep, config, **kwargs):
    """Perform wind motion correction on period dataframe

    Parameters
    ----------
    ec_prep : pandas.DataFrame
        Pandas DataFrame with prepared flux data.
    config : OrderedDict
        Dictionary with parsed configuration file.
    tilt_motion : numpy.array, optional keyword
        Passed to wind3D_correct.
    tilt_anemometer : numpy.array, optional keyword
        Passed to wind3D_correct.

    Returns
    -------
    pandas.DataFrame with wind corrected flux data.

    """
    # Extract all the config pieces
    imu2anem_pos = config["EC Motion Correction"]["imu2anemometer_pos"]
    sample_freq_hz = config["EC Inputs"]["sample_frequency"]
    Tcf = config["EC Motion Correction"]["complementary_filter_period"]
    Ta = config["EC Motion Correction"]["accel_highpass_cutoff"]
    # Extract needed components from ec_prepared -- Note we use .loc access
    # method to avoid confusion about working with copies. Unnecessary for
    # extracted components we won't be replacing.
    wind = ec_prep.loc[:, ["wind_speed_u", "wind_speed_v", "wind_speed_w"]]
    motion3d_names = ["acceleration_x", "acceleration_y", "acceleration_z",
                      "rate_x", "rate_y", "rate_z"]
    motion3d = ec_prep[motion3d_names]
    heading = ec_prep["heading"]
    sog = ec_prep["speed_over_ground"]
    # Pop keyword arguments, and default to the same defaults as
    # wind3D_correct.
    mot3d_phitheta = kwargs.pop("tilt_motion", np.array([0, 0]))
    wnd3d_phitheta = kwargs.pop("tilt_anemometer", np.array([0, 0]))

    # Save full tuple output and select later. Note that we the use the
    # interpolated, smoothed heading and speed over ground.
    UVW = wind3D_correct(wind.values,
                         motion3d.loc[:, :"acceleration_z"].values,
                         motion3d.loc[:, "rate_x":].values,
                         heading.values, sog.values, imu2anem_pos,
                         sample_freq_hz, Tcf, Ta)
    logger.debug("Motion corrected with unknown tilt angles")
    # Ship-referenced speeds
    UVW_ship = UVW[0]
    # Earth-referenced speeds
    UVW_earth = UVW[11]

    UVW_tilt = wind3D_correct(wind.values,
                              motion3d.loc[:, :"acceleration_z"].values,
                              motion3d.loc[:, "rate_x":].values,
                              heading.values, sog.values, imu2anem_pos,
                              sample_freq_hz, Tcf, Ta,
                              tilt_motion=mot3d_phitheta,
                              tilt_anemometer=wnd3d_phitheta)
    logger.debug("Motion corrected with calculated tilt angles")
    UVW_ship_tilt = UVW_tilt[0]
    UVW_earth_tilt = UVW_tilt[11]

    # Append corrected wind vectors to DataFrame
    wind_corr_names = ["wind_speed_u_ship_notilt",
                       "wind_speed_v_ship_notilt",
                       "wind_speed_w_ship_notilt",
                       "wind_speed_u_earth_notilt",
                       "wind_speed_v_earth_notilt",
                       "wind_speed_w_earth_notilt",
                       "wind_speed_u_ship_tilt",
                       "wind_speed_v_ship_tilt",
                       "wind_speed_w_ship_tilt",
                       "wind_speed_u_earth_tilt",
                       "wind_speed_v_earth_tilt",
                       "wind_speed_w_earth_tilt"]
    wind[wind_corr_names[0]] = pd.Series(UVW_ship[:, 0], index=wind.index)
    wind[wind_corr_names[1]] = pd.Series(UVW_ship[:, 1], index=wind.index)
    wind[wind_corr_names[2]] = pd.Series(UVW_ship[:, 2], index=wind.index)
    wind[wind_corr_names[3]] = pd.Series(UVW_earth[:, 0],
                                         index=wind.index)
    wind[wind_corr_names[4]] = pd.Series(UVW_earth[:, 1],
                                         index=wind.index)
    wind[wind_corr_names[5]] = pd.Series(UVW_earth[:, 2],
                                         index=wind.index)
    wind[wind_corr_names[6]] = pd.Series(UVW_ship_tilt[:, 0],
                                         index=wind.index)
    wind[wind_corr_names[7]] = pd.Series(UVW_ship_tilt[:, 1],
                                         index=wind.index)
    wind[wind_corr_names[8]] = pd.Series(UVW_ship_tilt[:, 2],
                                         index=wind.index)
    wind[wind_corr_names[9]] = pd.Series(UVW_earth_tilt[:, 0],
                                         index=wind.index)
    wind[wind_corr_names[10]] = pd.Series(UVW_earth_tilt[:, 1],
                                          index=wind.index)
    wind[wind_corr_names[11]] = pd.Series(UVW_earth_tilt[:, 2],
                                          index=wind.index)

    # Append results
    ec_wind_corr = pd.concat((ec_prep,
                              wind.loc[:, wind_corr_names[0]:]), axis=1)
    return ec_wind_corr


def main(config_file):
    """Perform flux analyses, given a configuration file

    Parameters
    ----------
    config_file : str
        Path to configuration file.

    Returns
    -------
    None

    Writes summary file and prints messages from process.

    """
    # Parse configuration file
    config = parse_config(config_file)
    ec_idir = config["EC Inputs"]["input_directory"]
    ec_files = config["EC Inputs"]["input_files"]
    colnames = config["EC Inputs"]["colnames"]
    ec_tilt_window_width = int(config["EC Motion Correction"]
                               ["tilt_window_width"])
    summary_file = config["EC Outputs"]["summary_file"]
    # Stop if we don't have any files
    if len(ec_files) < 1:
        raise FluxError("There are no input files")

    # Initialize tilt windows to work on
    ec_windows = TiltWindows(ec_files, ec_tilt_window_width)
    ec_windows.get_tilts_planarfit(config)
    tilts = ec_windows.tilts
    # Perhaps plot the tilt window data?
    tilt_figf = "ec_{0}min.png".format(ec_tilt_window_width)
    ec_windows.plot(tilt_figf)
    logger.info("Plot of tilt window calculations written to %s", tilt_figf)
    # Write tilt window calculations
    tilt_ofile = "tilts_{0}.csv".format(ec_tilt_window_width)
    ec_windows.tilts.to_csv(tilt_ofile, index_label="timestamp")
    logger.info("Tilt window calculations written to %s", tilt_ofile)
    # Summary DataFrame is also filled up and can be output
    ec_windows.prep_flags.to_csv(summary_file, index_label="input_file")
    logger.info("Pre-motion correction summary written to %s",
                summary_file)

    logger.setLevel(logging.INFO)
    for (w, l) in ec_windows.win_files:
        tilt_w = tilts.loc[w]
        if tilt_w.loc["failed_tilt_flag"]:
            continue
        mot3d_tilt = np.array([tilt_w.loc["phi_motion"],
                               tilt_w.loc["theta_motion"]])
        sonic_tilt = np.array([tilt_w.loc["phi_sonic"],
                               tilt_w.loc["theta_sonic"]])
        for ec_file in l:
            iname = osp.basename(ec_file)
            logger.info("Begin motion correction %s", iname)
            try:
                ec_prep, _ = prepare_period(ec_file, config)
            except (NavigationError, SonicError):
                continue
            else:
                iname_prefix = osp.splitext(iname)[0]
                ec_wind_corr = wind3D_correct_period(ec_prep, config,
                                                     tilt_motion=mot3d_tilt,
                                                     tilt_anemometer=sonic_tilt)  # noqa: E501
                # Save to file with suffix "_mc.csv"
                ec_wind_corr.to_csv(osp.join(ec_idir, iname_prefix +
                                             "_mc.csv"),
                                    index_label=colnames[1], na_rep="NaN")
            logger.info("End motion correction %s", iname)
            # TODO: Further flux processing
