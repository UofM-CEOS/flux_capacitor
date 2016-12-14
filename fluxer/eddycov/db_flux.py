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
from scipy.stats import circmean
from .settings import (FLUX_FLAGS)
from ..flux_config import parse_config
from .flux import (smooth_angle, wind3D_correct, despike_VickersMahrt,
                   planarfit, rotate_vectors)
from .tilt_windows import TiltWindows

__all__ = ["main", "prepare_period", "wind3D_correct_period"]

logger = logging.getLogger(__name__)
# Add the null handler if importing as library; whatever using this library
# should set up logging.basicConfig() as needed
logger.addHandler(logging.NullHandler())


# Exception classes for catching our conditions
class FluxError(Exception):
    """Base class for Exceptions in this module"""
    pass


class SonicError(FluxError):
    """Critical sonic anemometer Exception"""
    def __init__(self, message, flags):
        self.message = message
        self.flags = flags


class NavigationError(FluxError):
    """Critical navigation Exception"""
    def __init__(self, message, flags):
        self.message = message
        self.flags = flags


class MeteorologyError(FluxError):
    """Critical meteorology Exception"""
    def __init__(self, message, flags):
        self.message = message
        self.flags = flags


def prepare_period(period_file, config):
    """Parse input period file and set up data for flux analyses

    Parameters
    ----------
    period_file : str
        Path to file with pre-filtered, candidate, flux data.
    config : OrderedDict
        Dictionary with parsed configuration file.

    Returns
    -------

    Tuple of two pandas.DataFrame objects; first represents prepared flux
    data, and the second the flags encountered during tests.  The latter is
    also passed along with exception, if raised.

    """
    # Extract all the config pieces
    colnames = config["EC Inputs"]["colnames"]
    imu_xyz_idx = np.array(config["EC Motion Correction"]["imu_xyz_idx"],
                           dtype=int)
    imu2rhs_linmult = config["EC Motion Correction"]["imu2rhs_linaccel_mult"]
    imu2rhs_angmult = config["EC Motion Correction"]["imu2rhs_angaccel_mult"]
    dtypes = config["EC Inputs"]["dtypes"]
    # Read, specifying the options matching what we get in our database
    # output files
    ec = pd.read_csv(period_file, dtype=dtypes, header=1,
                     parse_dates=[0, 1], index_col=1, names=colnames,
                     na_values=["NAN"], true_values=["t"],
                     false_values=["f"])
    ec_nrows = len(ec.index)
    # Initial values for flags
    period_flags = dict.fromkeys(FLUX_FLAGS, False)
    # Put acceleration components in 3-column array.  Original comment:
    # read in angular rates in RH coordinate system, convert to rad/s.
    imu_linaccel_names = ["acceleration_x", "acceleration_y",
                          "acceleration_z"]
    imu_lin_xname = imu_linaccel_names[imu_xyz_idx[0]]
    imu_lin_yname = imu_linaccel_names[imu_xyz_idx[1]]
    imu_lin_zname = imu_linaccel_names[imu_xyz_idx[2]]
    imu_angaccel_names = ["rate_x", "rate_y", "rate_z"]
    imu_ang_xname = imu_angaccel_names[imu_xyz_idx[0]]
    imu_ang_yname = imu_angaccel_names[imu_xyz_idx[1]]
    imu_ang_zname = imu_angaccel_names[imu_xyz_idx[2]]
    # Now we can reorder input to RHS and scale units
    motion3d = pd.DataFrame({'acceleration_x': (imu2rhs_linmult[0] *
                                                ec[imu_lin_xname] * 9.81),
                             'acceleration_y': (imu2rhs_linmult[1] *
                                                ec[imu_lin_yname] * 9.81),
                             'acceleration_z': (imu2rhs_linmult[2] *
                                                ec[imu_lin_zname] * 9.81),
                             'rate_phi': (imu2rhs_angmult[0] *
                                          np.radians(ec[imu_ang_xname])),
                             'rate_theta': (imu2rhs_angmult[1] *
                                            np.radians(ec[imu_ang_yname])),
                             'rate_shi': (imu2rhs_angmult[2] *
                                          np.radians(ec[imu_ang_zname]))})
    wind = ec.loc[:, ["wind_speed_u", "wind_speed_v", "wind_speed_w"]]

    # [Original comment: check for any significant number of 'NAN's (not
    # worried about the odd one scattered here and there)].  [Original
    # comment: set open flag if gt 2% of records are 'NAN']
    win_width = np.int(config["EC Despiking"]["despike_win_width"])
    win_step = np.int(config["EC Despiking"]["despike_step"])
    nreps = np.int(config["EC Despiking"]["despike_nreps"])
    if (((ec.op_CO2_density.count() / float(ec_nrows)) < 0.98) or
        ((ec.op_H2O_density.count() / float(ec_nrows)) < 0.98) or
        ((ec.op_analyzer_status.count() / float(ec_nrows)) < 0.98)):
        period_flags["open_flag"] = True
        logger.debug("Setting open_flag -> too many missing records")
    # TODO: Here we need to prepare our check for the diagnostics from the
    # open path analyzers.  For now, keep using the rule of thumb
    elif (ec.op_analyzer_status.gt(249) |
          ec.op_analyzer_status.lt(240)).sum() > 0.02:
        period_flags["open_flag"] = True
        logger.debug("Setting open_flag -> too many bad analyzer status")
    else:
        op_CO2_VM = despike_VickersMahrt(ec.op_CO2_density,
                                         width=win_width, step=win_step,
                                         zscore_thr=3.5, nreps=nreps)
        ec.op_CO2_density = op_CO2_VM[0]
        op_H2O_VM = despike_VickersMahrt(ec.op_H2O_density,
                                         width=win_width, step=win_step,
                                         zscore_thr=3.5, nreps=nreps)
        ec.op_H2O_density = op_H2O_VM[0]
        op_Pr_VM = despike_VickersMahrt(ec.op_pressure,
                                        width=win_width, step=win_step,
                                        zscore_thr=3.5, nreps=nreps)
        ec.op_pressure = op_Pr_VM[0]
        if (((op_CO2_VM[1] / ec.op_CO2_density.count()) > 0.01) or
            ((op_H2O_VM[1] / ec.op_H2O_density.count()) > 0.01) or
            ((op_Pr_VM[1] / ec.op_pressure.count()) > 0.01)):
            period_flags["open_flag"] = True
            logger.debug("Setting open_flag -> too many despiked records")

    # [Original comment: set wind flag if gt 2% of records are 'NAN']
    if (((wind.wind_speed_u.count() / float(ec_nrows)) < 0.98) or
        ((wind.wind_speed_v.count() / float(ec_nrows)) < 0.98) or
        ((wind.wind_speed_w.count() / float(ec_nrows)) < 0.98) or
        ((ec.air_temperature_sonic.count() / float(ec_nrows)) < 0.98)):
        period_flags["sonic_flag"] = True
        logger.debug("Setting sonic_flag -> too many missing records")
    else:
        # [Original comment: now that we have looked for NANs, we may as
        # well fill in the NANs and any spikes using the shot filter].
        # [SPL: these changes are done outside the WIND array, which is the
        # one that is used later for motion correction, etc., so they are
        # lost.]
        wind_u_VM = despike_VickersMahrt(wind.wind_speed_u,
                                         width=win_width, step=win_step,
                                         zscore_thr=3.5, nreps=nreps)
        wind.wind_speed_u = wind_u_VM[0]
        wind_v_VM = despike_VickersMahrt(wind.wind_speed_v,
                                         width=win_width, step=win_step,
                                         zscore_thr=3.5, nreps=nreps)
        wind.wind_speed_v = wind_v_VM[0]
        wind_w_VM = despike_VickersMahrt(wind.wind_speed_w,
                                         width=win_width, step=win_step,
                                         zscore_thr=5.0, nreps=nreps)
        wind.wind_speed_w = wind_w_VM[0]
        # ec.air_temperature_sonic = shot_filter(ec.air_temperature_sonic)
        sonic_T_VM = despike_VickersMahrt(ec.air_temperature_sonic,
                                          width=win_width, step=win_step,
                                          zscore_thr=3.5, nreps=nreps)
        ec.air_temperature_sonic = sonic_T_VM[0]
        if (((wind_u_VM[1] / wind.wind_speed_u.count()) > 0.01) or
            ((wind_v_VM[1] / wind.wind_speed_v.count()) > 0.01) or
            ((wind_w_VM[1] / wind.wind_speed_w.count()) > 0.01) or
            ((sonic_T_VM[1] / ec.air_temperature_sonic.count()) > 0.01)):
            period_flags["sonic_flag"] = True
            logger.debug("Setting sonic_flag -> too many despiked records")

    # [Original comment: set motion flag if gt 2% of records are 'NAN']
    if (((motion3d.acceleration_x.count() / float(ec_nrows)) < 0.98) or
        ((motion3d.acceleration_y.count() / float(ec_nrows)) < 0.98) or
        ((motion3d.acceleration_z.count() / float(ec_nrows)) < 0.98) or
        ((motion3d.rate_phi.count() / float(ec_nrows)) < 0.98) or
        ((motion3d.rate_theta.count() / float(ec_nrows)) < 0.98) or
        ((motion3d.rate_shi.count() / float(ec_nrows)) < 0.98)):
        period_flags["motion_flag"] = True
        logger.debug("Setting motion_flag -> too many missing records")
    else:
        # [Original comment: shot filter the motion channels... this helps with
        # a problem where unreasonably high accelerations cause a 'NaN'
        # calculation]
        for col in motion3d.columns:
            motcol_VM = despike_VickersMahrt(motion3d[col], width=win_width,
                                             step=win_step, zscore_thr=3.5,
                                             nreps=nreps)
            motion3d[col] = motcol_VM[0]
    # Ready to assign to main dataframe
    ec[[imu_lin_xname, imu_lin_yname, imu_lin_zname]] = motion3d.iloc[:, 0:3]
    ec[[imu_ang_xname, imu_ang_yname, imu_ang_zname]] = motion3d.iloc[:, 3:]

    # Set closed flag is more than 2% of records are NaN.
    if (((ec.cp_CO2_fraction.count() / float(ec_nrows)) < 0.98) or
        ((ec.cp_H2O_fraction.count() / float(ec_nrows)) < 0.98) or
        ((ec.cp_pressure.count() / float(ec_nrows)) < 0.98)):
        period_flags["closed_flag"] = True
        logger.debug("Setting closed_flag -> too many missing records")
    else:
        cp_CO2_VM = despike_VickersMahrt(ec.cp_CO2_fraction,
                                         width=win_width, step=win_step,
                                         zscore_thr=3.5, nreps=nreps)
        ec.cp_CO2_fraction = cp_CO2_VM[0]
        cp_H2O_VM = despike_VickersMahrt(ec.cp_H2O_fraction,
                                         width=win_width, step=win_step,
                                         zscore_thr=3.5, nreps=nreps)
        ec.cp_H2O_fraction = cp_H2O_VM[0]
        cp_Pr_VM = despike_VickersMahrt(ec.cp_pressure,
                                        width=win_width, step=win_step,
                                        zscore_thr=3.5, nreps=nreps)
        ec.cp_pressure = cp_Pr_VM[0]
        if (((cp_CO2_VM[1] / ec.cp_CO2_fraction.count()) > 0.01) or
            ((cp_H2O_VM[1] / ec.cp_H2O_fraction.count()) > 0.01) or
            ((cp_Pr_VM[1] / ec.cp_pressure.count()) > 0.01)):
            period_flags["closed_flag"] = True
            logger.debug("Setting closed_flag -> too many despiked records")

    # [Original comment: now fill in the gaps by applying a moving
    # average... In this case, we use a 100 sample window (10 sec) moving
    # average... may need to tweak this value].  [SPL: perhaps a simple
    # linear interpolation is better; I don't know why this moving average
    # is used, where a window must be specified and may be introducing
    # bias.  Perhaps it doesn't matter.  Why aren't latitude and longitude
    # similarly interpolated?]
    cog, sog = smooth_angle(ec["course_over_ground"].values,
                            ec["speed_over_ground"].values, 50)
    ec["course_over_ground"] = cog
    ec["speed_over_ground"] = sog
    heading, _ = smooth_angle(ec["heading"].values, 1, 50)
    ec["heading"] = heading
    if ((ec["course_over_ground"].count() < len(cog)) or
        (ec["speed_over_ground"].count() < len(sog)) or
        (ec["heading"].count < len(heading))):
        period_flags["motion_flag"] = True
        logger.debug(("Setting motion_flag -> missing records"
                      "in COG, SOG, or heading"))

    # [Original comment: check for bad wind data: bad wind data can
    # usually be diagnosed by unusually high wind speeds.  This is
    # most obvious in the vertical wind where we wouldn't expect high
    # values bad sonic data can also turn up in the Tsonic before the
    # wind, check the deviation between Tsonic and mean air T.]
    air_temp_avg = ec["air_temperature"].mean()
    nbad_vertical_wind = abs(wind["wind_speed_w"]).gt(7).sum()
    nbad_air_temp_sonic = abs(ec["air_temperature_sonic"] -
                              air_temp_avg).gt(7).sum()
    if not (np.isfinite(air_temp_avg) or
            np.isfinite(ec.relative_humidity[0])):
        period_flags["bad_meteorology_flag"] = True
        logger.debug(("Setting bad_meteorology_flag -> missing mean "
                      "air temperature or relative_humidity"))

    # Now raise Exceptions to signal rest of analyses cannot continue.
    # Return flags so they can be caught and used later.

    # Set wind flag high if gt 0.5% of records are frost contaminated
    # [Original comment: check critical low frequency variables]
    if ((nbad_vertical_wind / float(ec_nrows)) > 0.05 or
        (nbad_air_temp_sonic / float(ec_nrows)) > 0.05):
        period_flags["sonic_flag"] = True
        logger.error(("Setting sonic_flag -> Too many "
                      "records where vertical wind was too high or "
                      "where sonic air temperature was aberrant "))
        raise SonicError("Bad sonic anemometer data", period_flags)
    # # Below will be needed at some point
    # sw_avg = ec.K_down[0]
    # lw_avg = ec.LW_down[0]
    # sog_avg = ec["speed_over_ground"].mean()

    # If we have no good COG, SOG, or heading, then we cannot continue.
    if (ec["course_over_ground"].count() < 1 or
        ec["speed_over_ground"].count() < 1 or
        ec["heading"].count() < 1):
        period_flags["bad_navigation_flag"] = True
        logger.error(("Setting bad_navigation_flag -> COG, "
                      "SOG, or heading are all missing "))
        raise NavigationError("Unusable COG, SOG, or heading records",
                              period_flags)

    return ec, period_flags


def wind3D_correct_period(ec_prep, config, **kwargs):
    """Perform wind motion correction on period dataframe

    Parameters
    ----------
    ec_prep : pandas.DataFrame
        Pandas DataFrame with prepared flux data.
    config : OrderedDict
        Dictionary with parsed configuration file.

    Keyword Parameters
    ------------------
    tilt_motion : numpy.array
        Passed to wind3D_correct.
    tilt_anemometer : numpy.array
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
    # [Original comment: create flags for the 6 possible sources of "bad"
    # data, flag=0 means data good]
    flags = dict.fromkeys(FLUX_FLAGS, False)
    # We set up a dataframe with all files to process as index, and all the
    # flags as columns.  This is the basis for our summary output file;
    # other columns (such as flux summary calculations for the period)
    # will be appended as we loop.
    osummary = pd.DataFrame(flags,
                            index=[osp.basename(x) for x in ec_files])

    # Iterate over the list of tuples (key=window time stamp, file-list)
    for (k, l) in ec_windows.win_files:
        logger.info("Begin window %s", k)
        ec_list = []
        ec_list_keys = []
        # Iterate over the file-list
        for ec_file in l:
            logger.info("Begin preparing %s", osp.basename(ec_file))
            # Get a file name prefix to be shared by the output files from this
            # period.  Note iname is THE SAME AS THE INDEX IN OSUMMARY
            iname = osp.basename(ec_file)
            iname_prefix = osp.splitext(iname)[0]
            # Try to prepare file; set flags, and populate the list of OK
            # files and summaries. Continue if preparation fails.
            try:
                ec_prep, prep_flags = prepare_period(ec_file, config)
            except (NavigationError, SonicError) as err:
                for flagname, flag in err.flags.items():
                    osummary.loc[iname, flagname] = flag
                    ec_windows.tilts.loc[k, "n" + flagname] += np.int(flag)
                logger.info("Skip %s", ec_file)
                continue
            else:
                for flagname, flag in prep_flags.iteritems():
                    osummary.loc[iname, flagname] = flag
                    ec_windows.tilts.loc[k, "n" + flagname] += np.int(flag)
                ec_list.append(ec_prep)
                ec_list_keys.append(ec_file)
            logger.info("End preparing %s", osp.basename(ec_file))

        n_ok = len(ec_list)
        ec_windows.tilts.loc[k, "nfiles_ok"] = n_ok
        if n_ok > 0:
            ec_win = pd.concat(ec_list, keys=ec_list_keys)
            logger.info("Combining %s periods in window %s", n_ok, k)
            wind_names = ["wind_speed_u", "wind_speed_v", "wind_speed_w"]
            wnd = ec_win.loc[:, wind_names]
            mot3d_names = ["acceleration_x", "acceleration_y",
                           "acceleration_z", "rate_x", "rate_y", "rate_z"]
            mot3d = ec_win[mot3d_names]
            # Calculate tilt angles for motion sensor
            mot3d_k, __ = planarfit(mot3d.values[:, 0:3])
            __, mot3d_tilt = rotate_vectors(mot3d.values[:, 0:3],
                                            k_vector=mot3d_k)
            ec_windows.tilts.loc[k, "phi_motion"] = mot3d_tilt[0]
            ec_windows.tilts.loc[k, "theta_motion"] = mot3d_tilt[1]
            # Calculate tilt angles for sonic anemometer
            wnd3d_k, __ = planarfit(wnd.values[:, 0:3])
            __, wnd3d_tilt = rotate_vectors(wnd.values[:, 0:3],
                                            k_vector=wnd3d_k)
            ec_windows.tilts.loc[k, "phi_sonic"] = wnd3d_tilt[0]
            ec_windows.tilts.loc[k, "theta_sonic"] = wnd3d_tilt[1]
            wind_direction_mean = circmean(ec_win["wind_direction"].values,
                                           high=np.degrees(2 * np.pi))
            ec_windows.tilts.loc[k, "wind_direction"] = wind_direction_mean

            for ec_file in ec_list_keys:
                logger.info("Begin processing %s", iname)
                ec_prep = ec_win.ix[ec_file]
                iname = osp.basename(ec_file)
                iname_prefix = osp.splitext(iname)[0]
                ec_wind_corr = wind3D_correct_period(ec_prep, config,
                                                     tilt_motion=mot3d_tilt,
                                                     tilt_anemometer=mot3d_tilt)  # noqa: E501
                # Save to file with suffix "_mc.csv"
                ec_wind_corr.to_csv(osp.join(ec_idir, iname_prefix +
                                             "_mc.csv"),
                                    index_label=colnames[1])

                # TODO: Further flux processing

                logger.info("End processing %s", iname)

        logger.info("End window %s", k)

    # Perhaps plot the tilt window data?
    tilt_figf = "ec_{0}min.png".format(ec_tilt_window_width)
    ec_windows.plot(tilt_figf)
    logger.info("Plot of tilt window calculations written to %s", tilt_figf)
    # Write tilt window calculations
    tilt_ofile = "tilts_{0}.csv".format(ec_tilt_window_width)
    ec_windows.tilts.to_csv(tilt_ofile, index_label="timestamp")
    logger.info("Tilt window calculations written to %s", tilt_ofile)

    # Now we have the summary DataFrame filled up and can work with it.
    osummary.to_csv(summary_file, index_label="input_file")
    logger.info("Summary of fluxes written to %s", summary_file)
