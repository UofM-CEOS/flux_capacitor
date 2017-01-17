# pylint: disable=too-many-locals,invalid-name,no-member

"""Parse and flag input eddy covariance file

"""

import logging
import numpy as np
import pandas as pd
from .settings import (FLUX_FLAGS)
from .flux import (smooth_angle, despike_VickersMahrt,
                   rotate_coordinates)

__all__ = ["prepare_period", "FluxError", "SonicError",
           "NavigationError", "MeteorologyError"]

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
    sonic_xoffset = config["EC Motion Correction"]["sonic_xoffset"]
    dtypes = config["EC Inputs"]["dtypes"]
    # Read, specifying the options matching what we get in our database
    # output files
    ec = pd.read_csv(period_file, dtype=dtypes, header=None,
                     parse_dates=[0, 1], index_col=1, names=colnames,
                     na_values=["NAN"], true_values=["t"],
                     false_values=["f"], low_memory=False)
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
    # Rotate the sonic vectors to align the instrument's coordinate frame
    # with the ship's
    wind.loc[:] = rotate_coordinates(wind.values, theta=sonic_xoffset,
                                     axis=2, rotate_vectors=False)

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

    # Set wind flag high if gt 0.5% of records are froth contaminated
    # [Original comment: check critical low frequency variables]
    if ((nbad_vertical_wind / float(ec_nrows)) > 0.05 or
        (nbad_air_temp_sonic / float(ec_nrows)) > 0.05):
        period_flags["sonic_flag"] = True
        logger.debug(("Setting sonic_flag -> Too many "
                      "records where vertical wind was too high or "
                      "where sonic air temperature was aberrant "))
    # # Below will be needed at some point
    # sw_avg = ec.K_down[0]
    # lw_avg = ec.LW_down[0]
    # sog_avg = ec["speed_over_ground"].mean()

    # Now raise Exceptions to signal rest of analyses cannot continue.
    # Return flags so they can be caught and used later.

    # If we have no good COG, SOG, or heading, then we cannot continue.
    if (ec["course_over_ground"].count() < 1 or
        ec["speed_over_ground"].count() < 1 or
        ec["heading"].count() < 1):
        period_flags["bad_navigation_flag"] = True
        period_flags["failed_prep_flag"] = True
        logger.error(("Setting bad_navigation_flag -> COG, "
                      "SOG, or heading are all missing "))
        raise NavigationError("Unusable COG, SOG, or heading records",
                              period_flags)

    return ec, period_flags
