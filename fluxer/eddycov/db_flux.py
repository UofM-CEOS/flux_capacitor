#! /usr/bin/env python
# pylint: disable=too-many-locals,invalid-name,no-member

"""Steps towards CO2 flux analyses, using data files output from PostgreSQL
database.

It takes a single argument, which should be the path to a configuration
file containing necessary set up information such as location of input
files and variables.
"""

import argparse
import os.path as osp
import numpy as np
import pandas as pd
# import psycopg2 as pg
# from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from fluxer.flux_config import parse_config
from fluxer.eddycov.flux import (smooth_angle, wind3D_correct,
                                 despike_VickersMahrt, planarfit_coef,
                                 rotate_vectors)

__all__ = ["main", "flux_period"]

plt.style.use("ggplot")


# Exception classes for catching our problems
class FluxError(Exception):
    """Base class for Exceptions in this module"""
    pass


def flux_period(period_file, config):
    """Perform required calculations on period."""
    # Extract all the config pieces
    colnames = config["EC Inputs"]["colnames"]
    imu2anem_pos = config["EC Motion Correction"]["imu2anemometer_pos"]
    imu_xyz_idx = np.array(config["EC Motion Correction"]["imu_xyz_idx"],
                           dtype=int)
    imu2rhs_linmult = config["EC Motion Correction"]["imu2rhs_linaccel_mult"]
    imu2rhs_angmult = config["EC Motion Correction"]["imu2rhs_angaccel_mult"]
    sample_freq_hz = config["EC Inputs"]["sample_frequency"]
    Tcf = config["EC Motion Correction"]["complementary_filter_period"]
    Ta = config["EC Motion Correction"]["accel_highpass_cutoff"]
    dtypes = config["EC Inputs"]["dtypes"]
    # Read, specifying the options matching what we get in our database
    # output files
    ec = pd.read_csv(period_file, dtype=dtypes, header=1,
                     parse_dates=[0, 1], index_col=1, names=colnames,
                     na_values=["NAN"], true_values=["t"],
                     false_values=["f"])
    ec_nrows = len(ec.index)
    # Initial values for flags
    open_flag, closed_flag = False, False
    sonic_flag, motion_flag = False, False
    bad_navigation_flag, bad_meteorology_flag = False, False
    # Put acceleration components in 3-column array and make copy to keep
    # uncorrected data.  Original comment: read in angular rates in RH
    # coordinate system, convert to rad/s.
    imu_linaccel_names = ["acceleration_x", "acceleration_y",
                          "acceleration_z"]
    imu_lin_xname = imu_linaccel_names[imu_xyz_idx[0]]
    imu_lin_yname = imu_linaccel_names[imu_xyz_idx[1]]
    imu_lin_zname = imu_linaccel_names[imu_xyz_idx[2]]
    imu_angaccel_names = ["rate_x", "rate_y", "rate_z"]
    imu_ang_xname = imu_angaccel_names[imu_xyz_idx[0]]
    imu_ang_yname = imu_angaccel_names[imu_xyz_idx[1]]
    imu_ang_zname = imu_angaccel_names[imu_xyz_idx[2]]
    # Preliminary data frame with input ordering
    motion3d_pre = pd.DataFrame({'acceleration_x': ec["acceleration_x"],
                                 'acceleration_y': ec["acceleration_y"],
                                 'acceleration_z': ec["acceleration_y"],
                                 'rate_phi': ec["rate_z"],
                                 'rate_theta': ec["rate_x"],
                                 'rate_shi': ec["rate_y"]})
    # Now we can reorder in to RHS and scale units
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
    del motion3d_pre
    wind = ec[["wind_speed_u", "wind_speed_v", "wind_speed_w"]].copy()
    # [Original comment: check for any significant number of 'NAN's (not
    # worried about the odd one scattered here and there)].  [Original
    # comment: set open flag if gt 2% of records are 'NAN']
    if (((ec.op_CO2_density.count() / float(ec_nrows)) < 0.98) or
        ((ec.op_H2O_density.count() / float(ec_nrows)) < 0.98) or
        ((ec.op_analyzer_status.count() / float(ec_nrows)) < 0.98)):
        open_flag = True
    # [Original comment: set wind flag if gt 2% of records are 'NAN']
    if (((wind.wind_speed_u.count() / float(ec_nrows)) < 0.98) or
        ((wind.wind_speed_v.count() / float(ec_nrows)) < 0.98) or
        ((wind.wind_speed_w.count() / float(ec_nrows)) < 0.98) or
        ((ec.air_temperature_sonic.count() / float(ec_nrows)) < 0.98)):
        sonic_flag = True
    # [Original comment: set motion flag if gt 2% of records are 'NAN']
    if (((motion3d.acceleration_x.count() / float(ec_nrows)) < 0.98) or
        ((motion3d.acceleration_y.count() / float(ec_nrows)) < 0.98) or
        ((motion3d.acceleration_z.count() / float(ec_nrows)) < 0.98) or
        ((motion3d.rate_phi.count() / float(ec_nrows)) < 0.98) or
        ((motion3d.rate_theta.count() / float(ec_nrows)) < 0.98) or
        ((motion3d.rate_shi.count() / float(ec_nrows)) < 0.98)):
        motion_flag = True
    # Set closed flag is more than 2% of records are NaN.
    if (((ec.cp_CO2_fraction.count() / float(ec_nrows)) < 0.98) or
        ((ec.cp_H2O_fraction.count() / float(ec_nrows)) < 0.98) or
        ((ec.cp_pressure.count() / float(ec_nrows)) < 0.98)):
        closed_flag = True

    # [Original comment: now that we have looked for NANs, we may as
    # well fill in the NANs and any spikes using the shot filter].
    # [SPL: these changes are done outside the WIND array, which is
    # the one that is used later for motion correction, etc., so they
    # are lost.]
    win_width = np.int(config["EC Despiking"]["despike_win_width"])
    win_step = np.int(config["EC Despiking"]["despike_step"])
    nreps = np.int(config["EC Despiking"]["despike_nreps"])

    if not sonic_flag:
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
            sonic_flag = True

    if not open_flag:
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
            open_flag = True

    if not closed_flag:
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
            closed_flag = True

    # TODO: Here we need to prepare our check for the diagnostics from the
    # open path analyzers.  For now, keep using the rule of thumb
    if (ec.op_analyzer_status.gt(249) |
        ec.op_analyzer_status.lt(240)).sum() > 0.02:
        open_flag = True

    # [Original comment: check for bad wind data: bad wind data can
    # usually be diagnosed by unusually high wind speeds.  This is
    # most obvious in the vertical wind where we wouldn't expect high
    # values bad sonic data can also turn up in the Tsonic before the
    # wind, check the deviation between Tsonic and mean air T.]
    air_temp_avg = ec["air_temperature"].mean()
    nbad_vertical_wind = abs(wind["wind_speed_w"]).gt(7).sum()
    nbad_air_temp_sonic = abs(ec["air_temperature_sonic"] -
                              air_temp_avg).gt(7).sum()
    # Set wind flag high if gt 0.5% of records are frost contaminated
    if ((nbad_vertical_wind / float(ec_nrows)) > 0.5 or
        (nbad_air_temp_sonic / float(ec_nrows)) > 0.5):
        sonic_flag = True
        raise FluxError("Bad sonic anemometer data")
    # [Original comment: check critical low frequency variabiles]
    if not (np.isfinite(air_temp_avg) or
            np.isfinite(ec.relative_humidity[0])):
        bad_meteorology_flag = True
        raise FluxError("RH or average air temperature unavailable")
    # # Below will be needed at some point
    # sw_avg = ec.K_down[0]
    # lw_avg = ec.LW_down[0]
    # sog_avg = ec["speed_over_ground"].mean()

    # [Original comment: now fill in the gaps by applying a moving
    # average... In this case, we use a 100 sample window (10 sec) moving
    # average... may need to tweak this value].  [SPL: perhaps a simple
    # linear interpolation is better; I don't know why this moving average
    # is used, where a window must be specified and may be introducing
    # bias.  Perhaps it doesn't matter.  Why aren't latitude and longitude
    # not similarly interpolated?]
    cog, sog = smooth_angle(ec["course_over_ground"].values,
                            ec["speed_over_ground"].values, 21)
    cog = pd.Series(cog, index=ec.index)
    sog = pd.Series(sog, index=ec.index)
    heading, _ = smooth_angle(ec["heading"].values, 1, 21)
    heading = pd.Series(heading, index=ec.index)

    if ((cog.count() < len(cog)) or (sog.count() < len(sog)) or
        (heading.count < len(heading))):
        motion_flag = True
    # If we have no good COG, SOG, or heading, then we cannot continue.
    if cog.count() < 1 or sog.count() < 1 or heading.count() < 1:
        bad_navigation_flag = True
        raise FluxError("Unusable COG, SOG, or heading records")
    # [Original comment: shot filter the motion channels... this helps with
    # a problem where unreasonably high accelerations cause a 'NaN'
    # calculation]
    for col in motion3d.columns:
        motcol_VM = despike_VickersMahrt(motion3d[col],
                                         width=win_width, step=win_step,
                                         zscore_thr=3.5, nreps=nreps)
        motion3d[col] = motcol_VM[0]

    # Tilt angles for sonic and IMU
    wnd_k, wnd_kcoefs = planarfit_coef(wind.values)
    mot3d_k, mot3d_kcoefs = planarfit_coef(motion3d.values)
    wndrot_dummy, wnd_phitheta = rotate_vectors(wind.values,
                                                k_vector=wnd_k)
    mot3d_dummy, mot3d_phitheta = rotate_vectors(motion3d.values[:, 0:3],
                                                 k_vector=mot3d_k)

    # # Output to Octave for debugging
    # import scipy.io as sio
    # sio.savemat(iname_prefix + "_wind_motion.mat",
    #             {'wind_speed': wind.values,
    #              'acceleration': motion3d.loc[:, :"acceleration_z"].values,
    #              'angular_rate': motion3d.loc[:, "rate_phi":].values,
    #              'heading': np.reshape(heading, (len(heading), 1)),
    #              'sog': np.reshape(sog, (len(sog), 1))})
    # # OR do a round-trip via Oct2Py!
    # from oct2py import octave
    # uvw_ship, uvw_earth = octave.motion_octave(wind.values,
    #                                            motion3d.values[:, :3],
    #                                            motion3d.values[:, 3:],
    #                                            np.reshape(heading,
    #                                                       (len(heading), 1)),
    #                                            np.reshape(sog, (len(sog), 1)),
    #                                            imu2anem_pos,
    #                                            sample_freq_hz, Tcf, Ta,
    #                                            [0.0, 0.0], [0.0, 0.0],
    #                                            {"uearth"})

    # Save full tuple output and select later. Note that we the use the
    # interpolated, smoothed heading and speed over ground.
    UVW = wind3D_correct(wind.values,
                         motion3d.loc[:, :"acceleration_z"].values,
                         motion3d.loc[:, "rate_phi":].values,
                         heading.values, sog.values, imu2anem_pos,
                         sample_freq_hz, Tcf, Ta, [0.0, 0.0], [0.0, 0.0])
    # Ship-referenced speeds
    UVW_ship = UVW[0]
    # Earth-referenced speeds
    UVW_earth = UVW[11]
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
    # Repeat applying instrument tilt angle correction, assuming both sonic
    # and IMU have the same angles
    UVW_tilt = wind3D_correct(wind.values,
                              motion3d.loc[:, :"acceleration_z"].values,
                              motion3d.loc[:, "rate_phi":].values,
                              heading.values, sog.values, imu2anem_pos,
                              sample_freq_hz, Tcf, Ta, mot3d_phitheta,
                              mot3d_phitheta)
    UVW_ship_tilt = UVW_tilt[0]
    UVW_earth_tilt = UVW_tilt[11]
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

    # # Plots comparing with Octave output
    # fig, axs = plt.subplots(3, 2, sharex = True)
    # axs[0, 0].plot(uvw_ship[:, 0], UVW_ship[:, 0])
    # # etc, etc.
    # # Plot smoothed and corrected data for each vector. Turned off for
    # # production of output files.
    # fig, axs = plt.subplots(3, 1, sharex=True)
    # fig.set_size_inches((11, 9))
    # wind[["wind_speed_u", "wind_speed_u_corr"]].plot(ax=axs[0],
    #                                                  legend=False)
    # axs[0].set_title("U"); axs[0].set_xlabel("")
    # wind[["wind_speed_v", "wind_speed_v_corr"]].plot(ax=axs[1],
    #                                                  legend=False)
    # axs[1].set_title("V")
    # axs[1].set_ylabel("Wind speed [m/s/s]"); axs[1].set_xlabel("")
    # wind[["wind_speed_w", "wind_speed_w_corr"]].plot(ax=axs[2],
    #                                                  rot=0,
    #                                                  legend=False)
    # axs[2].set_title("W"); axs[2].set_xlabel("")
    # leg = axs[2].legend(loc=9, bbox_to_anchor=(0.5, -0.1),
    #                     title=iname_prefix, frameon=False,
    #                     borderaxespad=0, ncol=3)
    # leg.get_texts()[0].set_text("Measured (de-spiked)")
    # leg.get_texts()[1].set_text("Corrected (ship reference)")
    # # axs[1].set_xticklabels(wind.index, rotation=0, ha="center")
    # fig.tight_layout()
    # fig.savefig(iname_prefix + ".png", bbox_extra_artists=(leg,),
    #             bbox_inches="tight")

    # Append results
    ec_wind_corr = pd.concat((ec, wind.loc[:, wind_corr_names[0]:]),
                             axis=1)
    return ec_wind_corr, dict(open_flag=open_flag, closed_flag=closed_flag,
                              sonic_flag=sonic_flag, motion_flag=motion_flag,
                              bad_navigation_flag=bad_navigation_flag,
                              bad_meteorology_flag=bad_meteorology_flag)


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
    summary_file = config["EC Outputs"]["summary_file"]
    # Stop if we don't have any files
    if len(ec_files) < 1:
        raise FluxError("There are no input files")

    # [Original comment: create flags for the 4 possible sources of "bad"
    # data, flag=0 means data good]
    flags = dict.fromkeys(["open_flag", "closed_flag", "sonic_flag",
                           "motion_flag", "bad_navigation_flag",
                           "bad_meteorology_flag"], False)
    # We set up a dataframe with all files to process as index, and all the
    # flags as columns.  This is the basis for our summary output file;
    # other columns (such as flux summary calculations for the period)
    # will be appended as we loop.
    osummary = pd.DataFrame(flags,
                            index=[osp.basename(x) for x in ec_files])
    for ec_file in ec_files:
        # print(ec_file)             # REMOVE FOR PRODUCTION
        # Get a file name prefix to be shared by the output files from this
        # period.  Note iname is THE SAME AS THE INDEX IN OSUMMARY
        iname = osp.basename(ec_file)
        iname_prefix = osp.splitext(iname)[0]

        try:
            ec_wind_corr, ec_flags = flux_period(ec_file, config)
        except FluxError as e:
            print("{0}: {1}".format(ec_file, e.message))
        else:
            # Save to file with suffix "_mc.csv"
            ec_wind_corr.to_csv(osp.join(ec_idir, iname_prefix + "_mc.csv"),
                                index_label=colnames[1])
            for k, v in ec_flags.iteritems():
                osummary.loc[iname, k] = v

    # Now we have the summary DataFrame filled up and can work with it.
    osummary.to_csv(summary_file, index_label="input_file")
    print "Summary of fluxes written to " + summary_file


if __name__ == "__main__":
    description = "Perform flux analyses, given a configuration file"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("config_file", type=str,
                        help="Path to configuration file")
    args = parser.parse_args()
    main(args.config_file)
