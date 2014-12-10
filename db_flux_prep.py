#! /usr/bin/env python

"""Steps towards CO2 flux analyses, using data files output from PostgreSQL
database.

"""

import numpy as np
from flux import shot_filter, smooth_angle, wind3D_correct
import pandas as pd
import os.path as osp
import glob
# import psycopg2 as pg
# from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

plt.style.use("ggplot")

ec_files = glob.glob("/home/sluque/Data/ArcticNet/2011/FromDB/EC_*.csv")
# ec_files = glob.glob("/home/sluque/tmp/EC_*.csv")
ec_files.sort()
# colnames = ["time_20min", "time_study", "longitude", "latitude",
#             "speed_over_ground", "course_over_ground", "heading",
#             "pitch", "roll", "heave", "atmospheric_pressure", 
#             "air_temperature", "relative_humidity", "surface_temperature",
#             "wind_speed", "wind_direction", "true_wind_speed",
#             "true_wind_direction", "PAR", "K_down",
#             "temperature_thermopile", "temperature_case",
#             "temperature_dome", "LW_down", "UV_sensor_temperature",
#             "UV_b", "UV_broad", "acceleration_x", "acceleration_y",
#             "acceleration_z", "rate_x", "rate_y", "rate_z", "wind_speed_u", 
#             "wind_speed_v", "wind_speed_w", "air_temperature_sonic",
#             "sound_speed", "anemometer_status", "op_CO2_fraction",
#             "op_CO2_density", "op_CO2_absorptance", "op_H2O_fraction",
#             "op_H2O_density", "op_H2O_absorptance", "op_pressure", 
#             "op_temperature", "op_temperature_base", "op_temperature_spar",
#             "op_temperature_bulb", "op_cooler_voltage", "op_bandwidth",
#             "op_delay_interval", "bad_chopper_wheel_temperature_flag", 
#             "bad_detector_temperature_flag",
#             "bad_optical_wheel_rotation_rate_flag", 
#             "bad_sync_flag", "op_AGC", "cp_analyzer_status",
#             "cp_CO2_fraction", "cp_CO2_density", "cp_CO2_dry_fraction",
#             "cp_CO2_absorptance", "cp_H2O_fraction", "cp_H2O_density",
#             "cp_H2O_dry_fraction", "cp_H2O_absorptance", 
#             "cp_pressure", "cp_temperature", "cp_temperature_in",
#             "cp_temperature_out", "cp_temperature_block",
#             "cp_temperature_cell", "cp_CO2_signal_strength", 
#             "cp_H2O_signal_strength"]
colnames = ["time_20min", "time_study", "longitude", "latitude",
            "speed_over_ground", "course_over_ground", "heading",
            "pitch", "roll", "heave", "atmospheric_pressure", 
            "air_temperature", "relative_humidity", "surface_temperature",
            "wind_speed", "wind_direction", "true_wind_speed",
            "true_wind_direction", "PAR", "K_down",
            "temperature_thermopile", "temperature_case",
            "temperature_dome", "LW_down", "UV_sensor_temperature",
            "UV_b", "UV_broad", "acceleration_x", "acceleration_y",
            "acceleration_z", "rate_x", "rate_y", "rate_z", "wind_speed_u", 
            "wind_speed_v", "wind_speed_w", "air_temperature_sonic",
            "sound_speed", "anemometer_status", "op_analyzer_status",
            "op_CO2_fraction", "op_CO2_density", "op_CO2_absorptance",
            "op_H2O_fraction", "op_H2O_density", "op_H2O_absorptance",
            "op_pressure", "op_temperature", "op_temperature_base",
            "op_temperature_spar", "op_temperature_bulb",
            "op_cooler_voltage", "op_CO2_signal_strength", "op_bandwidth",
            "op_delay_interval", "cp_analyzer_status", "cp_CO2_fraction",
            "cp_CO2_density", "cp_CO2_dry_fraction", "cp_CO2_absorptance",
            "cp_H2O_fraction", "cp_H2O_density", "cp_H2O_dry_fraction",
            "cp_H2O_absorptance", "cp_pressure", "cp_temperature",
            "cp_temperature_in", "cp_temperature_out",
            "cp_temperature_block", "cp_temperature_cell",
            "cp_CO2_signal_strength", "cp_H2O_signal_strength"]

# ppdf = PdfPages("winds_2011.pdf")
for ec_period in ec_files[0:3]:
    print ec_period             # REMOVE FOR PRODUCTION
    ec = pd.read_csv(ec_period, dtype=np.float, parse_dates=[0, 1],
                     index_col=1, names=colnames)
    ec_nrows = len(ec.index)
    # Get a file name prefix to be shared by the output files from this
    # period
    iname = osp.basename(ec_period)
    iname_prefix = osp.splitext(iname)[0]
    # Put acceleration components in 3-column array and make copy to keep
    # uncorrected data.  Original comment: read in angular rates in RH
    # coordinate system, convert to rad/s.
    motion3d = pd.DataFrame({"acceleration_x" : ec["acceleration_z"] * 9.81,
                             "acceleration_y" : -ec["acceleration_x"] * 9.81,
                             "acceleration_z" : -ec["acceleration_y"] * 9.81,
                             "rate_phi" : np.radians(ec["rate_z"]),
                             "rate_theta" : np.radians(ec["rate_x"]),
                             "rate_shi" : np.radians(ec["rate_y"])})
    wind = ec[["wind_speed_u", "wind_speed_v", "wind_speed_w"]].copy()
    # [Original comment: create flags for the 4 possible sources of "bad"
    # data, flag=0 means data good]
    flags = dict.fromkeys(["open_flag", "closed_flag",
                           "sonic_flag", "motion_flag"], False)
    # [Original comment: check for any significant number of 'NAN's (not
    # worried about the odd one scattered here and there)].  [Original
    # comment: set open flag if gt 2% of records are 'NAN']
    if (((ec.op_CO2_density.count() / ec_nrows) < 0.98) or
        ((ec.op_H2O_density.count() / ec_nrows) < 0.98) or
        ((ec.op_analyzer_status.count() / ec_nrows) < 0.98)):
        flags["open_flag"] = True
    # [Original comment: set wind flag if gt 2% of records are 'NAN']
    if (((wind.wind_speed_u.count() / ec_nrows) < 0.98) or
        ((wind.wind_speed_v.count() / ec_nrows) < 0.98) or
        ((wind.wind_speed_w.count() / ec_nrows) < 0.98) or
        ((ec.air_temperature_sonic.count() / ec_nrows) < 0.98)):
        flags["sonic_flag"] = True
    # [Original comment: set motion flag if gt 2% of records are 'NAN']
    if (((motion3d.acceleration_x.count() / ec_nrows) < 0.98) or
        ((motion3d.acceleration_y.count() / ec_nrows) < 0.98) or
        ((motion3d.acceleration_z.count() / ec_nrows) < 0.98) or
        ((motion3d.rate_phi.count() / ec_nrows) < 0.98) or
        ((motion3d.rate_theta.count() / ec_nrows) < 0.98) or
        ((motion3d.rate_shi.count() / ec_nrows) < 0.98)):
        flags["motion_flag"] = True

    # [Original comment: now that we have looked for NANs, we may as
    # well fill in the NANs and any spikes using the shot filter].
    # [SPL: these changes are done outside the WIND array, which is
    # the one that is used later for motion correction, etc., so they
    # are lost.]
    if not flags["sonic_flag"]:
        wind = wind.apply(shot_filter)
        ec.air_temperature_sonic = shot_filter(ec.air_temperature_sonic)

    if not flags["open_flag"]:
        ec["op_CO2_density"] = shot_filter(ec["op_CO2_density"])
        ec["op_H2O_density"] = shot_filter(ec["op_H2O_density"])
        ec["op_pressure"] = shot_filter(ec["op_pressure"])
    # [Original comment: this is necessary to check if there is still ugly
    # shot noise... if there is, we need to skip this]
    if any((abs(ec["op_CO2_density"] - np.mean(ec["op_CO2_density"]))) >
           (6 * np.std(ec["op_CO2_density"]))):
        flags["open_flag"] = True

    if not flags["closed_flag"]:
        ec["cp_CO2_fraction"] = shot_filter(ec["cp_CO2_fraction"])
        ec["cp_H2O_fraction"] = shot_filter(ec["cp_H2O_fraction"])
        ec["cp_pressure"] = shot_filter(ec["cp_pressure"])

    if (any((abs(ec["cp_CO2_fraction"] - np.mean(ec["cp_CO2_fraction"]))) >
            (6 * np.std(ec["cp_CO2_fraction"]))) or
        any((abs(ec["cp_H2O_fraction"] - np.mean(ec["cp_H2O_fraction"]))) >
            (6 * np.std(ec["cp_H2O_fraction"])))):
        flags["closed_flag"] = True

    # TODO: Here we need to prepare our check for the diagnostics from the
    # open path analyzers.  For now, keep using the rule of thumb
    # if ((ec.op_analyzer_status > 249) or (ec.op_analyzer_status < 240)):
    #     flags["open_flag"] = True
    
    # [Original comment: check for bad wind data: bad wind data can
    # usually be diagnosed by unusually high wind speeds.  this is
    # most obvious in the vertical wind where we wouldn't expect high
    # values bad sonic data can also turn up in the Tsonic before the
    # wind, check the deviation between Tsonic and mean air T.]
    air_temp_avg = np.mean(ec.air_temperature)
    nbad_vertical_wind = np.sum(wind["wind_speed_w"] > 7)
    nbad_air_temp_sonic = np.sum(abs(ec["air_temperature_sonic"] -
                                     air_temp_avg) > 7)
    if ((nbad_vertical_wind / ec_nrows) > 0.5 or
        (nbad_air_temp_sonic / ec_nrows) > 0.5):
        flags["sonic_flag"] = True
    # Set wind flag high if gt 0.5% of records are frost contaminated
    if flags["sonic_flag"]:
        print "Bad sonic anemometer data. Skipping."
        continue
    # [Original comment: check critical low frequency variabiles]
    if not (np.isfinite(air_temp_avg) or
            np.isfinite(ec.relative_humidity[0])):
        print "RH or average air temperature unavailable. Skipping."
        continue
    sw_avg = ec.K_down[0]
    lw_avg = ec.LW_down[0]
    sog_avg = np.mean(ec.speed_over_ground)

    # [Original comment: now fill in the gaps by applying a moving
    # average... In this case, we use a 100 sample window (10 sec) moving
    # average... may need to tweak this value].  [SPL: perhaps a simple
    # linear interpolation is better; I don't know why this moving average
    # is used, where a window must be specified and may be introducing
    # bias.  Perhaps it doesn't matter.  Why aren't latitude and longitude
    # not similarly interpolated?]
    cog, sog = smooth_angle(ec.course_over_ground.values,
                            ec.speed_over_ground.values, 21)
    cog = pd.Series(cog, index=ec.index)
    sog = pd.Series(sog, index=ec.index)
    heading, _ = smooth_angle(ec.heading.values, 1, 21)
    heading = pd.Series(heading, index=ec.index)

    if ((cog.count() < len(cog)) or (sog.count() < len(sog)) or
        (heading.count < len(heading))):
        flags["motion_flag"] = True
    # If we have no good COG, SOG, or heading, then we should skip
    # processing entirely.
    if cog.count() < 1 or sog.count() < 1 or heading.count() < 1:
        print "Unusable COG, SOG, or heading records. Skipping."
        continue
    # [Original comment: shot filter the motion channels... this helps with
    # a problem where unreasonably high accelerations cause a 'NaN'
    # calculation]
    motion3d = motion3d.apply(shot_filter)

    # # Output to Octave for debugging
    # import scipy.io as sio
    # sio.savemat(iname_prefix + '_wind_motion.mat',
    #             {"wind_speed": wind.values,
    #              "acceleration": motion3d.iloc[:, :3].values,
    #              "angular_rate": motion3d.iloc[:, 3:].values,
    #              "heading": heading.values.T, "sog": sog.values.T})

    # Save full tuple output and select later
    UVW = wind3D_correct(wind.values, motion3d.iloc[:, :3].values,
                         motion3d.iloc[:, 3:].values,
                         heading.values, sog.values, [1.7, 0, 2.725],
                         10.0, 10.0, 20.0, [0.0, 0.0], [0.0, 0.0])
    # Ship-referenced speeds
    uvw_ship = UVW[0]
    # Append corrected wind vectors to DataFrame
    wind_corr_names = ["wind_speed_u_corr", "wind_speed_v_corr",
                       "wind_speed_w_corr"]
    wind[wind_corr_names[0]] = pd.Series(UVW[11][:, 0], index=wind.index)
    wind[wind_corr_names[1]] = pd.Series(UVW[11][:, 1], index=wind.index)
    wind[wind_corr_names[2]] = pd.Series(UVW[11][:, 2], index=wind.index)

    # Plot smoothed and corrected data for each vector
    fig, axs = plt.subplots(3, 1, sharex=True)
    fig.set_size_inches((11, 9))
    wind[["wind_speed_u", "wind_speed_u_corr"]].plot(ax=axs[0],
                                                     legend=False)
    axs[0].set_title("U"); axs[0].set_xlabel('')
    wind[["wind_speed_v", "wind_speed_v_corr"]].plot(ax=axs[1],
                                                     legend=False)
    axs[1].set_title("V")
    axs[1].set_ylabel("Wind speed [m/s/s]"); axs[1].set_xlabel('')
    wind[["wind_speed_w", "wind_speed_w_corr"]].plot(ax=axs[2],
                                                     rot=0,
                                                     legend=False)
    axs[2].set_title("W"); axs[2].set_xlabel('')
    leg = axs[2].legend(loc=9, bbox_to_anchor=(0.5, -0.1),
                        title=iname_prefix, frameon=False,
                        borderaxespad=0, ncol=3)
    # axs[1].set_xticklabels(wind.index, rotation=0, ha="center")
    fig.tight_layout()
    fig.savefig(iname_prefix + ".png", bbox_extra_artists=(leg,),
                bbox_inches="tight")

# ppdf.close()



## TESTS ------------------------------------------------------------------

# angdeg = np.linspace(0, 360, 13)
# angdc = decompose(angdeg, 1)
# print recompose(angdc["x"], angdc["y"]) # all good

# angdeg = np.linspace(0, 365, 13)
# angdc = decompose(angdeg, 1)
# print recompose(angdc["x"], angdc["y"]) # wrapping taken care of

# cog_xy = smooth_angle(ec.course_over_ground.values,
#                       ec.speed_over_ground.values, 21) # catch at least 2 s
# plt.plot(ec.course_over_ground, 'b.')
# plt.plot(cog_xy["angle"])

# heading = smooth_angle(ec.heading.values, 1, 21)
# plt.plot(ec.heading, 'b.')
# plt.plot(heading["angle"])
