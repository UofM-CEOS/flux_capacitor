#! /usr/bin/env python
# $Id$

"""Steps towards CO2 flux analyses, using data files output from PostgreSQL
database.

"""

import numpy as np
from scipy.stats import zscore
from flux import shot_filter, smooth_angle, wind3D_correct
from flux_config import parse_config
import pandas as pd
import os.path as osp
# import psycopg2 as pg
# from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

plt.style.use('ggplot')

def do_flux(period_file, config):
    """Perform required calculations on period."""
    # Extract all the config pieces
    colnames = config['Inputs']['colnames']
    mot2anem_pos = config['Motion Correction']['motion2anemometer_pos']
    sample_freq_hz = config['Inputs']['sample_freq']
    Tcf = config['Motion Correction']['complementary_filter_period']
    Ta = config['Motion Correction']['accel_highpass_cutoff']
    ec = pd.read_csv(period_file, dtype=np.float, parse_dates=[0, 1],
                     index_col=1, names=colnames)
    ec_nrows = len(ec.index)
    # Put acceleration components in 3-column array and make copy to keep
    # uncorrected data.  Original comment: read in angular rates in RH
    # coordinate system, convert to rad/s.
    motion3d = pd.DataFrame({'acceleration_x' : ec['acceleration_z'] * 9.81,
                             'acceleration_y' : -ec['acceleration_x'] * 9.81,
                             'acceleration_z' : -ec['acceleration_y'] * 9.81,
                             'rate_phi' : np.radians(ec['rate_z']),
                             'rate_theta' : np.radians(ec['rate_x']),
                             'rate_shi' : np.radians(ec['rate_y'])})
    wind = ec[['wind_speed_u', 'wind_speed_v', 'wind_speed_w']].copy()
    # [Original comment: check for any significant number of 'NAN's (not
    # worried about the odd one scattered here and there)].  [Original
    # comment: set open flag if gt 2% of records are 'NAN']
    if (((ec.op_CO2_density.count() / ec_nrows) < 0.98) or
        ((ec.op_H2O_density.count() / ec_nrows) < 0.98) or
        ((ec.op_analyzer_status.count() / ec_nrows) < 0.98)):
        open_flag = True
    # [Original comment: set wind flag if gt 2% of records are 'NAN']
    if (((wind.wind_speed_u.count() / ec_nrows) < 0.98) or
        ((wind.wind_speed_v.count() / ec_nrows) < 0.98) or
        ((wind.wind_speed_w.count() / ec_nrows) < 0.98) or
        ((ec.air_temperature_sonic.count() / ec_nrows) < 0.98)):
        sonic_flag = True
    # [Original comment: set motion flag if gt 2% of records are 'NAN']
    if (((motion3d.acceleration_x.count() / ec_nrows) < 0.98) or
        ((motion3d.acceleration_y.count() / ec_nrows) < 0.98) or
        ((motion3d.acceleration_z.count() / ec_nrows) < 0.98) or
        ((motion3d.rate_phi.count() / ec_nrows) < 0.98) or
        ((motion3d.rate_theta.count() / ec_nrows) < 0.98) or
        ((motion3d.rate_shi.count() / ec_nrows) < 0.98)):
        motion_flag = True

    # [Original comment: now that we have looked for NANs, we may as
    # well fill in the NANs and any spikes using the shot filter].
    # [SPL: these changes are done outside the WIND array, which is
    # the one that is used later for motion correction, etc., so they
    # are lost.]
    if not sonic_flag:
        wind = wind.apply(shot_filter)
        ec.air_temperature_sonic = shot_filter(ec.air_temperature_sonic)

    if not open_flag:
        ec['op_CO2_density'] = shot_filter(ec['op_CO2_density'])
        ec['op_H2O_density'] = shot_filter(ec['op_H2O_density'])
        ec['op_pressure'] = shot_filter(ec['op_pressure'])
    # [Original comment: this is necessary to check if there is still ugly
    # shot noise... if there is, we need to skip this]
    if any(abs(zscore(ec['op_CO2_density'])) > 6):
        open_flag = True

    if not closed_flag:
        ec['cp_CO2_fraction'] = shot_filter(ec['cp_CO2_fraction'])
        ec['cp_H2O_fraction'] = shot_filter(ec['cp_H2O_fraction'])
        ec['cp_pressure'] = shot_filter(ec['cp_pressure'])

    if (any(abs(zscore(ec['cp_CO2_fraction'])) > 6) or
        any(abs(zscore(ec['cp_H2O_fraction'])) > 6)):
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
    air_temp_avg = ec['air_temperature'].mean()
    nbad_vertical_wind = wind['wind_speed_w'].gt(7).sum()
    nbad_air_temp_sonic = abs(ec['air_temperature_sonic'] -
                              air_temp_avg).gt(7).sum()
    # Set wind flag high if gt 0.5% of records are frost contaminated
    if ((nbad_vertical_wind / ec_nrows) > 0.5 or
        (nbad_air_temp_sonic / ec_nrows) > 0.5):
        sonic_flag = True
        print "Bad sonic anemometer data. Skipping."
        continue
    # [Original comment: check critical low frequency variabiles]
    if not (np.isfinite(air_temp_avg) or
            np.isfinite(ec.relative_humidity[0])):
        print "RH or average air temperature unavailable. Skipping."
        bad_meteorology_flag = True
        continue
    sw_avg = ec.K_down[0]
    lw_avg = ec.LW_down[0]
    sog_avg = ec['speed_over_ground'].mean()

    # [Original comment: now fill in the gaps by applying a moving
    # average... In this case, we use a 100 sample window (10 sec) moving
    # average... may need to tweak this value].  [SPL: perhaps a simple
    # linear interpolation is better; I don't know why this moving average
    # is used, where a window must be specified and may be introducing
    # bias.  Perhaps it doesn't matter.  Why aren't latitude and longitude
    # not similarly interpolated?]
    cog, sog = smooth_angle(ec['course_over_ground'].values,
                            ec['speed_over_ground'].values, 21)
    cog = pd.Series(cog, index=ec.index)
    sog = pd.Series(sog, index=ec.index)
    heading, _ = smooth_angle(ec['heading'].values, 1, 21)
    heading = pd.Series(heading, index=ec.index)

    if ((cog.count() < len(cog)) or (sog.count() < len(sog)) or
        (heading.count < len(heading))):
        motion_flag = True
    # If we have no good COG, SOG, or heading, then we should skip
    # processing entirely.
    if cog.count() < 1 or sog.count() < 1 or heading.count() < 1:
        print "Unusable COG, SOG, or heading records. Skipping."
        bad_navigation_flag = True
        continue
    # [Original comment: shot filter the motion channels... this helps with
    # a problem where unreasonably high accelerations cause a 'NaN'
    # calculation]
    motion3d = motion3d.apply(shot_filter)

    # # Output to Octave for debugging
    # import scipy.io as sio
    # sio.savemat(iname_prefix + '_wind_motion.mat',
    #             {'wind_speed': wind.values,
    #              'acceleration': motion3d.loc[:, :'acceleration_z'].values,
    #              'angular_rate': motion3d.loc[:, 'rate_phi':].values,
    #              'heading': np.reshape(heading, (len(heading), 1)),
    #              'sog': np.reshape(sog, (len(sog), 1))})

    # Save full tuple output and select later. Note that we the use the
    # interpolated, smoothed heading and speed over ground.
    UVW = wind3D_correct(wind.values,
                         motion3d.loc[:, :'acceleration_z'].values,
                         motion3d.loc[:, 'rate_phi':].values,
                         heading.values, sog.values, mot2anem_pos,
                         sample_freq_hz, Tcf, Ta, [0.0, 0.0], [0.0, 0.0])
    # Ship-referenced speeds
    uvw_ship = UVW[0]
    # Earth-referenced speeds
    uvw_earth = UVW[11]
    # Append corrected wind vectors to DataFrame
    wind_corr_names = ['wind_speed_u_ship', 'wind_speed_v_ship',
                       'wind_speed_w_ship', 'wind_speed_u_earth',
                       'wind_speed_v_earth', 'wind_speed_w_earth']
    wind[wind_corr_names[0]] = pd.Series(uvw_ship[:, 0], index=wind.index)
    wind[wind_corr_names[1]] = pd.Series(uvw_ship[:, 1], index=wind.index)
    wind[wind_corr_names[2]] = pd.Series(uvw_ship[:, 2], index=wind.index)
    wind[wind_corr_names[3]] = pd.Series(uvw_earth[:, 0],
                                         index=wind.index)
    wind[wind_corr_names[4]] = pd.Series(uvw_earth[:, 1],
                                         index=wind.index)
    wind[wind_corr_names[5]] = pd.Series(uvw_earth[:, 2],
                                         index=wind.index)

    # # Plot smoothed and corrected data for each vector. Turned off for
    # # production of output files.
    # fig, axs = plt.subplots(3, 1, sharex=True)
    # fig.set_size_inches((11, 9))
    # wind[['wind_speed_u', 'wind_speed_u_corr']].plot(ax=axs[0],
    #                                                  legend=False)
    # axs[0].set_title("U"); axs[0].set_xlabel('')
    # wind[['wind_speed_v', 'wind_speed_v_corr']].plot(ax=axs[1],
    #                                                  legend=False)
    # axs[1].set_title("V")
    # axs[1].set_ylabel("Wind speed [m/s/s]"); axs[1].set_xlabel('')
    # wind[['wind_speed_w', 'wind_speed_w_corr']].plot(ax=axs[2],
    #                                                  rot=0,
    #                                                  legend=False)
    # axs[2].set_title("W"); axs[2].set_xlabel('')
    # leg = axs[2].legend(loc=9, bbox_to_anchor=(0.5, -0.1),
    #                     title=iname_prefix, frameon=False,
    #                     borderaxespad=0, ncol=3)
    # leg.get_texts()[0].set_text("Measured (de-spiked)")
    # leg.get_texts()[1].set_text("Corrected (ship reference)")
    # # axs[1].set_xticklabels(wind.index, rotation=0, ha='center')
    # fig.tight_layout()
    # fig.savefig(iname_prefix + '.png', bbox_extra_artists=(leg,),
    #             bbox_inches='tight')

    # Append results
    ec_wind_corr = pd.concat((ec, wind.loc[:, 'wind_speed_u_ship':]),
                             axis=1)
    return ec_wind_corr, dict(open_flag=open_flag, closed_flag=closed_flag,
                              sonic_flag=sonic_flag, motion_flag,
                              bad_navigation_flag=bad_navigation_flag,
                              bad_meteorology_flag=bad_meteorology_flag)
    

# Parse configuration file
config = parse_config(config_file)
# Stop if we don't have any files
if (len(config["Inputs"]["input_files"]) < 1):
    raise Exception("There are no input files")

# [Original comment: create flags for the 4 possible sources of "bad"
# data, flag=0 means data good]
flags = dict.fromkeys(['open_flag', 'closed_flag', 'sonic_flag',
                       'motion_flag', 'bad_navigation_flag',
                       'bad_meteorology_flag'], False)
# We set up a dataframe with all files to process as index, and all the
# flags as columns.  This is the basis for our summary output file; other
# columns (such as flux summary calculations for the period) will be
# appended as we loop.
osummary = pd.DataFrame(flags, index=[osp.basename(x) for x in ec_files])

for ec_file in ec_files[0:5]:
    print ec_file             # REMOVE FOR PRODUCTION
    # Get a file name prefix to be shared by the output files from this
    # period.  Note iname is THE SAME AS THE INDEX IN OSUMMARY
    iname = osp.basename(ec_file)
    iname_prefix = osp.splitext(iname)[0]
    ec_wind_corr, ec_flags = do_flux(ec_file, config)
    # Save to file with suffix '_mc.csv'
    ec_wind_corr.to_csv(osp.join(ec_idir, iname_prefix + '_mc.csv'),
                        index_label=colnames[1])


# Now we have the summary file is filled up and can work with it.
osummary.to_csv(osummary_fname, index_label="input_file")
print "Summary of fluxes written to " + osummary_fname

## TESTS ------------------------------------------------------------------

# angdeg = np.linspace(0, 360, 13)
# angdc = decompose(angdeg, 1)
# print recompose(angdc['x'], angdc['y']) # all good

# angdeg = np.linspace(0, 365, 13)
# angdc = decompose(angdeg, 1)
# print recompose(angdc['x'], angdc['y']) # wrapping taken care of

# cog_xy = smooth_angle(ec.course_over_ground.values,
#                       ec.speed_over_ground.values, 21) # catch at least 2 s
# plt.plot(ec.course_over_ground, 'b.')
# plt.plot(cog_xy['angle'])

# heading = smooth_angle(ec.heading.values, 1, 21)
# plt.plot(ec.heading, 'b.')
# plt.plot(heading['angle'])

# # Replicate Miller's test code
# import scipy.io as sio
# miller_objs = sio.loadmat('motiondata.mat')
# UVW = wind3D_correct(miller_objs['uvwm'], miller_objs['acc'],
#                      miller_objs['rate'], np.squeeze(miller_objs['heading']),
#                      np.squeeze(miller_objs['speed']),
#                      np.squeeze(miller_objs['r']),
#                      np.squeeze(miller_objs['sf']), 10, 20, [0, 0], [0, 0])

# plt.plot(miller_objs['uvwm'][:, 0])
# plt.plot(UVW[0][:, 0])
