#! /usr/bin/env python
# pylint: disable=too-many-locals,invalid-name,no-member

"""Calculate tilt angles in flux files using non-overlapping moving windows

Usage
-----

For help using this script, type:

tilt_search.py -h

at command line.

"""

import sys
import itertools as itrt
import datetime
import logging
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import circmean
from fluxer import eddycov
from fluxer.eddycov.flux import (planarfit_coef, rotate_vectors)
from fluxer.eddycov.db_flux import (NavigationError, SonicError)
from fluxer.flux_config import parse_config


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def fname2time(fname):
    """Construct a time stamp from a file name"""
    return datetime.datetime.strptime(fname[-18:-4], "%Y%m%d%H%M%S")


def make_key(delta_t):
    """Return a function with attribute to track current upper bound"""
    def key(fname):
        """Set upper bound for current time window given a file name"""
        if key.upper_bound is None:
            key.upper_bound = fname2time(fname) + delta_t
        else:
            tstamp = fname2time(fname)
            while tstamp >= key.upper_bound:
                key.upper_bound += delta_t
        return key.upper_bound
    key.upper_bound = None
    return key


def main(config_file, win_hours, win_minutes, outfile):
    """Perform calculation of tilt angles and write output to file"""
    config = parse_config(config_file)
    ec_files = config["EC Inputs"]["input_files"]

    # Time delta
    tdelta = datetime.timedelta(hours=win_hours, minutes=win_minutes)
    # Build a list of tuples (key, file-list)
    win_files = [(pd.to_datetime(k), list(grp)) for k, grp in
                 itrt.groupby(ec_files, key=make_key(tdelta))]
    # Set up dataframe to place tilt angles and other details
    cols = (["nfiles", "nfiles_ok"] +
            ["n{0}".format(i) for i in eddycov.db_flux._FLUX_FLAGS] +
            ["phi_motion", "theta_motion", "phi_sonic", "theta_sonic"] +
            ["wind_direction"])
    tilts = pd.DataFrame(index=[x[0] for x in win_files], columns=cols)
    tilts.nfiles = [len(x[1]) for x in win_files]
    # Zero all the columns except the nfiles and parameters
    tilts.loc[:, cols[1:-5]] = 0

    # Iterate over the list of tuples (key, file-list)
    for (k, l) in win_files:
        logger.info("Begin window %s", k)
        ec_list = []
        for f in l:
            logger.info("Begin preparing %s", f)
            try:
                ec_prep, prep_flags = eddycov.prepare_period(f, config)
            except (NavigationError, SonicError) as err:
                for flagname, flag in err.flags.items():
                    tilts.loc[k, "n" + flagname] += np.int(flag)
                logger.info("Skip %s", f)
                continue
            else:
                for flagname, flag in prep_flags.iteritems():
                    tilts.loc[k, "n" + flagname] += np.int(flag)
                ec_list.append(ec_prep)
            logger.info("End preparing %s", f)

        n_ok = len(ec_list)
        tilts.loc[k, "nfiles_ok"] = n_ok
        if n_ok > 0:
            ec_win = pd.concat(ec_list)
            logger.info("Combining %s periods in window %s", n_ok, k)
            wind_names = ["wind_speed_u", "wind_speed_v", "wind_speed_w"]
            wnd = ec_win.loc[:, wind_names]
            mot3d_names = ["acceleration_x", "acceleration_y",
                           "acceleration_z", "rate_x", "rate_y", "rate_z"]
            mot3d = ec_win[mot3d_names]
            # Calculate tilt angles for motion sensor
            mot3d_k, _ = planarfit_coef(mot3d.values[:, 0:3])
            _, mot3d_phitheta = rotate_vectors(mot3d.values[:, 0:3],
                                               k_vector=mot3d_k)
            tilts.loc[k, "phi_motion"] = mot3d_phitheta[0]
            tilts.loc[k, "theta_motion"] = mot3d_phitheta[1]
            # Calculate tilt angles for sonic anemometer
            wnd3d_k, _ = planarfit_coef(wnd.values[:, 0:3])
            _, wnd3d_phitheta = rotate_vectors(wnd.values[:, 0:3],
                                               k_vector=wnd3d_k)
            tilts.loc[k, "phi_sonic"] = wnd3d_phitheta[0]
            tilts.loc[k, "theta_sonic"] = wnd3d_phitheta[1]
            wind_direction_mean = circmean(ec_win["wind_direction"].values,
                                           high=np.degrees(2 * np.pi))
            tilts.loc[k, "wind_direction"] = wind_direction_mean
        logger.info("End window %s", k)

    tilts.to_csv(outfile, index_label="timestamp")
    return tilts


def plot_tilt(tilt_file, fig_file, title=None):
    """Generate a plot of tilt angles and diagnostics

    Parameters
    ----------
    tilt_file : str
        Path to file to read Dataframe from
    fig_file : str
        Path to file to write plot to
    title : str
        Title for plot
    """
    plt.style.use('ggplot')
    # Workaround for stupid problem in Debian now:
    plt.rcParams['mathtext.fontset'] = 'stix'

    tilts = pd.read_csv(tilt_file, index_col=0, parse_dates=True)
    # # Get a scaled nfiles_ok
    # nfiles_scaled = ((tilts["nfiles_ok"] - tilts["nfiles_ok"].min()) /
    #                  tilts["nfiles_ok"].ptp())
    # Example color bar
    # cbar = axs[0].scatter(tilts.index, tilts["wind_direction"],
    #                       c=tilts["nfiles_ok"], s=30, cmap=plt.cm.coolwarm)
    # clb = fig.colorbar(cbar, ax=axs[0], orientation="horizontal", shrink=0.8,
    #                    aspect=70, pad=0.07, ticks=np.arange(11), format="%d")
    # clb.ax.set_xlabel("N periods in window", size=10)

    # Exclude missing data from plots so lines work
    # Set up wind direction vector
    def plot_quiver(tseries, ax, step=1, yoffset=1, **kwargs):
        """Return an axis with quiver plot for a time series

        Parameters
        ----------
        tseries : pandas.Series
            Time series with angles (radians) to plot.
        ax : matplotlib.axes
            Axis to plot in.
        step : int
            Step (stride) to subsample the time series by.
        yoffset : int
            Location for quiver in the y-coordinate.
        **kwargs : optional keyword arguments
            Arguments passed to pyplot.quiver.

        Returns
        -------
        Matplotlib axis
        """
        pivot = kwargs.pop("pivot", "mid")
        width = kwargs.pop("width", 0.015)
        scale = kwargs.pop("scale", 1 / 0.2)
        units = kwargs.pop("units", "inches")
        ts = tseries.dropna()[::step]
        tstamps = mpl.dates.date2num(ts.index.to_pydatetime())
        ts_cos = np.cos(ts).values.flatten()
        ts_sin = np.sin(ts).values.flatten()
        qvr = ax.quiver(tstamps, np.ones(ts.shape[0]) * yoffset,
                        ts_sin, ts_cos, pivot=pivot, units=units,
                        width=width, scale=scale, **kwargs)
        return qvr

    # Set up phi_motion vector
    mot_phi = tilts[["phi_motion"]]
    # Set up theta_motion vector
    mot_theta = tilts[["theta_motion"]]
    # Set up phi_sonic vector
    wnd3d_phi = tilts[["phi_sonic"]]
    # Set up theta_sonic vector
    wnd3d_theta = tilts[["theta_sonic"]]

    fig = plt.figure(figsize=(15, 12))
    # Main plotting grid
    grd = mpl.gridspec.GridSpec(2, 1)
    grd.update(hspace=0.1)
    # Top plots
    grd1 = mpl.gridspec.GridSpecFromSubplotSpec(5, 1, subplot_spec=grd[0],
                                                hspace=0.015)
    # Top line plot
    ax1 = plt.subplot(grd1[1:, :])
    # Corresponding quiver plot above
    ax0 = plt.subplot(grd1[0, :], sharex=ax1)
    # Bottom plots
    grd2 = mpl.gridspec.GridSpecFromSubplotSpec(5, 1, subplot_spec=grd[1],
                                                hspace=0.015)
    # Bottom line plot
    ax3 = plt.subplot(grd2[1:, :], sharex=ax1)
    # Corresponding quiver plot above
    ax2 = plt.subplot(grd2[0, :], sharex=ax1)
    # Plotting
    # Quiver colors as list of dictionaries to match the rest
    q_colors = [x for x in plt.rcParams["axes.prop_cycle"]]
    q_yticks = np.array([1, 1.1])  # quiver y ticks
    q_ylims = [np.min(q_yticks) - 0.05, np.max(q_yticks) + 0.05]
    # Top
    plot_quiver(mot_phi, ax0, step=3, yoffset=q_yticks[0],
                **q_colors[0])
    plot_quiver(wnd3d_phi, ax0, step=3, yoffset=q_yticks[1],
                **q_colors[1])
    ax0.yaxis.set_visible(False)
    ax0.set_ylim(q_ylims)
    ax0.tick_params(axis="x", bottom="off", labelbottom="off")
    ax1.plot("phi_motion", data=mot_phi.dropna())
    ax1.plot("phi_sonic", data=wnd3d_phi.dropna())
    ax1.set_ylabel(r"Roll $\phi$")
    ax1.xaxis.set_visible(False)
    # Bottom
    plot_quiver(mot_theta, ax2, step=3, yoffset=q_yticks[0],
                **q_colors[0])
    plot_quiver(wnd3d_theta, ax2, step=3, yoffset=q_yticks[1],
                **q_colors[1])
    ax2.yaxis.set_visible(False)
    ax2.set_ylim(q_ylims)
    ax2.tick_params(axis="x", bottom="off", labelbottom="off")
    ax3.plot("theta_motion", data=mot_theta.dropna())
    ax3.plot("theta_sonic", data=wnd3d_theta.dropna())
    ax3.set_ylabel(r"Pitch $\theta$")
    ax3.set_xlabel("")
    # Legend
    leg = ax3.legend(loc=9, bbox_to_anchor=(0.5, -0.05), frameon=False,
                     borderaxespad=0, ncol=2)
    leg.get_texts()[0].set_text("motion sensor")
    leg.get_texts()[1].set_text("sonic anemometer")
    fig.savefig(fig_file, bbox_inches="tight")


if __name__ == "__main__":
    _DESCRIPTION = ("Calculate tilt angles \\phi and \\theta, " +
                    "given a configuration file and temporal " +
                    "window to search.")
    _FORMATERCLASS = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(description=_DESCRIPTION,
                                     formatter_class=_FORMATERCLASS)
    parser.add_argument("config_file", type=str,
                        help="Path to configuration file")
    parser.add_argument("window_width", nargs=2, default=[2, 0],
                        help=("Integer hours and minutes for the moving " +
                              "temporal window width."))
    parser.add_argument("outfile", type=argparse.FileType("w"),
                        default=sys.stdout,
                        help=("Output file to write results to."))
    parser.add_argument("--log-file", type=str, default="tilt_search.log",
                        help="Path to log file")
    args = parser.parse_args()
    main(args.config_file, win_hours=args.window_width[0],
         win_minutes=args.window_width[1], outfile=args.outfile)
