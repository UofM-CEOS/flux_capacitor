# pylint: disable=too-many-locals,invalid-name,no-member

"""Calculate tilt angles in flux files using non-overlapping moving windows

"""

import itertools as itrt
import datetime
import os.path as osp
import logging
import numpy as np
import pandas as pd
from scipy.stats import circmean
import matplotlib as mpl
import matplotlib.pyplot as plt
from fluxer.eddycov.settings import FLUX_FLAGS
from fluxer.eddycov.parse_ecfile import (prepare_period, NavigationError,
                                         SonicError)
from fluxer.eddycov.flux import planarfit

__all__ = ["TiltWindows"]

logger = logging.getLogger(__name__)
# Add the null handler if importing as library; whatever using this library
# should set up logging.basicConfig() as needed
logger.addHandler(logging.NullHandler())


def _fname2time(fname):
    """Construct a time stamp from a file name"""
    return datetime.datetime.strptime(fname[-18:-4], "%Y%m%d%H%M%S")


def _make_key(delta_t):
    """Return a function with attribute to track current upper bound"""
    def key(fname):
        """Set upper bound for current time window given a file name"""
        if key.upper_bound is None:
            key.upper_bound = _fname2time(fname) + delta_t
        else:
            tstamp = _fname2time(fname)
            while tstamp >= key.upper_bound:
                key.upper_bound += delta_t
        return key.upper_bound
    key.upper_bound = None
    return key


def _make_windows(ec_files, win_minutes):
    """Build list of tuples linking list of files to a time window"""
    # Time delta
    tdelta = datetime.timedelta(minutes=win_minutes)
    # Build a list of tuples (key=window start timestamp, file-list)
    win_files = [(pd.to_datetime(k), list(grp)) for k, grp in
                 itrt.groupby(ec_files, key=_make_key(tdelta))]
    return win_files


class TiltWindows:
    """Well-defined pandas DataFrame"""
    def __init__(self, ec_files=[], win_minutes=0):
        """Set up attributes for handling tilt window data

        Parameters
        ----------
        ec_files : list, optional
            List of file names.  Default is empty list.
        win_minutes : int, optional
            Number of minutes to build windows, based on file names.
            Default is zero.

        """
        win_files = _make_windows(ec_files, win_minutes)
        # Set up dataframe to place tilt angles and other details
        cols = (["failed_tilt_flag", "nfiles", "nfiles_ok"] +
                ["n{0}".format(i) for i in FLUX_FLAGS] +
                ["phi_motion", "theta_motion", "phi_sonic", "theta_sonic"] +
                ["wind_direction"])
        tilts = pd.DataFrame(index=[x[0] for x in win_files], columns=cols)
        tilts.loc[:, cols[1]] = [len(x[1]) for x in win_files]
        tilts.loc[:, cols[-5:]] = tilts[cols[-5:]].astype(np.float).values
        # Zero all the columns except nfiles and parameters
        tilts.loc[:, cols[0]] = False
        tilts.loc[:, cols[2:-5]] = 0
        self.width = win_minutes
        self.tilts = tilts
        self.win_files = win_files
        # [Original comment: create flags for the 6 possible sources of
        # "bad" data, flag=0 means data good] in each input file
        flags = dict.fromkeys(FLUX_FLAGS, False)
        # We set up a dataframe with all files to process as index, and all
        # the flags as columns.  This is the basis for our summary output
        # file; other columns (e.g. flux summary calculations for the
        # period) will be appended as we loop.
        prep_flags = pd.DataFrame(flags,
                                  index=[osp.basename(x) for x in ec_files])
        self.prep_flags = prep_flags

    def __str__(self):
        # Show tilts attribute first, followed by the number of windows
        # found, followed by the width (window width) attribute
        msg = ("TiltWindows({0.tilts}, nwindows={1}, "
               "width={0.width!r})")
        return msg.format(self, len(self.win_files))

    __repr__ = __str__

    def get_tilts_planarfit(self, config):
        """Compute tilt angles in all windows using planar fit method

        Parameters
        ----------
        config : OrderedDict
            Dictionary with parsed configuration file.

        Returns
        -------

        pandas.DataFrame with tilt angles for each window.

        """
        wind_names = ["wind_speed_u", "wind_speed_v", "wind_speed_w"]
        mot3d_names = ["acceleration_x", "acceleration_y",
                       "acceleration_z", "rate_x", "rate_y", "rate_z"]

        # Iterate over the list of tuples (key=window time stamp, file-list)
        msg_suffix = " tilt angles via planar fit for window %s"
        for (k, l) in self.win_files:
            logger.info("Begin" + msg_suffix, k)
            ec_list = []
            ec_ok_files = []
            # Iterate over the file-list
            for ec_file in l:
                logger.info("Begin preparing %s", osp.basename(ec_file))
                # Get a file name prefix to be shared by the output files
                # from this period.  Note iname is THE SAME AS THE INDEX IN
                # prep_flags
                iname = osp.basename(ec_file)
                # Try to prepare file; set flags, and populate the list of
                # OK files and summaries. Continue if preparation fails.
                try:
                    ec_prep, prep_flags = prepare_period(ec_file, config)
                except (NavigationError, SonicError) as err:
                    self.prep_flags.loc[iname, "failed_prep_flag"] = True
                    for flagname, flag in err.flags.items():
                        self.prep_flags.loc[iname, flagname] = flag
                        self.tilts.loc[k, "n" + flagname] += np.int(flag)
                    logger.info("Skip %s", ec_file)
                    continue
                else:
                    for flagname, flag in prep_flags.iteritems():
                        self.prep_flags.loc[iname, flagname] = flag
                        self.tilts.loc[k, "n" + flagname] += np.int(flag)
                ec_list.append(pd.DataFrame(ec_prep.mean()).transpose())
                ec_ok_files.append(ec_file)
                logger.info("End preparing %s", osp.basename(ec_file))

            l = ec_ok_files
            n_ok = len(ec_list)
            self.tilts.loc[k, "nfiles_ok"] = n_ok
            if n_ok < 3:
                logger.info("Aborting window without enough adequate files %s",
                            k)
                for ec_file in l:
                    self.tilts.loc[k, "failed_tilt_flag"] = True
                continue
            else:
                ec_win = pd.concat(ec_list, keys=ec_ok_files)
                logger.info("Combining %s periods in window %s", n_ok, k)
                # Extract means from each period file for planar fitting
                wnd = ec_win.loc[:, wind_names]
                mot3d = ec_win.loc[:, mot3d_names[0:3]]
                # Tilt angles for motion sensor
                mot3d_pfit = planarfit(mot3d.values)
                self.tilts.loc[k, "phi_motion"] = mot3d_pfit.phi
                self.tilts.loc[k, "theta_motion"] = mot3d_pfit.theta
                # Tilt angles for sonic anemometer
                sonic_pfit = planarfit(wnd.values)
                self.tilts.loc[k, "phi_sonic"] = sonic_pfit.phi
                self.tilts.loc[k, "theta_sonic"] = sonic_pfit.theta
                wind_direction_mean = circmean(ec_win["wind_direction"].values,
                                               high=np.degrees(2 * np.pi))
                self.tilts.loc[k, "wind_direction"] = wind_direction_mean
            logger.info("End" + msg_suffix, k)

    def plot(self, fig_file, step=1, title=None):
        """Generate a plot of tilt angles and diagnostics

        Parameters
        ----------
        fig_file : str
            Path to file to write plot to
        step : int, optional
            Step (stride) to subsample the time series by.  Default is 1.
        title : str
            Title for plot.  Default is None.

        """
        # Workaround for stupid problem in Debian now:
        plt.rcParams['mathtext.fontset'] = 'stix'

        tilts = self.tilts

        # Exclude missing data from plots so lines work
        # Set up wind direction vector
        def plot_quiver(tseries, ax, step, yoffset=1, **kwargs):
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
        plot_quiver(mot_phi, ax0, step=step, yoffset=q_yticks[0],
                    **q_colors[0])
        plot_quiver(wnd3d_phi, ax0, step=step, yoffset=q_yticks[1],
                    **q_colors[1])
        ax0.yaxis.set_visible(False)
        ax0.set_ylim(q_ylims)
        ax0.tick_params(axis="x", bottom="off", labelbottom="off")
        ax1.plot("phi_motion", data=mot_phi.dropna())
        ax1.plot("phi_sonic", data=wnd3d_phi.dropna())
        ax1.set_ylabel(r"Roll $\phi$")
        ax1.xaxis.set_visible(False)
        # Bottom
        plot_quiver(mot_theta, ax2, step=step, yoffset=q_yticks[0],
                    **q_colors[0])
        plot_quiver(wnd3d_theta, ax2, step=step, yoffset=q_yticks[1],
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
    # Test defaults
    tlt = TiltWindows()
    # Generate 6 hours (20 min x 30) of time stamps for filenames
    tstamps = pd.date_range("2016-01-01 00:30:00", periods=30,
                            freq="20min")
    files = tstamps.strftime("EC_%Y%m%d%H%M%S.csv")
    tlt = TiltWindows(files, 120)
    tlt
