# pylint: disable=too-many-locals,invalid-name,no-member

"""Calculate tilt angles in flux files using non-overlapping moving windows

"""

import itertools as itrt
import datetime
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from .settings import FLUX_FLAGS

__all__ = ["TiltWindows"]


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
        ec_files : list
            List of file names.
        win_minutes : int
            Number of minutes to build windows, based on file names.

        """
        win_files = _make_windows(ec_files, win_minutes)
        # Set up dataframe to place tilt angles and other details
        cols = (["nfiles", "nfiles_ok"] +
                ["n{0}".format(i) for i in FLUX_FLAGS] +
                ["phi_motion", "theta_motion", "phi_sonic", "theta_sonic"] +
                ["wind_direction"])
        tilts = pd.DataFrame(index=[x[0] for x in win_files], columns=cols)
        tilts.loc[:, cols[0]] = [len(x[1]) for x in win_files]
        tilts.loc[:, cols[-5:]] = tilts[cols[-5:]].astype(np.float).values
        # Zero all the columns except the nfiles and parameters
        tilts.loc[:, cols[1:-5]] = 0
        self.width = win_minutes
        self.tilts = tilts
        self.win_files = win_files

    def __str__(self):
        # Show tilts attribute first, followed by the number of windows
        # found, followed by the width (window width) attribute
        msg = ("TiltWindows({0.tilts}, nwindows={1}, "
               "width={0.width!r})")
        return msg.format(self, len(self.win_files))

    def plot(self, fig_file, title=None):
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

        tilts = self.tilts
        # # Get a scaled nfiles_ok
        # nfiles_scaled = ((tilts["nfiles_ok"] - tilts["nfiles_ok"].min()) /
        #                  tilts["nfiles_ok"].ptp())
        # Example color bar
        # cbar = axs[0].scatter(tilts.index, tilts["wind_direction"],
        #                       c=tilts["nfiles_ok"], s=30,
        #                       cmap=plt.cm.coolwarm)
        # clb = fig.colorbar(cbar, ax=axs[0], orientation="horizontal",
        #                    shrink=0.8, aspect=70, pad=0.07,
        #                    ticks=np.arange(11), format="%d")
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
    # Test defaults
    tlt = TiltWindows()
    # Generate 6 hours (20 min x 30) of time stamps for filenames
    tstamps = pd.date_range("2016-01-01 00:30:00", periods=30,
                            freq="20min")
    files = tstamps.strftime("EC_%Y%m%d%H%M%S.csv")
    tlt = TiltWindows(files, 120)
    print tlt
