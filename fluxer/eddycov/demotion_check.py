#! /usr/bin/env python

"""Compare measured and motion-corrected spectra/cospectra of 3D wind"""

import argparse
import numpy as np
from scipy import signal
from scipy.stats import binned_statistic
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import fluxer.eddycov as ec
from fluxer.eddycov.flux import _cumtrapz as cumtrapz


_WIND3D_AXES_NAMES = ["U", "V", "W"]
_IMU_AXES_NAMES = ["x", "y", "z"]
# Math bracket assumed
_EULER_ANGLE_NAMES = [r"{{\phi}}", r"{{\theta}}", r"{{\psi}}"]
# Constants for legend placement
_LEG3X1_ANCHOR = (0.5, -0.23)
_LEG2X1_ANCHOR = (0.5, -0.2)


def _bin_spectra(f, pxx, nbins):
    """Compute exponentially binned (co)spectra

    This procedure follows EddyPro's option to reduce (co)spectra in order
    to reduce noise.  The exponential nature of the size increase of bins
    assures that increasingly more spectral values are averaged as it moves
    toward higher frequencies, thereby increasingly reducing noise.

    Parameters
    ----------
    f : array_like
        Frequencies of the spectral density `pxx`.
    pxx : array_like
        Spectral density estimates at each element in `f`.
    nbins : int
        Number of bins to generate compute averages at.

    Returns
    -------
    bin_mid : ndarray
        Array of length ``nbins`` with the average frequency in each bin.
    bin_means : ndarray
        Array of length ``nbins`` with the average spectral density at each
        bin.

    """
    f_bins = np.logspace(np.floor(np.log10(f[1])), np.log10(f[-1]), nbins)
    bin_means, bin_edges, binnum = binned_statistic(f, pxx, bins=f_bins)
    bin_mid = np.mean(np.column_stack((bin_edges[:-1], bin_edges[1:])), -1)
    # Remove empty bins
    bin_mid = bin_mid[~np.isnan(bin_means)]
    bin_means = bin_means[~np.isnan(bin_means)]
    return bin_mid, bin_means


def psd_fun(x, nbins=100, **kwargs):
    """Estimate spectral density using Welch's method

    This uses the Welch method to compute the PSD of `x`, using a Tukey
    window (1% tapered) of 1 min length, and an overlap of 15 s, by
    default.  By default also, the signal is detrended linearly before
    computing the DFT.

    In addition to the full resolution PSD, given the window, a reduced
    spectra is also provided with exponentially binned averages.

    Parameters
    ----------
    x : ndarray
        time series of measurement values.
    nbins : int, optional
        Number of bins to generate compute averages at.
    fs : float, optional
        Sampling frequency of the `x` time series.  Defaults to 10.0.
    window : tuple, optional
        Window to apply to the time series.  Defaults to ("tukey", 0.1).
    nperseg : int, optional
        Length of each segment for building averages.
        Defaults to 60 * 5 * `fs`.
    noverlap : int, optional
        Number of overlapping samples between segments.
        Defaults to `nperseg` // 4.
    detrend : str, optional
        Type of detrending to apply to time series.  Defaults to "linear".
    **kwargs
        Keyword arguments passed to signal.welch.

    Returns
    -------
    f : ndarray
        Array of sample frequencies.
    Pxx : ndarray
        Power spectral density of `x`.
    bin_mid : ndarray
        Array of length ``nbins`` with the average frequency in each bin.
    bin_means : ndarray
        Array of length ``nbins`` with the average spectral density at each
        bin.

    See Also
    --------
    scipy.signal.welch

    """
    if "window" not in kwargs:
        kwargs.update({"window": ("tukey", 0.1)})
    if "fs" not in kwargs:
        kwargs.update({"fs": 10.0})
    if "nperseg" not in kwargs:
        kwargs.update({"nperseg": 60 * 5 * kwargs.get("fs")})
    if "noverlap" not in kwargs:
        kwargs.update({"noverlap": kwargs.get("nperseg") // 4})
    if "detrend" not in kwargs:
        kwargs.update({"detrend": "linear"})
    kwargs.update({"scaling": "density"})
    f, Pxx = signal.welch(x, **kwargs)
    bin_mid, bin_means = _bin_spectra(f, Pxx, nbins=nbins)
    return f, Pxx, bin_mid, bin_means


# New parameter values for cospectra
def csd_fun(x, y, nbins=100, **kwargs):
    """Estimate cospectral density, Pxy, using Welch's method

    This uses the Welch method to compute the cross power spectral density
    of `x` and `y`, using a Tukey window (1% tapered) of 1 min length, and
    an overlap of 15 s.  The signal is detrended linearly before computing
    the DFT.

    In addition to the full resolution cospectra, given the window, a
    reduced cospectra is also provided with exponentially binned averages.

    Parameters
    ----------
    x : ndarray
        time series of measurement values.
    y : ndarray
        time series of measurement values.
    nbins : int, optional
        Number of bins to generate compute averages at.
    fs : float, optional
        Sampling frequency of the `x` time series.  Defaults to 10.0.
    window : tuple, optional
        Window to apply to the time series.  Defaults to ("tukey", 0.1).
    nperseg : int, optional
        Length of each segment for building averages.
        Defaults to 60 * 5 * `fs`.
    noverlap : int, optional
        Number of overlapping samples between segments.
        Defaults to `nperseg` // 4.
    detrend : str, optional
        Type of detrending to apply to time series.  Defaults to "linear".
    **kwargs
        Keyword arguments passed to signal.welch.

    Returns
    -------
    f : ndarray
        Array of sample frequencies.
    Pxy : ndarray
        Cross spectral density or cross power spectrum of `x`, `y`.
    bin_mid : ndarray
        Array of length ``nbins`` with the average frequency in each bin.
    bin_means : ndarray
        Array of length ``nbins`` with the average spectral density at each
        bin.

    See Also
    --------
    scipy.signal.welch

    """
    if "win" not in kwargs:
        kwargs.update({"win": ("tukey", 0.1)})
    if "fs" not in kwargs:
        kwargs.update({"fs": 10.0})
    if "nperseg" not in kwargs:
        kwargs.update({"nperseg": 60 * 5 * kwargs.get("fs")})
    if "noverlap" not in kwargs:
        kwargs.update({"noverlap": kwargs.get("nperseg") // 4})
    if "detrend" not in kwargs:
        kwargs.update({"detrend": "linear"})
    kwargs.update({"scaling": "density"})
    f, Pxy = signal.csd(x, y, scaling="density")
    Pxy = np.abs(Pxy)
    bin_mid, bin_means = _bin_spectra(f, Pxy, nbins=nbins)
    return f, Pxy, bin_mid, bin_means


def plot_wind3D_spectra(uvw_measured, uvw_corrected, **kwargs):
    """Compare wind spectra of measured and corrected 3D wind

    Parameters
    ----------
    uvw_measured : Nx3 array
        Array of measured U, V, and W wind velocity components.
    uvw_corrected : Nx3 array
        Array of motion-corrected U, V, and W wind velocity components.
    **kwargs
        Keyword arguments passed to psd_fun.

    Returns
    -------
    tuple
        Figure, Axes, and Legend instances consisting of 3x1 matrix
        of plots of measured vs motion-corrected spectra.

    """
    # Constants for spectra figures
    fig3x1 = (4, 10)                # fig size
    # Plot measured and motion-corrected spectra
    fig, axs = plt.subplots(3, 1, sharex=True, sharey=True)
    fig.set_size_inches(fig3x1)
    axs[0].set_title("Power spectra")
    axs[2].set_xlabel("Frequency (Hz)")

    # Iterate through wind components
    for idx, col in enumerate(_WIND3D_AXES_NAMES):
        msrd = uvw_measured[:, idx]
        corr = uvw_corrected[:, idx]
        f_msrd, pxx_msrd, f_bin_msrd, pxx_bin_msrd = psd_fun(msrd, **kwargs)
        f_corr, pxx_corr, f_bin_corr, pxx_bin_corr = psd_fun(corr, **kwargs)
        axs[idx].loglog(f_msrd, f_msrd * pxx_msrd, "C0.",
                        markersize=0.2, label="measured")
        axs[idx].loglog(f_bin_msrd, f_bin_msrd * pxx_bin_msrd, "C0",
                        label="measured (binned)")
        axs[idx].loglog(f_corr, f_corr * pxx_corr, "C1.",
                        markersize=0.2, label="corrected")
        axs[idx].loglog(f_bin_corr, f_bin_corr * pxx_bin_corr, "C1",
                        label="corrected (binned)")
        ylabel = r"$f \cdot S_{0}(f)\quad [m^2 \cdot s^{{-2}}]$"
        axs[idx].set_ylabel(ylabel.format(col.lower()))

    leg = axs[2].legend(loc=9, bbox_to_anchor=_LEG3X1_ANCHOR,
                        frameon=False, borderaxespad=0, ncol=2)
    return (fig, axs, leg)


def plot_wind3D_cospectra(uvw_measured, uvw_corrected, **kwargs):
    """Compare wind cospectra of measured and corrected 3D wind

    Parameters
    ----------
    uvw_measured : Nx3 array
        Array of measured U, V, and W wind velocity components.
    uvw_corrected : Nx3 array
        Array of motion-corrected U, V, and W wind velocity components.
    **kwargs
        Keyword arguments passed to csd_fun.

    Returns
    -------
    tuple
        Figure, Axes, and Legend instances consisting of 2x1 matrix
        of plots of measured vs motion-corrected cospectra.

    """
    # Constants for cospectra figure
    fig2x1 = (4, 7)                # fig size
    # Plot measured and motion-corrected u-w and v-w cospectra
    fig, axs = plt.subplots(2, 1, sharex=True, sharey=True)
    fig.set_size_inches(fig2x1)
    axs[0].set_title("Power cospectra")

    # Iterate through wind components
    for idx, col in enumerate(_WIND3D_AXES_NAMES[0:2]):
        msrd = uvw_measured
        corr = uvw_corrected
        f_msrd, pxx_msrd, f_bin_msrd, pxx_bin_msrd = csd_fun(msrd[:, idx],
                                                             msrd[:, 2],
                                                             **kwargs)
        f_corr, pxx_corr, f_bin_corr, pxx_bin_corr = csd_fun(corr[:, idx],
                                                             corr[:, 2],
                                                             **kwargs)
        axs[idx].semilogx(f_msrd, f_msrd * pxx_msrd, "C0.",
                          markersize=0.2, label="measured")
        axs[idx].semilogx(f_bin_msrd, f_bin_msrd * pxx_bin_msrd,
                          "C0", label="measured (binned)")
        axs[idx].semilogx(f_corr, f_corr * pxx_corr, "C1.",
                          markersize=0.2, label="corrected")
        axs[idx].semilogx(f_bin_corr, f_bin_corr * pxx_bin_corr,
                          "C1", label="corrected (binned)")
        ylabel = r"$f \cdot C_{{{0}w}}(f)$"
        axs[idx].set_ylabel(ylabel.format(col.lower()))

    axs[1].set_xlabel("Frequency (Hz)")
    leg = axs[1].legend(loc=9, bbox_to_anchor=_LEG2X1_ANCHOR,
                        frameon=False, borderaxespad=0, ncol=2)
    return (fig, axs, leg)


def _spectra_subplot_init(ncols, nrows, title, **kwargs):
    """Initialize figure and axes with requested subplots"""
    if ncols == 3 and nrows == 1:
        fig_size = (6, 10)
    elif ncols == 2 and nrows == 1:
        fig_size = (5, 8)
    fig, axs = plt.subplots(ncols, nrows, **kwargs)
    fig.set_size_inches(fig_size)
    axs[0].set_title(title)
    axs[ncols - 1].set_xlabel("Frequency (Hz)")
    return (fig, axs)


def _plot1_IMU_complementary_spectra(motion_corrected_obj, axs,
                                     sample_freq_hz, nperseg, **kwargs):
    """Plot spectral densities of Euler angles in motion corrected object"""
    euler_angles = motion_corrected_obj.euler_angles
    for i, angle in enumerate(_EULER_ANGLE_NAMES[:-1]):
        f, pxx, _, _ = psd_fun(euler_angles[:, i], fs=sample_freq_hz,
                               nperseg=nperseg)
        axs[i].loglog(f, f * pxx,
                      label="Complementary filter", **kwargs)
        ylabel = r"$f \cdot S_{0}(f)\quad [rad^2]$"
        axs[i].set_ylabel(ylabel.format(angle.lower()))


def _plot1_IMU_euler_angles_spectra(rates, accelerations, axs,
                                    sample_freq_hz, nperseg, **kwargs):
    """Plot spectra of IMU's angular rates and linear accelerations"""
    # Integrate angular rates (straight from wind3D_correct)
    rate0 = signal.detrend(rates, 0, "constant")
    EA_rate = cumtrapz(rate0, sample_freq_hz)
    # Normalize x and y acceleration for gravity
    g = np.sqrt(np.sum(np.mean(accelerations, 0) ** 2))
    EA_acc = np.column_stack((np.arcsin(accelerations[:, 1] / g),   # phi
                              np.arcsin(-accelerations[:, 0] / g)))  # theta
    for i, angle in enumerate(_EULER_ANGLE_NAMES[:-1]):
        # angular rates
        f, pxx, _, _ = psd_fun(EA_rate[:, i], fs=sample_freq_hz,
                               nperseg=nperseg)
        axs[i].loglog(f, f * pxx, "C0",
                      label="integrated angular rate", **kwargs)
        # acceleration
        f, pxx, _, _ = psd_fun(EA_acc[:, i], fs=sample_freq_hz,
                               nperseg=nperseg)
        axs[i].loglog(f, f * pxx, "C1",
                      label="normalized acceleration", **kwargs)
        ylabel = r"$f \cdot S_{0}(f)\quad [rad^2]$"
        axs[i].set_ylabel(ylabel.format(angle.lower()))


def _plot1_IMU_spectra(rates, accelerations, axs_rate, axs_acc,
                       sample_freq_hz, nperseg, **kwargs):
    """Plot IMU's angular rate and acceleration spectra independently"""
    for i, col in enumerate(_IMU_AXES_NAMES):
        f, pxx, _, _ = psd_fun(rates[:, i], fs=sample_freq_hz,
                               nperseg=nperseg)
        axs_rate[i].semilogy(f, pxx, "C0", **kwargs)
        rate_label = r"$S_{0}(f)\quad [rad^2 \cdot s^{{-2}}] / Hz$"
        axs_rate[i].set_ylabel(rate_label.format(col.lower()))
        # Now for acceleration, use default nbins
        f, pxx, _, _ = psd_fun(accelerations[:, i], fs=sample_freq_hz,
                               nperseg=nperseg)
        axs_acc[i].semilogy(f, pxx, "C0", **kwargs)
        acc_label = r"$S_{0}(f)\quad [m \cdot s^{{-2}}]^2 / Hz$"
        axs_acc[i].set_ylabel(acc_label.format(col.lower()))


def plot_IMU_euler_angles_spectra(config_file, nperseg=60 * 5 * 10,
                                  nfiles=None, **kwargs):
    """Compute and plot spectra for motion sensor (IMU) fusion

    The plot shows the spectral density for all individual runs from input
    files specified in the configuration file.

    Parameters
    ----------
    config_file : string
        Path to file with configuration settings for processing Eddy
        covariance files.
    nperseg : int, optional
        Number of samples per segment for computing the spectral density
        using Welch's method (see scipy.signal.welch).
    nfiles : int, optional
        Number of files to randomly sample from the configuration file list.
    **kwargs
        Keyword arguments passed to Axes.loglog.

    Returns
    -------
    tuple
        Figure, Axes, and Legend instances consisting of 2x1 matrix
        of plots of spectral densities for non-filtered :math:`\\phi` and
        :math:`\\theta` Euler angles estimated directly from integration of
        angular velocities and normalized linear accelerations from IMU.

    """
    if "linewidth" not in kwargs:
        kwargs.update({"linewidth": 0.1})
    config = ec.db_flux.parse_config(config_file)
    input_files = config["EC Inputs"]["input_files"]
    if nfiles is not None and nfiles < len(input_files):
        input_files = np.random.choice(input_files, nfiles, replace=False)
    sample_freq_hz = config["EC Inputs"]["sample_frequency"]
    # Set up subplots
    title = "IMU naive Euler roll and pitch spectra"
    fig, axs = _spectra_subplot_init(2, 1, title,
                                     sharex=True, sharey=True)
    for fname in input_files:
        ec_prep, dbflags = ec.db_flux.prepare_period(fname, config)
        acceleration = ec_prep[["acceleration_x",
                                "acceleration_y",
                                "acceleration_z"]].values
        rate = ec_prep[["rate_x", "rate_y", "rate_z"]].values
        _plot1_IMU_euler_angles_spectra(rate, acceleration, axs,
                                        sample_freq_hz, nperseg, **kwargs)
    _, labels = axs[1].get_legend_handles_labels()
    leg = axs[1].legend(labels[:2], loc=9, bbox_to_anchor=_LEG2X1_ANCHOR,
                        frameon=False, borderaxespad=0, ncol=2)
    return (fig, axs, leg)


def plot_IMU_spectra(config_file, nperseg=60 * 10, nfiles=None, **kwargs):
    """Compute and plot spectra for all signals in motion sensor (IMU)

    The plot shows the spectral density for each individual run from input
    files specified in the configuration file.

    Parameters
    ----------
    config_file : string
        Path to file with configuration settings for processing Eddy
        covariance files.
    nperseg : int, optional
        Number of samples per segment for computing the spectral density
        using Welch's method (see scipy.signal.welch).
    nfiles : int, optional
        Number of files to randomly sample from the configuration file list.
    **kwargs
        Optional keyword arguments passed to Axes.loglog.

    Returns
    -------
    tuple
        Two sets of Figure and Axes instances for angular velocities and
        linear acceleration, each consisting of 3x1 matrix of semi-log (y)
        plots of spectral densities for the three IMU angular velocity
        sensors (first set, with order: Figure and Axes instances), and for
        the three IMU linear acceleration sensors (second set, also with
        order Figure, and Axes instances).

    """
    if "linewidth" not in kwargs:
        kwargs.update({"linewidth": 0.1})
    config = ec.db_flux.parse_config(config_file)
    input_files = config["EC Inputs"]["input_files"]
    if nfiles is not None and nfiles < len(input_files):
        input_files = np.random.choice(input_files, nfiles, replace=False)
    sample_freq_hz = config["EC Inputs"]["sample_frequency"]
    # Prepare plots for angular rate
    fig_rate, axs_rate = _spectra_subplot_init(3, 1,
                                               "Angular rate spectra",
                                               sharex=True,
                                               sharey=True)
    # For linear acceleration
    fig_acc, axs_acc = _spectra_subplot_init(3, 1,
                                             "Acceleration spectra",
                                             sharex=True, sharey=True)
    for fname in input_files:
        ec_prep, dbflags = ec.db_flux.prepare_period(fname, config)
        acceleration = ec_prep[["acceleration_x",
                                "acceleration_y",
                                "acceleration_z"]].values
        rate = ec_prep[["rate_x", "rate_y", "rate_z"]].values
        _plot1_IMU_spectra(rate, acceleration, axs_rate, axs_acc,
                           sample_freq_hz, nperseg, **kwargs)
    return (fig_rate, axs_rate, fig_acc, axs_acc)


def plot1_IMU_complementary_spectra(config_file, file_idx, Tcf, Ta,
                                    nperseg=60 * 5 * 10, **kwargs):
    """Plot spectral densities used in complementary-filtered IMU data

    The plot shows the spectral densities for an individual run from input
    files specified in the configuration file.

    Parameters
    ----------
    config_file : string
        Path to file with configuration settings for processing Eddy
        covariance files.
    file_idx : int
        Index of the run file to plot in the data set specified in
        config_file.
    Tcf : float
        Complementary filter cutoff period
        (see fluxer.eddycov.flux.wind3D_correct)
    Ta : float
        Filter cutoff period (see fluxer.eddycov.flux.wind3D_correct)
    nperseg : int, optional
        Number of samples per segment for computing the spectral density
        using Welch's method (see scipy.signal.welch).
    **kwargs
        Optional keyword arguments passed to Axes.loglog.

    Returns
    -------
    tuple
        Figure, Axes, and Legend instances consisting of 2x1 matrix
        of plots of spectral densities for non-filtered :math:`\\phi` and
        :math:`\\theta` Euler angles estimated directly from integration of
        angular velocities, normalized linear accelerations from IMU, and
        complementary-filtered signals.

    """
    config = ec.db_flux.parse_config(config_file)
    input_files = config["EC Inputs"]["input_files"]
    imu2anem_pos = config["EC Motion Correction"]["imu2anemometer_pos"]
    sample_freq_hz = config["EC Inputs"]["sample_frequency"]
    ifile = input_files[file_idx]
    ec_prep, _ = ec.db_flux.prepare_period(ifile, config)
    wind = ec_prep[["wind_speed_u", "wind_speed_v", "wind_speed_w"]]
    acceleration = ec_prep[["acceleration_x",
                            "acceleration_y",
                            "acceleration_z"]]
    rate = ec_prep[["rate_x", "rate_y", "rate_z"]]
    ec_corr = ec.wind3D_correct(wind.values, acceleration.values,
                                rate.values, ec_prep.heading.values,
                                ec_prep.speed_over_ground.values,
                                imu2anem_pos, sample_freq_hz, Tcf, Ta)
    fig, axs = _spectra_subplot_init(2, 1,
                                     "IMU Euler angles spectra",
                                     sharex=True, sharey=True)
    _plot1_IMU_euler_angles_spectra(rate.values, acceleration.values,
                                    axs, sample_freq_hz, nperseg,
                                    **kwargs)
    _plot1_IMU_complementary_spectra(ec_corr, axs, sample_freq_hz,
                                     nperseg, color="C2", **kwargs)
    ylim = axs[1].get_ylim()
    for ax in axs:
        ax.vlines(1.0 / Tcf, ylim[0], ylim[1],
                  linestyles="dashed", linewidth=0.5)
        ax.text(1.0 / Tcf, ylim[0], r"$\tau_c = {0}s$".format(Tcf))
    leg = axs[1].legend(loc=9, bbox_to_anchor=_LEG2X1_ANCHOR,
                        frameon=False, borderaxespad=0, ncol=2)
    return (fig, axs, leg)


if __name__ == '__main__':
    from oct2py import octave
    from scipy import io as sio
    _DESCRIPTION = ("Plot (co)spectra of observed and motion corrected "
                    "data in Matlab (*.mat) file, using Matlab function "
                    "developed in Miller et al (2008) and ported to Octave. "
                    "The Octave function motion_correct_oursorig is assumed "
                    "to be in Octave's search path.")
    parser = argparse.ArgumentParser(description=_DESCRIPTION)
    parser.add_argument("mat_file", metavar="mat-file",
                        help=("Path to *mat file to pass to "
                              "motion_correct_oursorig"))
    parser.add_argument("--ofigure-file",
                        default="motion_correct_compare.pdf",
                        type=argparse.FileType("w"),
                        help="Path to output figure file.")
    args = parser.parse_args()
    ours, orig, _ = octave.motion_correct_oursorig(args.mat_file,
                                                   nout=3)
    miller = sio.loadmat(args.mat_file,
                         variable_names=["wind_speed"])
    uvw = miller["wind_speed"]

    with PdfPages(args.ofigure_file) as pdf:
        fig, axs, leg = plot_wind3D_spectra(uvw, ours["uvwc"],
                                            nperseg=60 * 10)
        pdf.savefig(fig, bbox_extra_artists=(leg,), bbox_inches="tight")
        plt.close(fig)
        fig, axs, leg = plot_wind3D_cospectra(uvw, ours["uvwc"],
                                              nperseg=60 * 10)
        pdf.savefig(fig, bbox_extra_artists=(leg,), bbox_inches="tight")
        plt.close(fig)


# -------------------------------------------------------------------------
# Other approaches (assuming data above)

# from scipy import signal

# # Try periodogram
# import pandas as pd
# f, pxx_den = signal.periodogram(signal.hamming(uvw.shape[0]) * uvw[:, 2],
#                                 10, window=("tukey", 0.1),
#                                 scaling="density")
# spectra = pd.DataFrame({"psd": pxx_den}, index=f)
# bins2 = np.logspace(10e-4, np.log2(np.max(f)), num=100, base=2)
# spectra["freq_bins"] = pd.cut(f, bins2, include_lowest=True,
#                               right=False).dropna()
# plt.loglog(spectra.index, spectra.index * spectra.psd, '.')
# spectra_binned = spectra.groupby("freq_bins").mean()
# spectra_binned.plot()
# plt.loglog(spectra_binned.freq,
#            spectra_binned.freq * spectra_binned.psd)

# # Try Welch
# f, pxx_den = signal.welch(uvw[:, 2], 10, nperseg=1800, scaling="density")
# plt.loglog(f, f * pxx_den)
# for i in range(256, 256 * 2, 128):
#     f, pxx_den = signal.welch(uvw[:, 2], 10,
#                               nperseg=i, noverlap=0,
#                               # noverlap=i // 3,
#                               # nfft=nfft,
#                               scaling="density")
#     plt.loglog(f, f * pxx_den)

# # Try multi-taper technique
# import nitime.algorithms as tsa
# f_msrd, psd_msrd, nu = tsa.multi_taper_psd(uvw[:, 2], 10, NW=15,
#                                            adaptive=True, jackknife=True)
# plt.loglog(f_msrd, f_msrd * psd_msrd)
# f_ours, psd_ours, nu = tsa.multi_taper_psd(ours["uvwc"][:, 2], 10, NW=15,
#                                            adaptive=False, jackknife=True)
# plt.loglog(f_ours, f_ours * psd_ours)
# plt.ylim((10e-5, 10e-1))
