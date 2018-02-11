#! /usr/bin/env python
# pylint: disable=too-many-locals,invalid-name,no-member

"""Generate spectral density plots of key motion correction variables

Usage
-----

For help using this script, type:

check_spectra.py -h

at command line.

"""

import argparse
import numpy as np
import fluxer.eddycov.demotion_check as check
import fluxer.eddycov.db_flux as db_flux
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def main(config_file, pre_nfiles=None, post_nfiles=None,
         cutoffs=[12.0, 60.0]):
    """Produce (co)spectral density plots"""
    imu_complement_fig = (config_file[:-4] +
                          "_IMU_complementary_spectra.pdf")
    imu_spectra_fig = (config_file[:-4] + "_IMU_spectra.pdf")
    nperseg1 = 60 * 5 * 10
    nperseg2 = 60 * 10
    fig, axs, leg = check.plot_IMU_euler_angles_spectra(config_file,
                                                        nperseg=nperseg2,
                                                        nfiles=pre_nfiles)
    fig.savefig(imu_complement_fig, bbox_extra_artists=(leg,),
                bbox_inches="tight")
    plt.close(fig)
    fig_rate, _, fig_acc, _ = check.plot_IMU_spectra(config_file,
                                                     nperseg=nperseg2,
                                                     nfiles=pre_nfiles)
    with PdfPages(imu_spectra_fig) as pdf:
        pdf.savefig(fig_rate, bbox_inches="tight")
        plt.close(fig_rate)
        pdf.savefig(fig_acc, bbox_inches="tight")
        plt.close(fig_acc)

    # Look at individual periods
    config = db_flux.parse_config(config_file)
    input_files = config["EC Inputs"]["input_files"]
    if post_nfiles is not None and post_nfiles < len(input_files):
        input_files = np.random.choice(input_files, post_nfiles,
                                       replace=False)
    imu2anem_pos = config["EC Motion Correction"]["imu2anemometer_pos"]
    sample_freq_hz = config["EC Inputs"]["sample_frequency"]
    Tcf, Ta = cutoffs
    for ifile in input_files:
        file_idx = ifile.index(ifile)
        ec_prep, _ = db_flux.prepare_period(ifile, config)
        wind = ec_prep[["wind_speed_u", "wind_speed_v", "wind_speed_w"]]
        acceleration = ec_prep[["acceleration_x",
                                "acceleration_y",
                                "acceleration_z"]]
        rate = ec_prep[["rate_x", "rate_y", "rate_z"]]
        ec_corr = db_flux.wind3D_correct(wind.values, acceleration.values,
                                         rate.values,
                                         ec_prep.heading.values,
                                         ec_prep.speed_over_ground.values,
                                         imu2anem_pos, sample_freq_hz,
                                         Tcf, Ta)
        fig, _, leg = check.plot1_IMU_complementary_spectra(config_file,
                                                            file_idx,
                                                            Tcf=Tcf,
                                                            Ta=Ta,
                                                            nperseg=nperseg1)
        fig.savefig((ifile[:-4] + "_euler_angle_spectra.pdf"),
                    bbox_extra_artists=(leg,), bbox_inches="tight")
        plt.close()
        fig, axs, leg = check.plot_wind3D_spectra(wind.values,
                                                  ec_corr.uvw_ship,
                                                  nperseg=nperseg2)
        fig.savefig((ifile[:-4] + "_uvw_spectra.pdf"),
                    bbox_extra_artists=(leg,), bbox_inches="tight")
        plt.close()
        # Plot corrected wind cospectra
        fig, _, leg = check.plot_wind3D_cospectra(wind.values,
                                                  ec_corr.uvw_ship,
                                                  nperseg=nperseg2)
        fig.savefig((ifile[:-4] + "_uvw_cospectra.pdf"),
                    bbox_extra_artists=(leg,), bbox_inches="tight")
        plt.close()


if __name__ == "__main__":
    _DESCRIPTION = ("Given a configuration file, optionally subset a number "
                    "of input files, and plot key spectral densities for "
                    "aiding in choosing motion correction parameters.")
    _FORMATERCLASS = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(description=_DESCRIPTION,
                                     formatter_class=_FORMATERCLASS)
    parser.add_argument("config_file", metavar="config-file",
                        type=argparse.FileType("r"),
                        help="Path to configuration file")
    parser.add_argument("--pre-nfiles", type=int, default=None,
                        help="Number of input files to subsample for "
                        "plots prior to complementary filter.  Defaults "
                        "to all files in configuration.")
    parser.add_argument("--post-nfiles", type=int, default=None,
                        help="Number of input files to subsample for "
                        "plots showing complementary filter.  Defaults "
                        "to all files in configuration.")
    parser.add_argument("--cutoffs", nargs=2, type=float,
                        metavar=("Tcf", "Ta"), default=[12.0, 60.0],
                        help="Cutoff period for complementary filter, and "
                        "for high-pass filtering accelerations.")
    args = parser.parse_args()
    # Note we just pass the config file name, "type" only used for checking
    # useability
    main(args.config_file.name, pre_nfiles=args.pre_nfiles,
         post_nfiles=args.post_nfiles, cutoffs=args.cutoffs)
