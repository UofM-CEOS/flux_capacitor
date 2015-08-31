#! /usr/bin/env python

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from subprocess import Popen, PIPE
import glob
from StringIO import StringIO
plt.style.use('ggplot')

def main(subset_prog, match_prog, bottle_files, uw_files, ofigure_file,
         target_depth=5, tol_diff=1, depth_fld=6, temperature_fld=7,
         salinity_fld=9, max_time_diff=600, min_flow=2, date_fld=3,
         time_fld=4, flow_fld=16):
    """Perform subsetting and matching of rosette and pCO2 data.

    Generate output figure file with regression of external temperature on
    equilibrator temperature.

    See parser help message for description of arguments.  All arguments
    are coerced to string during execution.

    """
    subset_cmd = (["gawk", "-f", subset_prog,
                   "-v", "target_depth=" + str(target_depth),
                   "-v", "tol_diff=" + str(tol_diff),
                   "-v", "depth_fld=" + str(depth_fld),
                   "-v", "temperature_fld=" + str(temperature_fld),
                   "-v", "salinity_fld=" + str(salinity_fld)] +
                  glob.glob(bottle_files))
    match_cmd = (["gawk", "-f", match_prog,
                  "-v", "max_time_diff=" + str(max_time_diff),
                  "-v", "min_flow=" + str(min_flow),
                  "-v", "date_fld=" + str(date_fld),
                  "-v", "time_fld=" + str(time_fld),
                  "-v", "flow_fld=" + str(flow_fld), "--", "-"] +
                 glob.glob(uw_files))
    bottles = Popen(subset_cmd, stdout=PIPE)
    bottle_matches = Popen(match_cmd, stdin=bottles.stdout, stdout=PIPE)
    rosette = pd.read_csv(StringIO(bottle_matches.communicate()[0]))
    # Fit model
    temp_fit = np.polyfit(rosette.equ_temperature, rosette.temperature, 1)
    temp_predict = np.poly1d(temp_fit)

    plt.figure(figsize=(7, 6))
    x_new = np.linspace(np.min(rosette.equ_temperature),
                        np.max(rosette.equ_temperature))
    plt.plot(x_new, temp_predict(x_new))
    plt.scatter(rosette.equ_temperature, rosette.temperature)
    plt.xlabel("Equilibrator Temperature (C)")
    plt.ylabel("Rosette Temperature (C)")
    plt.title("y={0[1]:.4f} + {0[0]:.4f}x".format(temp_fit))
    plt.tight_layout()
    plt.savefig(ofigure_file, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    import argparse
    import os.path as osp
    _DESCRIPTION = ("Subset Rosette bottles around specified depth, and " +
                    "select underway pCO2 data with a matching time stamp " +
                    "plus/minus specified time difference.")
    # Scripts' directory -- we set AWKPATH to directory where this script
    # is located, and then specify the programs loosely so that awk finds
    # them via its normal search path.  This ensures we find the programs
    # in the most natural way, and avoid hard-coding any location.
    _SCRIPTS_DIR = osp.dirname(osp.realpath(__file__))
    _SUBSET_PROG = osp.join(_SCRIPTS_DIR, "subset_bottles.awk")
    _MATCH_PROG = osp.join(_SCRIPTS_DIR, "pCO2_bottle_match.awk")
    parser = argparse.ArgumentParser(description=_DESCRIPTION)
    # We have to specify every argument to the underlying programs.
    parser.add_argument("ofigure_file", type=argparse.FileType("w"),
                        help="Path to output figure file.")
    parser.add_argument("bottle_files",
                        help="Glob pattern for location of bottle files.")
    parser.add_argument("uw_files",
                    help=("Glob pattern for location of underway pCO2 " +
                          "files (including suffix)."))
    parser.add_argument("--target_depth", default=5, type=float,
                        help="Depth (m) of bottles to extract.")
    parser.add_argument("--tol_diff", default=1, type=float,
                        help=("Maximum difference (m) allowed for " +
                              "deviation from target depth."))
    parser.add_argument("--depth_fld", default=6, type=int,
                        help=("Field where depth is located in average " +
                              "rows of bottle files."))
    parser.add_argument("--temperature_fld", default=7, type=int,
                        help=("Field where temperature is located in " +
                              "average rows of bottle files."))
    parser.add_argument("--salinity_fld", default=9, type=int,
                        help=("Field where salinity is located in " +
                              "average rows of bottle files."))
    parser.add_argument("--max_time_diff", default=600, type=float,
                        help=("Maximum time difference (seconds) allowed " +
                              "between bottle and underway pCO2 data."))
    parser.add_argument("--min_flow", default=2, type=float,
                        help=("Minimum water flow rate allowed in " +
                              "underway pCO2 data (units as in input)."))
    parser.add_argument("--date_fld", default=3, type=int,
                        help=("Field where date is located in underway " +
                              "pCO2 files."))
    parser.add_argument("--time_fld", default=4, type=int,
                        help=("Field where time is located in underway " +
                              " pCO2 files."))
    parser.add_argument("--flow_fld", default=16, type=int,
                        help=("Field where flow rate is located in " +
                              "underway pCO2 files."))
    args = parser.parse_args()
    main(subset_prog=_SUBSET_PROG, target_depth=args.target_depth,
         tol_diff=args.tol_diff, depth_fld=args.depth_fld,
         temperature_fld=args.temperature_fld,
         salinity_fld=args.salinity_fld, bottle_files=args.bottle_files,
         match_prog=_MATCH_PROG, max_time_diff=args.max_time_diff,
         min_flow=args.min_flow, date_fld=args.date_fld,
         time_fld=args.time_fld, flow_fld=args.flow_fld,
         uw_files=args.uw_files, ofigure_file=args.ofigure_file)
