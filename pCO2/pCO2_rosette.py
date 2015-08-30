#! /usr/bin/env python

import argparse
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from subprocess import Popen, PIPE
import glob
from StringIO import StringIO
plt.style.use('ggplot')

# Programs
subset_prog = "subset_bottles.awk"
match_prog = "pCO2_bottle_match.awk"

_DESCRIPTION = ("Subset Rosette bottles around specified depth, and " +
                "select underway pCO2 data with a matching time stamp " +
                "plus/minus specified time difference.")
parser = argparse.ArgumentParser(description=_DESCRIPTION)
parser.add_argument("output_file", type=argparse.FileType("w"),
                    help="Path to output figure file.")
parser.add_argument("bottle_files",
                    help="Glob pattern for location of bottle files.")
parser.add_argument("uw_files",
                    help=("Glob pattern for location of underway pCO2 " +
                          "files (including suffix)."))
parser.add_argument("--target_depth", default=5, type=float,
                    help="Depth (m) of bottles to extract.")
parser.add_argument("--tol_diff", default=1, type=float,
                    help=("Maximum difference (m) allowed for deviation " +
                          "from target depth."))
parser.add_argument("--depth_fld", default=6, type=int,
                    help=("Field where depth is located in average rows " +
                          "of bottle files."))
parser.add_argument("--temperature_fld", default=7, type=int,
                    help=("Field where temperature is located in average " +
                          "rows of bottle files."))
parser.add_argument("--salinity_fld", default=9, type=int,
                    help=("Field where salinity is located in average " +
                          "rows of bottle files."))
parser.add_argument("--max_time_diff", default=600, type=float,
                    help=("Maximum time difference (seconds) allowed " +
                          "between bottle and underway pCO2 data."))
parser.add_argument("--min_flow", default=2, type=float,
                    help=("Minimum water flow rate allowed in underway " +
                          "pCO2 data (units as in input)."))
parser.add_argument("--date_fld", default=3, type=int,
                    help=("Field where date is located in underway pCO2 " +
                          "files."))
parser.add_argument("--time_fld", default=4, type=int,
                    help=("Field where time is located in underway pCO2 " +
                          "files."))
parser.add_argument("--flow_fld", default=16, type=int,
                    help=("Field where flow rate is located in underway " +
                          "pCO2 files."))
args = parser.parse_args()

subset_cmd = (["gawk", "-f", subset_prog,
               "-v", "target_depth=" + str(args.target_depth),
               "-v", "tol_diff=" + str(args.tol_diff),
               "-v", "depth_fld=" + str(args.depth_fld),
               "-v", "temperature_fld=" + str(args.temperature_fld),
               "-v", "salinity_fld=" + str(args.salinity_fld)] +
              glob.glob(args.bottle_files))
match_cmd = (["gawk", "-f", match_prog,
              "-v", "max_time_diff=" + str(args.max_time_diff),
              "-v", "min_flow=" + str(args.min_flow),
              "-v", "date_fld=" + str(args.date_fld),
              "-v", "time_fld=" + str(args.time_fld),
              "-v", "flow_fld=" + str(args.flow_fld), "--", "-"] +
             glob.glob(args.uw_files))
bottles = Popen(subset_cmd, stdout=PIPE)
bottle_matches = Popen(match_cmd, stdin=bottles.stdout, stdout=PIPE)
rosette = pd.read_csv(StringIO(bottle_matches.communicate()[0]))

plt.figure(figsize=(7, 6))
temp_fit = np.polyfit(rosette.equ_temperature, rosette.temperature, 1)
temp_predict = np.poly1d(temp_fit)
x = np.linspace(np.min(rosette.equ_temperature),
                np.max(rosette.equ_temperature))
plt.plot(x, temp_predict(x))
plt.scatter(rosette.equ_temperature, rosette.temperature)
plt.xlabel("Equilibrator Temperature (C)")
plt.ylabel("Rosette Temperature (C)")
plt.title("y={0[1]:.4f} + {0[0]:.4f}x".format(temp_fit))
# leg = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), frameon=False,
#                  borderaxespad=0, ncol=2)
# leg.get_texts()[0].set_text("Tim")
# leg.get_texts()[1].set_text("Database")
plt.tight_layout()
plt.savefig(args.output_file, bbox_inches='tight')
# plt.savefig("external_temperature_2010.png", bbox_extra_artists=(leg,),
#             bbox_inches='tight')
plt.close()

