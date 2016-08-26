#! /usr/bin/env python
# pylint: disable=too-many-locals,invalid-name,no-member

"""Replace data in one set of files with those from another set

This is a wrapper layer for batch processing files in a directory set,
performing the operations in underlying AWK script of same name. Check the
assumptions in that script before using this.

Usage
-----

For help using this script, type:

match_replace_cols.py -h

at command line.

"""

import argparse
import os.path as osp
from subprocess import Popen
import re
import glob
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
plt.style.use('ggplot')

__version__ = "0.1.0"

def main(base, source, **kwargs):
    """Replace data in base files with data from source files

    See parser help message for description of arguments.
    """
    match_prog = kwargs.get("match_prog")
    timecols_base = ",".join(str(v) for v in kwargs.get("timecols_base"))
    timecols_source = ",".join(str(v) for v in kwargs.get("timecols_source"))
    skip_base = str(kwargs.get("skip_base"))
    skip_source = str(kwargs.get("skip_source"))
    cols_base = ",".join(str(v) for v in kwargs.get("cols_base"))
    cols_source = ",".join(str(v) for v in kwargs.get("cols_source"))
    match_cmd_pre = ["gawk", "-f", match_prog,
                     "-v", "timecols_a=" + timecols_base,
                     "-v", "timecols_b=" + timecols_source,
                     "-v", "skip_a=" + skip_base,
                     "-v", "skip_b=" + skip_source,
                     "-v", "cols_a=" + cols_base,
                     "-v", "cols_b=" + cols_source, "--"]
    src_files = glob.glob(source)
    for basef in glob.glob(base):
        # Get a file name prefix to be shared by the output file
        iname = osp.basename(basef)
        iname_prefix = osp.splitext(iname)[0]
        # Construct a name for file with replacements
        oname = (osp.join(osp.dirname(basef), iname_prefix) + "_alt" +
                 osp.splitext(iname)[1])
        # Extract date part of file name
        basef_date = iname_prefix[:10]
        # Find source files from current base's file date
        src_m = [f for f in src_files if re.match(basef_date, osp.basename(f))]
        # If no files match, then skip to next source file.  Or do
        # something else?
        if not src_m:
            # match_cmd = match_cmd_pre + [basef] + [basef]
            continue
        else:                   # we have matching files so do replacements
            match_cmd = match_cmd_pre + [basef] + src_m

        # Write stdout to file
        with open(oname, "w") as rpl_file:
            rpl = Popen(match_cmd, stdout=rpl_file)
            rpl.communicate()


if __name__ == "__main__":
    _DESCRIPTION = ("Based on matching time stamp, substitute data in "
                    "selected columns in a set of base files with data "
                    "from a set of source files having comparable structure")
    # Scripts' directory -- we set AWKPATH to directory where this script
    # is located, and then specify the programs loosely so that awk finds
    # them via its normal search path.  This ensures we find the programs
    # in the most natural way, and avoid hard-coding any location
    _SCRIPTS_DIR = osp.dirname(osp.realpath(__file__))
    _MATCH_PROG = osp.join(_SCRIPTS_DIR, "match_replace_cols.awk")
    parser = argparse.ArgumentParser(description=_DESCRIPTION,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    group = parser.add_argument_group("required arguments")
    # We have to specify every argument to the underlying program
    parser.add_argument("base",
                        help="Glob pattern for location of base files")
    parser.add_argument("source",
                        help=("Glob pattern for location of source files "
                              "supplying data to substitute in base files"))
    # group.add_argument("--ofigure-file", required=True,
    #                    type=argparse.FileType("w"),
    #                    help="Path to output figure file.")
    parser.add_argument("--timecols-base", nargs='+', type=int, default=[7, 8],
                        help="Date and time column(s) in base set of files")
    parser.add_argument("--timecols-source", nargs='+', type=int,
                        default=[6, 7],
                        help="Date and time column(s) in source set of files")
    parser.add_argument("--skip-base", type=int, default=7,
                        help=("Number of lines to skip from top of base "
                              "set of files"))
    parser.add_argument("--skip-source", type=int, default=7,
                        help=("Number of lines to skip from top of source "
                              "set of files"))
    parser.add_argument("--cols-base", nargs='+', type=int,
                        default=[19, 20, 21, 22],
                        help="Columns to replace in base set of files")
    parser.add_argument("--cols-source", nargs='+', type=int,
                        default=[14, 15, 16, 17],
                        help="Columns to replace in source set of files")
    parser.add_argument("--version", action="version",
                        version="%(prog)s {}".format(__version__))
    args = parser.parse_args()
    main(args.base, args.source, match_prog=_MATCH_PROG,
         timecols_base=args.timecols_base,
         timecols_source=args.timecols_source,
         skip_base=args.skip_base, skip_source=args.skip_source,
         cols_base=args.cols_base, cols_source=args.cols_source)
