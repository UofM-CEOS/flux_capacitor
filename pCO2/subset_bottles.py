#! /usr/bin/env python
# -*- coding: utf-8
# pylint: disable=too-many-locals,invalid-name,no-member

"""Subset data from Rosette bottle files meeting conditions

This is a port of the AWK utility of the same name, but is approximately
five times slower.  This is the cost of avoiding the inability of
setuptools to adequately handle the UTF-8 encoding required to install the
equivalent AWK tool.

Usage
-----

For help using this script, type:

subset_bottles.py -h

at command line.

"""

import re
from datetime import datetime
import argparse
import sys
import csv


def _french_month(month):
    """Convert French month string to month number"""
    mois = "janvfévrmarsavr-mai-juinjuilaoûtseptoct-nov-déc-"
    mois_loc = re.search(month.lower(), mois.lower())
    if mois_loc:
        mois_no = (mois_loc.start() + 4) / 4
        return "0{}".format(mois_no)


def main(files, target_depth, tol_diff, depth_fld, temperature_fld,
         salinity_fld):
    """Perform bottle data extraction with given conditions"""
    ymd_fmt = "{0}-{1}-{2}"
    # Note change of line terminator to make it compatible with AWK
    writer = csv.writer(sys.stdout, lineterminator='\n')
    writer.writerow(["time", "bottle_number", "depth",
                     "temperature", "salinity"])
    for bottle_f in files:
        with bottle_f:
            for lineno, line in enumerate(bottle_f, start=1):
                if re.search('^(\*|#)', line):
                    continue
                if re.search('^ *bottle', line.lower()):
                    databeg = lineno + 2
                    continue
                if (lineno - databeg) % 2 == 0:
                    parsedl = line.split()
                    mm = _french_month(parsedl[1])
                    if mm:
                        fmt = '%Y-%m-%d %H:%M:%S'
                        yyyymmdd = ymd_fmt.format(parsedl[3], mm, parsedl[2])
                    else:
                        fmt = '%Y-%b-%d %H:%M:%S'
                        yyyymmdd = ymd_fmt.format(parsedl[3], parsedl[1],
                                                  parsedl[2])
                    bottle_no = parsedl[0]
                    depth = parsedl[depth_fld - 1]
                    temperature = parsedl[temperature_fld - 1]
                    salinity = parsedl[salinity_fld - 1]
                else:
                    if lineno <= databeg:
                        continue
                    hhmmss = line.split()[0]
                    ymd = datetime.strptime("{0} {1}".format(yyyymmdd,
                                                             hhmmss), fmt)
                    depth_dev = abs(float(depth) - target_depth)
                    if (depth_dev < tol_diff):
                        line_out = (ymd.strftime('%Y-%m-%d %H:%M:%S'),
                                    bottle_no, depth, temperature,
                                    salinity)
                        writer.writerow(line_out)


if __name__ == "__main__":
    """Main script for subsetting Rosette bottle data from a depth range"""
    _DESCRIPTION = "Subset data from Rosette bottle files meeting conditions"
    _FORMATERCLASS = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(description=_DESCRIPTION,
                                     formatter_class=_FORMATERCLASS)
    parser.add_argument("file", nargs="+", type=argparse.FileType(),
                        help="Path to rosette bottle file(s)")
    parser.add_argument("--target_depth", type=float, default=5,
                        help="Depth (m) of bottles to extract")
    parser.add_argument("--tol_diff", type=float, default=2,
                        help=("Number of meters for maximum deviation from" +
                              "target_depth"))
    parser.add_argument("--depth_fld", type=int, default=6,
                        help=("Field (integer) where depth is located in" +
                              "average rows"))
    parser.add_argument("--temperature_fld", type=int, default=7,
                        help=("Field (integer) where temperature is located" +
                              "in average rows"))
    parser.add_argument("--salinity_fld", type=int, default=8,
                        help=("Field (integer) where salinity is located in" +
                              "average rows"))
    args = parser.parse_args()
    main(args.file, target_depth=args.target_depth, tol_diff=args.tol_diff,
         depth_fld=args.depth_fld, temperature_fld=args.temperature_fld,
         salinity_fld=args.salinity_fld)
