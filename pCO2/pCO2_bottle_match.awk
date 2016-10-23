#! /usr/bin/gawk -f
# Author: Sebastian Luque
# Created: 2013-11-21T13:40:36+0000
# Last-Updated: 2016-10-23T13:16:10+0000
#           By: Sebastian P. Luque
# copyright (c) 2013-2016 Sebastian P. Luque
# -------------------------------------------------------------------------
# This program is Free Software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with GNU Emacs; see the file COPYING. If not, write to the
# Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
# -------------------------------------------------------------------------
# Commentary:
#
# Watch shebang and make sure it matches system!
#
# This script takes a bottle file, like those returned by
# subset_bottles.awk, as first argument and any number of pCO2 data files.
# For each bottle record, it selects the pCO2 data that are within a given
# maximum time difference from it, and that are higher than a given flow
# threshold.
#
# Variables required:
#
#     max_time_diff:   Maximum time difference (seconds).
#     min_flow:        Minimum water flow allowed.
#     date_fld:        Date field index in pCO2 files.  We assume format
#                      "DD/MM/YY"...
#     time_fld:        Time field in pCO2 files.  We assume format
#                      "HH:MM:SS".
#     flow_fld:        Water flow field in pCO2 files.
#
# We assume:
#
#     o Bottle and pCO2 files have a single header line.
#     o Full time stamp is available in bottle file, as single field (1st).
#     o The first field in pCO2 files is an identifier, and we select only
#       records matching the string "EQU" (case irrelavant).
#
# Example:
#
# ./pCO2_bottle_match.awk -v max_time_diff=900 -v min_flow=1.5 \
#     -v date_fld=3 -v time_fld=4 -v flow_fld=16 bottle_2to7m.csv *
#
# -------------------------------------------------------------------------
# Code:

BEGIN { FS="[,\t]"; OFS="," }	# commas OR tabs separate fields

NR == FNR && FNR > 1 {		# process bottle file (skip header line)
    split($1, a, /[ :-]/)
    tstr_bottle=sprintf("%s %s %s %s %s %s",
			a[1], a[2], a[3], a[4], a[5], a[6])
    t_bottle=mktime(tstr_bottle)
    # Save everything from bottle (bottle_number, depth, temperature,
    # salinity)
    bottle[t_bottle]=sprintf("%s,%s,%s,%s", $2, $3, $4, $5)
    next
}

# Work on the pCO2 files (skip header line)
FNR > 1 && (tolower($1) == "equ") {   # matching "EQU" rows exactly
    d_pCO2=fix_date_string($date_fld)	# date part
    split($time_fld, a, /[ :]/)
    tstr_pCO2=sprintf("%s %s %s %s", d_pCO2, a[1], a[2], a[3])
    t_pCO2=mktime(tstr_pCO2)
    if (t_pCO2 < 0) {
	printf "Faulty data %s, %s\n", FILENAME, FNR > "/dev/stderr"
	next
    }
    # Save all pCO2 data.  Just equilibration temperature for now.
    pCO2[t_pCO2]=sprintf("%s,%s", $5, $flow_fld)
}

END {
    printf "%s,%s\n", "bottle_time,pCO2_time,bottle_number,depth",
	"temperature,salinity,equ_temperature,H2O_flow"
    # Loop through pCO2, compare time stamps with bottle.  Sort first.
    n=asorti(pCO2, pCO2_srt)	# now pCO2_srt contains sorted indices
    m=asorti(bottle, bottle_srt) # bottle_srt contains sorted indices
    for (i=1; i <= n; i++) {	# scan each sorted index in pCO2
    	split(pCO2[pCO2_srt[i]], p, ",") # extract saved data
    	for (j=1; j <= m; j++) {	 # scan each sorted index in bottle
    	    tdiff=bottle_srt[j] - pCO2_srt[i]
    	    tdiff=tdiff >= 0 ? tdiff : -tdiff
	    # If current difference is smaller than previous, or we're at
	    # first iteration, record a minimum
	    if ((!tdiff_min) || tdiff < tdiff_min) {
		tdiff_min=tdiff	     # time difference
		bottle_t=bottle_srt[j]   # bottle time
		bottle_data=bottle[bottle_srt[j]]
	    }
    	}
	# We've gone through all bottle records, and we can determine
	# whether minimum is smaller than the maximum difference allowed,
	# and if flow is greater than minimum, then print record
	if ((tdiff_min < max_time_diff) && (p[2] > min_flow))
	    printf "%s,%s,%s,%s\n",
		strftime("%F %T", bottle_t),
	    	strftime("%F %T", pCO2_srt[i]), bottle_data,
		pCO2[pCO2_srt[i]]
	tdiff_min=tdiff		# reset the difference
    }
}


# FUNCTIONS ---------------------------------------------------------------

# Fix stupid input date format
function fix_date_string(s,	a) {
  split(s, a, /[ \/]/)
  # Assuming current century!
  return sprintf("%s %s %s", 2000+a[3], a[2], a[1])
}




#_ + Emacs local variables
# Local variables:
# allout-layout: (1 + : 0)
# End:
