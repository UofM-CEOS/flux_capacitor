#! /usr/bin/awk -f
# -*- coding: utf-8
# Author: Sebastian P. Luque
# Created: 2013-11-16T19:47:33+0000
# copyright (c) 2013-2017 Sebastian P. Luque
# ------------------------------------------------------------------------
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
# ------------------------------------------------------------------------
# Commentary:
#
# Variables to pass to this script:
#
#     target_depth:    depth of bottles to search
#     tol_diff:        maximum deviation from target depth to allow
#     depth_fld:       depth field index.
#     temperature_fld: temperature field index.
#     salinity_fld:    salinity field index.
#
# We assume data consist of alternating average and standard deviation
# rows, in that order, and:
#
#     bottle# -> 1st field
#     year    -> 4th field
#     month   -> 2nd field
#     day     -> 3rd field
#
## -----------------------------------------------------------------------
# Code:

BEGIN {
    OFS=","
    pn="subset_bottles.awk"
    for (i=1; i < ARGC && ARGV[i] != "-h"; i++){
    	if (ARGV[i] ~ /target_depth/) target_depth=target_depth
    	if (ARGV[i] ~ /tol_diff/) tol_diff=tol_diff
    	if (ARGV[i] ~ /depth_fld/) depth_fld=depth_fld
    	if (ARGV[i] ~ /temperature_fld/) temperature_fld=temperature_fld
    	if (ARGV[i] ~ /salinity_fld/) salinity_fld=salinity_fld
    }
    if (i != ARGC) {usage(pn); exit}
    if (!target_depth) target_depth = 5
    if (!tol_diff) tol_diff = 2
    if (!depth_fld) depth_fld = 6
    if (!temperature_fld) temperature_fld = 7
    if (!salinity_fld) salinity_fld = 8
    print "time,bottle_number,depth,temperature,salinity"
}

/^(\*|#)/ { data_beg=0; next }	# reset the start line, if on header lines
tolower($1) == "bottle" { data_beg = FNR + 2 ; next } # find start line
FNR < data_beg { next }		# skip header lines

(FNR - data_beg) % 2 == 0 {
    mon=fix_month($2)
    date_str=sprintf("%s-%s-%s", $4, mon, $3) # year, month, day
    bottle_no=$1
    depth=$depth_fld
    temperature=$temperature_fld
    salinity=$salinity_fld
    next
}

{
    depth_dev=depth - target_depth
    abs_depth_dev=depth_dev < 0 ? -depth_dev : depth_dev
    if (abs_depth_dev < tol_diff) {
    	printf "%s %s,%s,%s,%s,%s\n", date_str, $1, bottle_no,
    	    depth, temperature, salinity
    }
}


# FUNCTIONS ---------------------------------------------------------------
# Help message
function usage(pn) {
    printf "USAGE:\n\t%s [variables to initialize] file ...\n\n", pn
    printf "VARIABLES:\n\t%-10s %s\n\t%-10s %s\n\t%-10s %s\n\t%-10s %-10s\n\n",
	"target_depth=", "Depth (m) of bottles to extract (Default: 5)",
	"tol_diff=", "Number of meters for maximum deviation from target_depth (Default: 2)",
	"depth_fld=", "Field (integer) where depth is located in average rows (Default: 6)",
	"temperature_fld=", "Field (integer) where temperature is located in average rows (Default: 7)",
	"salinity_fld=", "Field (integer) where salinity is located in average rows (Default: 8)"
    printf "%s\n", "DETAILS:"
    printf "\t%s\n", "Extract date, bottle number, depth, temperature, and salinity of bottles close to a target depth."
}

function fix_month(date) {	# string containing month name substring
    # Translate english/french months to number
    if (match("janfebmaraprmayjunjulaugsepoctnovdec",
	      substr(tolower(date), 1, 3))) {
	return sprintf("%02d", (RSTART + 2) / 3)
    }
    if (match("janvfévrmarsavr-mai-juinjuilaoûtseptoct-nov-déc-",
	      substr(tolower(date), 1, 4))) {
	return sprintf("%02d", (RSTART + 3) / 4)
    }
}



#_ + Emacs local variables
# Local variables:
# allout-layout: (1 + : 0)
# End:
