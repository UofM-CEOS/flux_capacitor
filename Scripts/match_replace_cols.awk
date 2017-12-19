#! /usr/bin/gawk -f
# Author: Sebastian Luque
# Created: 2016-08-21T13:40:36+0000
# copyright (c) 2016, 2017 Sebastian P. Luque
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
# Watch shebang and make sure it matches system's gawk!
#
# Given two sets of files, this script takes a file from the first set as
# its first argument and a (sub)set of files from the second set as the
# next arguments.  It replaces data from one of the files in the first set
# with those from the second set of files based on time stamp.  It returns
# lines from the file given as first argument, with the requested
# replacement data.
#
# Variables required (passed via -v switches):
#
#     timecols_a:   Date and time column(s) in first set of files
#     timecols_b:   Date and time column(s) in second set of files
#     skip_a:       Number of initial lines to skip in first set of files
#     skip_b:       Number of initial lines to skip in second set of files
#     cols_a:       Columns to replace in first set of files
#     cols_b:       Columns to replace with from second set of files
#     nomatch_null: Nullify cols_a columns if no match is found
#     nomatch_val:  Integer value to use when no match is found
#
# We assume:
#
#     o The first characters of the file name in both sets of files have
#       the following convention: YYYY-MM-DDTHHMMSS
#     o Time stamps referred to by timecols_a and timecols_b are assumed to
#       have the same format, and can be compared as given without any
#       manipulation
#     o The line following the skipped lines in both sets of files is a
#       column header line.
#
# Example:
#
# Below replaces data in the Metek folder with those from Gill folder:
#
# ./match_replace_cols.awk -v timecols_a='7,8' -v timecols_b='6,7' \
#     -v skip_a=7 skip_b=7 -v cols_a='19,20,21,22' -v cols_b='14,15,16,17' \
#     Metek/2015-05-12T* Gill/2015-05-12T234646_AIU-0723.data
#
# -------------------------------------------------------------------------
# Code:

BEGIN { FS="[,\t]"; OFS=","	# commas OR tabs separate fields
    pn="match_replace_cols.awk"
    if (length(ARGV) < 2) {usage(pn); exit}
    if (!timecols_a) {
	split("7,8", tcolsarr_a, ",")
    } else split(timecols_a, tcolsarr_a, ",")
    if (!timecols_b) {
	split("6,7", tcolsarr_b, ",")
    } else split(timecols_b, tcolsarr_b, ",")
    if (!skip_a) skip_a=7
    if (!skip_b) skip_b=7
    if (!cols_a) {
	split("19,20,21,22", colsarr_a, ",")
    } else split(cols_a, colsarr_a, ",")
    if (!cols_b) {
	split("14,15,16,17", colsarr_b, ",")
    } else split(cols_b, colsarr_b, ",")
}

NR == FNR && FNR > skip_a {	# process first file (skip header lines)
    if (FNR == (skip_a + 1)) {	# capture header with column names
	hdr_a=$0
	next
    }
    if (length(tcolsarr_a) > 1) { # build time stamp
	time_a=sprintf("%s %s", $tcolsarr_a[1], $tcolsarr_a[2])
    } else time_a=$tcolsarr_a[1]
    file_a[time_a]=$0		# save and index all data
    # Keep track of row number to recover order at the end. TODO: beware of
    # duplicates; index with time stamp and record number of records?
    order[FNR]=time_a
    next
}

# Work on second set of files (skip header line)
FNR > (skip_b + 1) {		      # skip header and non-data lines
    if (length(tcolsarr_b) > 1) { # build time stamp
	time_b=sprintf("%s %s", $tcolsarr_b[1], $tcolsarr_b[2])
    } else time_b=$tcolsarr_b[1]
    # Save and index all data
    file_b[time_b]=$0
}

END {
    # Print column names header
    split(hdr_a, hdr_arr)
    for (i=1; i <= length(hdr_arr); i++) {
	printf "%s,", hdr_arr[i]
    }
    print "is_replaced_flag"
    asorti(order, order_arr, "@ind_num_asc") # sort row numbers in first file
    # Loop first file array in order
    for (i in order_arr) {
	split(file_a[order[order_arr[i]]], irow_a) # parse row data
	# If time stamp matches row in second set of files
	if (order[order_arr[i]] in file_b) {
	    # Parse matched row data
	    split(file_b[order[order_arr[i]]], irow_b)
	    # Loop columns to replace and replace accordingly
	    for (col in colsarr_a) {
		irow_a[colsarr_a[col]]=irow_b[colsarr_b[col]]
	    }
	    replaced_flag=1
	} else {
	    for (col in colsarr_a) {
		if (nomatch_null) irow_a[colsarr_a[col]]=""
		if (nomatch_val) irow_a[colsarr_a[col]]=nomatch_val
	    }
	    replaced_flag=0
	}
	for (j=1; j <= length(irow_a); j++) { # print all columns
	    printf "%s,", irow_a[j]
	}
	print replaced_flag
    }
}


# FUNCTIONS ---------------------------------------------------------------

# Help message
function usage(pn) {
    printf "USAGE:\n\t%s [variables to initialize] file ...\n\n", pn
    printf "VARIABLES:\n\t%-12s %s\n\t%-12s %s\n\t%-12s %s\n\t%-12s %-12s\n",
	"timecols_a", "Date and time column(s) in first set of files (Default: 7,8)",
	"timecols_b", "Date and time column(s) in second set of files (Default: 6,7)",
	"skip_a", "Number of initial lines to skip in first set of files (Default: 7)",
	"skip_b", "Number of initial lines to skip in second set of files (Default: 7)"
    printf "\t%-12s %s\n\t%-12s %s\n\t%-12s %s\n\t%-12s %s\n\n",
	"cols_a", "Columns to replace in first set of files (Default: 19,20,21,22)",
	"cols_b", "Columns to replace in second set of files (Default: 14,15,16,17)",
	"nomatch_null", "Nullify cols_a columns if no match is found (Default: false)",
	"nomatch_val", "Integer value to use when no match is found (Default: original)"
    printf "%s\n", "DETAILS:"
    printf "\t%s\n", "Replace data in first file with matching data from rest of files."
}



#_ + Emacs local variables
# Local variables:
# allout-layout: (1 + : 0)
# End:
