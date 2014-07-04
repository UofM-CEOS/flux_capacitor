;; $Id$
;; Author: Sebastian P. Luque
;; Created: 2013-11-26T23:41:42+0000
;; Last-Updated: 2014-05-04T13:59:50+0000
;;           By: Sebastian Luque
;; ------------------------------------------------------------------------
;; Commentary:
;;
;; All the matching and lining up of data, selection criteria for periods
;; that are adequate for flux analyses, and computations for low frequency
;; computations are now implemented in the database via a number of queries
;; of the data as collected and sanitized.  If you want access to that,
;; download LibreOffice (you will use this as a client for the PostgreSQL
;; server hosting the database):
;;
;; www.libreoffice.org
;;
;; Do this:
;;
;; 1. Start the "Base" (for "Database") LibreOffice application.
;;
;; 2. Select "Connect to an existing database", and select "PostgreSQL"
;;    from the drop-down menu. Click "Next".
;;
;; 3. Type: "host=net82.ceos.umanitoba.ca port=5433 dbname=gases" (without
;;    the quotes) at the "Data Source URL" box. Click "Next".
;;
;; 4. Type: "ceos" for the user name. Click "Test Connection", entering
;;    password: "Gases2014". Click "Next".
;;
;; 5. Click to choose to register the database so that you have a file to
;;    save your work into.
;;
;; You will see a collection of tables.  Refer to the data model I sent you
;; and Tonya to understand the relationships between them.  LibreOffice
;; does not differentiate between tables and "views".  The latter do not
;; have any data per se, but rather are queries that combine data in
;; regular tables to produce the presentation we want.  These are the ones
;; you will find most useful.  The relevant queries for 2013 are those with
;; name suffix "_2013".  They progressively build the data we need from the
;; underlying data tables, ending with "flux_10Hz_2013", which gives us the
;; data we want for flux calculations in IDL.  It contains only the 20
;; minute periods that have passed the different selection criteria.
;;
;; Selecting the complex queries ("views") from the server takes some time
;; (up 5.5 h for "flux_10Hz_2013"), and your connection to the server may
;; time out.  This is something I need to tweak in the server for
;; connections from other computers in the network; this doesn't happen for
;; local connections.  So I have output this query onto a number of files,
;; one per valid 20 minute period.
;;
;; I have collected the files into the compressed archive "EC_2013.zip"
;; (you should have gotten an email with the link):
;;
;; https://www.dropbox.com/s/p9gozhym0eei1vj/EC_2013.zip
;;
;; It contains the following:
;;
;; fluxes.csv              --> Final flux summary for each period.
;; FromDB                  --> Folder containing full data for each period.
;; FromDB/Motion_Corrected --> Subfolder with motion corrected files.
;;
;; The IDL code needed to run the analyses again (using the files in the
;; FromDB folder) is in the archive "flux_capacitor_2013.zip".  Make sure
;; you have the AstroLib and Coyote libraries (see the IDL website for
;; downloading them) *and* are available in your IDL installation.  Edit
;; the control_2013.inc file to reflect where you have unpacked the
;; archive.  Then call "@flux_startup.pro" from the command line, followed
;; by the commands in this script.
;; ------------------------------------------------------------------------
;;; Code:

@control_2013.inc

;; Calculate fluxes

DB_FLUX, ec_period_dir, ec_std_template, 1, ec_daily_rate, $
         ec_period, motpak_offset, sog_thr, lfreq_thr, hfreq_thr, $
         xover_freq_thr, ec_motcorr_dir, ec_footprint_dir, ec_fluxes_file, $
         /overwrite



;;;_ + Emacs Local Variables
;; Local variables:
;; allout-layout: (-2 + : 0)
;; End:
;;
;;; flux_capacitor.pro ends here
