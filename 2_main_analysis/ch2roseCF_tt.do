****************************************************
* ch2roseCF_tt.do
*
* Main analysis of the Rotasiil (ROSE) trial: composite and primary effects.
*   AG model (stcox, robust SE clustered on child)  -> composite effect
*   CF model (strmcure, gamma shared frailty,
*             strata by event order)                -> primary effect, theta
* Time-varying effect: vactime = group * log(t), evaluated on data split at
* every observed time (ROSE_split_tt.csv, 1,368,708 rows; built by
* 1_data_preparation/prepare_data.R). Stata imports the R column
* "vac.time" as "vactime".
* NOTE: without lastpool(), strmcure keeps every value of event_count as its
* own stratum (0, 1, 2, 3; only 2 failures in stratum 2). The Methods state
* that episodes of order >= 2 were pooled into one stratum, which would be
* strmcure ..., strata(event_count) lastpool(1). Decide which to report.
*
* Expected: AG b0 = -.549, b1 = .186;
*           CF b0 = -.5715, b1 = .1776, theta = .8236 (95% CI .7507-.8965)
*
* Usage: set the data folder below (it must contain ROSE_split_tt.csv),
*        run from this folder; the log is written next to this file.
* Requires: strmcure (Xu & Cheung, Stata J 2015; see ../strmcure_ado).
****************************************************

capture log close
log using ch2roseCF_tt, text replace

global datadir "../data"     // <-- folder containing ROSE_split_tt.csv

clear
set more off
set seed 123

import delimited "$datadir/ROSE_split_tt.csv"



stset end, fail(event) id(id) time0(start) exit(time .)

* AG model
stcox group vactime, cluster(id) nohr
   
* CF model
strmcure group vactime, shared(id) strata(event_count) dots ///
   reps(100) repspart(5) obspart(30) rspart(15)

log close
