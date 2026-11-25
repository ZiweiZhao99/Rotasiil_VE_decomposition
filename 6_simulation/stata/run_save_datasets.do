****************************************************
* run_save_datasets.do
*
* Generate and SAVE 100 simulated datasets at theta=1 and 100 at theta=0.2,
* using the SAME data-generating process and seeds as the estimation runs
* (run_parallel_tt_full.do + rep_one_tt_full.do).
*
* Seeds: seed = seed0 + r, with seed0 = 100000 and r = 1..100, matching the
* original driver. Provided the DGP and seed are identical, the saved dataset
* for replicate r is bit-for-bit the same data that was analyzed in that
* replicate. That is what makes the Stata-vs-R comparison a same-data
* comparison.
*
* Output:
*   out_tt_full_1\simdata_001.dta ... simdata_100.dta   (theta = 1)
*   out_tt_full_0.2\simdata_001.dta ... simdata_100.dta (theta = 0.2)
* plus matching .csv for easy reading into R.
*
* NOTE: This runs serially in a single Stata session (no shell/parallel),
* so it needs no external Stata path. 100 datasets are cheap to simulate;
* the cost in the original runs was the estimation, not the simulation.
* Edit `repdo' and the two `outdir' paths to your machine.
****************************************************

version 15.1
clear all
set more off

*--------------------------
* Settings — EDIT THESE PATHS
*--------------------------
local R      = 100
local seed0  = 100000
local repdo  "rep_one_tt_save_data.do"

* Two output folders, one per theta
local outdir1   "out_tt_full_1"
local outdir02  "out_tt_full_0.2"

cap mkdir "`outdir1'"
cap mkdir "`outdir02'"

*--------------------------
* theta = 1
*--------------------------
di as text _n "==== Saving 100 datasets at theta = 1 ===="
forvalues r = 1/`R' {
    local seed = `seed0' + `r'
    do "`repdo'" `r' `seed' "`outdir1'" 1
}

*--------------------------
* theta = 0.2
*--------------------------
di as text _n "==== Saving 100 datasets at theta = 0.2 ===="
forvalues r = 1/`R' {
    local seed = `seed0' + `r'
    do "`repdo'" `r' `seed' "`outdir02'" 0.2
}

di as result _n "======================================================"
di as result "Done. Datasets in:"
di as result "  `outdir1'\simdata_001.dta ... simdata_100.dta   (theta=1)"
di as result "  `outdir02'\simdata_001.dta ... simdata_100.dta  (theta=0.2)"
di as result "======================================================"

exit 0
