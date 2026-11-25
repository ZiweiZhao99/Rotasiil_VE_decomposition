****************************************************
* summarize_partial.do
*
* Check which replicates have finished (done flag exists),
* load and combine only those, compute summary stats.
*
* Can be run at any time while simulation is in progress.
*
* Usage:
*   do summarize_partial.do [outdir]
****************************************************

version 15.1
clear all
set more off

args outdir
if "`outdir'"=="" local outdir = "C:\Users\e1139775\NUS Dropbox\Ziwei Zhao\STATA simulation tt splitting\out_tt_full"

local R = 100

local b0_true = -0.57
local b1_true = 0.17
scalar auc_true = 1 - exp(`b0_true') * (20^`b1_true') / (`b1_true' + 1)
local auc_true = scalar(auc_true)

****************************************************
* Find completed replicates
****************************************************
local done_list ""
local n_done = 0
local n_missing = 0

forvalues r = 1/`R' {
    capture confirm file "`outdir'\done_`=string(`r',"%03.0f")'.txt"
    if (!_rc) {
        capture confirm file "`outdir'\rep_`=string(`r',"%03.0f")'.dta"
        if (!_rc) {
            local done_list "`done_list' `r'"
            local n_done = `n_done' + 1
        }
        else {
            local n_missing = `n_missing' + 1
        }
    }
}

di as text "Completed replicates: `n_done' out of `R'"
di as text "Missing dta despite done flag: `n_missing'"

if (`n_done' == 0) {
    di as error "No completed replicates found in `outdir'"
    exit 601
}

****************************************************
* Load and combine completed reps
****************************************************
local first = 1
foreach r of local done_list {
    if (`first') {
        use "`outdir'\rep_`=string(`r',"%03.0f")'.dta", clear
        local first = 0
    }
    else {
        append using "`outdir'\rep_`=string(`r',"%03.0f")'.dta"
    }
}

di as text "Loaded `n_done' replicates."
di as text ""

****************************************************
* Basic descriptives
****************************************************
di as text "======================================================"
di as text "PARTIAL SUMMARY (`n_done' / `R' reps completed)"
di as text "True: beta0=`b0_true'  beta1=`b1_true'  AUC=" %6.4f `auc_true'
di as text "======================================================"

summarize failures fe_n_events, detail

* Event timing (if available)
capture confirm variable med_first
if (!_rc) {
    quietly summarize med_first, meanonly
    di as text _n "--- Event timing (mean across completed reps) ---"
    di as text "First-event median: " %5.1f r(mean) "m"
    capture quietly summarize med_later, meanonly
    if (!_rc) di as text "Later-event median: " %5.1f r(mean) "m"
    capture quietly summarize n_later, meanonly
    if (!_rc) di as text "Mean later events:  " %5.1f r(mean)
    capture quietly summarize pct_first_b8, meanonly
    if (!_rc) di as text "First % <8m:        " %5.1f r(mean) "%"
    capture quietly summarize pct_later_b8, meanonly
    if (!_rc) di as text "Later % <8m:        " %5.1f r(mean) "%"
}

****************************************************
* Empirical SE
****************************************************

* --- FE ---
quietly summarize fe_b0hat if fe_b0hat < .
scalar fe_b0_mean  = r(mean)
scalar fe_b0_empSE = r(sd)
scalar fe_b0_N     = r(N)

quietly summarize fe_b1hat if fe_b1hat < .
scalar fe_b1_mean  = r(mean)
scalar fe_b1_empSE = r(sd)

quietly summarize fe_auc if fe_auc < .
scalar fe_auc_mean  = r(mean)
scalar fe_auc_empSE = r(sd)

* --- AG ---
quietly summarize ag_b0hat if ag_b0hat < .
scalar ag_b0_mean  = r(mean)
scalar ag_b0_empSE = r(sd)
scalar ag_b0_N     = r(N)

quietly summarize ag_b1hat if ag_b1hat < .
scalar ag_b1_mean  = r(mean)
scalar ag_b1_empSE = r(sd)

quietly summarize ag_auc if ag_auc < .
scalar ag_auc_mean  = r(mean)
scalar ag_auc_empSE = r(sd)

* --- CF (if present) ---
capture confirm variable cf_b0hat
local has_cf = (!_rc)

if (`has_cf') {
    capture confirm variable cf_failed
    if (!_rc) {
        quietly summarize cf_b0hat if cf_b0hat < . & cf_failed == 0
    }
    else {
        quietly summarize cf_b0hat if cf_b0hat < .
    }
    scalar cf_b0_mean  = r(mean)
    scalar cf_b0_empSE = r(sd)
    scalar cf_b0_N     = r(N)

    capture {
        if (!missing("`cf_failed'")) {
            quietly summarize cf_b1hat if cf_b1hat < . & cf_failed == 0
        }
        else {
            quietly summarize cf_b1hat if cf_b1hat < .
        }
    }
    scalar cf_b1_mean  = r(mean)
    scalar cf_b1_empSE = r(sd)

    capture {
        if (!missing("`cf_failed'")) {
            quietly summarize cf_auc if cf_auc < . & cf_failed == 0
        }
        else {
            quietly summarize cf_auc if cf_auc < .
        }
    }
    scalar cf_auc_mean  = r(mean)
    scalar cf_auc_empSE = r(sd)

    capture {
        if (!missing("`cf_failed'")) {
            quietly summarize theta if theta < . & cf_failed == 0
        }
        else {
            quietly summarize theta if theta < .
        }
    }
    scalar theta_mean  = r(mean)
    scalar theta_empSE = r(sd)
}

****************************************************
* Coverage
****************************************************

* --- FE ---
gen fe_cover_b0 = (fe_b0hat - 1.96*scalar(fe_b0_empSE) <= `b0_true') & ///
                  (`b0_true' <= fe_b0hat + 1.96*scalar(fe_b0_empSE)) ///
                  if fe_b0hat < .

gen fe_cover_b1 = (fe_b1hat - 1.96*scalar(fe_b1_empSE) <= `b1_true') & ///
                  (`b1_true' <= fe_b1hat + 1.96*scalar(fe_b1_empSE)) ///
                  if fe_b1hat < .

gen fe_cover_auc = (fe_auc - 1.96*scalar(fe_auc_empSE) <= `auc_true') & ///
                   (`auc_true' <= fe_auc + 1.96*scalar(fe_auc_empSE)) ///
                   if fe_auc < .

* --- AG ---
gen ag_cover_b0 = (ag_b0hat - 1.96*scalar(ag_b0_empSE) <= `b0_true') & ///
                  (`b0_true' <= ag_b0hat + 1.96*scalar(ag_b0_empSE)) ///
                  if ag_b0hat < .

gen ag_cover_b1 = (ag_b1hat - 1.96*scalar(ag_b1_empSE) <= `b1_true') & ///
                  (`b1_true' <= ag_b1hat + 1.96*scalar(ag_b1_empSE)) ///
                  if ag_b1hat < .

gen ag_cover_auc = (ag_auc - 1.96*scalar(ag_auc_empSE) <= `auc_true') & ///
                   (`auc_true' <= ag_auc + 1.96*scalar(ag_auc_empSE)) ///
                   if ag_auc < .

* --- CF (if present) ---
if (`has_cf') {
    capture confirm variable cf_failed
    local cf_cond = cond(!_rc, "& cf_failed == 0", "")

    gen cf_cover_b0 = (cf_b0hat - 1.96*scalar(cf_b0_empSE) <= `b0_true') & ///
                      (`b0_true' <= cf_b0hat + 1.96*scalar(cf_b0_empSE)) ///
                      if cf_b0hat < . `cf_cond'

    gen cf_cover_b1 = (cf_b1hat - 1.96*scalar(cf_b1_empSE) <= `b1_true') & ///
                      (`b1_true' <= cf_b1hat + 1.96*scalar(cf_b1_empSE)) ///
                      if cf_b1hat < . `cf_cond'

    gen cf_cover_auc = (cf_auc - 1.96*scalar(cf_auc_empSE) <= `auc_true') & ///
                       (`auc_true' <= cf_auc + 1.96*scalar(cf_auc_empSE)) ///
                       if cf_auc < . `cf_cond'
}

****************************************************
* Print results
****************************************************
di as text _n "--- First-Event (N=" scalar(fe_b0_N) ") ---"
di as text "beta0:  mean=" %9.4f scalar(fe_b0_mean) "  bias=" %9.4f scalar(fe_b0_mean)-`b0_true' "  empSE=" %9.4f scalar(fe_b0_empSE)
di as text "beta1:  mean=" %9.4f scalar(fe_b1_mean) "  bias=" %9.4f scalar(fe_b1_mean)-`b1_true' "  empSE=" %9.4f scalar(fe_b1_empSE)
di as text "AUC:    mean=" %9.4f scalar(fe_auc_mean) "  bias=" %9.4f scalar(fe_auc_mean)-`auc_true' "  empSE=" %9.4f scalar(fe_auc_empSE)

quietly summarize fe_cover_b0, meanonly
di as text "Coverage beta0: " %5.1f 100*r(mean) "%"
quietly summarize fe_cover_b1, meanonly
di as text "Coverage beta1: " %5.1f 100*r(mean) "%"
quietly summarize fe_cover_auc, meanonly
di as text "Coverage AUC:   " %5.1f 100*r(mean) "%"

di as text _n "--- AG (N=" scalar(ag_b0_N) ") ---"
di as text "beta0:  mean=" %9.4f scalar(ag_b0_mean) "  bias=" %9.4f scalar(ag_b0_mean)-`b0_true' "  empSE=" %9.4f scalar(ag_b0_empSE)
di as text "beta1:  mean=" %9.4f scalar(ag_b1_mean) "  bias=" %9.4f scalar(ag_b1_mean)-`b1_true' "  empSE=" %9.4f scalar(ag_b1_empSE)
di as text "AUC:    mean=" %9.4f scalar(ag_auc_mean) "  bias=" %9.4f scalar(ag_auc_mean)-`auc_true' "  empSE=" %9.4f scalar(ag_auc_empSE)

quietly summarize ag_cover_b0, meanonly
di as text "Coverage beta0: " %5.1f 100*r(mean) "%"
quietly summarize ag_cover_b1, meanonly
di as text "Coverage beta1: " %5.1f 100*r(mean) "%"
quietly summarize ag_cover_auc, meanonly
di as text "Coverage AUC:   " %5.1f 100*r(mean) "%"

if (`has_cf') {
    di as text _n "--- CF (N=" scalar(cf_b0_N) ") ---"
    di as text "beta0:  mean=" %9.4f scalar(cf_b0_mean) "  bias=" %9.4f scalar(cf_b0_mean)-`b0_true' "  empSE=" %9.4f scalar(cf_b0_empSE)
    di as text "beta1:  mean=" %9.4f scalar(cf_b1_mean) "  bias=" %9.4f scalar(cf_b1_mean)-`b1_true' "  empSE=" %9.4f scalar(cf_b1_empSE)
    di as text "theta:  mean=" %9.4f scalar(theta_mean)  "  bias=" %9.4f scalar(theta_mean)-1 "  empSE=" %9.4f scalar(theta_empSE)
    di as text "AUC:    mean=" %9.4f scalar(cf_auc_mean) "  bias=" %9.4f scalar(cf_auc_mean)-`auc_true' "  empSE=" %9.4f scalar(cf_auc_empSE)

    quietly summarize cf_cover_b0, meanonly
    di as text "Coverage beta0: " %5.1f 100*r(mean) "%"
    quietly summarize cf_cover_b1, meanonly
    di as text "Coverage beta1: " %5.1f 100*r(mean) "%"
    quietly summarize cf_cover_auc, meanonly
    di as text "Coverage AUC:   " %5.1f 100*r(mean) "%"

    capture confirm variable cf_failed
    if (!_rc) {
        quietly summarize cf_failed, meanonly
        di as text "CF failed: " %3.0f r(sum) " out of " %3.0f r(N)
    }
}

di as text _n "--- Key comparison ---"
di as text "FE beta1 = " %7.4f scalar(fe_b1_mean)
di as text "AG beta1 = " %7.4f scalar(ag_b1_mean)
di as text "Gap      = " %7.4f scalar(ag_b1_mean) - scalar(fe_b1_mean)
if (`has_cf') {
    di as text "CF beta1 = " %7.4f scalar(cf_b1_mean)
}

* Runtimes (if available)
capture quietly summarize rt_sim_sec, meanonly
if (!_rc) di as text _n "Mean runtime (sec): sim=" %5.1f r(mean)
capture quietly summarize rt_fe_sec, meanonly
if (!_rc) di as text "                    FE=" %5.1f r(mean)
capture quietly summarize rt_ag_sec, meanonly
if (!_rc) di as text "                    AG=" %5.1f r(mean)
capture quietly summarize rt_cf_sec, meanonly
if (!_rc) di as text "                    CF=" %5.1f r(mean)

di as result _n "======================================================"
di as result "Partial summary complete (`n_done'/`R' reps)."
di as result "Re-run this script later to include more completed reps."
di as result "======================================================"

exit 0
