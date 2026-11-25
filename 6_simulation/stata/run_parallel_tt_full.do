****************************************************
* run_parallel_tt_full.do
*
* Launch 100 reps of rep_one_tt_full.do in parallel,
* combine results, compute empirical SE and coverage.
****************************************************

version 15.1
clear all
set more off

*--------------------------
* Settings
*--------------------------
local R      = 50
local cores  = 25
local seed0  = 100000

local stataexe C:\Program Files (x86)\Stata15\StataMP-64.exe
local repdo   C:\Users\e1139775\NUS Dropbox\Ziwei Zhao\STATA simulation tt splitting\rep_one_tt_full.do
local outdir  C:\Users\e1139775\NUS Dropbox\Ziwei Zhao\STATA simulation tt splitting\out_tt_full_0.2

cap mkdir "`outdir'"

*--------------------------
* Clean old files
*--------------------------
forvalues r = 1/`R' {
    cap erase "`outdir'\done_`=string(`r',"%03.0f")'.txt"
    cap erase "`outdir'\rep_`=string(`r',"%03.0f")'.dta"
    cap erase "`outdir'\log_`=string(`r',"%03.0f")'.smcl"
}
cap erase "`outdir'\results_all.dta"
cap erase "`outdir'\results_all.csv"
cap erase "`outdir'\summary_stats.dta"
cap erase "`outdir'\summary_stats.csv"

*--------------------------
* Launch in batches
*--------------------------
local nbatch = ceil(`R'/`cores')

forvalues b = 0/`=`nbatch'-1' {
    local start = `b'*`cores' + 1
    local end   = min(`R', `start' + `cores' - 1)

    di as text "Launching batch `=`b'+1' of `nbatch' : reps `start' to `end'"

    forvalues r = `start'/`end' {
        local seed = `seed0' + `r'
        shell start "" /min "`stataexe'" /e do "`repdo'" `r' `seed' "`outdir'"
    }

    local pending = 1
    while (`pending') {
        local pending = 0
        forvalues r = `start'/`end' {
            capture confirm file "`outdir'\done_`=string(`r',"%03.0f")'.txt"
            if (_rc) local pending = 1
        }
        if (`pending') sleep 10000
    }

    di as result "Batch `=`b'+1' complete."
}

****************************************************
* Combine all replicate results
****************************************************
capture confirm file "`outdir'\rep_001.dta"
if (_rc) {
    di as error "rep_001.dta not found."
    exit 601
}

use "`outdir'\rep_001.dta", clear

forvalues r = 2/`R' {
    capture confirm file "`outdir'\rep_`=string(`r',"%03.0f")'.dta"
    if (_rc) {
        di as error "Missing: rep_`=string(`r',"%03.0f")'.dta"
        continue
    }
    append using "`outdir'\rep_`=string(`r',"%03.0f")'.dta"
}

save "`outdir'\results_all.dta", replace
export delimited using "`outdir'\results_all.csv", replace

****************************************************
* True values
****************************************************
quietly summarize beta0_true, meanonly
local b0_true = r(mean)
quietly summarize beta1_true, meanonly
local b1_true = r(mean)
quietly summarize theta_true, meanonly
local th_true = r(mean)
quietly summarize auc_true, meanonly
local auc_true = r(mean)

****************************************************
* Empirical SE
****************************************************

* --- FE ---
quietly summarize fe_b0hat if fe_b0hat < .
scalar fe_b0_mean   = r(mean)
scalar fe_b0_empSE  = r(sd)
scalar fe_b0_N      = r(N)

quietly summarize fe_b1hat if fe_b1hat < .
scalar fe_b1_mean   = r(mean)
scalar fe_b1_empSE  = r(sd)

quietly summarize fe_auc if fe_auc < .
scalar fe_auc_mean  = r(mean)
scalar fe_auc_empSE = r(sd)

* --- AG ---
quietly summarize ag_b0hat if ag_b0hat < .
scalar ag_b0_mean   = r(mean)
scalar ag_b0_empSE  = r(sd)
scalar ag_b0_N      = r(N)

quietly summarize ag_b1hat if ag_b1hat < .
scalar ag_b1_mean   = r(mean)
scalar ag_b1_empSE  = r(sd)

quietly summarize ag_auc if ag_auc < .
scalar ag_auc_mean  = r(mean)
scalar ag_auc_empSE = r(sd)

* --- CF (exclude failed) ---
quietly summarize cf_b0hat if cf_b0hat < . & cf_failed == 0
scalar cf_b0_mean   = r(mean)
scalar cf_b0_empSE  = r(sd)
scalar cf_b0_N      = r(N)

quietly summarize cf_b1hat if cf_b1hat < . & cf_failed == 0
scalar cf_b1_mean   = r(mean)
scalar cf_b1_empSE  = r(sd)

quietly summarize cf_auc if cf_auc < . & cf_failed == 0
scalar cf_auc_mean  = r(mean)
scalar cf_auc_empSE = r(sd)

quietly summarize theta if theta < . & cf_failed == 0
scalar theta_mean   = r(mean)
scalar theta_empSE  = r(sd)

****************************************************
* Coverage using empirical SE
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

* --- CF ---
gen cf_cover_b0 = (cf_b0hat - 1.96*scalar(cf_b0_empSE) <= `b0_true') & ///
                  (`b0_true' <= cf_b0hat + 1.96*scalar(cf_b0_empSE)) ///
                  if cf_b0hat < . & cf_failed == 0

gen cf_cover_b1 = (cf_b1hat - 1.96*scalar(cf_b1_empSE) <= `b1_true') & ///
                  (`b1_true' <= cf_b1hat + 1.96*scalar(cf_b1_empSE)) ///
                  if cf_b1hat < . & cf_failed == 0

gen cf_cover_auc = (cf_auc - 1.96*scalar(cf_auc_empSE) <= `auc_true') & ///
                   (`auc_true' <= cf_auc + 1.96*scalar(cf_auc_empSE)) ///
                   if cf_auc < . & cf_failed == 0

gen cf_cover_theta = (theta - 1.96*scalar(theta_empSE) <= `th_true') & ///
                     (`th_true' <= theta + 1.96*scalar(theta_empSE)) ///
                     if theta < . & cf_failed == 0

save "`outdir'\results_all.dta", replace

****************************************************
* Print summary
****************************************************
di as text _n "======================================================"
di as text "SUMMARY: Event-time split (tt), declining baseline"
di as text "True: beta0=`b0_true'  beta1=`b1_true'  theta=`th_true'  AUC=" %6.4f `auc_true'
di as text "======================================================"

quietly summarize cf_failed, meanonly
di as text "CF failed: " %3.0f r(sum) " out of " %3.0f r(N) _n

di as text "--- First-Event ---"
di as text "N=" scalar(fe_b0_N)
di as text "beta0:  mean=" %9.4f scalar(fe_b0_mean) "  bias=" %9.4f scalar(fe_b0_mean)-`b0_true' "  empSE=" %9.4f scalar(fe_b0_empSE)
di as text "beta1:  mean=" %9.4f scalar(fe_b1_mean) "  bias=" %9.4f scalar(fe_b1_mean)-`b1_true' "  empSE=" %9.4f scalar(fe_b1_empSE)
di as text "AUC:    mean=" %9.4f scalar(fe_auc_mean) "  bias=" %9.4f scalar(fe_auc_mean)-`auc_true' "  empSE=" %9.4f scalar(fe_auc_empSE)

quietly summarize fe_cover_b0, meanonly
di as text "Coverage beta0: " %5.1f 100*r(mean) "%"
quietly summarize fe_cover_b1, meanonly
di as text "Coverage beta1: " %5.1f 100*r(mean) "%"
quietly summarize fe_cover_auc, meanonly
di as text "Coverage AUC:   " %5.1f 100*r(mean) "%"

di as text _n "--- AG ---"
di as text "N=" scalar(ag_b0_N)
di as text "beta0:  mean=" %9.4f scalar(ag_b0_mean) "  bias=" %9.4f scalar(ag_b0_mean)-`b0_true' "  empSE=" %9.4f scalar(ag_b0_empSE)
di as text "beta1:  mean=" %9.4f scalar(ag_b1_mean) "  bias=" %9.4f scalar(ag_b1_mean)-`b1_true' "  empSE=" %9.4f scalar(ag_b1_empSE)
di as text "AUC:    mean=" %9.4f scalar(ag_auc_mean) "  bias=" %9.4f scalar(ag_auc_mean)-`auc_true' "  empSE=" %9.4f scalar(ag_auc_empSE)

quietly summarize ag_cover_b0, meanonly
di as text "Coverage beta0: " %5.1f 100*r(mean) "%"
quietly summarize ag_cover_b1, meanonly
di as text "Coverage beta1: " %5.1f 100*r(mean) "%"
quietly summarize ag_cover_auc, meanonly
di as text "Coverage AUC:   " %5.1f 100*r(mean) "%"

di as text _n "--- CF (strmcure, converged only) ---"
di as text "N=" scalar(cf_b0_N)
di as text "beta0:  mean=" %9.4f scalar(cf_b0_mean) "  bias=" %9.4f scalar(cf_b0_mean)-`b0_true' "  empSE=" %9.4f scalar(cf_b0_empSE)
di as text "beta1:  mean=" %9.4f scalar(cf_b1_mean) "  bias=" %9.4f scalar(cf_b1_mean)-`b1_true' "  empSE=" %9.4f scalar(cf_b1_empSE)
di as text "theta:  mean=" %9.4f scalar(theta_mean)  "  bias=" %9.4f scalar(theta_mean)-`th_true' "  empSE=" %9.4f scalar(theta_empSE)
di as text "AUC:    mean=" %9.4f scalar(cf_auc_mean) "  bias=" %9.4f scalar(cf_auc_mean)-`auc_true' "  empSE=" %9.4f scalar(cf_auc_empSE)

quietly summarize cf_cover_b0, meanonly
di as text "Coverage beta0:  " %5.1f 100*r(mean) "%"
quietly summarize cf_cover_b1, meanonly
di as text "Coverage beta1:  " %5.1f 100*r(mean) "%"
quietly summarize cf_cover_theta, meanonly
di as text "Coverage theta:  " %5.1f 100*r(mean) "%"
quietly summarize cf_cover_auc, meanonly
di as text "Coverage AUC:    " %5.1f 100*r(mean) "%"

****************************************************
* Save summary
****************************************************
* Grab coverage rates
quietly summarize fe_cover_b0, meanonly
local fe_cov_b0 = r(mean)
quietly summarize fe_cover_b1, meanonly
local fe_cov_b1 = r(mean)
quietly summarize fe_cover_auc, meanonly
local fe_cov_auc = r(mean)

quietly summarize ag_cover_b0, meanonly
local ag_cov_b0 = r(mean)
quietly summarize ag_cover_b1, meanonly
local ag_cov_b1 = r(mean)
quietly summarize ag_cover_auc, meanonly
local ag_cov_auc = r(mean)

quietly summarize cf_cover_b0, meanonly
local cf_cov_b0 = r(mean)
quietly summarize cf_cover_b1, meanonly
local cf_cov_b1 = r(mean)
quietly summarize cf_cover_auc, meanonly
local cf_cov_auc = r(mean)
quietly summarize cf_cover_theta, meanonly
local cf_cov_theta = r(mean)

* Mean runtimes and event counts
quietly summarize rt_sim_sec, meanonly
local rt_sim = r(mean)
quietly summarize rt_fe_sec, meanonly
local rt_fe = r(mean)
quietly summarize rt_ag_sec, meanonly
local rt_ag = r(mean)
quietly summarize rt_cf_sec, meanonly
local rt_cf = r(mean)
quietly summarize failures, meanonly
local mf = r(mean)
quietly summarize fe_n_events, meanonly
local mfe = r(mean)

preserve
    clear
    set obs 1

    gen beta0_true     = `b0_true'
    gen beta1_true     = `b1_true'
    gen theta_true     = `th_true'
    gen auc_true       = `auc_true'
    gen mean_total_ev  = `mf'
    gen mean_first_ev  = `mfe'

    * FE
    gen fe_N           = scalar(fe_b0_N)
    gen fe_b0_mean     = scalar(fe_b0_mean)
    gen fe_b0_bias     = scalar(fe_b0_mean) - `b0_true'
    gen fe_b0_empSE    = scalar(fe_b0_empSE)
    gen fe_b0_cover    = `fe_cov_b0'
    gen fe_b1_mean     = scalar(fe_b1_mean)
    gen fe_b1_bias     = scalar(fe_b1_mean) - `b1_true'
    gen fe_b1_empSE    = scalar(fe_b1_empSE)
    gen fe_b1_cover    = `fe_cov_b1'
    gen fe_auc_mean    = scalar(fe_auc_mean)
    gen fe_auc_bias    = scalar(fe_auc_mean) - `auc_true'
    gen fe_auc_empSE   = scalar(fe_auc_empSE)
    gen fe_auc_cover   = `fe_cov_auc'

    * AG
    gen ag_N           = scalar(ag_b0_N)
    gen ag_b0_mean     = scalar(ag_b0_mean)
    gen ag_b0_bias     = scalar(ag_b0_mean) - `b0_true'
    gen ag_b0_empSE    = scalar(ag_b0_empSE)
    gen ag_b0_cover    = `ag_cov_b0'
    gen ag_b1_mean     = scalar(ag_b1_mean)
    gen ag_b1_bias     = scalar(ag_b1_mean) - `b1_true'
    gen ag_b1_empSE    = scalar(ag_b1_empSE)
    gen ag_b1_cover    = `ag_cov_b1'
    gen ag_auc_mean    = scalar(ag_auc_mean)
    gen ag_auc_bias    = scalar(ag_auc_mean) - `auc_true'
    gen ag_auc_empSE   = scalar(ag_auc_empSE)
    gen ag_auc_cover   = `ag_cov_auc'

    * CF
    gen cf_N           = scalar(cf_b0_N)
    gen cf_b0_mean     = scalar(cf_b0_mean)
    gen cf_b0_bias     = scalar(cf_b0_mean) - `b0_true'
    gen cf_b0_empSE    = scalar(cf_b0_empSE)
    gen cf_b0_cover    = `cf_cov_b0'
    gen cf_b1_mean     = scalar(cf_b1_mean)
    gen cf_b1_bias     = scalar(cf_b1_mean) - `b1_true'
    gen cf_b1_empSE    = scalar(cf_b1_empSE)
    gen cf_b1_cover    = `cf_cov_b1'
    gen theta_mean     = scalar(theta_mean)
    gen theta_bias     = scalar(theta_mean) - `th_true'
    gen theta_empSE    = scalar(theta_empSE)
    gen theta_cover    = `cf_cov_theta'
    gen cf_auc_mean    = scalar(cf_auc_mean)
    gen cf_auc_bias    = scalar(cf_auc_mean) - `auc_true'
    gen cf_auc_empSE   = scalar(cf_auc_empSE)
    gen cf_auc_cover   = `cf_cov_auc'

    * Runtimes
    gen rt_sim_sec     = `rt_sim'
    gen rt_fe_sec      = `rt_fe'
    gen rt_ag_sec      = `rt_ag'
    gen rt_cf_sec      = `rt_cf'

    save "`outdir'\summary_stats.dta", replace
    export delimited using "`outdir'\summary_stats.csv", replace
restore

di as result _n "======================================================"
di as result "Finished."
di as result "Results: `outdir'\results_all.dta"
di as result "Summary: `outdir'\summary_stats.dta"
di as result "======================================================"

exit 0
