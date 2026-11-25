****************************************************
* rep_one_tt_full.do
*
* One simulation replicate with declining baseline hazard.
* Split at unique event times (mimics R's tt()).
* Fit: FE, AG, CF (strmcure reps=0).
*
* Same DGP and seeds as sim_save_data.do / rep_one_timevar.do
* but using event-time splitting instead of weekly.
*
* Usage:
*   do rep_one_tt_full.do <rep_id> <seed> <outdir>
****************************************************

version 15.1
clear all
set more off
set rmsg on
set processors 1

*--------------------------
* Parse arguments
*--------------------------
args rep_id seed outdir
if "`rep_id'"=="" local rep_id = 1
if "`seed'"==""   local seed   = 100001
if "`outdir'"=="" local outdir = "C:\Users\e1139775\NUS Dropbox\Ziwei Zhao\STATA simulation tt splitting\out_tt_full_0.2"

cap mkdir "`outdir'"

capture log close
log using "`outdir'\log_`=string(`rep_id',"%03.0f")'.smcl", replace

di as text "Running replicate `rep_id' with seed `seed'"

*--------------------------
* User settings (IDENTICAL to sim_save_data.do)
*--------------------------
local n           = 3400
local Tau         = 20
local beta0_true  = -0.57
local beta1_true  = 0.17
local theta_true  = 0.2
local event_dep   = -1

local baseline_peak = 0.037
local decay         = 0.18

local tail_frac   = 0.05

set seed `seed'

tempfile simdata
tempname posth

****************************************************
* Simulation program (identical to sim_save_data.do)
****************************************************
capture program drop sim_one_replicate
program define sim_one_replicate
    syntax using/ , N(integer) Tau(real) Peak(real) Decay(real) ///
    B0(real) B1(real) Theta(real) Eventdep(real) Posth(name)

    tempname maxgamma
    scalar `maxgamma' = .

    preserve
        clear
        set obs `n'
        gen long id = _n
        gen byte Z  = (id <= `n'/2)
        gen double frailty = rgamma(1/`theta', `theta')
        quietly summarize frailty, meanonly
        scalar `maxgamma' = r(max)
        tempfile subj
        save `subj', replace
    restore

    scalar lambda_upper = 1.5 * `peak' * scalar(`maxgamma')

    postfile `posth' long id double start double end byte event byte Z ///
        using "`using'", replace

    use `subj', clear
    sort id

    quietly {
        forvalues i = 1/`n' {
            quietly summarize Z if id==`i', meanonly
            local Zi = r(mean)
            quietly summarize frailty if id==`i', meanonly
            local gi = r(mean)

            local Tstar = 0
            local prev  = 0
            local j     = 0

            while (`Tstar' < `tau') {
                local R = -ln(runiform()) / scalar(lambda_upper)
                local Tstar = `Tstar' + `R'
                if (`Tstar' >= `tau') {
                    continue, break
                }
                local v = runiform()
                local lam0 = `peak' * exp(-`decay' * `Tstar')
                local lamt  = `lam0' * exp((`b0' + `b1'*ln(`Tstar')) * `Zi')
                local inten = ((`j'+1)^(`eventdep')) * `gi' * `lamt'
                if (`inten' >= `v' * scalar(lambda_upper)) {
                    local j = `j' + 1
                    post `posth' (`i') (`prev') (`Tstar') (1) (`Zi')
                    local prev = `Tstar'
                }
            }
            post `posth' (`i') (`prev') (`tau') (0) (`Zi')
        }
    }

    postclose `posth'
    use "`using'", clear
end

****************************************************
* Simulate
****************************************************
timer clear 99
timer on 99
sim_one_replicate using "`simdata'", n(`n') tau(`Tau') ///
    peak(`baseline_peak') decay(`decay') ///
    b0(`beta0_true') b1(`beta1_true') theta(`theta_true') ///
    eventdep(`event_dep') posth(`posth')
timer off 99

****************************************************
* stset + event order
****************************************************
stset end, failure(event) id(id) time0(start) exit(time .)

stgen event_order = nfailures()
replace event_order = event_order + 1

count if _d==1
scalar nfail_total = r(N)

****************************************************
* Collect unique event times for splitting
****************************************************
preserve
    keep if _d == 1
    keep _t
    rename _t etime
    duplicates drop
    sort etime
    local cuts_et ""
    local ncuts = _N
    forvalues k = 1/`ncuts' {
        local cuts_et "`cuts_et' `=etime[`k']'"
    }
    di as text "Unique event times: `ncuts'"
restore

****************************************************
* Choose lastpool(J)
****************************************************
tempname J
scalar `J' = .

preserve
    keep if _d==1
    keep event_order
    contract event_order, freq(fails)
    quietly summarize fails, meanonly
    local totalfails = r(sum)
    gsort -event_order
    gen cum = sum(fails)
    gen frac = cum / `totalfails'
    quietly summarize event_order if frac >= `tail_frac', meanonly
    scalar `J' = r(max)
restore

****************************************************
* Truth
****************************************************
scalar auc_true = 1 - exp(`beta0_true') * (`Tau'^`beta1_true') / (`beta1_true' + 1)

****************************************************
* 1. FIRST-EVENT (event-time split)
****************************************************
scalar fe_b0hat    = .
scalar fe_b1hat    = .
scalar fe_auc      = .
scalar fe_n_events = .

timer clear 3
timer on 3

preserve
    sort id start
    by id: gen _cumev = sum(event)
    keep if _cumev == 0 | (_cumev == 1 & event == 1)
    drop _cumev

    count if event == 1
    scalar fe_n_events = r(N)

    stset end, failure(event) id(id) time0(start) exit(time .)

    stsplit _ts, at(`cuts_et')
    gen double ln_t   = ln(_t)
    gen double Z_ln_t = Z * ln_t

    capture noisily stcox Z Z_ln_t, nohr

    if (!_rc) {
        scalar fe_b0hat = _b[Z]
        scalar fe_b1hat = _b[Z_ln_t]
        if (fe_b0hat < . & fe_b1hat < .) {
            scalar fe_auc = 1 - exp(fe_b0hat) * (`Tau'^fe_b1hat) / (fe_b1hat + 1)
        }
    }
restore

timer off 3

****************************************************
* 2. AG MODEL (event-time split on full recurrent data)
****************************************************
scalar ag_b0hat = .
scalar ag_b1hat = .
scalar ag_auc   = .

timer clear 2
timer on 2

stsplit _tsplit, at(`cuts_et')
gen double ln_t   = ln(_t)
gen double Z_ln_t = Z * ln_t

stcox Z Z_ln_t, vce(cluster id) nohr

scalar ag_b0hat = _b[Z]
scalar ag_b1hat = _b[Z_ln_t]
if (ag_b0hat < . & ag_b1hat < .) {
    scalar ag_auc = 1 - exp(ag_b0hat) * (`Tau'^ag_b1hat) / (ag_b1hat + 1)
}

timer off 2

****************************************************
* 3. CF MODEL — strmcure reps(0), event-time split
****************************************************
scalar cf_failed = 0
scalar cf_rc     = .
scalar cf_b0hat  = .
scalar cf_b1hat  = .
scalar theta     = .
scalar cf_auc    = .

timer clear 1
timer on 1

capture noisily strmcure Z Z_ln_t, shared(id) strata(event_order) ///
    lastpool(`=scalar(`J')') reps(0) dots

scalar cf_rc = _rc

timer off 1

if (_rc) {
    di as error "CF model failed for rep `rep_id' with return code " _rc
    scalar cf_failed = 1
}
else {
    matrix b_cf = e(b)

    local pZ_cf   = colnumb(b_cf, "Z")
    local pZln_cf = colnumb(b_cf, "Z_ln_t")
    local pth_cf  = colnumb(b_cf, "_cons")

    if (`pZ_cf' < .) {
        scalar cf_b0hat = b_cf[1,`pZ_cf']
    }
    if (`pZln_cf' < .) {
        scalar cf_b1hat = b_cf[1,`pZln_cf']
    }
    if (`pth_cf' < .) {
        scalar theta = b_cf[1,`pth_cf']
    }
    if (`pZ_cf' < . & `pZln_cf' < .) {
        scalar cf_auc = 1 - exp(cf_b0hat) * (`Tau'^cf_b1hat) / (cf_b1hat + 1)
    }
}

****************************************************
* Runtimes
****************************************************
timer list 99
scalar runtime_sim = r(t99)
timer list 3
scalar runtime_fe = r(t3)
timer list 2
scalar runtime_ag = r(t2)
timer list 1
scalar runtime_cf = r(t1)

****************************************************
* Save one-row result
****************************************************
clear
set obs 1

gen rep          = `rep_id'
gen seed         = `seed'
gen failures     = scalar(nfail_total)
gen fe_n_events  = scalar(fe_n_events)
gen lastpoolJ    = scalar(`J')

gen beta0_true   = `beta0_true'
gen beta1_true   = `beta1_true'
gen theta_true   = `theta_true'
gen auc_true     = scalar(auc_true)

gen fe_b0hat     = scalar(fe_b0hat)
gen fe_b1hat     = scalar(fe_b1hat)
gen fe_auc       = scalar(fe_auc)

gen ag_b0hat     = scalar(ag_b0hat)
gen ag_b1hat     = scalar(ag_b1hat)
gen ag_auc       = scalar(ag_auc)

gen cf_failed    = scalar(cf_failed)
gen cf_rc        = scalar(cf_rc)
gen cf_b0hat     = scalar(cf_b0hat)
gen cf_b1hat     = scalar(cf_b1hat)
gen theta        = scalar(theta)
gen cf_auc       = scalar(cf_auc)

gen rt_sim_sec   = scalar(runtime_sim)
gen rt_fe_sec    = scalar(runtime_fe)
gen rt_ag_sec    = scalar(runtime_ag)
gen rt_cf_sec    = scalar(runtime_cf)

save "`outdir'\rep_`=string(`rep_id',"%03.0f")'.dta", replace

file close _all
file open fd using "`outdir'\done_`=string(`rep_id',"%03.0f")'.txt", write replace text
file write fd "done" _n
file close fd

di as text "Replicate `rep_id' complete."

log close
exit 0
