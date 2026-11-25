****************************************************
* rep_one_tt_save_data.do
*
* One simulation replicate — SAVES THE SIMULATED DATASET.
*
* Identical data-generating process and seeds to rep_one_tt_full.do.
* The only change is that after simulating the recurrent-event data,
* the raw counting-process dataset (id start end event Z) is saved to
*   <outdir>\simdata_<rep>.dta
* BEFORE any stsplit, so the same dataset can be re-analyzed in R with
* R's own tt()/splitting. This makes the Stata-vs-R comparison use
* identical data.
*
* theta_true is set below; run once with theta_true = 1 (into out_tt_full_1)
* and once with theta_true = 0.2 (into out_tt_full_0.2).
*
* Usage:
*   do rep_one_tt_save_data.do <rep_id> <seed> <outdir> <theta_true>
****************************************************

version 15.1
clear all
set more off
set rmsg on
set processors 1

*--------------------------
* Parse arguments
*--------------------------
args rep_id seed outdir theta_arg
if "`rep_id'"==""    local rep_id = 1
if "`seed'"==""      local seed   = 100001
if "`outdir'"==""    local outdir = "out_tt_full_save"
if "`theta_arg'"=="" local theta_arg = 0.2

cap mkdir "`outdir'"

capture log close
log using "`outdir'\log_save_`=string(`rep_id',"%03.0f")'.smcl", replace

di as text "Running replicate `rep_id' with seed `seed', theta=`theta_arg'"

*--------------------------
* User settings (IDENTICAL to rep_one_tt_full.do)
*--------------------------
local n           = 3400
local Tau         = 20
local beta0_true  = -0.57
local beta1_true  = 0.17
local theta_true  = `theta_arg'      // <-- set by argument: 1 or 0.2
local event_dep   = -1

local baseline_peak = 0.037
local decay         = 0.18

local tail_frac   = 0.05

set seed `seed'

tempfile simdata
tempname posth

****************************************************
* Simulation program (IDENTICAL to rep_one_tt_full.do)
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
sim_one_replicate using "`simdata'", n(`n') tau(`Tau') ///
    peak(`baseline_peak') decay(`decay') ///
    b0(`beta0_true') b1(`beta1_true') theta(`theta_true') ///
    eventdep(`event_dep') posth(`posth')

****************************************************
* Build event_count (number of PREVIOUS events per id), to match
* the R analysis dataset (event_count = 0 for first interval, etc.)
****************************************************
sort id start
by id: gen long event_count = sum(event) - event

****************************************************
* Label and save the raw counting-process dataset
****************************************************
gen int    rep        = `rep_id'
gen long   seed_used  = `seed'
gen double theta_true = `theta_true'
gen double beta0_true = `beta0_true'
gen double beta1_true = `beta1_true'

order rep seed_used id start end event event_count Z theta_true beta0_true beta1_true
label data "Simulated recurrent-event data, rep `rep_id', theta=`theta_true'"

save "`outdir'\simdata_`=string(`rep_id',"%03.0f")'.dta", replace
export delimited using "`outdir'\simdata_`=string(`rep_id',"%03.0f")'.csv", replace

di as text "Saved simdata_`=string(`rep_id',"%03.0f")'.dta (and .csv)"

file close _all
file open fd using "`outdir'\donesave_`=string(`rep_id',"%03.0f")'.txt", write replace text
file write fd "done" _n
file close fd

log close
exit 0
