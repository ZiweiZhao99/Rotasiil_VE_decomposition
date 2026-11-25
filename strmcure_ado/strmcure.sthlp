{smcl}
{* 23Jan2013}{...}
{cmd:help strmcure}{right: ({browse "http://www.stata-journal.com/article.html?article=st0374":SJ15-1: st0374})}
{hline}

{title:Title}

{p2colset 5 17 19 2}{...}
{p2col :{hi:strmcure} {hline 2}}Frailty models and frailty-mixture
models for recurrent event times {p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 11 2}{cmd:strmcure} {it:varlist}
{ifin}{cmd:,}
{cmdab:sh:ared(}{it:varname}{cmd:)}
[{cmdab:st:rata(}{it:varname}{cmd:)}
{cmdab:xoff:set(}{it:varname}{cmd:)}
{cmdab:zlist:(}{it:varlist}{cmd:)}
{cmdab:zoff:set(}{it:varname}{cmd:)}
{opt znocons:tant}
{cmdab:link(logistic}|{cmd:probit}|{cmd:cloglog)}
{cmdab:lastp:ool(}{it:#}{cmd:)}
{cmdab:tailcm(zerotail}|{cmd:exponential}|{cmd:weibull)}
{cmdab:phi0(}{it:#}{cmd:)}
{cmdab:iter:ate(}{it:#}{cmd:)}
{cmdab:tol:erance(}{it:#}{cmd:)}
{opt log}
{cmdab:baseh0(}{it:varname}{cmd:)}
{opt saving}
{cmdab:seed(}{it:#}{cmd:)}
{cmdab:reps(}{it:#}{cmd:)}
{cmdab:repsp:art(}{it:#}{cmd:)}
{opt dots}
{cmdab:l:evel(}{it:#}{cmd:)}
{cmdab:obspart(}{it:#}{cmd:)}
{cmdab:rspart(}{it:#}{cmd:)}]

{p 4 4 2} You must {cmd:stset} your data before using {cmd:strmcure};
see {manhelp stset ST}.

{p 4 4 2} {it:varlist} may contain factor variables; see {help fvvarlist}.

{p 4 4 2} {cmd: bootstrap}, {cmd:by}, {cmd: jackknife}, {cmd:statsby},
{cmd: stepwise}, and {cmd: xi} are allowed; see {help prefix}.


{title:Description}

{p 4 4 2}
{cmd:strmcure} fits frailty models and frailty-mixture models (Xu et al.
2012) on {helpb st} data.  All models include a gamma-distributed
frailty term.  Users specify whether there is stratification for event
order or a cured fraction.  The parameter estimation is performed by the
expectation-maximization (EM) algorithm, and the variances are then
estimated by the Louis formula (Louis 1982) for the EM algorithm.


{title:Options}

{p 4 8 2}
{cmd:shared(}{it:varname}{cmd:)} specifies a variable in the dataset
that identifies the subjects.  It overrides the
{cmd:id(}{it:varname}{cmd:)} specification in {cmd:stset} if it has been
specified in {cmd:stset}.  {cmd:shared()} is required.

{p 4 8 2}
{cmd:strata(}{it:varname}{cmd:)} stratifies the recurrent event times
according to the specified event-order variable {it: varname}.
Observations with the same value belong to the same stratum, and the
baseline hazard function is unique to each stratum.  Not specifying
{cmd:strata()} means fitting a model without stratification when the
model does not involve a cured fraction.  This option must be specified
for the frailty-mixture model.

{p 4 8 2}
{cmd:xoffset(}{it:varname}{cmd:)} specifies an offset variable for
modeling the log(hazard).  The regression coefficient for the specified
variable {it:varname} is constrained to one.

{p 4 8 2}
{cmd:zlist(}{it:varlist}{cmd:)} specifies the list of variable names to
be included in modeling the probability of being noncured in the
frailty-mixture model.

{p 4 8 2}
{cmd:zoffset(}{it:varname}{cmd:)} specifies an offset variable for
modeling the probability of being noncured in the frailty-mixture model.
The regression coefficient is constrained to one.

{p 4 8 2}
{cmd:znoconstant} suppresses the constant term (intercept) in the
modeling of the probability of being noncured in the frailty-mixture
model.

{p 4 8 2}
{cmd:link(logistic}|{cmd:probit}|{cmd:cloglog)} specifies the link
function for modeling the probability of being noncured in the
frailty-mixture model.

{p 4 8 2}
{cmd:lastpool(}{it:#}{cmd:)} specifies an integer j such that the jth
and later events as specified by the {it:varname} in
{cmd:strata(}{it:varname}{cmd:)} will be pooled to form one stratum.

{p 4 8 2}
{cmd:tailcm(zerotail}|{cmd:exponential}|{cmd:weibull)} specifies the
tail-completion method for the frailty-mixture model.  The default is
{cmd:tailcm(zerotail)}.

{p 4 8 2}
{cmd:phi0(}{it:#}{cmd:)} specifies the starting value for the parameter
(that is, the variance of the gamma-distributed frailty).  The default
value is {cmd:phi0(2)}.  If a negative value is specified, then the
specification will be ignored and replaced with the default value of
{cmd:phi0(2)}.

{p 4 8 2}
{cmd:iterate(}{it:#}{cmd:)} specifies the maximum number of iterations
for the EM estimation.  The default is {cmd:iterate(5000)}.

{p 4 8 2}
{cmd:tolerance(}{it:#}{cmd:)} specifies the minimum tolerance level for
the regression coefficient estimation.  The default is
{cmd:tolerance(10^(-5))}.  The {cmd:tolerance(}{it:#}{cmd:)} works
together with {cmd:iterate(}{it:#}{cmd:)}: the EM estimation procedure
is regarded as converged if either criterion is met.

{p 4 8 2}
{cmd:log} requests that the parameter estimates from each iteration of
the EM estimation procedure be displayed.

{p 4 8 2}
{cmd:baseh0(}{it:varname}{cmd:)} requests that the estimated baseline
hazard (using the Nelson-Aalen estimator) be saved.  The variable
{it:varname} is replaced if it already exists.

{p 4 8 2}
{cmd:saving} requests that the iterative parameter estimates from the EM
estimation procedure be saved.  The saved result is in the matrix
{cmd:e(EMiter_par)}.  If the column or row of the matrix
{cmd:e(EMiter_par)} exceeds the matsize limit, then the user may not be
able to read the matrix into Stata using the command {cmd:matlist
e(EMiter_par)}. In this case, the user can read the matrix into Mata or
by converting the matrix to variables in a Stata dataset.

{p 4 8 2}
{cmd:seed(}{it:#}{cmd:)} sets the random-number seed when initiating the
process of variance estimation using Louis's formula.  The default is
{cmd:seed(0)}.

{p 4 8 2}
{cmd:reps(}{it:#}{cmd:)} specifies the number of realizations for the
Monte Carlo simulation when estimating the variance using Louis's
formula.  The specified number should be a multiple of 10.  If not, it
will be rounded up.  The default is {cmd:reps(1000)}.  A specification
of {cmd:reps(0)} indicates no estimation of variance using Louis's
formula.  Users can specify the {cmd:reps(0)} option if they use the
{cmd:bootstrap} or {cmd:jackknife} methods for variance estimation.

{p 4 8 2}
{cmd:repspart(}{it:#}{cmd:)} determines how many blocks to split the
Monte Carlo replications into for the variance estimation using Louis's
formula.  For each block, K=(reps/repspart) realizations of the latent
variables are generated.  The default is {cmd:repspart(10)}.  A larger K
requires less computation time but a larger matrix size.

{p 4 8 2}
{cmd:dots} displays the replication dots in the estimation of variance
using Louis's formula.  One dot is displayed for each successful
completion with one block of K Monte Carlo realizations.

{p 4 8 2}
{cmd:level(}{it:#}{cmd:)} specifies the confidence level for estimating
the confidence interval.  The default is {cmd:level(95)}.

{p 4 8 2}
{cmd:obspart(}{it:#}{cmd:)} specifies the number of partitions for the
total number of observations.  This is used for partitioning the at-risk
set in variance estimation using Louis's formula.  The default is
{cmd:obspart(N1)}, where {cmd:N1} is (total number of observations)/500,
rounded up to the next integer.  A smaller number requires less
computation time but a larger matrix size.

{p 4 8 2}
{cmd:rspart(}{it:#}{cmd:)} specifies the number of partitions for the
total number of events in variance estimation.  The default is
{cmd:rspart(N2)}, where {cmd:N2} is (total number of events)/500, rounded
up to the next integer.  A smaller number requires less computation time
but a larger matrix size.


{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. use data_rhdnase}

{phang}{cmd:{ul:Example of fitting a frailty model without stratification}}{p_end}

{pstd}Declare data to be survival-time data using the counting-process
time scale.  The variables {cmd:id} and {cmd:status} are the subject
identifier and event indicator, respectively.  The two variables
{cmd:time1} and {cmd:time2} together define the observation time
interval of {cmd:(time1, time2]} for each data record.  The variable
{cmd:trt} is a binary variable indicating the treatment allocation (for
example, 1 for intervention, 0 for placebo).{p_end}

{phang2}{cmd:. stset time2, fail(status) id(id) time0(time1) exit(time .)}{p_end}

{pstd}Fit a frailty model  {p_end}
{phang2}{cmd:. strmcure trt, shared(id) dots}

{phang}{cmd:{ul:Example of fitting a frailty model, with stratification for event order}}{p_end}

{pstd}Generate a variable corresponding to the event order {p_end}
{phang2}{cmd:. stgen event_order=nfailures() }{p_end}
{phang2}{cmd:. replace event_order=event_order+1}{p_end}
 
{pstd}Fit a frailty model with stratification for event order (note: the
fourth and fifth events are pooled in the same stratum) {p_end}
{phang2}{cmd:. strmcure trt, shared(id) strata(event_order) lastpool(4) dots }{p_end}

{phang}{cmd:{ul:Example of fitting a frailty-mixture model, with stratification for event order}}{p_end}

{pstd}Declare data to be survival-time data using the gap time {p_end}
{phang2}{cmd:. gen time=time2-time1}{p_end}
{phang2}{cmd:. stset time, fail(status) exit(time .)}{p_end}
{phang2}{cmd:. strmcure trt, shared(id) strata(event_order) lastpool(4) zlist(trt) link(logistic) tailcm(ze) dots }{p_end}
{phang2}{cmd:. strmcure trt fev, shared(id) strata(event_order) lastpool(4) zlist(trt fev) link(logistic) tailcm(ze) dots }{p_end}


{title:Stored results}

{pstd}
{cmd:strmcure} stores the following in {cmd:e()}:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_sub)}}number of subjects{p_end}
{synopt:{cmd:e(N_fail)}}number of failures{p_end}
{synopt:{cmd:e(N_g)}}number of groups{p_end}
{synopt:{cmd:e(g_min)}}smallest group size{p_end}
{synopt:{cmd:e(g_max)}}largest group size{p_end}
{synopt:{cmd:e(g_avg)}}average group size{p_end}
{synopt:{cmd:e(total_time)}}total analysis time{p_end}
{synopt:{cmd:e(CurePresent)}}whether a cured fraction is specified (0 no and 1 yes){p_end}
{synopt:{cmd:e(phi_initial)}}initial value for the parameter phi (that is, the variance of the gamma-distributed frailty) {p_end}
{synopt:{cmd:e(tolv_cvg)}}the square root of the distance between the two parameter estimate vectors at the last iteration and its preceding iteration of the EM estimation procedure{p_end}
{synopt:{cmd:e(iter_cvg)}}the total number of iterations at the convergence of the EM estimation procedure{p_end}
{synopt:{cmd:e(Tail)}}code number of the tail-completion method for the frailty-mixture model{p_end}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:strmcure}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(title)}}frailty model or frailty-mixture model{p_end}
{synopt:{cmd:e(offsetX)}}name of the offset variable for modeling the log(hazard){p_end}
{synopt:{cmd:e(offsetZ)}}name of the offset variable for modeling the probability of being noncured for the frailty-mixture model{p_end}
{synopt:{cmd:e(Tailname)}}tail-completion method for the frailty-mixture model{p_end}
{synopt:{cmd:e(link)}}name of the link function for the frailty-mixture model{p_end}
{synopt:{cmd:e(ties)}}method used for handling ties (Breslow method is always used){p_end}
{synopt:{cmd:e(shared)}}frailty grouping variable{p_end}
{synopt:{cmd:e(t0)}}{cmd:_t0}{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synopt:{cmd:e(strata)}}strata variable{p_end}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}
{synopt:{cmd:e(EMiter_par)}}matrix storing the history of the parameter estimates during the iterative EM estimation procedure, if {cmd:saving} option was specified in the {cmd:strmcure} command{p_end}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}


{title:References}

{p 4 8 2} Louis, T. A.  1982.  Finding the observed information matrix
when using the EM algorithm. {it: Journal of the Royal Statistical Society, Series B} 44: 226-233.

{p 4 8 2} Xu, Y., Y. B. Cheung, K. F. Lam, and P. Milligan.  2012.
Estimation of summary protective efficacy using a frailty mixture model
for recurrent event time data. {it:Statistics in Medicine} 31:
4023-4039.


{title:Authors}

{pstd}Ying Xu{p_end}
{pstd}Center for Quantitative Medicine{p_end}
{pstd}Office of Clinical Sciences{p_end}
{pstd}Duke-NUS Graduate Medical School{p_end}
{pstd}Singapore{p_end}
{pstd}and Department of Biostatistics{p_end}
{pstd}Singapore Clinical Research Institute{p_end}
{pstd}Singapore{p_end}
{pstd}{browse "mailto:tina.xuying@duke-nus.edu.sg":tina.xuying@duke-nus.edu.sg}

{pstd}Yin Bun Cheung{p_end}
{pstd}Center for Quantitative Medicine{p_end}
{pstd}Office of Clinical Sciences{p_end}
{pstd}Duke-NUS Graduate Medical School{p_end}
{pstd}Singapore{p_end}
{pstd}and Department of International Health{p_end}
{pstd}University of Tampere{p_end}
{pstd}Finland{p_end}
{pstd}{browse "mailto:yinbun.cheung@duke-nus.edu.sg":yinbun.cheung@duke-nus.edu.sg}


{title:Also see}

{p 4 14 2}Article:  {it:Stata Journal}, volume 15, number 1: {browse "http://www.stata-journal.com/article.html?article=st0374":st0374}

{p 7 14 2}
Help:  {manhelp st ST}, {manhelp stset ST}
{p_end}
