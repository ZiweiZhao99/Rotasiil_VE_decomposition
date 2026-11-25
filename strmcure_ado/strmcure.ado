*! Version 3.0 16-Oct-2014
capture program drop strmcure
program strmcure, eclass sortpreserve byable(onecall) prop(st)
version 11.0
	if _by() {
		local by "by `_byvars'`_byrc0':"
	}
	if replay() {
	      if "`e(cmd)'" != "strmcure" exit 301
		if _by()    error 190
	      Replay `0'
		exit
	}
	`by' Estimate `0'
end


capture program drop Replay
program Replay
	syntax [, Level(cilevel)]
	ereturn display, level(`level') 
end


capture program drop Estimate
program Estimate, eclass byable(recall)
  st_is 2 analysis 
  syntax  varlist(numeric fv) [if] [in],  SHared(varname)        ///     specify the sequence of events occured in one subject
    [                                                             /// 
	STrata(varname numeric)                             											///
    XOFFset(varname numeric)                                                                    ///      specify offset term in modeling the event hazard
    ZLIST(varlist numeric fv)    ZOFFset(varname numeric)   ZNOCONStant    LINK(string)            ///     sepcify covariates and offset term in modeling the susceptability probability  
    ITERate(numlist max=1)  TOLerance(numlist max=1)              ///     specify the link function for the probability of being susceptible
    LOG       DOTS     LASTPool(numlist max=1 integer>0)          ///    maxmum event stratum. Strata >=lastpool will be pooled together assuming same baseline hazard functions
    TAILCM(string)              ///           "zero" "weibull"  "exponential"  
    SEED(numlist max=1)         ///          seed used to generate the random variates for K and W during the variance estimation using Louis's formula
    SAVING                      ///             specify to save all the parameter estimates until EM convergence
    RSPART(numlist max=1 )      ///         partition risk set (in wake of computing resournces and matsize limitations)
    OBSPART(numlist max=1)      ///           partition subjects (in wake of computing resournces and matsize limitations)
    REPS(numlist max=1)         ///            replications for simulating unobserved frailty terms K and W
    REPSPart(numlist max=1)     ///            replications be partitioned
    BASEH0(string) Level(cilevel) ///
    PHI0(numlist max=1)           ///         specify the initial value for phi-variance of the gamma frailty term
	]   
   			
  						
qui{
   cap ereturn clear
   **** display and maximization options ****
   marksample touse
   qui replace `touse' = 0 if _st==0

   *mlopts mlopts, `options'
   local display_options level(`level') 

   **** check subject identifier correctly specified or not ****
   if("`shared'" != "") {
             local idcluster `shared'
   }
   else {
   }                     

   **** Specify analysis set ****
   markout `touse' `strata' `varlist' `zlist' `zoffset' `xoffset' 
   markout `touse' `idcluster', strok 

   **** check empty data set ****
   qui count if `touse'
   if(r(N)==0) {
			noi: di as error "no observations"
			exit
   }

   **** check offset options for X and Z ****
   local xoffopt
   local zoffopt
   local constant
   if("`xoffset'" != "")   	local xoffopt "offset(`xoffset')"  
   if("`zoffset'" != "")   	local zoffopt "offset(`zoffset')"  
   if("`znoconstant'" != "")   local constant noconstant   

   *** check strata() and lastpool() options ***
   if("`strata'" != "") {
       	qui summ `strata', meanonly
            if "`lastpool'" == ""  {
 				scalar maxTE=r(max)
		}
		else {
                        scalar maxTE=min(`lastpool', r(max))
                        if `lastpool' < r(max)   noi: dis as text "strata " maxTE " to " r(max) " are pooled together for hazard estimation"
            } 
   }
   else {
        scalar maxTE=1
        if "`lastpool'" != ""  noi: dis as text "No stratification: lastpool() option ignored"   
   }  
 
  ****check covariates for link function***
     local z_index
     local z_name
     local x_index
     local x_name
     if("`zlist'" != "") {
		     _rmcoll `zlist' if `touse', expand
             local z_index `r(varlist)'            
             foreach vari of local z_index {
                local z_name `z_name' Theta:`vari'
             }
     }

     if("`znoconstant'" == "")  local z_name `z_name' Theta:_cons  
        _rmcoll `varlist' if `touse', expand

     local x_index `r(varlist)' 
     foreach vari of local x_index {
          local x_name `x_name' Beta:`vari'
     }

 ***check link function****
  local linkarg
  scalar cureI=1
  tokenize `link'
  local temp0 : word count `link'
  if (`temp0' == 1) {
     	  	if substr(lower("`link'"),1,4) == "prob" {
            		 local linkarg = "probit"
			}
			else if substr(lower("`link'"),1,4) == "clog" {
					local linkarg ="cloglog"
			}
			else if substr(lower("`link'"),1,4) == "logi" {
					local linkarg = "logit"
			}
			else {
					noi: di as error "option link() incorrectly specified" 
					noi: di as error "recognized: (logistic | complementary log-log | probit)"
					error 198
					exit
			}
   }
   else if (`temp0' == 0) {
			if("`zlist'" != "" | "`znoconstant'" != ""  | "`zoffset'" !="") {
						noi: di as error "link function must be specified"
						error 197
						exit
			}
			else {
                   	scalar cureI=0
					local z_index 
					local z_name  
                         
			}
   }
   else if (`temp0' > 1) {
       			noi: di as error "too many specifications for option link()"
       			noi: di as error "recognized: (logistic | complementary log-log | probit)"
       			noi: error 198
 				exit
   }
   else {
                noi: error 199
				exit
   }

   if ("`linkarg'" != "" & "`zlist'" == "" & "`znoconstant'" != "")  {
       			noi: di as error "too few variables specified for the link function"
 				exit
   }

 
 ***check whether strata() option is provided if a cure model is specified  ****
   if(cureI==1 & "`strata'" == "")  {
       			noi: di as error "strata() option must be specified for estimation of frailty mixture models"
 				exit       
   }
 
 
 ***check tail completion method (default: tailmet=1, zero-tail completion method)*
	tokenize `tailcm'
  	local temp0 : word count `tailcm'
 	scalar tailmet=1
    local tailname "zero-tail"
 	if(`temp0' == 1) {
     	  	if substr(lower("`tailcm'"),1,2) == "ze" {
            	 scalar tailmet = 1
      		 	 local tailname "zero-tail"
       	}
      	else if substr(lower("`tailcm'"),1,2) == "ex" {
      		 	 scalar tailmet = 2
      		 	 local tailname "exponential"
       	}
      	else if substr(lower("`tailcm'"),1,2) == "we" {
      		 	 scalar tailmet = 3
      		 	 local tailname "weibull"
  		}
  		else {
       			 noi: di as error "option tailcm() incorrectly specified" 
       	 		 noi: di as error "recognized: (zerotail | exponential | weibull)"
       	 		 noi: error 198
				 exit
  	  	}
   	}
   	else if ( `temp0' > 1) {
       			noi: di as error "too many specifications for tailcm() option"	 
				noi: di as error "recognized: (zerotail | exponential | weibull)"
       	 		noi: error 198
				exit
   	}

    if(cureI==0) {
 		scalar tailmet=0
      	local tailname
    }

 ****verifying other parameter settings, set to default options if needed
 local iterEM=5000
 if("`iterate'" != "")  local iterEM=`iterate'

 local tolvEM=10^(-5)
 if("`tolerance'" != "") local tolvEM=`tolerance'
 
 if("`reps'" == "") local reps=1000
 local iseed=0
 if("`seed'" != "") local iseed=`seed'
 
 local phi_initial=2
 if("`phi0'" != "") local phi_initial=`phi0'
 if(`phi_initial'<0) local phi_initial=2
 
 noi: st_show
 noi: di _n
****************************************************************************************
************   		parameter estimation via EM algorithm
****************************************************************************************
	tempname temp new_est old_est aa result_save

      local drop_var_list "subject grp_count Di temp0 old_ev event_order first_id start_1 end_1 ind_Constr"
      foreach drop_var of local drop_var_list {
			    cap drop ``drop_var''
				tempvar `drop_var'
	}

	egen `subject'=group(`idcluster')     if `touse'  
	bysort `subject': egen `Di'=sum(_d)   if `touse' 
    bysort `subject': gen `grp_count'=_N  if `touse'
	
	****
	****specify event order according to sequential order of strata(variable)****
	****
	if ("`strata'" != "") {
			clonevar `temp0'=`strata'  				if  `touse'
			replace `temp0'=maxTE  					if `strata'>=maxTE  &  `touse'
			egen `old_ev'=group(`strata')           if `touse' 
			egen `event_order'=group(`temp0')       if `touse' 
	}
	else {
			gen byte `old_ev'=1    		            if `touse' 
			gen byte `event_order'=1                if `touse' 
      }
	sort `subject' `old_ev' _t0
	by `subject': gen byte `first_id'=cond(_n==1,1,0)  if `touse' 

    ****
	****check data, if all subjects had first events occurred, then no cured fraction in the population
	****
	sum _d if `first_id'==1
	if (r(min)==1 & cureI==1) {
  		noi: di as error _newline(1) "All subjects have failures: mixture cure model can not be fitted." 
		exit
	}
	if (r(min)==r(max))  {
  		noi: di as error _newline(1) "Outcome does not vary." 
		exit
	}
	
	sort `subject' `event_order'  _t0  
	by `subject': egen `start_1'=min(_t0) if  `event_order'==1  &  `touse'
	by `subject': egen `end_1'=max(_t)    if  `event_order'==1   &  `touse'
	by `subject': replace `start_1'=`start_1'[1] if  _n~=1  &  `touse'
	by `subject': replace `end_1'=`end_1'[1]     if  _n~=1  &  `touse'

	sum _t if _d==1 & `event_order'==1 & `touse'
	scalar maxT=r(max)

	if (cureI==1) {
            gen byte `ind_Constr'=cond(`event_order'==1 & `end_1'>maxT  & `Di'==0, 1, 0)  if `touse'  
	}
	else {
            gen byte `ind_Constr'=0  if `touse'  
	}

	************
	summ `subject'  if `touse'  
	local N_sub=r(max)

	summ `grp_count' if `first_id'==1
	local N_g=r(N)
	local g_min=r(min)
	local g_avg=r(mean)
	local g_max=r(max)
	
	summ _d if `touse'
	local N_fail=r(sum)
	summ _t0 if `touse'
	local time_ar0=r(sum)
 
 	summ _t if `touse'
	local time_ar1=r(sum)
	local time_ar=`time_ar1'-`time_ar0'
	
 scalar iter=0
 scalar tolv=`tolvEM'+1
 local maxTE=maxTE
 
 ******************************************************
 ****************zero-tail constraints or no constraint
 ******************************************************
 if(tailmet==0 | tailmet==1) {     
	local drop_var_list "temp_pi temp_expB temp_marg temp_marg_in"
	foreach drop_var of local drop_var_list {
					   cap drop ``drop_var''
					   tempvar `drop_var'
	}

	******************
	if (cureI==1) {
			 `linkarg' _d  `z_index' if `first_id'==1, `zoffopt'  `constant'
			 predict `temp_pi' if `touse', pr 
			 mat `new_est'=e(b)
	}
	else {
		 	 gen byte `temp_pi'=1  if `touse'
	}

	******************
    stcox `x_index' if  `ind_Constr'==0, strata(`event_order') mgale(`temp_marg_in') `xoffopt'  
	predict `temp_expB'  if  `ind_Constr'==0, hr nooffset 
	if(cureI==1) {
              mat `temp'=`new_est', e(b), `phi_initial'
   		 	  mat `new_est'=`temp'
	}
	else {
  		  	  mat `new_est'=e(b), `phi_initial'
	}      
	gen `temp_marg'=_d-`temp_marg_in'        if  `ind_Constr'==0 

  	******************
  	while (tolv>`tolvEM'  &  iter<`iterEM')  {
     	 scalar iter=iter+1
     	 mat `old_est'=`new_est'
     	 scalar Phi=`old_est'[1,colsof(`old_est')] 
     	 scalar IPhi=1/Phi
     	 
		 local drop_var_list "temp_cumH est_q est_p temp_Ai temp_Ci temp_offset temp_diff_g12"
     	 foreach drop_var of local drop_var_list {
			cap drop ``drop_var''
			tempvar `drop_var'
     	 }
     	 bysort `subject': egen `temp_cumH'=sum(`temp_marg')    if `touse'

     	 if(cureI==1) {
     			gen `est_q'=cond(`Di'==0, `temp_pi'/((1+Phi*`temp_cumH')^IPhi), 1)               if `touse'
     			gen `est_p'=cond(`Di'==0, `est_q'/(`est_q'+1-`temp_pi')*(1-`ind_Constr'), 1)     if `touse'
     	 }
     	 else{
     			gen byte `est_q'=1    if `touse'
     			gen byte `est_p'=1    if `touse'
		 }
     	 gen `temp_Ai'=IPhi+`Di'         if `ind_Constr'==0 
     	 gen `temp_Ci'=IPhi+`temp_cumH'  if `ind_Constr'==0 
     	 gen `temp_diff_g12'=digamma(`temp_Ai')-ln(`temp_Ci')-`temp_Ai'/`temp_Ci'   if `ind_Constr'==0 

     	 if ("`xoffset'" != "") {
            		gen `temp_offset'=ln(`est_p'*`temp_Ai'/`temp_Ci')+`xoffset'   if `ind_Constr'==0    
     	 }
     	 else {
				    gen `temp_offset'=ln(`est_p'*`temp_Ai'/`temp_Ci')             if `ind_Constr'==0 
     	 }
     	 ******************
     	 local drop_var_list "temp_pi temp_expB temp_marg temp_marg_in"
     	 foreach drop_var of local drop_var_list {
							  cap drop ``drop_var''
							  tempvar `drop_var'
     	 }

         ******************        
     	 if(cureI == 1) {
                 cap drop `temp_zt'
                 tempvar temp_zt 
                 if("`constant'" == "") {                   
                    	sum `est_p' if `first_id'==1, meanonly
                    	if("`linkarg'" == "logit") {
			     		    local b0=logit(r(mean))                     
                        }
                   	    else if ("`linkarg'" == "probit") {
			     		    local b0=invnormal(r(mean))                  
                    	}
                    	else if ("`linkarg'" == "cloglog") {
			     		    local b0=ln(-ln(1-r(mean)))  
                    	}  
                    	else {
                    	}
                        local initopt init(_cons=`b0') search(quietly)
                 }

				ml model d2 strmcure_my`linkarg'_d2 (mu: `est_p' = `z_index' , `zoffopt'  `constant')   if `first_id'==1, `initopt'  maximize   
				predict `temp_zt' if `touse', xb
				mat `new_est'=e(b)

				if("`linkarg'" == "logit") {
					gen `temp_pi'=invlogit(`temp_zt')      if `touse'                     
				}
				else if ("`linkarg'" == "probit") {
					gen `temp_pi'=normal(`temp_zt')        if `touse'                      
				}
				else if ("`linkarg'" == "cloglog") {
					gen `temp_pi'=1-exp(-exp(`temp_zt'))   if `touse'     
				}  
				else {
				}                 
     	 }
     	 else {
			    gen byte `temp_pi'=1     if `touse'
     	 }

     	 ******************
    	 stcox `x_index' if `ind_Constr'==0, offset(`temp_offset') strata(`event_order') mgale(`temp_marg_in')  
     	 predict `temp_expB'  if `ind_Constr'==0, hr nooffset 
     	 if(cureI==1) {
			  mat `temp'=`new_est', e(b)
			  mat `new_est'=`temp'
     	 }
     	 else {
			  mat `new_est'=e(b)
     	 }	

         if "`xoffset'" != ""  {
     			  gen `temp_marg'=(_d-`temp_marg_in')/exp(`temp_offset')*exp(`xoffset')   if `ind_Constr'==0 
         }
         else {
     			  gen `temp_marg'=(_d-`temp_marg_in')/exp(`temp_offset')  		          if `ind_Constr'==0 
         }

 		******************   
    	ml model d2 strmcure_myGamma_d2 (phi: `temp_diff_g12' `est_p'=)  if `first_id'==1 & `ind_Constr'==0, search(on) maximize 
		mat `temp'=`new_est', [phi]_b[_cons]
     	mat `new_est'=`temp'

 		******************
  		matrix `aa'=(`new_est'-`old_est')*(`new_est'-`old_est')'
  		scalar tolv=sqrt(`aa'[1,1])

  		matrix `aa'=iter, `new_est'
 		if(iter==1)  {
						mata: result_savem=st_matrix("`aa'")
		}
  		else {
						mata:result_savem=result_savem\st_matrix("`aa'")
		}
        matname `aa' iter `z_name' `x_name' frailty:_cons, columns(.) explicit
        ******************
      	if ("`log'" !="") {
              	if (iter==1)  {
 			       noi: display "{txt}{hline 60}"
                  	 noi: di as text "......EM Iterative Parameter Estimation Process Starts......"
                  	 noi: matlist `aa', names(columns)
              	}
		  	    else {
				     noi: matlist `aa', nonames
                }
      	}

  	}
  	******************
   	if("`log'" != "") {
              noi: di as text "......EM Iterative Parameter Estimation Process Ends......"  
 	          noi: display "{txt}{hline 60}"
    }
 
 }
 else if (tailmet==2) {
 ******************************************************
 ****************exponential tail completion
 ******************************************************
     local drop_var_list "temp_pi temp_expB temp_marg temp_marg_in temp_H0"
     foreach drop_var of local drop_var_list {
	                   cap drop ``drop_var''
					   tempvar `drop_var'
     }

     ******************
     `linkarg' _d  `z_index' if `first_id'==1, `zoffopt'  `constant'
     predict `temp_pi' if `touse', pr 
     mat `new_est'=e(b)
     	
     ******************
     stcox `x_index' if `touse', strata(`event_order') mgale(`temp_marg_in') basec(`temp_H0') `xoffopt'  
	 predict `temp_expB' if `touse', hr nooffset 
     mat `temp'=`new_est', e(b), `phi_initial'
   	 mat `new_est'=`temp'
     gen `temp_marg'=_d-`temp_marg_in'              if `touse'
	 
		
    ******************
  	while (tolv>`tolvEM'  &  iter<`iterEM')  {
    		scalar iter=iter+1
     		mat `old_est'=`new_est'
     		scalar Phi=`old_est'[1,colsof(`old_est')] 
     		scalar IPhi=1/Phi

     		local drop_var_list "temp_cumH est_q est_p temp_Ai temp_Ci temp_offset temp_diff_g12"
     		foreach drop_var of local drop_var_list {
					cap drop ``drop_var''
					tempvar `drop_var'
     		}

     		bysort `subject': egen `temp_cumH'=sum(`temp_marg')  if `touse'  
			
     		gen `est_q'=cond(`Di'==0, `temp_pi'/((1+Phi*`temp_cumH')^IPhi), 1)   if `touse'
     		gen `est_p'=cond(`Di'==0, `est_q'/(`est_q'+1-`temp_pi'), 1)          if `touse'
			
    		gen `temp_Ai'=IPhi+`Di'             if `touse'
     		gen `temp_Ci'=IPhi+`temp_cumH'      if `touse'
     		gen `temp_diff_g12'=digamma(`temp_Ai')-ln(`temp_Ci')-`temp_Ai'/`temp_Ci'   if `touse'

            if ("`xoffset'" != "") {
            		gen `temp_offset'=ln(`est_p'*`temp_Ai'/`temp_Ci')+`xoffset'   if `touse'
		    }
            else {
				    gen `temp_offset'=ln(`est_p'*`temp_Ai'/`temp_Ci')             if `touse'
            }

     		******************
     		local drop_var_list "temp_pi temp_expB temp_marg temp_marg_in temp_H0"
     		foreach drop_var of local drop_var_list {
							  cap drop ``drop_var''
							  tempvar `drop_var'
	        }

     		******************
     		cap drop `temp_zt'
     		tempvar temp_zt 
     		if("`constant'" == "") {                   
                    	sum `est_p' if `first_id'==1, meanonly
                    	if("`linkarg'" == "logit") {
			     		    local b0=logit(r(mean))                     
                        }
                   	    else if ("`linkarg'" == "probit") {
			     		    local b0=invnormal(r(mean))                  
                    	}
                    	else if ("`linkarg'" == "cloglog") {
			     		    local b0=ln(-ln(1-r(mean)))  
                    	}  
                    	else {
                    	}
                        local initopt init(_cons=`b0') search(quietly)
     		}
     		ml model d2 strmcure_my`linkarg'_d2 (mu: `est_p' = `z_index' , `zoffopt'  `constant')   if `first_id'==1, `initopt'  maximize  
     		predict `temp_zt' if `touse', xb
     		mat `new_est'=e(b)

     		if("`linkarg'" == "logit") {
					gen `temp_pi'=invlogit(`temp_zt')      if `touse'                     
     		}
     		else if ("`linkarg'" == "probit") {
					gen `temp_pi'=normal(`temp_zt')        if `touse'                      
     		}
     		else if ("`linkarg'" == "cloglog") {
					gen `temp_pi'=1-exp(-exp(`temp_zt'))   if `touse'     
     		}  
     		else {
     		}                            

        ******************
    	stcox `x_index' if `touse', offset(`temp_offset') strata(`event_order') mgale(`temp_marg_in')  basec(`temp_H0') 
     	predict `temp_expB' if `touse', hr nooffset

        mat `temp'=`new_est', e(b)
   		mat `new_est'=`temp'

		gen `temp_marg'=(_d-`temp_marg_in')/exp(`temp_offset')       if `touse'
        
	    forv i=1/1 {
			sum `temp_H0' if `event_order'==`i' & _d==1
			local Etail_H=r(max)
			sum _t if `event_order'==`i' & _d==1
			local maxTi=r(max)
			replace `temp_marg'=`temp_marg'+(`Etail_H'/`maxTi')*(_t-`maxTi')*`temp_expB'    if `event_order'==`i'  &  _d==0  &  _t>=`maxTi'  & _t0<`maxTi'
		    replace `temp_marg'=(`Etail_H'/`maxTi')*(_t-_t0)*`temp_expB'        if `event_order'==`i'  &  _d==0  &  _t>=`maxTi'  & _t0>`maxTi'	
	    }

	    if "`xoffset'" != "" {
     		  	replace `temp_marg'=`temp_marg'*exp(`xoffset')    	 if `touse'	
        }	
      	******************
        ml model d2 strmcure_myGamma_d2 (phi: `temp_diff_g12' `est_p'=)  if `first_id'==1,  search(on) maximize 

     	mat `temp'=`new_est', [phi]_b[_cons]
     	mat `new_est'=`temp'

        ******************
  		matrix `aa'=(`new_est'-`old_est')*(`new_est'-`old_est')'
  		scalar tolv=sqrt(`aa'[1,1])

  		matrix `aa'=iter, `new_est'
 		if(iter==1)  {
						mata: result_savem=st_matrix("`aa'")
		}
  		else {
						mata: result_savem=result_savem\st_matrix("`aa'")
		}
 
        matname `aa' iter `z_name' `x_name' frailty:_cons, columns(.) explicit
   
        ******************
      	if "`log'" !="" {
              	if iter==1  {
 			         noi: display "{txt}{hline 60}"
                  	 noi: di as text "......EM Iterative Parameter Estimation Process Starts......"
                  	 noi: matlist `aa', names(columns)
              	}
		  	    else {
				     noi: matlist `aa', nonames
		  	    }
      	}

  	}
  	******************
   	if("`log'" != "") {
              noi: di as text "......EM Iterative Parameter Estimation Process Ends......"  
 	          noi: display "{txt}{hline 60}"
   	}

 }
 else if (tailmet==3) {
 ******************************************************
 ****************weibull tail completion
 ******************************************************

	local drop_var_list "temp_pi temp_expB temp_marg temp_marg_in temp_H0 temp_wbl_t temp_wbl_t0 temp_wbl_coef"
	foreach drop_var of local drop_var_list {
	                           cap drop ``drop_var''
							   tempvar `drop_var'
	}

	******************
	`linkarg' _d `z_index' if `first_id'==1, `zoffopt'  `constant'
	predict `temp_pi' if `touse', pr 
	mat `new_est'=e(b)

	******************
	stcox `x_index' if `touse', strata(`event_order') mgale(`temp_marg_in') basec(`temp_H0') `xoffopt' 
	predict `temp_expB' if `touse', hr 

	mat `temp'=`new_est', e(b), `phi_initial'
   	mat `new_est'=`temp'

	gen `temp_marg'=_d-`temp_marg_in'     	     if  `touse'
	
  	******************
  	while (tolv>`tolvEM'  &  iter<`iterEM')  {
    		scalar iter=iter+1
     		mat `old_est'=`new_est'
     		scalar Phi=`old_est'[1,colsof(`old_est')] 
     		scalar IPhi=1/Phi

     		local drop_var_list "temp_cumH est_q est_p temp_Ai temp_Ci temp_offset temp_diff_g12"
     		foreach drop_var of local drop_var_list {
					    cap drop ``drop_var''
		    			tempvar `drop_var'
     		}
     		bysort `subject': egen `temp_cumH'=sum(`temp_marg')  if `touse'
       
     		gen `est_q'=cond(`Di'==0, `temp_pi'/((1+Phi*`temp_cumH')^IPhi), 1)   if `touse'
     		gen `est_p'=cond(`Di'==0, `est_q'/(`est_q'+1-`temp_pi'), 1)          if `touse'

    		gen `temp_Ai'=IPhi+`Di'             if `touse'
     		gen `temp_Ci'=IPhi+`temp_cumH'      if `touse'
     		gen `temp_diff_g12'=digamma(`temp_Ai')-ln(`temp_Ci')-`temp_Ai'/`temp_Ci'   if `touse'
            if ("`xoffset'" != "") {
            		gen `temp_offset'=ln(`est_p'*`temp_Ai'/`temp_Ci')+`xoffset'   if `touse'  
            }
            else {
				gen `temp_offset'=ln(`est_p'*`temp_Ai'/`temp_Ci')              if `touse'
            }

     		******************
     		local drop_var_list "temp_pi temp_expB temp_marg temp_marg_in temp_H0 temp_wbl_t temp_wbl_t0 temp_wbl_coef"
     		foreach drop_var of local drop_var_list {
							  cap drop ``drop_var''
							  tempvar `drop_var'
            }

     		******************
     		cap drop `temp_zt'
     		tempvar temp_zt 
     		if("`constant'" == "") {                   
                    	sum `est_p' if `first_id'==1, meanonly
                    	if("`linkarg'" == "logit") {
			     		    local b0=logit(r(mean))                     
                        }
                   	    else if ("`linkarg'" == "probit") {
			     		    local b0=invnormal(r(mean))                  
                    	}
                    	else if ("`linkarg'" == "cloglog") {
			     		    local b0=ln(-ln(1-r(mean)))  
                    	}  
                    	else {
                    	}
                        local initopt init(_cons=`b0') search(quietly)
     		}
     		ml model d2 strmcure_my`linkarg'_d2 (mu: `est_p' = `z_index' , `zoffopt'  `constant')   if `first_id'==1, `initopt'  maximize  
     		predict `temp_zt' if `touse', xb
     		mat `new_est'=e(b)

     		if("`linkarg'" == "logit") {
					gen `temp_pi'=invlogit(`temp_zt')      if `touse'                     
     		}
     		else if ("`linkarg'" == "probit") {
					gen `temp_pi'=normal(`temp_zt')        if `touse'                      
     		}
     		else if ("`linkarg'" == "cloglog") {
					gen `temp_pi'=1-exp(-exp(`temp_zt'))   if `touse'     
     		}  
     		else {
     		}  
			
     		******************
    		stcox `x_index' if `touse', offset(`temp_offset') strata(`event_order') mgale(`temp_marg_in')  basec(`temp_H0')
     		predict `temp_expB' if `touse', hr nooffset

            mat `temp'=`new_est', e(b)
   		    mat `new_est'=`temp'

			gen `temp_marg'=(_d-`temp_marg_in')/exp(`temp_offset')     	     if  `touse'
	
			forv i=1/1 {
				sum `temp_H0' if `event_order'==`i' & _d==1
				local Wtail_H=r(max)
				sum _t if `event_order'==`i' & _d==1
				local maxTi=r(max)
				gen  `temp_wbl_t'=_t/`maxTi'                     if `event_order'==`i'
				gen  `temp_wbl_t0'=_t0/`maxTi'                   if `event_order'==`i'
				gen  `temp_wbl_coef'=`temp_expB'*`Wtail_H'*exp(`temp_offset')      if `event_order'==`i'
			    ml model d2 strmcure_myWeibL_d2 (WbL: _d  `temp_wbl_t'  `temp_wbl_coef'  `temp_wbl_t0'= )  if `event_order'==`i',  maximize  search(on)

				local Wtail_p=[WbL]_b[_cons]
				local Wtail_lam=(`Wtail_H')^(1/`Wtail_p')/`maxTi'
	
				replace `temp_marg'=`temp_marg'+(`Wtail_lam'*_t)^(`Wtail_p')*`temp_expB'-(`Wtail_lam'*`maxTi')^(`Wtail_p')*`temp_expB'  if `event_order'==`i'  &  _d==0  &  _t>`maxTi'  & _t0<`maxTi'
				replace `temp_marg'=(`Wtail_lam'*_t)^(`Wtail_p')*`temp_expB'-(`Wtail_lam'*_t0)^(`Wtail_p')*`temp_expB'      if `event_order'==`i'  &  _d==0  &  _t>`maxTi'  & _t0>`maxTi'		
		    }			
			
            if "`xoffset'" != "" {
     		  	replace `temp_marg'=`temp_marg'*exp(`xoffset')    	 if `touse'	
            }
	
     		******************
			ml model d2 strmcure_myGamma_d2 (phi: `temp_diff_g12' `est_p'=)  if `first_id'==1,  search(on) maximize 

     		mat `temp'=`new_est', [phi]_b[_cons]
     		mat `new_est'=`temp'

     		******************
      	    matrix `aa'=(`new_est'-`old_est')*(`new_est'-`old_est')'
      	    scalar tolv=sqrt(`aa'[1,1])

      	    matrix `aa'=iter, `new_est'
			if(iter==1)  {
							mata: result_savem=st_matrix("`aa'")
			}
			else {
							mata:result_savem=result_savem\st_matrix("`aa'")
			}

            matname `aa' iter `z_name' `x_name' frailty:_cons, columns(.) explicit

      	    ****************** 
      	    if "`log'" !="" {
             	 if iter==1  {
 			         noi: display "{txt}{hline 60}"
                  	 noi: di as text "......EM Iterative Parameter Estimation Process Starts......"
                  	 noi: matlist `aa', names(columns)
              	 }
		  	     else {
				     noi: matlist `aa', nonames
      	         }
      	   }

  	}
  	******************
   	if("`log'" != "") {
                noi: di as text "......EM Iterative Parameter Estimation Process Ends......"  
 	            noi: display "{txt}{hline 60}"
    }
	
}
****************************************************************************************  
else {
}
****************************************************************************************

****************************************************************************************
   cap mat drop `temp'
   cap mat drop `aa'
   cap mat drop `old_est'

   tempname b V
   mat `b'=`new_est'
   mat `V'=J(colsof(`new_est'), colsof(`new_est'), 0)

   *noi di _newline(1)
****************************************************************************************

****************************************************************************************
************   			Baseline hazard estimate 
****************************************************************************************
if("`baseh0'" !="") {
    
	   scalar Phi=`new_est'[1,colsof(`new_est')] 
       scalar IPhi=1/Phi

       local drop_var_list "temp_cumH est_q est_p temp_Ai temp_Ci temp_offset"
       foreach drop_var of local drop_var_list {
					    cap drop ``drop_var''
		    			tempvar `drop_var'
     	}
       bysort `subject': egen `temp_cumH'=sum(`temp_marg')  if `touse'    
      
       if(cureI==1) {       
     			gen `est_q'=cond(`Di'==0, `temp_pi'/((1+Phi*`temp_cumH')^IPhi), 1)   if `touse'
     			gen `est_p'=cond(`Di'==0, `est_q'/(`est_q'+1-`temp_pi'), 1)          if `touse'
            	if(tailmet==1) {
                  	           replace `est_p'=0          	if `ind_Constr'==1
			    }
	   }
	   else {
     			gen byte `est_q'=1		    if `touse'
     			gen byte `est_p'=1  		 if `touse'
	   }

       gen `temp_Ai'=IPhi+`Di'             if `touse'
       gen `temp_Ci'=IPhi+`temp_cumH'      if `touse'
	  
       if ("`xoffset'" != "") {
            	 gen `temp_offset'=ln(`est_p'*`temp_Ai'/`temp_Ci')+`xoffset'    if `touse'   
       }
       else {
				 gen `temp_offset'=ln(`est_p'*`temp_Ai'/`temp_Ci')              if `touse'
       }
  
       cap drop `est_h0'
       cap drop `rs_id'
       cap drop `est_h0_con'
       tempvar est_h0 rs_id est_h0_con
   
       gen `est_h0'=.	   
       gen `est_h0_con'=`temp_expB'*exp(`temp_offset')  if `touse'

       sort `event_order' _t 
	   egen `rs_id'=group(`event_order' _t)      if _d==1  &  `touse' 
	   
	   tempname total_D 

	   sum `rs_id' if `touse'
       local `total_D'=r(max)

       forv j=1/``total_D'' {						
							sum _t  if  `rs_id'==`j' 
							local rs_time=r(min)
							local rs_N=r(N)
							
							sum `event_order' if  `rs_id'==`j' 
							local rs_evr=r(min)
							
							cap drop `RS_ind'
							tempvar RS_ind							
							gen `RS_ind'=cond(`event_order'==`rs_evr' & _t0<`rs_time' & _t>=`rs_time', 1, 0) 

							cap drop `temp1'
							tempvar temp1
							egen `temp1'=sum(`est_h0_con')   if `RS_ind'==1
							replace `est_h0'=`rs_N'/`temp1'  if `rs_id'==`j' 
		  
	   }

	 cap confirm new variable `baseh0', exact
	 if(_rc==0) {
	 	 clonevar `baseh0'=`est_h0' if  `touse' 
		 noi: di as text "baseline hazard estimate stored in variable: `baseh0'"
	 }
	 else {    
	    	  noi: di as text "`baseh0' already exists. Replaced with baseline hazard estimates"
              cap drop `baseh0'
	 	      clonevar `baseh0'=`est_h0' if  `touse'   
	 }
}
 
 ereturn clear	
***************************************************
************variance estimation by Louis' Formula
***************************************************				 
if(`reps' == 0) {
       noi di _n as text "no variance estimation required."
}
else { 
**********************************begin of else***** 	
	  scalar Phi=`new_est'[1,colsof(`new_est')] 
      scalar IPhi=1/Phi
      scalar WM_p1=digamma(IPhi)+ln(Phi)-1
      scalar WM_p2=(IPhi-IPhi^2*trigamma(IPhi))*IPhi^2
     
      local drop_var_list "temp_cumH est_q est_p temp_Ai temp_Ci "
      foreach drop_var of local drop_var_list {
						    cap drop ``drop_var''
		    				tempvar `drop_var'
      }
      bysort `subject': egen `temp_cumH'=sum(`temp_marg')  if `touse'    
      
      if(cureI==1) {       
     			gen `est_q'=cond(`Di'==0, `temp_pi'/((1+Phi*`temp_cumH')^IPhi), 1)   if `touse'
     			gen `est_p'=cond(`Di'==0, `est_q'/(`est_q'+1-`temp_pi'), 1)          if `touse'
            	if(tailmet==1) {
                  	           replace `est_p'=0          	if `ind_Constr'==1
			    }
	  }
	  else {
     			gen byte `est_q'=1		    if `touse'
     			gen byte `est_p'=1  		 if `touse'
	  }

      gen `temp_Ai'=IPhi+`Di'             if `touse'
      gen `temp_Ci'=IPhi+`temp_cumH'      if `touse'
	 
	********matrix invariant to the unknown ki and wi
    cap mat drop  L1_beta_part1
	mat vecaccum  L1_beta_part1=_d `x_index' if `touse', noc

	if(cureI==1) {
      	  cap mat drop L1_theta_part1
	      cap mat drop EL2_theta
		  tempvar temp1 temp2
            
            if("`linkarg'" == "logit") {
 			       gen `temp1'=`temp_pi'   		  if `first_id'==1  
			       gen `temp2'=-`temp_pi'*(1-`temp_pi')  if `first_id'==1            
            }
            else if ("`linkarg'" == "probit") {
			       gen `temp1'=normalden(`temp_zt')/(1-`temp_pi')         if `first_id'==1   
                   gen `temp2'=normalden(`temp_zt')/(1-`temp_pi')*(`temp_zt'-normalden(`temp_zt')/(1-`temp_pi'))   if `first_id'==1 
            }
            else {
			       gen `temp1'=exp(`temp_zt')       if `first_id'==1 
                   gen `temp2'=-`temp1'             if `first_id'==1 
            }

 	      mat vecaccum L1_theta_part1=`temp1' `z_index'     if `first_id'==1,  `constant' 
	      mat accum EL2_theta_part1=`z_index' [iw=`temp2']   if `first_id'==1,  `constant'    
  	}
	**********************************
	tempfile temp_save0
    save `temp_save0', replace
 
**********************************begin of capture***************************************************
capture { 
    keep if `touse'
    if(tailmet==1) {
              keep if `ind_Constr'==0  
    } 	
    *******
    cap drop subject 
    sort `subject' _t _t0
    egen subject=group(`subject') 

    local keep_var_list "first_id  temp_expB est_p temp_Ai temp_Ci event_order old_ev"
    foreach keep_var of local keep_var_list {
									 cap drop `keep_var'
									 rename ``keep_var'' `keep_var'
    }  

    if "`xoffset'" != "" {
         	 replace temp_expB=temp_expB*exp(`xoffset')  if `touse'
    }

    if(cureI==1) {   
	      cap drop tmp_piL1
		  cap drop tmp_piL2
            if("`linkarg'" == "logit") {
 			    gen tmp_piL1=first_id       
			    gen tmp_piL2=0                          
            }
            else if ("`linkarg'" == "probit") {
				gen tmp_piL1=normalden(`temp_zt')*(1/(1-`temp_pi')+1/`temp_pi')*first_id      				
				gen tmp_piL2=(-normalden(`temp_zt')/`temp_pi'*(`temp_zt'+normalden(`temp_zt')/`temp_pi')-`temp2')*first_id     
                replace  tmp_piL2=0      if   missing(tmp_piL2)
            }
            else {
				gen tmp_piL1=exp(`temp_zt')/`temp_pi'*first_id     		
                gen tmp_piL2=tmp_piL1*(1-tmp_piL1+exp(`temp_zt'))          
            } 
    }    

    sort event_order _t _t0 subject old_ev 
    cap drop seqN
    gen seqN=_n 
 
    cap drop Enum
    gen Enum=sum(_d)
    replace Enum=. if _d==0
		  
    summ _d
    scalar maxRS=r(sum)
  
    sum seqN
    scalar Nobs=r(max)

    cap drop start_t
    cap drop end_t
    cap drop status
    rename _t0 start_t
    rename _t end_t
    rename _d status 
    
    if("`rspart'" == "")  {
  				scalar RS_part=ceil(maxRS/500)
    }
    else {
			  	scalar RS_part=`rspart'
    }

    if("`obspart'" == "") {
  			      scalar obs_part=ceil(Nobs/500)
    }
    else {
		  	      scalar obs_part=`obspart'
    }

    cap drop DS
    cap matrix drop `RS_list'
   
    sort Enum	
    mata: strmcure_FillRS_MATA("Enum event_order end_t", st_numscalar("RS_part"), st_numscalar("obs_part"), st_numscalar("maxRS"), st_numscalar("Nobs")) 
    cap drop DS

    **************
    local Reps=ceil(`reps'/10)*10
    if("`repspart'" == "")   {
        local Rep_W=`Reps'/10
        local max_W=`Reps'/`Rep_W'
    } 
    else {
        local max_W=ceil(`Reps'/`repspart')
        local Rep_W=`Reps'/`max_W'      	
    }	  
    scalar Rep_W=`Rep_W'      
    scalar max_W=`max_W'
 
    ***********************************************    
   if(cureI==1) {  
	scalar zcons=1
	if("`znoconstant'" != "") scalar zcons=0
     	cap matrix drop EL1_par_total
    	cap matrix drop EL2_total
    	cap drop Sim_k_*
    	cap drop Sim_w_*
    	forv i=1/`max_W' {
			 gen byte Sim_k_`i'=.
    	}
    	unab Simk: Sim_k_*
    	forv i=1/`max_W' {
			 gen Sim_w_`i'=.
    	}
    	unab Simw: Sim_w_*

    	set seed `iseed'

    	forv Rep_i=1/`Rep_W' {
           scalar Rep_i_comp=`Rep_i'       
           gsort subject -first_id			  
           forv i=1/`max_W' {							  
				replace Sim_k_`i'=runiform()<=est_p  if first_id==1  
				replace Sim_w_`i'=cond(Sim_k_`i'==1, rgamma(temp_Ai, 1/temp_Ci),0)  if first_id==1
				replace Sim_w_`i'=1     if  missing(Sim_w_`i')  &  first_id==1
				by subject: replace Sim_k_`i'=Sim_k_`i'[1] if _n>1 & first_id~=1
				by subject: replace Sim_w_`i'=Sim_w_`i'[1] if _n>1 & first_id~=1
           }
           if("`dots'" != "") {
           		if(`Rep_i'==1) nois _dots 0, title(Estimating Variance (`reps' realizations in `Rep_W' loops). Loop running...) reps(`Rep_W')
           		nois _dots `Rep_i' 0
           } 

           sort seqN
           mata: strmcure_VARFILLC_MATA(st_numscalar("Rep_W"), st_numscalar("Rep_i_comp"), "first_id", "tmp_piL1", "tmp_piL2", "`x_index'", "temp_expB", "`Simk'", "`Simw'",  st_matrix("L1_beta_part1"), st_numscalar("RS_part"), st_numscalar("obs_part"), st_numscalar("zcons"), "`z_index'", st_matrix("L1_theta_part1"), st_numscalar("WM_p1"), st_numscalar("WM_p2"), st_numscalar("IPhi"), st_matrix("EL2_theta_part1"))				   
	}
      noi di _newline(1)
   }
   else {
     	cap matrix drop EL2_total
    	cap drop Sim_w_*
    	forv i=1/`max_W' {
      	    gen Sim_w_`i'=.
    	}
    	unab Simw: Sim_w_*	
    	set seed `iseed'

    	forv Rep_i=1/`Rep_W' {
          	 scalar Rep_i_comp=`Rep_i'       
          	 gsort subject -first_id					  
          	 forv i=1/`max_W' {								  
				replace Sim_w_`i'=rgamma(temp_Ai, 1/temp_Ci)     if first_id==1
				replace Sim_w_`i'=1                              if missing(Sim_w_`i') & first_id==1
				by subject: replace Sim_w_`i'=Sim_w_`i'[1]       if _n>1 & first_id~=1
           	 }

           	 if("`dots'" != "") {
           		if(`Rep_i'==1) nois _dots 0, title(Estimating Variance (`reps' realizations in `Rep_W' loops). Loop running...) reps(`Rep_W')
           		nois _dots `Rep_i' 0
             }

             sort seqN
             mata: strmcure_VARFILL_MATA(st_numscalar("Rep_W"), st_numscalar("Rep_i_comp"), "first_id", "`x_index'", "temp_expB", "`Simw'", st_matrix("L1_beta_part1"),  st_numscalar("RS_part"), st_numscalar("obs_part"),  st_numscalar("WM_p1"), st_numscalar("WM_p2"), st_numscalar("IPhi"))					  
     } 
      noi di _newline(1)      
  }		
}
**********************************end of capture***************************************************
    local rc1=_rc
    use `temp_save0',clear

    cap matrix drop `RS_list'
    cap matrix drop EL2_theta
    cap matrix drop L1_theta_part1
    cap matrix drop EL2_theta_part1
    cap matrix drop L1_beta_part1
    cap matrix drop EL1_par_total
    cap matrix drop EL2_total
    cap drop matrix VAR_EST

    if(`rc1'==0) {
		   tempname nT nA
		   scalar `nT'=colsof(`new_est')
               capture {
       	          mata: st_matrix("VAR_EST",invsym(st_matrix("L2_ad"))[1..st_numscalar("`nT'"),1..st_numscalar("`nT'")])  
               }
               if(_rc==0)    mat `V'=VAR_EST
    }
**********************************end of else***** 	
}

   cap matrix drop VAR_EST
   cap matrix drop  L2_ad
   ereturn clear

   matname `b' `z_name' `x_name' frailty:_cons, columns(.) explicit
   matname `V' `z_name' `x_name' frailty:_cons, explicit

   count if `touse'
   local rN=r(N)
   ereturn post `b' `V', esample(`touse')  obs(`rN')

   if("`strata'" !="") ereturn local strata "`strata'"  
   ereturn local depvar "_t"
   ereturn local failure "_d"
   ereturn local t0 "_t0"
   ereturn local shared "`idcluster'"
   ereturn local ties "breslow"
   

   ereturn scalar N_sub=`N_sub'
   ereturn scalar N_fail=`N_fail'
   ereturn scalar N_g=`N_g'
   ereturn scalar g_min=`g_min'
   ereturn scalar g_max =`g_max'
   ereturn scalar g_avg=`g_avg'
   ereturn scalar total_time=`time_ar'
   
   if("`saving'"~="") {
   mata: st_matrix("`result_save'", result_savem)
   matname `result_save' iter `z_name'  `x_name'  frailty:_cons, columns(.) explicit
   ereturn matrix EMiter_par `result_save'
   }
   
   ereturn scalar CurePresent=cureI
   ereturn scalar tolv_cvg=tolv
   ereturn scalar iter_cvg=iter
   ereturn scalar phi_initial=`phi_initial'
   
   if(cureI==1) {
        ereturn scalar Tail=tailmet
        ereturn local link `linkarg'
        ereturn local Tailname "`tailname'"
        if ("`zoffset'" != "") ereturn local offsetZ "`zoffopt'"
   }
   if ("`xoffset'" != "") ereturn local offsetX "`xoffopt'"

   if(cureI==1) 	ereturn local title "frailty mixture model"
   if(cureI==0) 	ereturn local title "frailty model"
   ereturn local cmd "strmcure" 
   ereturn local cmdline "strmcure `0'"
   
   local h = cond(`"`strata'"'=="", "No stratification", "Stratified analysis")
   local h1="Breslow method for ties" 
   local crlen=21	
   local title=upper(substr(`"`e(title)'"',1,1)) + substr(`"`e(title)'"',2,.)

 noi{ 
        di as txt "`title' -- `h' "
        di _col(10) as txt "`h1'" _col(49) as txt /*
		*/ "Number of obs" _col(68) "=" _col(70) as res %9.0g e(N)
        di _col(10) as txt "Gamma shared frailty" _col(49) as txt /*
		*/ "Number of groups" _col(68) "=" _col(70) as res %9.0g e(N_g)
        di as txt "Group variable: " as res abbrev("`e(shared)'",12)    

        di _n as txt %-`crlen's "No. of subjects" " = " as res %12.0g  e(N_sub) _col(49)  /*
                */ as txt "Obs per group: min = " as res %9.0g   e(g_min)

        di as txt %-`crlen's "No. of failures" " = " as res %12.0g  e(N_fail)   /*
                   */ _col(64) as txt "avg" " = " as res %9.0g    e(g_avg)
        di as txt %-`crlen's "Total analysis time" " = " as res %12.0g   e(total_time)  /*
                 */  _col(64) as txt "max" " = " as res %9.0g    e(g_max)
        di
    
        Replay, `display_options'	
   }
   
   if("`rc1'" !="") {
      if(`rc1' !=0)  noi error `rc1'
   }
}	   
end






**********************MATA files
cap mata: mata drop strmcure_VARFILLC_MATA()
mata void strmcure_VARFILLC_MATA(real scalar Rep_W, real scalar Rep_i, string scalar varname0, string scalar varname1, string scalar varname2, string scalar x_index, string scalar eb, string scalar klist, string scalar wlist, real matrix L1_b1, real scalar RS_part, real scalar obs_part, real scalar zcons, string scalar z_index, real matrix L1_t1, real scalar wp1, real scalar wp2, real scalar iphi, real matrix EL2_t1)
{ 
  matrix x, z, z1, k, w, KZ, expB, kw, kwx
  matrix temp0, temp1, temp2, temp_L2_1, temp_L2_2, tempNs1, tempNs, tempNq, tempNq1, xsxq
  matrix temp2_0, L1_beta_F, EL1_par, EL2_beta, L2_par, a, EL_MEAN
  matrix L1_theta, L1_phi, EL1_theta, EL2_theta, EL1_phi, EL2_phi, L1_par, EL1L1_par
  matrix tempRS, tempN, temp_L1_1, L1_beta, EL1_beta
  real scalar mm, counter, i, j, mk, s, q, colsk
  string scalar name
  matrix tmp_pi1, tmp_pi2

  x=.
  st_view(x,.,tokens(x_index),0)
  z=.
  st_view(z,.,tokens(z_index),0)
  k=.
  st_view(k,.,tokens(klist),0)
  w=.
  st_view(w,.,tokens(wlist),0)
  expB=.
  st_view(expB,.,tokens(eb),0)
  tmp_pi1=.
  st_view(tmp_pi1,.,tokens(varname1),0)  
  tmp_pi2=.
  st_view(tmp_pi2,.,tokens(varname2),0)  
  first_id=.
  st_view(first_id,.,tokens(varname0),0) 
 
  colsk=cols(k)
  if(zcons==0) 	z1=z
  if(zcons==1) 	z1=z,J(rows(z),1,1)
  temp0=z1:*tmp_pi1
  temp1=temp0'*k
  L1_theta=temp1:-L1_t1'
  EL1_theta=rowsum(L1_theta)/colsk
  
  temp0=.
  temp1=.
  temp_L2_1=.
  EL2_theta=.
  temp0=k:*tmp_pi2
  temp1=rowsum(temp0)/colsk
  temp_L2_1=cross(z1,temp1,z1)  
  EL2_theta=temp_L2_1+EL2_t1
  
  temp0=editmissing(log(w),0)
  temp1=w-temp0
  temp2=(k:*(temp1:+wp1))*(iphi^2)
  L1_phi=colsum(temp2:*first_id)
  EL1_phi=sum(L1_phi)/colsk

  
  temp2=k:*(-2*(temp1:+wp1)*(iphi^3):+wp2)
  EL2_phi=sum(temp2:*first_id)/colsk

  kw=k:*w
  kwx=expB:*kw

  mm=cols(x)
  temp_L1_1=J(mm, colsk,0)
  temp_L2_1=J(mm,mm,0)
  temp_L2_2=J(mm,mm,0)
  counter=0
 
 for (i=1; i<=RS_part; i++)
  { 
    for (j=1; j<=obs_part; j++)
    {
      name="RS_"+strofreal(i)+strofreal(j)
      temp0=st_matrix(name)
      if(j==1) tempRS=temp0
      if(j>1)  tempRS=tempRS,temp0
     } 
     
     temp0=cross(tempRS', kwx)

	 for(s=1; s<=mm; s++)
	 {
	     tempNs1=cross(tempRS', x[|.,s|], kwx)
	     tempNs=tempNs1:/temp0
         tempNs1=.

	   for(q=s; q<=mm; q++)
	   {
	      if(q==s)   tempNq=tempNs
		  if(q>s) 
		  {
		       tempNq1=cross(tempRS', x[|.,q|], kwx)	
               tempNq=tempNq1:/temp0
	           tempNq1=.
		  }
		  
		   xsxq=x[|.,s|]:*x[|.,q|]
		   temp2_0=cross(tempRS', xsxq, kwx)	
		    
		   if(s==1) temp_L1_1[|q,.|]=colsum(tempNq)
		  
		   if(counter==0) 
		   {
		        temp_L2_1[s,q]=sum(tempNs:*tempNq)-sum(temp2_0:/temp0)
		    	if(q~=s) temp_L2_1[q,s]=temp_L2_1[s,q]
		   } 
	
		   if(counter>0) 
		   {
		    	temp_L2_1[s,q]=temp_L2_1[s,q]+sum(tempNs:*tempNq)-sum(temp2_0:/temp0)
			    if(q~=s) temp_L2_1[q,s]=temp_L2_1[s,q]
		   } 
		   counter=counter+1
		}
		   
     }
 	 if(i==1) L1_beta=temp_L1_1
	 if(i>1)  L1_beta=L1_beta+temp_L1_1
   }

  kwx=.
  temp0=.
  tempRS=.
  tempNs=.
  tempNq=.
  temp2_0=.

  L1_beta_F=L1_b1':-L1_beta
  EL1_beta=rowsum(L1_beta_F)/colsk
 
  EL1_par=EL1_theta\EL1_beta\EL1_phi
  if(Rep_i==1)  st_matrix("EL1_par_total", EL1_par)
  if(Rep_i>1)  
  {
   temp0=st_matrix("EL1_par_total")
   stata("cap matrix drop EL1_par_total")
   temp2_0=EL1_par+temp0
   st_matrix("EL1_par_total", temp2_0)
   temp0=.
   temp2_0=.
  }

  EL2_beta=(temp_L2_1)/colsk

  L2_par=.
  L2_par=blockdiag(blockdiag(EL2_theta, EL2_beta),EL2_phi)

  L1_par=L1_theta\L1_beta_F\L1_phi

  mm=rows(L1_par)
  mk=cols(L1_par)

  for (i=1; i<=mk; i++)
  {
    a=L1_par[|.,i|]
    L2_par=L2_par+a*a':/mk
  }


  if(Rep_i==1)  st_matrix("EL2_total", -L2_par)
  if(Rep_i>1)
  {
    temp0=st_matrix("EL2_total")
    stata("cap matrix drop EL2_total")
    temp2_0=-L2_par+temp0
    temp0=.
    st_matrix("EL2_total", temp2_0)
    temp2_0=.
  }
  if(Rep_i==Rep_W) 
  {
	temp0=st_matrix("EL1_par_total")/Rep_W
	temp2_0=st_matrix("EL2_total")/Rep_W

	temp1=temp2_0+temp0*temp0'
	
	stata("cap matrix drop EL1_par_total")
	stata("cap matrix drop EL2_total")  
    stata("cap matrix drop L2_ad")


	st_matrix("L2_ad", temp1)
      
	temp0=.
	temp2_0=.
	temp1=.
 
  }
} 


***************
cap mata: mata drop strmcure_VARFILL_MATA()
mata void strmcure_VARFILL_MATA(real scalar Rep_W, real scalar Rep_i, string scalar varname0, string scalar x_index, string scalar eb, string scalar wlist, real matrix L1_b1, real scalar RS_part, real scalar obs_part, real scalar wp1, real scalar wp2, real scalar iphi)
{ 
  matrix x, z, z1, k, w, KZ, expB, kw, kwx
  matrix temp0, temp1, temp2, temp_L2_1, temp_L2_2, tempNs1, tempNs, tempNq, tempNq1, xsxq
  matrix temp2_0, L1_beta_F, EL1_par, EL2_beta, L2_par, a, EL_MEAN
  matrix L1_theta, L1_phi, EL1_theta, EL2_theta, EL1_phi, EL2_phi, L1_par, EL1L1_par
  matrix tempRS, tempN, temp_L1_1, L1_beta, EL1_beta
  real scalar mm, counter, i, j, mk, s, q, colsk
  string scalar name

  x=.
  st_view(x,.,tokens(x_index),0)
  w=.
  st_view(w,.,tokens(wlist),0)
  expB=.
  st_view(expB,.,tokens(eb),0)
  first_id=.
  st_view(first_id,.,tokens(varname0),0) 
 
  colsk=cols(w)  
  temp0=editmissing(log(w),0)
  temp1=w-temp0
  temp2=(temp1:+wp1)*(iphi^2)
  L1_phi=colsum(temp2:*first_id)
  EL1_phi=sum(L1_phi)/colsk

  
  temp2=(-2*(temp1:+wp1)*(iphi^3):+wp2)
  EL2_phi=sum(temp2:*first_id)/colsk

  kwx=expB:*w

  mm=cols(x)
  temp_L1_1=J(mm, colsk,0)
  temp_L2_1=J(mm,mm,0)
  temp_L2_2=J(mm,mm,0)
  counter=0
 
 for (i=1; i<=RS_part; i++)
  { 
    for (j=1; j<=obs_part; j++)
    {
      name="RS_"+strofreal(i)+strofreal(j)
      temp0=st_matrix(name)
      if(j==1) tempRS=temp0
      if(j>1)  tempRS=tempRS,temp0
     } 
     
     temp0=cross(tempRS', kwx)

	 for(s=1; s<=mm; s++)
	 {
	     tempNs1=cross(tempRS', x[|.,s|], kwx)
	     tempNs=tempNs1:/temp0
         tempNs1=.

	   for(q=s; q<=mm; q++)
	    {
	      if(q==s)   tempNq=tempNs
		  if(q>s) 
		  {
		       tempNq1=cross(tempRS', x[|.,q|], kwx)	
               tempNq=tempNq1:/temp0
	           tempNq1=.
		  }
		  
		  xsxq=x[|.,s|]:*x[|.,q|]
		  temp2_0=cross(tempRS', xsxq, kwx)	
		    
		  if(s==1) temp_L1_1[|q,.|]=colsum(tempNq)
		  
		  if(counter==0) 
		  {
		      temp_L2_1[s,q]=sum(tempNs:*tempNq)-sum(temp2_0:/temp0)
		    	if(q~=s) temp_L2_1[q,s]=temp_L2_1[s,q]
		  } 
	
		  if(counter>0) 
		  {
		    	temp_L2_1[s,q]=temp_L2_1[s,q]+sum(tempNs:*tempNq)-sum(temp2_0:/temp0)
			    if(q~=s) temp_L2_1[q,s]=temp_L2_1[s,q]
		  } 
		  counter=counter+1
		}
		   
     }
 	 if(i==1) L1_beta=temp_L1_1
	 if(i>1)  L1_beta=L1_beta+temp_L1_1
   }

  kwx=.
  temp0=.
  tempRS=.
  tempNs=.
  tempNq=.
  temp2_0=.

  L1_beta_F=L1_b1':-L1_beta
  EL1_beta=rowsum(L1_beta_F)/colsk
 
  EL1_par=EL1_beta\EL1_phi
  if(Rep_i==1)  st_matrix("EL1_par_total", EL1_par)
  if(Rep_i>1)  
  {
   temp0=st_matrix("EL1_par_total")
   stata("cap matrix drop EL1_par_total")
   temp2_0=EL1_par+temp0
   st_matrix("EL1_par_total", temp2_0)
   temp0=.
   temp2_0=.
  }

  EL2_beta=(temp_L2_1)/colsk

  L2_par=.
  L2_par=blockdiag(EL2_beta,EL2_phi)

  L1_par=L1_beta_F\L1_phi

  mm=rows(L1_par)
  mk=cols(L1_par)

  for (i=1; i<=mk; i++)
  {
    a=L1_par[|.,i|]
    L2_par=L2_par+a*a':/mk
  }


  if(Rep_i==1)  st_matrix("EL2_total", -L2_par)
  if(Rep_i>1)
  {
    temp0=st_matrix("EL2_total")
    stata("cap matrix drop EL2_total")
    temp2_0=-L2_par+temp0
    temp0=.
    st_matrix("EL2_total", temp2_0)
    temp2_0=.
  }
  if(Rep_i==Rep_W) 
  {
	temp0=st_matrix("EL1_par_total")/Rep_W
	temp2_0=st_matrix("EL2_total")/Rep_W

	temp1=temp2_0+temp0*temp0'
	
	stata("cap matrix drop EL1_par_total")
	stata("cap matrix drop EL2_total")  
    stata("cap matrix drop L2_ad")


	st_matrix("L2_ad", temp1)
      
	temp0=.
	temp2_0=.
	temp1=.
 
  }
} 

*********
cap mata: mata drop strmcure_FillRS_MATA()
mata void strmcure_FillRS_MATA(string scalar varname, real scalar RS_part, real scalar obs_part, real scalar maxRS, real scalar obs)
{ 
  matrix x2, temp_RS, temp2, temp3
  real scalar blocki, blockj, i, j, starti, startj, endi, endj, totalrow
  string scalar RS_list, name
  temp2=.
  temp3=. 
  temp_RS=.
  x2=.

  st_view(x2,.,tokens(varname),0)

  blocki=floor(maxRS/RS_part)
  blockj=floor(obs/obs_part)
 
  RS_list=""
  for (i=1; i<=RS_part; i++)
  { 
    starti=1+(i-1)*blocki 
    endi=maxRS*(i==RS_part)+i*blocki*(i~=RS_part)
    temp_RS=.
	totalrow=endi-starti+1
	temp_RS=J(totalrow, obs, 0)

	for (j=starti; j<=endi; j++)
	{
    	st_numscalar("delta_str", x2[j,2])
	    st_numscalar("delta_time", x2[j,3])
        stata("qui cap drop DS")
	    stata("qui gen DS=1 if start_t<delta_time & end_t>=delta_time & event_order==delta_str")
        temp2=.
        st_view(temp2,.,tokens("seqN DS"),0)
        temp3=temp2[|.,1|]
	    temp_RS[j-starti+1,temp3']=J(1,rows(temp3),1)
   	}	  

	for (j=1; j<=obs_part; j++)
    {
       startj=1+(j-1)*blockj 
       endj=obs*(j==obs_part)+j*blockj*(j~=obs_part)

       name="RS_"+strofreal(i)+strofreal(j)
       RS_list=RS_list+ "  " +name
       st_matrix(name, temp_RS[|1,startj \ totalrow, endj|])
    } 
	
   }

  st_local("RS_list", RS_list)
}


