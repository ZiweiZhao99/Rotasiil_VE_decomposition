capture program drop strmcure_mycloglog_d2

program define strmcure_mycloglog_d2
 version 11.0
 args todo b lnf g H 
 tempvar xb lj1 
 mleval `xb'=`b'
 qui {
   gen double `lj1'=1-exp(-exp(`xb'))
   mlsum `lnf'=$ML_y1*ln(`lj1')+(1-$ML_y1)*ln(1-`lj1')

   if(`todo'==0 | `lnf'>=.) exit
   mlvecsum `lnf' `g'=($ML_y1-`lj1')*exp(`xb')/`lj1', eq(1)

   if(`todo'==1 | `lnf'>=.) exit
   tempvar nng1 
   gen double `nng1'=(1-`lj1')*exp(`xb')/`lj1'
   mlmatsum `lnf' `H'=-($ML_y1*`nng1'*(`nng1'-1+exp(`xb'))+(1-$ML_y1)*exp(`xb')), eq(1,1)
 }
end

