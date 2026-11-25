capture program drop strmcure_myprobit_d2

program define strmcure_myprobit_d2
 version 11.0
 args todo b lnf g H 
 tempvar xb lj1 
 mleval `xb'=`b'
 qui {
   gen double `lj1'=normal(`xb')
   mlsum `lnf'=$ML_y1*ln(`lj1')+(1-$ML_y1)*ln(1-`lj1')

   if(`todo'==0 | `lnf'>=.) exit
   mlvecsum `lnf' `g'=($ML_y1-`lj1')*normalden(`xb')/(`lj1'*(1-`lj1')), eq(1)

   if(`todo'==1 | `lnf'>=.) exit
   tempvar nng1 nng2
   gen double `nng1'=normalden(`xb')/`lj1'
   gen double `nng2'=normalden(`xb')/(1-`lj1')
   mlmatsum `lnf' `H'=-($ML_y1*`nng1'*(`nng1'+`xb')+(1-$ML_y1)*`nng2'*(`nng2'-`xb')), eq(1,1)
 }
end
