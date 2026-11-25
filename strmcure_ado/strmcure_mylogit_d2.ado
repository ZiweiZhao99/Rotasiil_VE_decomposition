capture program drop strmcure_mylogit_d2

program define strmcure_mylogit_d2
version 11.0
 args todo b lnf g H 
 tempvar xb lj1  
 qui mleval `xb'=`b'
 qui {
   gen double `lj1'=invlogit(`xb')
   mlsum `lnf'=$ML_y1*ln(`lj1')+(1-$ML_y1)*ln(1-`lj1')

   if(`todo'==0 | `lnf'>=.) exit
   mlvecsum `lnf' `g'=$ML_y1-`lj1', eq(1)

   if(`todo'==1 | `lnf'>=.) exit
   mlmatsum `lnf' `H'=-`lj1'*(1-`lj1'), eq(1,1)
 }
end

