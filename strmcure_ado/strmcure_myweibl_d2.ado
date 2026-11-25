capture program drop strmcure_myWeibL_d2

program define strmcure_myWeibL_d2
 version 11.0
 args todo b lnf g H
 tempvar xb 
 mleval `xb'=`b', scalar

 qui {
   mlsum `lnf'=$ML_y1*(ln(`xb')+`xb'*ln($ML_y2))-$ML_y3*(($ML_y2)^(`xb')-($ML_y4)^(`xb'))
      
  if(`todo'==0 | `lnf'>=.) exit
   tempvar g1
   generate `g1'=ln($ML_y4) 
   replace `g1'=0 if `g1'==.
   mlvecsum `lnf' `g'=$ML_y1*(1/`xb'+ln($ML_y2))-$ML_y3*((($ML_y2)^(`xb'))*ln($ML_y2)-(($ML_y4)^(`xb'))*`g1'), eq(1)
 
  if(`todo'==1 | `lnf'>=.) exit
   mlmatsum `lnf' `H'=-$ML_y1*(1/`xb')^2-$ML_y3*(($ML_y2)^(`xb')*(ln($ML_y2))^2 -($ML_y4)^(`xb')*(`g1')^2), eq(1,1)

 }

end

