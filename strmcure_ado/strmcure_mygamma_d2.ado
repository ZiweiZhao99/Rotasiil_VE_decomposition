capture program drop strmcure_myGamma_d2

program define strmcure_myGamma_d2
 version 11.0
 args todo b lnf g H 
 tempvar xb lj
 mleval `xb'=`b', scalar

 qui {
    gen double `lj'=($ML_y1/`xb'-lngamma(1/`xb')-ln(`xb')/`xb')*$ML_y2
    mlsum `lnf'=`lj'

    if(`todo'==0 | `lnf'>=.) exit
    tempvar g1
    generate double `g1'=( -$ML_y1/`xb'^2 + digamma(1/`xb')/`xb'^2 - 1/`xb'^2 + ln(`xb')/`xb'^2)*$ML_y2
    mlvecsum `lnf' `g'=`g1', eq(1)
 
    if(`todo'==1 | `lnf'>=.) exit
    mlmatsum `lnf' `H'=( 2*$ML_y1/`xb'^3 - 2*digamma(1/`xb')/`xb'^3 - trigamma(1/`xb')/`xb'^4 + 3/`xb'^3 - 2*ln(`xb')/`xb'^3 )*$ML_y2, eq(1,1)

 }

end
