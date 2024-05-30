****************************************************************
* Code for estimation programs
* Rafael Perez Ribas
* November 16, 2006
****************************************************************

* Univariate Probit for grouped data

cap prog drop bprobit2
program define bprobit2
   args lnf theta

   qui replace `lnf'=$ML_y1*ln(normprob(`theta'))+(1-$ML_y1)*ln(1-normprob(`theta'))

end


* Bivariate Switching Probit for grouped data

cap prog drop gbivprob
program define gbivprob
   args lnf theta1 theta2 theta3 theta4

   tempvar a1 a2 a3 a4 y1 y2 y3 y4
   qui {
      g double `a1'=binorm(`theta1',`theta2',`theta4')
      g double `a2'=binorm(`theta1',-`theta2',-`theta4')
      g double `a3'=binorm(-`theta1',`theta3',-`theta4')
      g double `a4'=binorm(-`theta1',-`theta3',`theta4')
      g double `y1'=min($ML_y1,$ML_y2)
      g double `y2'=max(0,$ML_y1-$ML_y2)
      g double `y3'=max(0,$ML_y2-$ML_y1)
      g double `y4'=1-max($ML_y1,$ML_y2)
      replace `lnf'=`y1'*ln(`a1')+`y2'*ln(`a2')+`y3'*ln(`a3')+`y4'*ln(`a4')
   }

end


* Bivariate Switching Probit without Selection

cap prog drop markovreg
program define markovreg
   args lf xb1 xb2 xb3 xb4

   local rho (exp(2*`xb4')-1) / (exp(2*`xb4')+1)

   qui {
      replace `lf' = binorm(`xb1',`xb2',`rho')   if $ML_y1 & $ML_y2
      replace `lf' = binorm(-`xb1',`xb3',-`rho') if !$ML_y1 & $ML_y2
      replace `lf' = binorm(`xb1',-`xb2',-`rho') if $ML_y1 & !$ML_y2
      replace `lf' = binorm(-`xb1',-`xb3',`rho') if !$ML_y1 & !$ML_y2
      replace `lf' = ln(`lf')
   }

end


* Trivariate Switching Probit with Selection

capture program drop markovreg2
program define markovreg2

   args lnf x1 x2 x3 x4 c21 c31 c32
   tempvar sp2 sp3 sp4 k1 k2 k3

qui {
   g double `k1' = 2*$ML_y1 - 1
   g double `k2' = 2*$ML_y2 - 1
   g double `k3' = 2*$ML_y3 - 1
   tempname cf21 cf22 cf31 cf32 cf33 C1 C2
   su `c21', meanonly
   sca `cf21' = r(mean)
   su `c31', meanonly
   sca `cf31' = r(mean)
   su `c32', meanonly
   sca `cf32' = r(mean)
   sca `cf22' = sqrt( 1 - `cf21'^2 )
   sca `cf33' = sqrt( 1 - `cf31'^2 - `cf32'^2 )
   mat `C1' = (1, 0 , 0 \ `cf21', `cf22', 0 \ `cf31' , `cf32' , `cf33')
   mat `C2' = (1, 0 \ `cf21', `cf22')
   egen `sp3' = mvnp(`x1' `x2' `x3') if $ML_y1==1 & $ML_y2==1, chol(`C1') ///
      dr($dr) prefix(z) signs(`k1' `k2' `k3')
   egen `sp4' = mvnp(`x1' `x2' `x4') if $ML_y1==1 & $ML_y2==0, chol(`C1') ///
      dr($dr) prefix(z) signs(`k1' `k2' `k3')
   egen `sp2' = mvnp(`x1' `x2' ) if $ML_y1==0, chol(`C2') dr($dr) prefix(z) ///
      signs(`k1' `k2')
   replace `lnf'= ln(`sp3') if $ML_y1==1 & $ML_y2==1
   replace `lnf'= ln(`sp4') if $ML_y1==1 & $ML_y2==0
   replace `lnf'= ln(`sp2') if $ML_y1==0
}

end

****************************************************************
* End of Do file
****************************************************************
