****************************************************************
* Regressions of Markov Models
* Rafael Perez Ribas
* January 19, 2007
****************************************************************

****************************************************************
* System parameters
****************************************************************

set mem 900m
set mat 5000
set seed 123456789
set mo off

cd "C:\PME"


****************************************************************
* Database
****************************************************************

u pme_set, clear


****************************************************************
* Survey data setting
****************************************************************

svyset psu [pw=v215], str(stratum)


****************************************************************
* Descriptive Analysis
****************************************************************

cd Log

log using description.log, replace


* Panel Attrition

svy: ta p0id1 p0id0 if p0id1!=. & (retent==0 | (retent==1 & p0id0!=.)) ///
   , miss for(%6.2f) per
svy: ta p0id1 p0id0 if p0id1!=. & (retent==0 | (retent==1 & p0id0!=.)) ///
   , miss for(%6.2f) per row


* True transitions with the retained sample

svy: ta p0id1 p0id0, for(%6.2f) per
svy: ta p0id1 p0id0, for(%6.2f) per row

log off


* Matrix of estimated transitions with pseudo-panel data

tempvar perm
g `perm' = min(p0gd0,p0gd1)
qui mean `perm' p0gd0 p0gd1 if repres==1 [aw=pesoa]

loc p1t = 100*_b[p0gd1]
loc p0t = 100 - (100*_b[p0gd1])
loc pt1 = 100*_b[p0gd0]
loc pt0 = 100 - (100*_b[p0gd0])

loc p11 = 100*_b[`perm']
loc p10 = `p1t' - (100*_b[`perm'])
loc p01 = `pt1' - (100*_b[`perm'])
loc p00 = `p0t' - `p01'

loc p111 = 100*`p11'/`p1t'
loc p101 = 100*`p10'/`p1t'
loc p010 = 100*`p01'/`p0t'
loc p000 = 100*`p00'/`p0t'

mat pseudo = `p00', `p01', `p0t' \ ///
		 `p10', `p11', `p1t' \ ///
		 `pt0', `pt1', 100
mat rown pseudo = "nonpoor t-1" "poor t-1" "total"
mat coln pseudo = "nonpoor t" "poor t" "total"

mat pseudor = `p000', `p010', 100 \ ///
		  `p101', `p111', 100 \ ///
		  `pt0' , `pt1' , 100
mat rown pseudor = "nonpoor t-1" "poor t-1" "total"
mat coln pseudor = "nonpoor t" "poor t" "total"

log on


* Estimated transitions with pseudo-panel data

mat l pseudo, f(%6.2f)
mat l pseudor, f(%6.2f)

log close

cd ..


****************************************************************
* Program Codes
****************************************************************

do models.do


****************************************************************
* Test of Instrumental Variables
****************************************************************

cd Log

log using regressions.log, replace

log off

ml model lf gbivprob ///
   (initial: p0gd1 = a2000-a2003 NE RJ_SP sexo age agee2 escol4-escol11 p0gd2) ///
   (permanence: p0gd0 = a2000-a2003 NE RJ_SP sexo age agee2 escol4-escol11) ///
   (transition: a2000-a2003 NE RJ_SP sexo age agee2 escol4-escol11) /rho ///
   if repres==1 [fw=pesof], tech(dfp nr) max dif search(off)
est store pseudo

log on

* For cohorts

ml display

log off

ml model lf markovreg ///
   (initial: p0id1 = a2000-a2003 NE RJ_SP sexo age agee2 escol4-escol11 p0gd2) ///
   (permanence: p0id0 = a2000-a2003 NE RJ_SP sexo age agee2 escol4-escol11) ///
   (transition: a2000-a2003 NE RJ_SP sexo age agee2 escol4-escol11) /athrho ///
   , tech(dfp nr) max dif search(off) svy
est store panel

log on

* For individuals without selection (Bivariate Probit)

ml display

log off


qui {
   probit retent a2000-a2003 NE RJ_SP sexo age agee2 escol4-escol11 ///
   		     entrev2-entrev4 [pw=v215]
   mat b1 = e(b)
   mat coleq b1 = retention

   probit p0id1 a2000-a2003 NE RJ_SP sexo age agee2 escol4-escol11 ///
   		    p0gd2 [pw=v215]
   mat b2 = e(b)
   mat coleq b2 = initial

   probit p0id0 a2000-a2003 NE RJ_SP sexo age agee2 escol4-escol11 ///
   		    if p0id1==1 [pw=v215]
   mat b3 = e(b)
   mat coleq b3 = permanence

   probit p0id0 a2000-a2003 NE RJ_SP sexo age agee2 escol4-escol11 ///
   		    if p0id1==0 [pw=v215]
   mat b4 = e(b)
   mat coleq b4 = transition

   mat b0 = b1, b2, b3, b4
}

capture drop z*
mdraws, dr(100) neq(3) prefix(z) antithetics random
global dr = r(n_draws)
   
ml model lf markovreg2 (retention: retent = a2000-a2003 NE RJ_SP ///
   sexo age agee2 escol4-escol11 entrev2-entrev4) ///
   (initial: p0id1 = a2000-a2003 NE RJ_SP sexo age agee2 escol4-escol11 p0gd2) ///
   (permanence: p0id0 = a2000-a2003 NE RJ_SP sexo age agee2 escol4-escol11) ///
   (transition: a2000-a2003 NE RJ_SP sexo age agee2 escol4-escol11) ///
   /c21 /c31 /c32 ///
   if p0id1!=. & ((p0id0==. & retent==0) | (p0id0!=. & retent==1)) ///
   , tech(dfp nr) missing max dif search(off) init(b0) svy
est store panel2
capture drop z*

log on

* For individuals with selection (Joint Trivariate Probit)

ml display

nlcom ( r21: [c21]_b[_cons] ) ///
	( r31: [c31]_b[_cons] ) ///
	( r32: [c21]_b[_cons]*[c31]_b[_cons] ///
	   + sqrt(1 - [c21]_b[_cons]^2)*[c32]_b[_cons] )


hausman panel2 panel1, force c a
hausman panel2 pseudo, force c a
hausman panel1 pseudo, force c a


* For individuals with selection (Two-step Trivariate Probit)

biprobit (p0id1 = a2000-a2003 NE RJ_SP sexo age agee2 escol4-escol11 p0gd2) ///
   (retent = a2000-a2003 NE RJ_SP sexo age agee2 escol4-escol11 entrev2-entrev4) ///
   [pw=v215]

qui {
   tempvar xb2 xb3 k2 k3 p2 p3
   predict `xb2', xb1
   predict `xb3', xb2

   g `k2' = 2*p0id1 - 1
   g `k3' = 2*retent - 1
   
   g `p2' = ((-normalden(-`k2'*`xb2')*(1 - normal((`k3'*(-`xb3' + e(rho)*`xb2'))/((1 - (e(rho)^2))^.5)))) ///
   - (`k2'*`k3'*e(rho)*normalden(-`k3'*`xb3')*(1 - normal((`k2'*(-`xb2' + e(rho)*`xb3'))/((1 - (e(rho)^2))^.5)))))/ ///
   binormal(`xb2', `xb3', e(rho))

   g `p3' = ((-normalden(-`k3'*`xb3')*(1 - normal((`k2'*(-`xb2' + e(rho)*`xb3'))/((1 - (e(rho)^2))^.5)))) ///
   - (`k2'*`k3'*e(rho)*normalden(-`k2'*`xb2')*(1 - normal((`k3'*(-`xb3' + e(rho)*`xb2'))/((1 - (e(rho)^2))^.5)))))/ ///
   binormal(`xb2', `xb3', e(rho))

   g r12 = (`p2' - e(rho)*`p3')/(1 - e(rho)^2)
   g r13 = (`p3' - e(rho)*`p2')/(1 - e(rho)^2)
}

svy: probit p0id0 p0id1 a2000-a2003 NE RJ_SP sexo age agee2 escol4-escol11 r12 r13
   		    
log close

cd ..

****************************************************************
* End of Do file
****************************************************************
