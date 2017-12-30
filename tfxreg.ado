*! tfxreg v1.5.0	01jan2018
* Calculates teacher value-added using matrix inversion identities and
* sum-to-zero constraints from felsdvredgm (Mihaly et al 2010).

*03.14.12 Added prediction (xb) capability.
*04.19.12 Fixed abbrevation error on variable parsing.
*08.03.12 Exploit block inversion identities for faster estimation.
*04.27.17 Removed clustering option and updated error reporting.
*12.30.17 Revised solution of standard errors to avoid constructing JxJ matrix.
*         Added limited clustering options.
*         Deprecated reference collection feature.

program tfxreg2, eclass sortpreserve
	version 11
	syntax varlist [if] [in], Teacher(varname) ///
    fe(name) se(name) [ xb(name) Weights(varname) ///
	ROBUST ]
	marksample touse

tempvar sumwt wt sample miss res tid const
tempname b V var_e var_jhat var_j var_e_g var_jhat_g var_j_g

tokenize `varlist'
local depvar `1'
macro shift
local indepvar `*'

// Generate internal teacher id and sort.
qui egen `tid'=group(`teacher')
sort `tid'

// Drop existing variables and generate regression sample
novarabbrev {
	cap drop `fe'
	if _rc ==0	di "note: `fe' already exists and will be replaced."
	cap drop `se'
	if _rc ==0	di "note: `se' already exists and will be replaced."
	cap drop `xb'
	if _rc==0 di "note: `xb' already exists and will be replaced."
}

if "`weights'"!=""	qui egen `miss'=rowmiss(`varlist' `teacher' `weights')
else qui egen `miss'=rowmiss(`varlist' `teacher')
qui gen `sample'=0
qui replace `sample'=1 if `touse' & `miss'==0

// Reference collection.
if "`reff'"!="" {
  di in red "WARNING: Reference collection feature deprecated."
  di in red "Check regression output to ensure all necessary"
  di in red "variables have been included."
  di in red "If you need this feature, consider using felsdvregdm."
}
qui gen `const'=1
local indepvar `indepvar' `const'


// Check for collinear regressors.
_rmcoll `indepvar', noconstant
local indepvar = r(varlist)

// Rescale so that weights sum to the number of observations.
if "`weights'"!="" {
	qui sum `weights' if `sample'==1
	qui gen `wt'=`weights'/r(mean)
}
else {
	qui gen `wt'=1
}

// Send to Mata
mata: partreg("`depvar'", "`indepvar'", "`tid'", "`wt'", "`robust'", "`sample'", "`fe'", "`se'")

// Process results
ereturn post `b' `V', depname("`depvar'") obs(`N')

di ""
ereturn display

end


mata:

void partreg(string scalar yvar, string rowvector xvar, string scalar teachid, string scalar weight, string scalar robust, string scalar sample, string scalar fevar, string scalar sevar)
{

  // Set up workspace.
  real scalar n, nj, nk, n0
  real colvector y, j, w, r, beta, gamma, fe, se, vx, vfx
  string colvector groupstrip
  real matrix X, PJ, XX, XF, FFFX, E

  jindex 	= st_varindex(teachid)
  xindex 	= st_varindex(tokens(xvar))
  yindex 	= st_varindex(yvar)
  windex 	= st_varindex(weight)
  st_view(j=.,.,jindex,sample)
  st_view(X=.,.,(xindex),sample)
  st_view(y=.,.,yindex,sample)
  st_view(w=.,.,windex,sample)

  PJ   = panelsetup(j,1)
  nk   = cols(X)
  nj   = rows(PJ)-1
  n    = rows(X)

  // Solve for coefficients.
  beta  = solvebeta(X,y,w,PJ,nk)
  r     = y - X[,(1::(nk-1))]*beta
  gamma = solvegamma(r,w,PJ,nj)
  rmean = mean(gamma)
  beta  = (beta \ rmean)
  gamma = gamma :- rmean
  
  feindex = st_addvar("double", fevar)
  st_view(fe=.,.,feindex,sample)
  makefx(fe,PJ,gamma)
  r = r - fe

  // Estimate variance of coefficients.
  sig2e= quadcross(r,w,r) / (n-nk-nj)
  n0   = quadsum(panelsubmatrix(w,1,PJ))
  vecn = matN(PJ,w,nj)
  XX   = quadcross(X,w,X)
  XF   = matXF(X,PJ,w,nj,nk)
  FFFX = matFFFX(XF,vecn,n0,nj,nk)
  E    = invsym(XX - XF*FFFX)
  if (robust=="") {
	vx   = sig2e*E
    vfx  = sig2e*varmatF(E,FFFX,vecn,n0,nj,nk)
  }
  	
  // Return coefficient vector and variance-covariance matrix to Stata
  st_matrix(st_macroexpand("`"+"b"+"'"),beta')
  st_matrixcolstripe(st_macroexpand("`"+"b"+"'"),(J(cols(X),1," "),((st_varname(xindex)))'))
  st_matrix(st_macroexpand("`"+"V"+"'"),vx)
  st_matrixcolstripe(st_macroexpand("`"+"V"+"'") ,(J(cols(X),1," "),((st_varname(xindex)))'))
  st_matrixrowstripe(st_macroexpand("`"+"V"+"'") ,(J(cols(X),1," "),((st_varname(xindex)))'))
  st_local("N", strofreal(n))

  // Create fixed effects and standard errors.
  seindex = st_addvar("double", sevar)
  st_view(se=.,.,seindex,sample)
  makese(se,PJ,vfx)
  
  //if (xb!="") {
  //  xbindex = st_addvar("double", xbvar)
  //  st_view(xb=.,.,xbindex,sample)
  //  makeresid(X,y,J,PJ,beta,gamma,u,xb)
  //}

}

real vector solvebeta(real matrix X, real vector y, real vector w, real matrix PJ, real scalar nk) {

  real matrix XXdm, xj
  real vector Xydm, yj, wj, beta
  
  XXdm = J(nk-1,nk-1,0)
  Xydm = J(nk-1,1,0)
  for(j=1; j<=rows(PJ); j++){
    xj = panelsubmatrix(X,j,PJ)[,(1::(nk-1))]
	yj = panelsubmatrix(y,j,PJ)
	wj = panelsubmatrix(w,j,PJ)
	XXdm = XXdm + quadcrossdev(xj,mean(xj,wj),wj,xj,mean(xj,wj))
	Xydm = Xydm + quadcrossdev(xj,mean(xj,wj),wj,yj,mean(yj,wj))
  }
  beta = invsym(XXdm)*Xydm
  return(beta)
}

real vector solvegamma(real vector r, real vector w, real matrix PJ, real scalar nj) {
  
  gamma = J(nj+1,1,0)
  
  for(j=1; j<=rows(PJ); j++){
    rj = panelsubmatrix(r,j,PJ)
    wj = panelsubmatrix(w,j,PJ)
    gamma[j] = mean(rj,wj)
  }
  return(gamma)
}

real vector matN(real matrix PJ, real colvector w, real scalar nj)
{
  real vector vecn
  vecn = J(nj,1,0)
  for (i=2; i<=nj+1; i++) {
    vecn[i-1]=quadsum(panelsubmatrix(w,i,PJ))
  }
	return(vecn)
}

real matrix matXF(real matrix X, real matrix PJ, real colvector w, real scalar nj, real scalar nk) {
  real vector xji, wji
  real matrix XF

  XF = J(nk,nj,0)
  for (ji=1; ji<=rows(PJ); ji++) {
    if (ji==1) {
	  xji = panelsubmatrix(X,ji,PJ)
	  wji = panelsubmatrix(w,ji,PJ)
	  XF[(1::nk),(1::nj)] = XF[(1::nk),(1::nj)] :- quadcross(xji,wji)
	}
	else {
	  xji = panelsubmatrix(X,ji,PJ)
	  wji = panelsubmatrix(w,ji,PJ)
	  XF[(1::nk),ji-1] = XF[(1::nk),ji-1] + quadcross(xji,wji)
	}
  }
  return(XF)
}

real matrix matFFFX(real matrix XF, real vector vecn, real scalar n0, real scalar nj, real scalar nk) {

  real scalar mult
  real matrix FFFX

  FFFX = J(nj,nk,0)
  for (i=1; i<=cols(XF); i++){
    FFFX[i,(1::nk)] = (vecn[i]^-1) * XF[(1::nk),i]'
  }
  mult = n0 / (1 + n0 * quadsum(vecn:^-1))
  FFFX = FFFX - mult :* quadcross((vecn:^-1)',quadcross((vecn:^-1),XF'))
  return(FFFX)
}

real vector varmatF(real matrix E, real matrix FFFX, real vector vecn, real scalar n0, real scalar nj, real scalar nk) {

  real vector vfx, rowk
  real scalar prodij, mult
  
  vfx = J(nj+1,1,0)
  mult = n0 / (1 + n0 * quadsum(vecn:^-1))
 
  for(i=1; i<=nj; i++){
    rowk=J(1,nk,0)
	for(k=1; k<=nk; k++){
	  rowk[k] = FFFX[i,(1::nk)] * E[(1::nk),k]
	}
    for(j=1; j<i; j++){
	  prodij = rowk * FFFX[j,(1::nk)]' - mult*(1/vecn[i])*(1/vecn[j])
	  vfx[1] = vfx[1] + 2*prodij
	}
	prodij = rowk * FFFX[i,(1::nk)]' - mult*(1/vecn[i]^2) + (1/vecn[i])
	vfx[i+1] = vfx[i+1] + prodij
	vfx[1] = vfx[1] + prodij
  }
  return(vfx)
}

void makefx(real colvector fe, real matrix PJ, real vector gamma)
{
  real vector feji
  panelsubview(feji=.,fe,1,PJ)
  feji[.] = J(PJ[1,2],1,gamma[1])
  for (j=2; j<=rows(PJ); j++) {
    panelsubview(feji=.,fe,j,PJ)
    feji[.] = J(PJ[j,2]-PJ[j-1,2],1,gamma[j])
  }
}

void makese(real colvector se, real matrix PJ, real vector varfx)
{
  real vector seji
  panelsubview(seji=.,se,1,PJ)
  seji[.] = J(PJ[1,2],1,sqrt(varfx[1]))
  for (j=2; j<=rows(PJ); j++) {
    panelsubview(seji=.,se,j,PJ)
    seji[.] = J(PJ[j,2]-PJ[j-1,2],1,sqrt(varfx[j]))
  }
}

end
