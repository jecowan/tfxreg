*! tfxreg v1.5	29jan2019
* Calculates teacher value-added using matrix inversion identities and
* sum-to-zero constraints from felsdvredgm (Mihaly et al 2010).

*03.14.12 Added prediction (xb) capability.
*04.19.12 Fixed abbrevation error on variable parsing.
*08.03.12 Exploit block inversion identities for faster estimation.
*04.27.17 Removed clustering option and updated error reporting.
*01.02.18 Revised solution of standard errors to avoid constructing JxJ matrix.
*         Deprecated reference collection feature.
*01.29.19 Added clustering subcommand based on Moulton correction.
*         Added fvvarlist support.
*         Switched internal coefficient estimation to areg.

program tfxreg, eclass sortpreserve
	syntax varlist(fv) [if] [in], Teacher(varname) ///
    fe(name) se(name) [ xb(name) Resid(name) cluster(varname) Weights(varname) ]
	marksample touse

tempvar sumwt wt sample miss res tid cid resid_dr resid_r cflag
tempname b V var_e var_jhat var_j var_e_g var_jhat_g var_j_g

tokenize `varlist'
local depvar `1'
macro shift
local indepvar `*'

// Expand factor variables.
local fvops = "`s(fvops)'" == "true" | _caller() >= 11
if `fvops' {
  fvexpand `indepvar'
  local `indepvar' r(varlist)
}

// Generate internal teacher id and sort.
qui egen `tid' = group(`teacher')
sort `tid'
if "`class'" != "" {
  egen `cid' = group(`teacher' `class')
  sort `tid' `class'
}

// Drop existing variables and generate regression sample
novarabbrev {
	cap drop `fe'
	if _rc == 0  di "note: `fe' already exists and will be replaced."
	cap drop `se'
	if _rc == 0  di "note: `se' already exists and will be replaced."
	cap drop `xb'
	if _rc == 0  di "note: `xb' already exists and will be replaced."
}

if "`weights'"!=""	qui egen `miss'=rowmiss(`varlist' `teacher' `weights')
else qui egen `miss'=rowmiss(`varlist' `teacher')
qui gen `sample'=0
qui replace `sample'=1 if `touse' & `miss'==0


// Rescale so that weights sum to the number of observations.
if "`weights'"!="" {
	qui sum `weights' if `sample'==1
	qui gen `wt'=`weights'/r(mean)
}
else {
	qui gen `wt'=1
}

// Obtain coefficients by areg.
if "`cluster'" == "" qui areg `depvar' `indepvar' [aw = `wt'], absorb(`tid')
if "`cluster'" != "" {
  qui areg `depvar' `indepvar' [aw = `wt'], ///
    absorb(`tid') vce(cluster `cluster')
  predict `resid_r', r
  qui icc `resid_r' `cid'
  local rho = r(icc i)
  egen `cflag' = tag(`cid')
  by `tid' `cid': egen `class_size' = total(`sample')
  by `tid': egen `m' = mean(cond(`cflag' == 1, `class_size', .))
  by `tid': egen `vm' = sd(cond(`cflag' == 1, `class_size', .))
  qui replace `vm' = `vm' ^ 2
  drop `resid_r' `cflag' `class_size'
}
predict `resid_dr', dr
if "`xb'" != ""  predict `xb', xb
if "`resid'" != ""  predict `resid', r

// Check for collinear regressors.
_rmcoll `indepvar', noconstant
local indepvar = r(varlist)

// Send to Mata
mata: partreg("`depvar'", "`indepvar'", "`tid'", "`cid'","`wt'", "`robust'", "`sample'", "`fe'", "`se'")

if `cluster' != "" {
  qui replace `se' = `se' * (1 + (`vm'/`m' + `m' - 1) * `rho')
}

// Process results
ereturn post `b' `V', depname("`depvar'") obs(`N')

di ""
di "Variance of classroom effects: `rho'"
di "Variance of teacher effects: `varfx'"
ereturn display

end


mata:

void partreg(string scalar yvar, string rowvector xvar, string scalar teachid, string scalar classid, string scalar weight, string scalar robust, string scalar sample, string scalar fevar, string scalar sevar)
{

  // Set up workspace.
  real scalar n, nj, nk, n0
  real colvector y, j, w, r, beta, gamma, fe, se, vx, vfx
  string rowvector xname
  string colvector groupstrip
  real matrix X, PJ, XX, XF, FFFX, E

  jindex 	= st_varindex(teachid)
  xindex 	= st_varindex(tokens(xvar))
  yindex 	= st_varindex(yvar)
  windex 	= st_varindex(weight)
  st_view(j = ., ., jindex, sample)
  st_view(X = ., ., (xindex), sample)
  st_view(y = ., ., yindex, sample)
  st_view(w = ., ., windex, sample)

  PJ   = panelsetup(j, 1)
  nk   = cols(X)
  nj   = rows(PJ) - 1
  n    = rows(X)

  // Solve for teacher effects.
  gamma = solvegamma(r, w, PJ, nj)
  rmean = mean(gamma)
  gamma = gamma :- rmean

  feindex = st_addvar("double", fevar)
  st_view(fe = ., ., feindex, sample)
  makefx(fe, PJ, gamma)
  r = r - fe :- rmean

  // Estimate variance of coefficients.
  n0   = quadsum(panelsubmatrix(w,1,PJ))
  vecn = matN(PJ, w, nj)
  XX   = quadcross(X, w, X)
  XF   = matXF(X, PJ, w, nj, nk)
  FFFX = matFFFX(XF, vecn, n0, nj, nk)
  E    = invsym(XX - XF * FFFX)
  sig2e= quadcross(r, w, r) / (n - nk - nj)
  vfx  = sig2e * varmatF(E, FFFX, vecn, n0, nj, nk)

  // Create fixed effects and standard errors.
  seindex = st_addvar("double", sevar)
  st_view(se = ., ., seindex, sample)
  makese(se, PJ, vfx)

  varfx = sum(gamma :^ 2) / (rows(gamma - 1))
  st_local("varfx", varfx)

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

  real vector vfx, rowi
  real scalar prodij, mult

  vfx = J(nj+1,1,0)
  mult = n0 / (1 + n0 * quadsum(vecn:^-1))
  for(i=1; i<=nj; i++){
	  rowi = FFFX[i,] * E
    for(j=1; j<i; j++){
	    vfx[1] = vfx[1] + 2*(rowi * FFFX[j,]' - mult*(1/vecn[i])*(1/vecn[j]))
	  }
	  prodij = rowi * FFFX[i,]' - mult*(1/vecn[i]^2) + (1/vecn[i])
	  vfx[i+1] = prodij
	  vfx[1] = vfx[1] + prodij
  }
  return(vfx)
}

real matrix matXgX(real matrix X, real colvector w, real matrix PC, real scalar sig2e, real scalar sig2c, real scalar nk) {
  real vector xl, wl
  real matrix XgX

  XgX = J(nk,nk,0)
  for(l=1; l<=rows(PC); l++){
    xl = panelsubmatrix(X,l,PC)
    wl = panelsubmatrix(w,l,PC) :^ 0.5 // check weights.
    XgX = XgX + quadcross(xl,wl)*quadcross(wl,xl)
  }
  XgX = sig2e*quadcross(X,w,X) + sig2c*XgX
  return(XgX)
}

real matrix matXgF(real matrix X, real matrix XF, real colvector c, real colvector w, real matrix PJ, real scalar sig2e, real scalar sig2c, real scalar nk) {

  real scalar ncl
  real vector xcl, wcl
  real matrix XgF

  XgF = J(nk,nj,0)
  for (j=1; j<=rows(PJ); j++) {
    if (j==1) {
      cj = panelsubmatrix(c,j,PJ)
      xj = panelsubmatrix(X,j,PJ)
      wj = panelsubmatrix(w,j,PJ)
      pcj = panelsetup(cj,1)
      for (l=1; l<=rows(pcj); l++){
        xcl = panelsubmatrix(xj,l,pcj)
        wcl = panelsubmatrix(wj,l,pcj)
        ncl = rows(wcl)
        XgF[(1::nk),(1::nj)] = XgF[(1::nk),(1::nj)] :- ncl*quadcross(xji,wji)
      }
    }
    else {
      cj = panelsubmatrix(c,j,PJ)
      xj = panelsubmatrix(X,j,PJ)
      wj = panelsubmatrix(w,j,PJ)
      pcj = panelsetup(cj,1)
      for (l=1; l<=rows(pcj); l++){
        xcl = panelsubmatrix(xj,l,pcj)
        wcl = panelsubmatrix(wj,l,pcj)
        ncl = rows(wcl)
        XgF[(1::nk),ji-1] = XgF[(1::nk),ji-1] + ncl*quadcross(xji,wji)
      }
    }
  }
  XgF = sig2e*XF + sig2c*XgF
  return(XgF)
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
