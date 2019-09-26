*! tfxreg v1.5	26sep2019
* Calculates teacher value-added using matrix inversion identities and
* sum-to-zero constraints from felsdvredgm (Mihaly et al 2010).

*03.14.12 Added prediction (xb) capability.
*04.19.12 Fixed abbrevation error on variable parsing.
*08.03.12 Exploit block inversion identities for faster estimation.
*04.27.17 Removed clustering option and updated error reporting.
*09.20.18 Revised solution of standard errors to avoid constructing JxJ matrix.
*         Deprecated reference collection feature.
*					Added fvvarlist support.
*					Dropped direct support for residuals and projections, but these are
*					now available using the predict command post-estimation.

program tfxreg, eclass sortpreserve
	syntax varlist(numeric fv) [if] [in], Teacher(varname) ///
    fe(name) se(name) [	Weights(varname) ]

marksample touse

tempvar wt sample tid resid

tokenize `varlist'
local depvar `1'
macro shift
local indepvar `*'


// Generate internal teacher id and sort.
qui egen `tid' = group(`teacher')
sort `tid'

// Drop existing variables.
novarabbrev {
	cap drop `fe'
	if _rc == 0  di "note: `fe' already exists and will be replaced."
	cap drop `se'
	if _rc == 0  di "note: `se' already exists and will be replaced."
}

// Parse variable list.
// Expand factor variables and remove base levels.
local fvops = "`s(fvops)'" == "true" | _caller() >= 11
if `fvops' {
  fvexpand `indepvar'
  local indepvar `r(varlist)'
}
// Drop collinear regressors.
qui _rmcoll `indepvar', noconstant
local indepvar `r(varlist)'

// Estimate regressions.
if "`weights'" == "" {
	qui areg `depvar' `indepvar' if `touse', absorb(`tid')
}
if "`weights'" != "" {
	qui areg `depvar' `indepvar' [aw = `weights'] if `touse', absorb(`tid')
}
qui predict `resid', dr
qui replace `touse' = e(sample)

// Rescale so that weights sum to the number of observations.
if "`weights'" != "" {
	qui sum `weights' if `touse' == 1
	qui gen `wt'=`weights' / r(mean)
}
else {
	qui gen `wt' = 1
}


// Set up Mata workspace.
mata: X = st_data(., "`indepvar'", "`touse'")
mata: X = (J(rows(X), 1, 1), X)
mata: weight = st_data(., "`wt'", "`touse'")
mata: teach = st_data(., "`tid'", "`touse'")
mata: resid = st_data(., "`resid'", "`touse'")
mata: PJ = panelsetup(teach, 1)
mata: nj = rows(PJ)
mata: nk = cols(X)
mata: dof = st_numscalar("e(df_r)")
mata: rss = st_numscalar("e(rss)")

// Solve for teacher effects.
mata: gamma = solve_fe(PJ, resid, weight, nj)

// Estimate standard errors.
mata: vecn  = count_obs(PJ, weight, nj)
mata: XX    = quadcross(X, weight, X)
mata: XF    = make_XF(PJ, X, weight, nj, nk)
mata: FFFX  = make_FFFX(XF, vecn, nj, nk)
mata: C     = invsym(XX - XF * FFFX)
mata: vfx   = make_var(FFFX, C, vecn, nj, rss, dof)

// Return fixed effects and standard errors.
mata: feindex = st_addvar("double", "`fe'")
mata: st_view(fe = ., ., feindex, "`touse'")
mata: seindex = st_addvar("double", "`se'")
mata: st_view(se = ., ., seindex, "`touse'")
mata: return_ests(PJ, gamma, vfx, nj, fe, se)

ereturn display
end


mata:

real matrix solve_fe(real matrix PJ, real colvector resid, real colvector weight, real scalar nj) {
	gamma = J(nj, 1, 0)
	for(j = 1; j <= nj; j++){
		rj = panelsubmatrix(resid, j, PJ)
		wj = panelsubmatrix(weight, j, PJ)
		gamma[j] = mean(rj, wj)
	}
	rmean = mean(gamma)
	gamma = gamma :- rmean
	return(gamma)
}

real matrix count_obs(real matrix PJ, real colvector weight, real scalar nj) {
	vecn = J(nj, 1, 0)
	for (i = 1; i <= nj; i++) {
		vecn[i] = quadsum(panelsubmatrix(weight, i, PJ))
	}
	return(vecn)
}

real matrix make_XF(real matrix PJ, real matrix X, real colvector weight, real scalar nj, real scalar nk) {
	XF = J(nk, nj - 1, 0)
	for (j = 2; j <= nj; j++) {
		XF[, j - 1] = XF[, j - 1] +
			quadcross(panelsubmatrix(X, j, PJ), panelsubmatrix(weight, j, PJ)) -
			quadcross(panelsubmatrix(X, 1, PJ), panelsubmatrix(weight, 1, PJ))
	}
	return(XF)
}

real matrix make_FFFX(real matrix XF, real colvector vecn, real scalar nj, real scalar nk) {
	FFFX = J(nj - 1, nk, 0)
	for (j = 2; j <= nj; j++) {
		FFFX[j - 1, ] = (vecn[j] ^ -1) * XF[, j - 1]'
	}
	mult = vecn[1] / (1 + vecn[1] * quadsum(vecn[2::nj] :^ -1))
	FFFX = FFFX - mult :* quadcross((vecn[2::nj] :^ -1)',
													quadcross((vecn[2::nj] :^ -1), XF'))
	return(FFFX)
}

real matrix make_var(real matrix FFFX, real matrix C, real colvector vecn, real scalar nj, real scalar rss, real scalar dof) {
	vfx = J(nj, 1, 0)
	mult = vecn[1] / (1 + vecn[1] * quadsum(vecn[2::nj] :^ -1))

	for(i = 2; i <= nj; i++){
		rowi = FFFX[(i - 1), ] * C
		for(j = 2; j < i; j++){
			vfx[1] = vfx[1] + 2 * (rowi * FFFX[(j - 1), ]' -
			mult * (1 / vecn[i]) * (1 / vecn[j]))
		}
		prodij = rowi * FFFX[(i - 1), ]' - mult * (1 / vecn[i]^2) + (1 / vecn[i])
		vfx[i] = prodij
		vfx[1] = vfx[1] + prodij
	}
	vfx  = (rss / dof) :* vfx
	return(vfx)
}

void return_ests(real matrix PJ, real colvector gamma, real colvector vfx, real scalar nj, fe, se) {
	panelsubview(feji = ., fe, 1, PJ)
	feji[.] = J(PJ[1, 2], 1, gamma[1])
	for (j = 2; j <= nj; j++) {
		panelsubview(feji = ., fe, j, PJ)
		feji[.] = J(PJ[j, 2] - PJ[j - 1, 2], 1, gamma[j])
	}

	panelsubview(seji = ., se, 1, PJ)
	seji[.] = J(PJ[1, 2], 1, sqrt(vfx[1]))
	for (j = 2; j <= nj; j++) {
		panelsubview(seji = ., se, j, PJ)
		seji[.] = J(PJ[j, 2] - PJ[j - 1, 2], 1, sqrt(vfx[j]))
	}
}

end
