*! tfxreg v1.4	27apr2017
* Calculates teacher value-added using matrix inversion identities and
* sum-to-zero constraints from felsdvredgm (Mihaly et al 2010).

*03.14.12 Added prediction (xb) capability.
*04.19.12 Fixed abbrevation error on variable parsing.
*08.03.12 Exploit block inversion identities for faster estimation.
*04.27.17 Removed clustering option and updated error reporting.

program tfxreg, eclass sortpreserve
	version 11
	syntax varlist [if] [in], Teacher(varname numeric) ///
    fe(name) se(name) xb(name) ///
		Reff(varname) [ Weights(varname) ROBUST ]
	marksample touse

tempvar sumwt wt sample miss sample  minreff maxreff numg res tempreff
tempname b V var_e var_jhat var_j var_e_g var_jhat_g var_j_g

tokenize `varlist'
local depvar `1'
macro shift
local indepvar `*'

//Verify that data is sorted by group and teacher id.
qui d, varl
local sortlist = r(sortlist)
local sort1: word 1 of `sortlist'
local sort2: word 2 of `sortlist'
if ("`sort1'"!="`reff'") | ("`sort2'"!="`teacher'") {
	sort `reff' `teacher'
}

//Drop existing variables and generate regression sample
novarabbrev {
	cap drop `fe'
	if _rc ==0	di "note: `fe' already exists and will be replaced."
	cap drop `se'
	if _rc ==0	di "note: `se' already exists and will be replaced."
	cap drop `xb'
	if _rc==0 di "note: `xb' already exists and will be replaced."
}

if "`weights'"!=""	qui egen `miss'=rowmiss(`varlist' `teacher' `student' `weights')
else				qui egen `miss'=rowmiss(`varlist' `teacher' `student')
qui gen `sample'=0
qui replace `sample'=1 if `touse' & `miss'==0


//Check that groups nest teacher id.
preserve
qui bysort `teacher': egen `minreff' = min(`reff') if `sample'==1
qui bysort `teacher': egen `maxreff' = max(`reff') if `sample'==1
qui count if `maxreff' != `minreff'
if r(N)>0 {
    di in red "The reference collections (`reff') do not nest the teacher id (`teacher')."
    di in red "Redefine teacher id, so that each value is unique to a reference collection."
    exit(197)
}
restore

//Check for collinear regressors.
_rmcoll `indepvar', noconstant
local indepvar = r(varlist)

//Verify that levels of refference group variable are included as regressors.
qui levelsof `reff', local(lreff)
qui _rmcoll ibn.`reff' `indepvar', noconstant
local nref: word count `lreff'
if r(k_omitted) != `nref' {
	qui _rmcollright ibn.`reff' `indepvar', noconstant
	di in red "Identifiers for each of the reference groups must be included"
	di in red "in the regression for proper interpretation of teacher effects."
  di in red "Often this error message results because you have attempted to"
  di in red "include sets of variables that are collinear with the reference"
  di in red "group."
	di in red "Included reference group indicators:", r(dropped)
	di in red "Stata found these reference group values:", r(block1)
	exit(197)
}

//Rescale so that weights sum to the number of observations.
if "`weights'"!="" {
	qui sum `weights' if `sample'==1
	qui gen `wt'=`weights'/r(mean)
}
else {
	qui gen `wt'=1
}

//Send to Mata
mata: partreg("`depvar'", "`indepvar'", "`teacher'", "`wt'", "`robust'", "`reff'", "`sample'", "`res'", "`fe'", "`se'", "`xb'", "`tcol1'", "`tcol2'")

//Process results
ereturn post `b' `V', depname("`depvar'") obs(`N')

di ""
ereturn display

end


mata:

void partreg(string scalar yvar, string rowvector xvar, string scalar teachid, string scalar weight, string scalar robust, string scalar reffid, string scalar sample, string scalar res, string scalar fevar, string scalar sevar, string scalar xbvar, string scalar tcol1, string scalar tcol2)
{

	real scalar numj, k
	real colvector y, j, w, u, coef, beta, gamma, fe, se, gr2
	string colvector groupstrip
	real matrix A, B, X, varbeta, vargamma, PG, PJ, J

	//Set up views
	jindex 	= st_varindex(teachid)
	gindex	= st_varindex(reffid)
	xindex 	= st_varindex(tokens(xvar))
	yindex 	= st_varindex(yvar)
	windex 	= st_varindex(weight)

	st_view(j=.,.,jindex,sample)
	st_view(g=.,.,gindex,sample)
	st_view(X=.,.,(xindex),sample)
	st_view(y=.,.,yindex,sample)
	st_view(w=.,.,windex,sample)

	//Set up group and teacher indicators
	PG = panelsetup(g,1)
	PJ = panelsetup(j,1)
	J  = maketeach(PG,PJ,j,w)
	G  = panelsetup(J,5)

	k	 = cols(X)
	numj = rows(J) - rows(G)
	n	 = rows(X)

	//Solve for coefficients.
	A 		= crossprodX(X,G,J,PJ,w,numj,k)
	B 		= crossprody(X,y,w,J,PJ,numj,k)
	coef 	= A*B
	beta	= coef[(numj+1::numj+k)]
	gamma	= coef[(1::numj)]


	//Generate residuals, predicted values and calculate standard errors.
	uindex 	= st_addvar("double", res)
	xbindex = st_addvar("double", xbvar)
	st_view(u=.,.,uindex,sample)
	st_view(xb=.,.,xbindex,sample)
	makeresid(X,y,J,PJ,beta,gamma,u,xb)
	if (robust=="") {
		var=olsvar(A,u,w,rows(X),k,numj)
	}
	else if (robust!="") {
		var=hcvar(X,A,w,u,J,PJ,rows(X),k,numj)
	}
	varbeta = var[|numj+1,numj+1 \ numj+k, numj+k|]
	vargamma = var[(1::numj),(1::numj)]

	//Return coefficient vector and variance-covariance matrix to Stata
	st_matrix(st_macroexpand("`"+"b"+"'"),beta[1..cols(X)]')
	st_matrixcolstripe(st_macroexpand("`"+"b"+"'") ,(J(cols(X),1," "),((st_varname(xindex)))'))

	st_matrix(st_macroexpand("`"+"V"+"'"),varbeta)
	st_matrixcolstripe(st_macroexpand("`"+"V"+"'") ,(J(cols(X),1," "),((st_varname(xindex)))'))
	st_matrixrowstripe(st_macroexpand("`"+"V"+"'") ,(J(cols(X),1," "),((st_varname(xindex)))'))

	st_local("N", strofreal(rows(X)))

	//Create fixed effects and standard errors.
	feindex = st_addvar("double", fevar)
	seindex = st_addvar("double", sevar)
	st_view(fe=.,.,feindex,sample)
	st_view(se=.,.,seindex,sample)
	makeeffects(fe,se,J,PJ,gamma,vargamma)

}


real matrix maketeach(real matrix PG, real matrix PJ, real colvector j, real colvector w)
{
	/*
	Links teacher id to teacher FE columns
	Input: panelsetup matrices for groups and teachers.
	Outputs (NJx5) matrix that contains indicator for omitted teacher, column of F matrix, column start & stop (for omitted teachers) and weighted student count.
	*/

	real scalar nf, nj, njg
	real matrix J

	J = J(rows(PJ),5,0)

	nj=1	// Counter for row of J matrix
	nf=1	// Counter for column of F matrix

	for (gi=1; gi<=rows(PG); gi++) {
		njg=rows(uniqrows(panelsubmatrix(j,gi,PG)))

		//Omitted teacher
		J[nj,1] = 1
		J[nj,2] = nf
		J[nj,3] = nf+njg-2
		J[nj,4] = quadsum(panelsubmatrix(w,nj,PJ))
		J[nj,5] = gi
		nj=nj+1

		//Other teachers
		for (ji=2; ji<=njg; ji++) {
			J[nj,2] = nf
			J[nj,4] = quadsum(panelsubmatrix(w,nj,PJ))
			J[nj,5] = gi
			nj=nj+1
			nf=nf+1
		}
	}
	return(J)

}


real matrix crossprodX(real matrix X, real matrix G, real matrix J, real matrix PJ, real colvector w, real scalar nj, real scalar nk) {

	real scalar jf, jmin, jmax
	real vector fjg, fdg, fdf, xji, wji
	real matrix A, E, XX, XF, FF, fg

	XX 	= quadcross(X,w,X)
	XF 	= J(nk,nj,0)
	FF = J(nj,nj,0)

	//Invert F'F
	for (gi=1; gi<=rows(G); gi++) {
		fjg = panelsubmatrix(J,gi,G)
		fdg = fjg[(2::rows(fjg)),4]:^(-1)
		fdf = -1*fjg[1,4]/(1+fjg[1,4]*quadsum(fdg))
		fg = fdf*(fdg*fdg')
		_diag(fg,fdg+diagonal(fg))
		FF[|fjg[1,2],fjg[1,2] \ fjg[1,3],fjg[1,3] |] = fg
	}

	//Generate X'F
	for (ji=1; ji<=rows(J); ji++) {
		if (J[ji,1]==0) {
			xji = panelsubmatrix(X,ji,PJ)
			wji = panelsubmatrix(w,ji,PJ)
			jf = J[ji,2]
			XF[(1::nk),jf] = XF[(1::nk),jf] + quadcross(xji,wji)
		}
		if (J[ji,1]==1) {
			xji = panelsubmatrix(X,ji,PJ)
			wji = panelsubmatrix(w,ji,PJ)
			jmin = J[ji,2]
			jmax = J[ji,3]
			XF[(1::nk),(jmin::jmax)] = XF[(1::nk),(jmin::jmax)] :- quadcross(xji,wji)
		}
	}

	E = invsym(XX - XF*FF*(XF'))
	A = ( FF + FF*(XF')*E*XF*FF, -FF*(XF')*E \ -E*XF*FF, E )

	return(A)
}

real colvector crossprody(real matrix X, real colvector y, real colvector w, real matrix J, real matrix PJ, real scalar nj, real scalar k)
{

	real scalar jf, jmin, jmax
	real vector yji, wji
	real colvector B

	B	= ( J(nj,1,0) \ quadcross(X,w,y) )

	for (ji=1; ji<=rows(J); ji++) {
		if (J[ji,1]==0) {
			yji = panelsubmatrix(y,ji,PJ)
			wji = panelsubmatrix(w,ji,PJ)
			jf = J[ji,2]
			B[jf,1] = B[jf,1] + quadcross(yji,wji)
		}
		if (J[ji,1]==1) {
			yji = panelsubmatrix(y,ji,PJ)
			wji = panelsubmatrix(w,ji,PJ)
			jmin = J[ji,2]
			jmax = J[ji,3]
			B[(jmin::jmax),1] = B[(jmin::jmax),1] :- quadcross(yji,wji)
		}
	}
	return(B)
}

void makeresid(real matrix X, real colvector y, real matrix J, real matrix PJ, real colvector beta, real colvector gamma, real colvector u, real colvector xb)
{

	real scalar f, f1, f2
	real vector yji, uji, xbji
	real matrix xji

	for (ji=1; ji<=rows(J); ji++) {
		yji = panelsubmatrix(y,ji,PJ)
		xji = panelsubmatrix(X,ji,PJ)
		panelsubview(uji=.,u,ji,PJ)
		panelsubview(xbji=.,xb,ji,PJ)
		if (J[ji,1]==0) {
			f = J[ji,2]
			uji[.] = yji - quadcross(xji',beta) :- gamma[f]
			xbji[.] = quadcross(xji',beta) :+ gamma[f]
		}
		if (J[ji,1]==1) {
			f1 = J[ji,2]
			f2 = J[ji,3]
			uji[.] = yji - quadcross(xji',beta) :+ quadsum(gamma[(f1::f2)])
			xbji[.] = quadcross(xji',beta) :- quadsum(gamma[(f1::f2)])
		}
	}
}

real matrix olsvar(real matrix A, real colvector u, real colvector w, real scalar n, real scalar k, real scalar j)
{
	/*
	Calculate OLS standard errors.
	A = inverted cross-product matrix, u = vector of residuals,  w = vector of weights, n = # obs, k = # vars, j = # teachers
	*/

	real scalar sigsq
	real matrix V

	sigsq 	= quadcross(u,w,u)/(n-k-j)
	V 		= sigsq:*A
	return(V)
}

void makeeffects(real colvector fe, real colvector se, real matrix J, real matrix PJ, real vector gamma, real matrix vargamma)
{

	real scalar n, f, f1, f2
	real vector feji, seji


	for (ji=1; ji<=rows(J); ji++) {
		panelsubview(feji=.,fe,ji,PJ)
		panelsubview(seji=.,se,ji,PJ)
		n = rows(feji)
		if (J[ji,1]==0) {
			f = J[ji,2]
			feji[.] = J(n,1,gamma[f])
			seji[.] = J(n,1,sqrt(vargamma[f,f]))
		}
		if (J[ji,1]==1) {
			f1 = J[ji,2]
			f2 = J[ji,3]
			feji[.] = -1*J(n,1,quadsum(gamma[(f1::f2)]))
			seji[.] = J(n,1,sqrt(quadsum(vargamma[|f1,f1 \ f2,f2|])))
		}
	}
}

real matrix hcvar(real matrix X, real matrix A, real colvector w, real colvector u, real matrix J, real matrix PJ, real scalar n, real scalar nk, real scalar nj)
{

	real scalar f, jmin, jmax
	real vector uj, wj
	real matrix xj, M, var

	M = J(nj+nk,nj+nk,0)
	for (i=1; i<=rows(X); i++) {
		M[|nj+1,nj+1 \ nj+nk,nj+nk|] = M[|nj+1,nj+1 \ nj+nk,nj+nk|] + (w[i]:^2):*(u[i]:^2):*quadcross(X[i,.],X[i,.])
	}

	for(j=1; j<=rows(J); j++) {
		uj=panelsubmatrix(u,j,PJ)
		wj=panelsubmatrix(w,j,PJ)
		xj=panelsubmatrix(X,j,PJ)
		if (J[j,1]==0) {
			f=J[j,2]
			M[f,f] = M[f,f] + quadsum((wj:^2):*(uj:^2))
			M[(nj+1::nj+nk),f] = M[(nj+1::nj+nk),f] + quadcolsum((wj:^2):*(uj:^2):*xj)'
			M[f,(nj+1::nj+nk)] = M[f,(nj+1::nj+nk)] + quadcolsum((wj:^2):*(uj:^2):*xj)
		}
		else if (J[j,1]!=0) {
			jmin=J[j,2]
			jmax=J[j,3]
			M[|jmin,jmin \ jmax,jmax|] = M[|jmin,jmin \ jmax,jmax|] :+ quadsum((wj:^2):*(uj:^2))
			M[(jmin::jmax),(nj+1::nj+nk)] = M[(jmin::jmax),(nj+1::nj+nk)] :- quadcolsum((wj:^2):*(uj:^2):*xj)
			M[(nj+1::nj+nk),(jmin::jmax)] = M[(nj+1::nj+nk),(jmin::jmax)] :- quadcolsum((wj:^2):*(uj:^2):*xj)'
		}
	}
	var = (n/(n-nk-nj)):*A*M*A
	return(var)
}



end
