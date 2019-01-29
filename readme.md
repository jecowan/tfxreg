**tfxreg** -- Efficiently estimate teacher fixed effects with standard errors.  
v1.5.2
1/29/2019


### Syntax

```stata
tfxreg varlist [if] [in], teacher(varname) fe(name) se(name) ///
[ xb(name) weights(varname) robust ]
```

### Description

tfxreg provides the capability to estimate teacher value-added measures.
tfxreg computes teacher fixed effects and standard errors with sum-to-zero constraints (Mihaly et al., 2010) and permits weights and a limited capacity for clustering in estimation. Estimation exploits block inversion identities to speed estimation and avoid construction of large matrices.

### Options

**teacher**(*varname*) specifies the teacher identifier.

**fe**(*name*) specifies the name of a new variable to store the teacher effects (caution: tfxreg will drop *name* if it already exists).

**se**(*name*) specifies the name of a new variable to store the standard errors of the teacher effects (caution: tfxreg will drop *name* if it already exists).

**xb**(*name*) specifies the optional name of a new variable to store the predicted student outcomes (caution: tfxreg will drop *name* if it already exists).

**weights**(*varname*) specifies the name of a variable that contains the weights attached to each observation. This might be the fraction of a student's time spent with the teacher (e.g., Hock & Isenberg, 2017).

**robust** estimates Hubert-White heteroskedasticity-robust standard errors.

**cluster** estimates cluster-adjusted standard errors using the Moulton correction (Moulton, 1986) suggested by Jacob & Lefgren (2008).


### Example

```stata
tfxreg y x1 x2, teacher(teachid) fe(fe) se(se)
```

### References

Cornelissen, T. (2008).  The Stata command felsdvreg to fit a linear model with two high-dimensional fixed effects. *Stata Journal, 8*, 170-189.

Henderson, H. and Searle, S. (1981). On deriving the inverse of a sum of matrices. *SIAM Review, 23*(1), 53-60.

Hock, H., & Isenberg, E. (2017). Methods for accounting for co-teaching in value-added models. *Statistics and Public Policy, 4*(1), 1–11.

Jacob, B. A., & Lefgren, L. (2008). Can principals identify effective teachers? Evidence on subjective performance evaluation in education. *Journal of Labor Economics, 26*(1), 101–135.

Mihaly, K. et. al. (2010). Centering and reference groups for estimates of fixed effects: Modifications to felsdvreg. *Stata Journal, 10*, 82-103.

Moulton, B. R. (1986). Random group effects and the precision of regression estimates. *Journal of Econometrics, 32*, 385–397.
