**tfxreg** -- Efficiently estimate teacher fixed effects with standard errors.  
v1.5
9/26/2019


### Syntax

```stata
tfxreg varlist [if] [in], teacher(varname) fe(name) se(name) [ weights(varname) ]
```

### Description

tfxreg provides the capability to estimate teacher value-added measures.
tfxreg computes teacher fixed effects and standard errors with sum-to-zero constraints (Mihaly et al., 2010). Estimation exploits block inversion identities to speed estimation and avoid construction of large matrices.

### Options

**teacher**(*varname*) specifies the teacher identifier.

**fe**(*name*) specifies the name of a new variable to store the teacher effects (caution: tfxreg will drop *name* if it already exists).

**se**(*name*) specifies the name of a new variable to store the standard errors of the teacher effects (caution: tfxreg will drop *name* if it already exists).

**weights**(*varname*) specifies the name of a variable that contains the weights attached to each observation. This might be the fraction of a student's time spent with the teacher (e.g., Hock & Isenberg, 2017).


### Example

```stata
tfxreg y x1 x2, teacher(teachid) fe(fe) se(se)
```

### References

Cornelissen, T. (2008).  The Stata command felsdvreg to fit a linear model with two high-dimensional fixed effects. *Stata Journal, 8*, 170-189.

Henderson, H. and Searle, S. (1981). On deriving the inverse of a sum of matrices. *SIAM Review, 23*(1), 53-60.

Hock, H., & Isenberg, E. (2017). Methods for accounting for co-teaching in value-added models. *Statistics and Public Policy, 4*(1), 1–11.

Mihaly, K. et. al. (2010). Centering and reference groups for estimates of fixed effects: Modifications to felsdvreg. *Stata Journal, 10*, 82-103.
