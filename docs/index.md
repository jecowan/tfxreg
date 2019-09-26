
<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>


`tfxreg` is a user-written Stata command to estimate teacher fixed effects and standard errors.

## Installation

`. net install github, from("https://haghish.github.io/github/") `

`. github install jecowan/tfxreg`

## Estimation

We assume a student achievement model of the form
$$ A_{ijt} = X_{ijt} \beta + \psi_{j} + \epsilon_{ijt} $$
where $$X_{ijt}$$ is a vector of observable characteristics of students (such as prior test scores), $$\psi_{j}$$, $$j= 0, 1, ..., J$$, is teacher value-added, and $$\epsilon_{ijt}$$ is a residual. McCaffrey et al. (2012) summarize a number of other routines in Stata that can estimate teacher value-added and standard errors with fixed effects models. Relative to those commands, `tfxreg` accomplishes several objectives:

* For models without a second high-dimensional fixed effect, it avoids constructing and inverting $$J \times J$$ matrices, saving on computation time and memory
* It imposes sum-to-zero constraints, so teacher value-added is interpreted as a deviation from the average teacher in the sample (Mihaly et al., 2010)
* It permits observation weights to support commonly used co-teaching models (Hock & Isenberg, 2017)


Estimation of the coefficients is straightforward using built-in Stata commands (`areg` and `reg`). After estimation, we recover the estimated teacher effects by summing residuals over the set of students assigned a particular teacher (and demeaning):

$$ \hat{\psi}_{j} = \frac{1}{N_{j}}\sum_{i=1}^{N_{j}}(y_{ij} - x_{ij}\hat{\beta}) -
\frac{1}{J+1}\sum_{l=0}^{J}(\bar{y}_{l} - \bar{x}_{l}\hat{\beta}) $$

For the calculation of standard errors, I use the parameterization in Mihaly et al. (2010), which replaces the $$N \times J$$ matrix of teacher indicators

$$ F = (F_{1}, F_{2}, \ldots, F_{J}) $$

with the matrix

$$ \tilde{F} = (F_{1} - F_{0}, F_{2} - F_{0}, \ldots, F_{J} - F_{0}).  $$

The OLS variance estimator is then

$$  
\begin{split}
\hat{\text{var}} \left( \begin{array}{c} \hat{\beta} \\ \hat{\psi} \end{array} \right) &=
\hat \sigma^{2}_{\epsilon} \left( \begin{array}{cc} X'X & X'\tilde{F} \\ \tilde{F}'X & \tilde{F}'\tilde{F} \end{array} \right)^{-1} \\
&=
\hat \sigma^{2}_{\epsilon} A \\
\end{split}
$$

The memory requirements to construct $$A^{-1}$$ can be substantial with the number of teachers typically included in state administrative datasets. To avoid building this matrix, we rely on matrix inversion identities and construct the final variance matrix in separate pieces. By the block inversion identity, the lower right entry of $$A$$ (i.e., the variance of the teacher effects) is

$$ A_{22} = (\tilde{F}'\tilde{F})^{-1} + (\tilde{F}'\tilde{F})^{-1}\tilde{F}'X(X'X - X'\tilde{F}(\tilde{F}'\tilde{F})^{-1}\tilde{F}'X)^{-1}X'\tilde{F}(\tilde{F}'\tilde{F})^{-1} $$

Note that we can write this as

$$ A_{22} = (\tilde{F}'\tilde{F})^{-1} + BCB' $$

where

$$ B = (\tilde{F}'\tilde{F})^{-1}\tilde{F}'X $$

and

$$
\begin{split}
C &= (X'X - X'\tilde{F}(\tilde{F}'\tilde{F})^{-1}\tilde{F}'X)^{-1} \\
	&= (X'X - X'\tilde{F}B)^{-1}
\end{split}
$$

With the exception of $$(\tilde{F}'\tilde{F})^{-1}$$, $$A_{22}$$ can be constructed from $$K \times K$$ and $$K \times J$$ matrices that feasibly can be stored in memory. However, $$\tilde{F}'\tilde{F}$$ has a simple form that facilitates inversion. By an identity in Henderson and Searle (1981),

$$
(\tilde{F}'\tilde{F})^{-1}_{ij} = \begin{cases}
N_{i}^{-1} - \frac{N_0}{N_{i}^{2}(1 + N_0 \sum_{l=1}^{J}N_l^{-1})}
  & \text{if} \; i=j \\
- \frac{N_0}{N_{i}N_{j}(1 + N_0 \sum_{l=1}^{J}N_l^{-1})}
  & \text{if} \; i \neq j
\end{cases}
$$

where $$N_{i}$$ is the number of observations associated with teacher $$i$$. Therefore,

$$ [(\tilde{F}'\tilde{F})^{-1} \tilde{F'}X]_{ij} = N_{i}^{-1} (F_{i} - F_{0})'X_{j} - \sum_{l=1}^{J} \frac{N_0}{N_{i}N_{l}(1 + N_0 \sum_{k=1}^{J} N_k^{-1})} (F_{l} - F_{0})'X_{j} $$

Finally, note that
$$ (BCB')_{ij} = B_{i}CB_{j}', $$
where $$B_{i}$$ is the $$i$$th row of $$B$$. We construct $$A_{22}$$ element-by-element, saving only the matrix sum (required for the variance of $$\psi_{0}$$) and the diagonal entries.

## Suggested Citation

Cowan, J. (2019). tfxreg.ado v1.5. [Software]. Available from https://jecowan.github.io/tfxreg/.

## References

Henderson, H. V., & Searle, S. R. (1981). On deriving the inverse of a sum of matrices. *SIAM Review, 23*(1), 53–60.

Hock, H., & Isenberg, E. (2017). Methods for accounting for co-teaching in value-added models. *Statistics and Public Policy, 4*(1), 1–11.

McCaffrey, D. F., Lockwood, J., Mihaly, K., & Sass, T. R. (2012). A review of Stata routines for fixed effects estimation in normal linear models. *The Stata Journal, 12*(3), 406–432.

Mihaly, K., McCaffrey, D. F., Lockwood, J., & Sass, T. R. (2010). Centering and reference groups for estimates of fixed effects: Modifications to felsdvreg. *The Stata Journal, 10*(1), 82–103.
