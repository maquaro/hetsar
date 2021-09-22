# hetsar: R package to estimate spatial autoregressive models with heterogeneous coefficient

Estimation and inference of heterogeneous spatial autoregressive (HSAR) panel
data models. These are models where the spatial lag coefficients are allowed to
differ over the cross-section units, in addition to the fixed effects generally
allowed for in the literature. The model also features weakly exogenous
regressors, such as lagged values of the dependent variable and heteroskedastic
error. The estimation is performed via quasi maximum-likelihood. See Aquaro,
Bailey and Pesaran (2021) for technical details.

## Installation

```r
# install release version from CRAN
install.packages("hetsar")

# install development version from GitHub
devtools::install_github("maquaro/hetsar")
```

## Usage

```r
library("hetsar")

data(hetsarDataDemo)
df_data <- hetsarDataDemo[["data"]]
m_C <- hetsarDataDemo[["network_matrix"]]

# row normalise
m_W <- m_C / rowSums(m_C)

# y = Wy + 1 + x + Wx + err
# other model specifications available, see `?hetsar`
out <- hetsar(
  data = df_data, 
  W = m_W,
  indices = c("id", "time"),
  dependent = "y",
  explanatory = "x",
  space_lags = "x")

summary(out)
summary(out, MG = TRUE)
```

## Vignettes

Coming soon.

## Other software

**hetsar** is also available in
MATLAB (coming soon),
[Python](https://pypi.org/project/hetsar/), and
[Stata](https://ideas.repec.org/c/boc/bocode/s458926.html).

## References

M. Aquaro, N. Bailey and M. H. Pesaran (2021).
"Estimation and inference for spatial models with heterogeneous coefficients: An application to US house prices".
Journal of Applied Econometrics 36(1): 18-44.
doi:https://doi.org/10.1002/jae.2792
