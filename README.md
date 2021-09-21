# hetsar

Estimation of heterogeneous spatial autoregressive models.

## Installation

```r
# install release version from CRAN
install.packages("hetsar")

# install development version from GitHub
devtools::install_github("")
```

## Usage

```r
data(hetsarDataDemo)
df_data <- hetsarDataDemo[["data"]]
m_C <- hetsarDataDemo[["network_matrix"]]
m_W <- m_C / rowSums(m_C) # row normalise

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

## References

M. Aquaro, N. Bailey and M. H. Pesaran (2021).
"Estimation and inference for spatial models with heterogeneous coefficients: An application to US house prices".
Journal of Applied Econometrics 36(1): 18-44.
doi:https://doi.org/10.1002/jae.2792
