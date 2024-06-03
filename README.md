
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bonmi

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Overview

Matrix completion has attracted attention in many fields, including statistics, applied mathematics, and electrical engineering. Most of the works focus on the independent sampling models under which the observed entries are sampled independently. Motivated by applications in the integration of knowledge graphs derived from multi-source biomedical data such as those from Electronic Health Records (EHR) and biomedical text, we propose the {\bf B}lock-wise {\bf O}verlapping {\bf N}oisy {\bf M}atrix {\bf I}ntegration (BONMI) to treat blockwise missingness of symmetric matrices representing relatedness between entity pairs. Our idea is to exploit the orthogonal Procrustes problem to align the eigenspace of the two sub-matrices, then complete the missing blocks by the inner product of the two low-rank components. 

:triangular_flag_on_post: **TO BE ADDED** :exclamation:

## Installation

Install the released version of bonmi from CRAN:
:triangular_flag_on_post: **TO BE ADDED** :exclamation:

<!-- ``` r -->
<!-- install.packages("bonmi") -->
<!-- ``` -->

Or install the development version from GitHub with:

``` r
install.packages("remotes")
remotes::install_github("celehs/bonmi")
```

## Usage

:triangular_flag_on_post: **TO BE ADDED** :exclamation:

**Example**

``` r
library(bonmi)
  set.seed(1)
  N = 3000
  r = 10
  m = 5
  p = 0.1
  X0 = matrix(rnorm(N*r),nrow=N)
  W0 = X0 %*% t(X0)
  rownames(W0) = colnames(W0) = paste0('code', 1:N)

  W = list()
  for(s in 1:m){
    ids = which(runif(N) < p)
    Ns = length(ids)
    Es = matrix(rnorm(Ns * Ns, sd = s * 0.01), nrow = Ns)
    Es = Es + t(Es)
    Ws = W0[ids,ids] + Es
    W[[s]] = Ws
  }

  Xhat <- bonmi(W, r, weights=NULL, rsvd.use=FALSE)
  codes = rownames(Xhat)
  What = Xhat %*% t(Xhat)

  id = match(codes,rownames(W0))
  Wstar = W0[id, id]

  # BONMI's result
  result = norm(What - Wstar, 'F')/norm(Wstar, 'F')
  
  # test
  expect_equal(result, 0.0037467037)
```

<!-- See the [getting started guide](https://celehs.github.io/bonmi/articles/main.html) -->
<!-- to learn how to use bonmi -->

## Citations

:triangular_flag_on_post: **TO BE ADDED** :exclamation:
