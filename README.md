
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lcc <img src="man/figures/logo.png" align="right" height = 150/>

[![Development](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/lcc)](https://cran.r-project.org/package=lcc)
[![Main
workflow](https://github.com/Prof-ThiagoOliveira/lcc/actions/workflows/main.yml/badge.svg?branch=main)](https://github.com/Prof-ThiagoOliveira/lcc/actions/workflows/main.yml)
[![codecov](https://codecov.io/github/Prof-ThiagoOliveira/lcc/branch/main/graph/badge.svg?token=AAufXfxoWH)](https://codecov.io/github/Prof-ThiagoOliveira/lcc)
[![downloads](https://cranlogs.r-pkg.org/badges/lcc)](https://cranlogs.r-pkg.org/badges/lcc)
[![total downloads](https://cranlogs.r-pkg.org/badges/grand-total/lcc)](https://cranlogs.r-pkg.org/badges/grand-total/lcc)

`lcc` provides estimation routines for longitudinal concordance
correlation (LCC), longitudinal Pearson correlation (LPC), and
longitudinal accuracy (LA) via polynomial mixed-effects regression.
With it you can quantify agreement profiles across methods, browse
numeric and graphical summaries, and extract confidence intervals for
each metric.

Key capabilities include:

- Unified interfaces for fitting agreement models and producing tidy
  summaries.
- Numerical and graphical diagnostics for fitted and sampled
  trajectories.
- Support for balanced or unbalanced designs, variance structures with
  time-dependent weights, and additional covariates in the fixed-effects
  predictor.

The project is maintained by Thiago de Paula Oliveira \[cre, aut\],
Rafael de Andrade Moral \[aut\], John Hinde \[aut\], Silvio Sandoval
Zocchi \[ctb\], and Clarice Garcia Borges Demétrio \[ctb\].

It has been available on CRAN since 2018
(<https://CRAN.R-project.org/package=lcc>). CRAN hosts the recommended
stable release; see [NEWS.md](NEWS.md) for the development change log.

Automated checks for the development branch run on Linux, macOS, and
Windows via GitHub Actions, with code coverage reported through
[Codecov](https://codecov.io/github/Prof-ThiagoOliveira/lcc). To verify
the test suite locally you can run:

``` r
devtools::test()
devtools::check()
```

This github page has its version under development. New functions will
be added as experimental work and, once it is done and running
correctly, we will synchronize the repositories and add it to the CRAN.

We worked hard to release a new stable version allowing users to analyze
data sets, where the objective is studied the extent of the agreement
profile among methods considering time as covariable.

`lcc` comprises a set of functions that help you build and summarise the
fitted model, compute point estimates and bootstrap confidence intervals
for LCC, LPC, and LA, and generate publication-ready plots. Some helper
functions remain internal and are not intended for direct use.

## Bootstrap enhancements

The latest development cycle introduced a richer bootstrap engine for
agreement metrics:

- Confidence intervals now support normal, percentile, and
  bias-corrected and accelerated (BCa) estimators for LCC, LPC, and LA.
- Parallel resampling respects reproducible seeding via `boot.seed` and
  works across every supported bootstrap scheme.
- Degenerate-response safeguards ensure intervals collapse gracefully
  when the observed variability is insufficient for inference.
- Multiple resampling schemes are available through `boot.scheme`,
  covering subject-level, residual-based, semi-parametric, and fully
  parametric workflows.

``` r
fit <- lcc(
  data              = your_data,
  resp              = "response",
  subject           = "subject_id",
  method            = "assay",
  time              = "visit",
  ci                = TRUE,
  ci.method         = "bca",
  boot.scheme       = "np_case_resid_gr",
  nboot             = 1000,
  boot.seed         = 2025,
  numCore           = 4,
  keep.boot.models  = FALSE,
  components        = TRUE
)

summary(fit)
plot(fit)
```

# Installation

## Installed from CRAN:

``` r
install.packages("lcc")
```

## Installed the development version from Github:

``` r
install.packages("devtools")
devtools::install_github("Prof-ThiagoOliveira/lcc")
```

If you use Windows, first install
[Rtools](https://CRAN.R-project.org/bin/windows/Rtools/). If you run into
issues, launch the installer with the *Run as Administrator* option. On
macOS, install Xcode from the App Store.

`lcc` can also be installed by downloading the appropriate files
directly at the CRAN web site and following the instructions given in
the section `6.3 Installing Packages` of the [R Installation and
Administration](https://CRAN.R-project.org/doc/manuals/R-admin.pdf)
manual.

# Longitudinal Concordance Correlation App

We hope you learn more about the LCC using the [LCC
App](https://prof-thiagooliveira.shinyapps.io/lccApp/). The application
helps illustrate how the model parameters influence the agreement
profile over time. Have fun!

# Tutorials

You can read `lcc` tutorials in our PeerJ article
(https://doi.org/10.7717/peerj.9850) or by using the link below:

[LCC paper](https://peerj.com/articles/9850/).
