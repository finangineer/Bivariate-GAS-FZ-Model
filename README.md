# Bivariate-GAS-FZ-Model

This repository contains the R code used to reproduce the empirical results in the manuscript of Bivariate GAS–FZ models for portfolio VaR/ES forecasting.

## Repository structure

- `R/` : core functions (data loading, loss functions, forecasting, runners, plotting)
- `examples/` : runnable replication scripts (end-to-end examples)
- `README.md` : instructions

## Data

Price data are downloaded from Yahoo Finance via `quantmod::getSymbols()` (no raw data are stored in this repository).
Because Yahoo Finance data can be revised, results may differ slightly over time.

## Requirements

### R packages
The code uses (at minimum):
- `quantmod`, `nloptr`, `zoo`, `rugarch`
- `PerformanceAnalytics` (for `VaRTest`, if not already installed)

Install packages in R:
```r
install.packages(c("quantmod","nloptr","zoo","rugarch","PerformanceAnalytics"))

