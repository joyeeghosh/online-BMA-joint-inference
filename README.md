# Online BMA Joint Inference

This repository contains the learning analytics data and R code used for the real data analysis in the paper ONLINE BAYESIAN MODEL AVERAGING WITH JOINT UNCERTAINTY QUANTIFICATION FOR MODELS AND REGRESSION COEFFICIENTS IN BINARY REGRESSION by Joyee Ghosh and Aixin Tan.

## Files

* `data.train.RData`: training data used in the learning analytics analysis.
* `data.pred.RData`: test data used for prediction in the learning analytics analysis.
* `allmodels.R`: R function for enumerating the model space.
* `junyi-analysis-github.R`: R code for running the online and offline BMA algorithms, followed by analyses including posterior prediction, credible intervals, AUC calculation, running time comparison, and top model analyses.

## Running the code

The main analysis file is

`junyi-analysis-github.R`

The default setting uses a short 10 batch run for a quick demonstration (the first batch can take a few minutes):

```r
n <- c(304, rep(1, 9))
```

To reproduce the full analysis reported in the paper, replace this by

```r
n <- c(304, rep(1, 1500))
```

which gives 1501 batches in total.

The section containing the top model analysis is disabled by default because it requires output from batch 1501. It should be run only after completing the full 1501 batch analysis.

## R packages

The analysis uses the following R packages:

* `MASS`
* `devtools`
* `RenewGLM`
* `pROC`
* `plotrix`

The `RenewGLM` package is available from GitHub at:

`https://github.com/luolsph/RenewGLM_pkg`

The main analysis script contains modified versions of selected `RenewGLM` functions to accommodate sequential updates with batch size 1.
