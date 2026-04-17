# opGMMassessment

Optimized Automated Gaussian Mixture Assessment

## Overview

`opGMMassessment` provides optimized automated evaluation of the number and parameters of Gaussian mixtures in one-dimensional data. The package implements various methods for parameter estimation and for determining the number of modes in Gaussian mixture models (GMMs).

## Motivation

Gaussian mixture models are probabilistic models commonly used in biomedical research to detect subgroup structures in data sets with one-dimensional information. Reliable model parameterization requires that the number of modes (i.e., states of the generating process) is known. However, this is rarely the case for empirically measured biomedical data.

## Key Features

- Multiple algorithms for automatic GMM parameterization
- Various methods for mode number determination including:
  - Likelihood Ratio (LR) test
  - Akaike Information Criterion (AIC)
  - Bayesian Information Criterion (BIC)
  - Gap statistic (GAP)
  - Silhouette Index (SI)
  - NbClust methods
  - Figueiredo-Jain (FM) method

- Multiple GMM parameter estimation methods:
  - Markov Chain Monte Carlo (MCMC)
  - ClusterR GMM
  - Density-based Model-based Clustering (densityMclust)
  - Distribution Optimization (DO)
  - EM algorithm (normalmixEM)

- Parallel processing support for efficient computation
- Automatic calculation of Bayes decision boundaries
- Visualization with ggplot2
- Kolmogorov-Smirnov test for goodness of fit

## Installation
### Install from CRAN

```r
install.packages("opGMMassessment")
```

## Or install development version from GitHub
### Newer version opGMMassessment_0.4.5

```r
devtools::install_github("username/opGMMassessment")
```

## Usage

```r
library(opGMMassessment)

# Basic usage with default parameters (MCMC + Likelihood Ratio)
result <- opGMMassessment(Data = your_data)

# With custom parameters
result <- opGMMassessment(
  Data = your_data,
  FitAlg = "MCMC",           # Algorithm for GMM fitting
  Criterion = "LR",           # Criterion for mode number determination
  MaxModes = 8,               # Maximum number of modes to test
  MaxCores = 2,               # Number of cores for parallel processing
  PlotIt = TRUE,              # Display the plot
  KS = TRUE                   # Perform Kolmogorov-Smirnov test
)

# Access results
classes <- result$Cls          # Cluster assignments
means <- result$Means          # Means of Gaussian components
sds <- result$SDs              # Standard deviations
weights <- result$Weights      # Mixture weights
boundaries <- result$Boundaries # Bayes decision boundaries
plot <- result$Plot            # ggplot2 object
ks_test <- result$KS           # KS test results
```

## Performance

Based on comparative evaluation on artificial and real-world data sets, the combination of the **Likelihood Ratio (LR) test** for mode number determination and **Markov Chain Monte Carlo (MCMC)** for parameter estimation consistently performs best and outperforms other available implementations.

## Citation

If you use this package in your research, please cite:

Lötsch, J., Malkusch, S., & Ultsch, A. (2022). Comparative assessment of automated algorithms for the separation of one-dimensional Gaussian mixtures. *Informatics in Medicine Unlocked*, 34, 101113. https://doi.org/10.1016/j.imu.2022.101113

## License

GPL-3

## References

Lötsch, J., Malkusch, S., & Ultsch, A. (2022). Comparative assessment of automated algorithms for the separation of one-dimensional Gaussian mixtures. *Informatics in Medicine Unlocked*, 34, 101113.
