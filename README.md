# Reduction-of-Small-Sample-Bias

# Reduction of Small-Sample Bias of GLM Parameter Estimates

A comprehensive study comparing three bias reduction methods for Generalized Linear Models (GLMs) in small-sample settings.

## Overview

This project examines the performance of three bias reduction techniques for Maximum Likelihood Estimation in GLMs:

1. **Asymptotic Bias Correction** - Adjusts MLE estimates using first-order bias approximation
2. **Firth's Method** - Modifies the score function to reduce bias while preserving invariance properties
3. **Log-F Prior Method** - A Bayesian-inspired technique that allows incorporation of prior information

## Authors

- **Annie Yao** - Log-F Prior Method implementation and analysis
- **Ravleen Bajaj**  - Firth's Method implementation and analysis  
- **Samir Arora** - Asymptotic Bias Correction implementation and analysis
- **All authors** - Collaborative simulation study and results analysis

## Repository Structure

```
├── README.md
├── src/
│   ├── asymptotic_bias_correction.R
│   ├── firth_method.R
│   ├── log_f_prior.R
│   └── simulation_study.R
├── data/
│   └── example_datasets.R
├── results/
│   ├── figures/
│   └── simulation_results.csv
├── docs/
│   └── STAT_851_Final_Project_Report.pdf
└── examples/
    ├── asymptotic_example.R
    ├── firth_example.R
    └── log_f_example.R
```

## Key Findings

- **Small samples (n < 20)**: Firth's method and Log-F prior show superior stability compared to asymptotic bias correction
- **Moderate samples (n = 20-50)**: Asymptotic bias correction performs well when MLE converges
- **Large samples (n > 50)**: All methods converge to similar performance
- **Extreme scenarios**: Penalized likelihood methods (Firth, Log-F) maintain robustness when standard MLE fails

## Installation and Dependencies

### Required R Packages

```r
install.packages(c("brglm2", "tidyverse", "ggplot2"))
```

### Load Required Libraries

```r
library(brglm2)
library(tidyverse)
library(ggplot2)
```

## Quick Start Examples

### Asymptotic Bias Correction
```r
source("src/asymptotic_bias_correction.R")

# Fit model with bias correction
fit_corrected <- glm(y ~ X1, 
                    family = binomial(logit), 
                    data = your_data,
                    method = brglm_fit, 
                    type = "correction")
```

### Firth's Method
```r
source("src/firth_method.R")

# Fit model with Firth's bias reduction
fit_firth <- glm(y ~ X1, 
                family = binomial(logit), 
                data = your_data,
                method = brglm_fit, 
                type = "AS_mean")
```

### Log-F Prior
```r
source("src/log_f_prior.R")

# Apply Log-F prior with penalty parameter m=1
fit_logf <- apply_log_f_prior(data = your_data, 
                             formula = y ~ X1, 
                             m_values = c(0, 1))  # m=0 for intercept, m=1 for X1
```

## Simulation Study

Run the complete simulation study:

```r
source("src/simulation_study.R")

# Run simulation for different scenarios
results_easy <- run_simulation_study(scenario = "easy", n_range = c(10, 20, 30, 50, 100))
results_moderate <- run_simulation_study(scenario = "moderate", n_range = c(10, 20, 30, 50, 100))
results_extreme <- run_simulation_study(scenario = "extreme", n_range = c(10, 20, 30, 50, 100))
```

## Method Comparison Summary

| Method | Advantages | Disadvantages | Best Use Case |
|--------|------------|---------------|---------------|
| Asymptotic Bias Correction | Second-order efficiency, straightforward implementation | Requires MLE convergence, loses invariance | Moderate-large samples |
| Firth's Method | Finite estimates, preserves invariance | Fixed penalty, may over-shrink | Small samples, separation issues |
| Log-F Prior | Adjustable penalties, incorporates prior knowledge | Requires prior specification | When prior information available |

## Performance by Sample Size

- **n ≤ 10**: Log-F (with good priors) > Firth > Asymptotic > MLE
- **10 < n ≤ 30**: Firth ≈ Log-F > Asymptotic > MLE  
- **n > 30**: Asymptotic ≈ Firth ≈ Log-F > MLE

## Contributing

This project was completed as part of STAT 851 coursework. Each contributor implemented their assigned method and collaborated on the simulation study and analysis.

## Citation

If you use this code or methodology, please cite this repository:

```
Yao, A., Bajaj, R., & Arora, S. (2025). Reduction of Small-Sample Bias of GLM Parameter Estimates. 
```

## License

This project is available for academic and research purposes.

## Contact

For questions about specific implementations:
- Asymptotic Bias Correction: Samir Arora
- Firth's Method: Ravleen Bajaj  
- Log-F Prior Method: Annie Yao
- Simulation Study: All authors

## Acknowledgments

Special thanks to Dr. Rachel Altman for guiding us through this project and the developers of the `brglm2` package which provided essential functionality for bias-reduced GLM fitting.