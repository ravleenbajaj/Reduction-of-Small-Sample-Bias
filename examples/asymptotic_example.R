# Installing the required libraries
install.packages("brglm2")

# Loading the required packages
library(brglm2)

# Fitting the asymptotic bias corrected binomial glm
fit_binomial_bias_corrected <- glm(cbind(y, n-y) ~ Grade + Sex + Participate,
binomial(logit), data = school_binomial,
method = brglm_fit, type = "correction")

# Getting the summary output
summary(fit_binomial_bias_corrected)
