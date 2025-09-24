# Loading the required libraries
library(brglm2)
library(tidyverse)

# Function to check coverage for a model
check_coverage <- function(fit, true_betas) {
ci <- confint.default(fit) # Wald CI
covers <- (true_betas >= ci[,1]) & (true_betas <= ci[,2])
return(as.numeric(covers))
}
# A function to run simulation
run_simulation <- function(N=1000, n=100, beta_0=-1, beta_1=0.5,
m0 = 0, m1 = 0){
# Initialize results data frames
MLE_results <- data.frame(
beta_0 = numeric(N), beta_1 = numeric(N),
cover_beta_0 = numeric(N), cover_beta_1 = numeric(N)
)
MC_results <- data.frame(
beta_0 = numeric(N), beta_1 = numeric(N),
cover_beta_0 = numeric(N), cover_beta_1 = numeric(N)
)
Firth_results <- data.frame(
beta_0 = numeric(N), beta_1 = numeric(N),
cover_beta_0 = numeric(N), cover_beta_1 = numeric(N)
)

log_F_results <- data.frame(
beta_0 = numeric(N), beta_1 = numeric(N),
cover_beta_0 = numeric(N), cover_beta_1 = numeric(N)
)
for (i in 1:N){
# Generate data
df <- generate_data(n=n, beta_0=beta_0, beta_1=beta_1)
# Fit models
fit_MLE <- glm(y ~ X1, family=binomial(link="logit"), data=df)
fit_MC <- glm(y ~ X1, family=binomial(link="logit"), data=df,
method=brglmFit, type="correction")
fit_Firth <- glm(y ~ X1, family=binomial(link="logit"), data=df,
method=brglmFit, type="AS_mean")
#log-F
df$intercept = 1
df$weight = 1
curr_coef = 0
for (m in c(m0, m1)){
df <- add_row(df,
y = 0.5,
intercept = ifelse(curr_coef==0,1,0),
X1 = ifelse(curr_coef==1,1,0),
weight = m)
curr_coef = curr_coef + 1
}
fit_log_F <- glm(y ~ -1 + intercept + X1,
family=binomial(link="logit"),
data=df,
weight = weight)
# Store estimates
MLE_results[i,1:2] <- coef(fit_MLE)
MC_results[i,1:2] <- coef(fit_MC)
Firth_results[i,1:2] <- coef(fit_Firth)
log_F_results[i,1:2] <- coef(fit_log_F)
# Compute Wald CIs and check coverage
true_betas <- c(beta_0, beta_1)
MLE_results[i,3:4] <- check_coverage(fit_MLE, true_betas)
MC_results[i,3:4] <- check_coverage(fit_MC, true_betas)
Firth_results[i,3:4] <- check_coverage(fit_Firth, true_betas)
log_F_results[i,3:4] <- check_coverage(fit_log_F, true_betas)
}
# Return results

return(list(
MLE = MLE_results,
Asymptotic_Bias_Correction = MC_results,
Firth = Firth_results,
log_F = log_F_results
))
}
# A function to compute estimated bias, its CI, and coverage probability
show_results <- function(results, true_betas){
results_MLE <- list(estimated_bias = vector(mode="numeric",
length=length(true_betas)),
bias_se = vector(mode="numeric",
length = length(true_betas)),
bias_CI = list(beta_0 = NA, beta_1 = NA),
coverage_prob = vector(mode="numeric",
length = length(true_betas)))
results_Asymptotic_Bias_Correction <- list(estimated_bias = vector(mode="numeric",
length=length(true_betas)),
bias_se = vector(mode="numeric",
length = length(true_betas)),
bias_CI = list(beta_0 = NA, beta_1 = NA),
coverage_prob = vector(mode="numeric",
length = length(true_betas)))
results_Firth <- list(estimated_bias = vector(mode="numeric",
length=length(true_betas)),
bias_se = vector(mode="numeric",
length = length(true_betas)),
bias_CI = list(beta_0 = NA, beta_1 = NA),
coverage_prob = vector(mode="numeric",
length = length(true_betas)))
results_log_F <- list(estimated_bias = vector(mode="numeric",
length=length(true_betas)),
bias_se = vector(mode="numeric",
length = length(true_betas)),
bias_CI = list(beta_0 = NA, beta_1 = NA),
coverage_prob = vector(mode="numeric",
length = length(true_betas)))
N = length(results$MLE[,1])
# Computing estimated bias for each method
for (i in 1:length(true_betas)){
results_MLE$estimated_bias[i] = mean(results$MLE[,i] - true_betas[i])
results_MLE$bias_se[i] = sd(results$MLE[,i] - true_betas[i])/sqrt(N)
results_MLE$bias_CI[[i]] = c(results_MLE$estimated_bias[i] -
1.96*results_MLE$bias_se[i],
results_MLE$estimated_bias[i] +
1.96*results_MLE$bias_se[i])
results_MLE$coverage_prob[i] = sum(results$MLE[,2+i])/N

results_Asymptotic_Bias_Correction$estimated_bias[i] =
mean(results$Asymptotic_Bias_Correction[,i] -
true_betas[i])
results_Asymptotic_Bias_Correction$bias_se[i] =
sd(results$Asymptotic_Bias_Correction[,i] -
true_betas[i])/sqrt(N)
results_Asymptotic_Bias_Correction$bias_CI[[i]] =
c(results_Asymptotic_Bias_Correction$estimated_bias[i] -
1.96*results_Asymptotic_Bias_Correction$bias_se[i],
results_Asymptotic_Bias_Correction$estimated_bias[i] +
1.96*results_Asymptotic_Bias_Correction$bias_se[i])
results_Asymptotic_Bias_Correction$coverage_prob[i] =
sum(results$Asymptotic_Bias_Correction[,2+i])/N
results_Firth$estimated_bias[i] = mean(results$Firth[,i] - true_betas[i])
results_Firth$bias_se[i] = sd(results$Firth[,i] - true_betas[i])/sqrt(N)
results_Firth$bias_CI[[i]] = c(results_Firth$estimated_bias[i] -
1.96*results_Firth$bias_se[i],
results_Firth$estimated_bias[i] +
1.96*results_Firth$bias_se[i])
results_Firth$coverage_prob[i] = sum(results$Firth[,2+i])/N
results_log_F$estimated_bias[i] = mean(results$log_F[,i] - true_betas[i])
results_log_F$bias_se[i] = sd(results$log_F[,i] - true_betas[i])/sqrt(N)
results_log_F$bias_CI[[i]] = c(results_log_F$estimated_bias[i] -
1.96*results_log_F$bias_se[i],
results_log_F$estimated_bias[i] +
1.96*results_log_F$bias_se[i])
results_log_F$coverage_prob[i] = sum(results$log_F[,2+i])/N
}
return(list(
MLE= results_MLE,
Asymptotic_Bias_Correction = results_Asymptotic_Bias_Correction,
Firth = results_Firth,
log_F = results_log_F
))
}
# A function to run experiments
run_exp <- function(beta_0 = 1, beta_1, m1, n_min, n_max){
n = seq(from=n_min, to=n_max, by=10)
true_betas = c(beta_0, beta_1)
output = list()
for (i in 1:length(n)){
cat("n = ", n[i], "\n")
results = run_simulation(N=1000, n=n[i], beta_0=beta_0, beta_1=beta_1,
m0 = 0, m1 = m1)
output[[i]] = show_results(results, true_betas)
}

return(output)
}
# Running experiments
data <- run_exp(beta_0 = -1, beta_1 = 0.5, m1 = 1, n_min = 20, n_max = 100)
# Extract results for beta_1
extract_results <- function(data, n_values) {
methods <- c("MLE", "Asymptotic_Bias_Correction", "Firth", "log_F")
results_list <- list()
for (i in seq_along(n_values)) {
for (method in methods) {
results_list <- append(results_list, list(
data.frame(
n = n_values[i],
method = method,
bias = data[[i]][[method]]$estimated_bias[2],
lower = data[[i]][[method]]$bias_CI[[2]][1],
upper = data[[i]][[method]]$bias_CI[[2]][2]
)
))
}
}
return(do.call(rbind, results_list))
}
# Generate the data for plotting
n_values <- seq(from = 20, to = 100, by = 10)
plot_data <- extract_results(data, n_values)
plot_data_log <- plot_data %>%
mutate(log_bias = sign(bias) * log1p(abs(bias)),
log_lower = sign(lower) * log1p(abs(lower)),
log_upper = sign(upper) * log1p(abs(upper)))
# plot with confidence bands
ggplot(plot_data, aes(x = n, y = bias, color = method, fill = method)) +
geom_line(size = 1) + # Line for estimated bias
geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
labs(
title = expression("Estimated bias of " * hat(beta)[1] *
"Confidence Bands with n >= 30, " * beta[0] *
" = -1, " * beta[1] * " = 0.5 and m1 = 1"),
x = "Sample Size (n)",
y = "Estimated Bias"
) +
theme_minimal() +
theme(legend.title = element_blank())

# Plot with confidence bands (log scale)
ggplot(plot_data_log, aes(x = n, y = log_bias, color = method, fill = method)) +
geom_line(size = 1) + # Line for estimated bias
geom_ribbon(aes(ymin = log_lower, ymax = log_upper), alpha = 0.2) +
labs(
title = expression("Log Transformed Estimated Bias of " * hat(beta)[1] *
" with Confidence Bands with n >= 30, " * beta[0] *
" = -1, " * beta[1] *
" = 0.5, and m1 = 1"),
x = "Sample Size (n)",
y = "Log Transformed Estimated Bias"
) +
theme_minimal() +
theme(legend.title = element_blank())
# Function to extract coverage probability for beta_1
extract_coverage <- function(data, n_values) {
methods <- c("MLE", "Asymptotic_Bias_Correction", "Firth", "log_F")
results_list <- list()
for (i in seq_along(n_values)) {
for (method in methods) {
results_list <- append(results_list, list(
data.frame(
n = n_values[i],
method = method,
coverage = data[[i]][[method]]$coverage_prob[2]
)
))
}
}
return(do.call(rbind, results_list))
}
# Generate the data for plotting
plot_data <- extract_coverage(data, n_values)
# Plot coverage probability for beta_1
ggplot(plot_data, aes(x = n, y = coverage, color = method)) +
geom_line(size = 1) +
geom_point(size = 2) +
labs(
title = expression("Coverage Probability of " * hat(beta)[1] *
" Across Sample Sizes when " * beta[0] *
" = -1, " * beta[1] *
" = 0.5 and m1 = 1"),
x = "Sample Size (n)",
y = "Coverage Probability",
color = "Method"
) +
geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +

theme_minimal()

