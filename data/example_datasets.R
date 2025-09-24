library(brglm2)
set.seed(2025)

generate_data <- function(n, beta_0, beta_1){
X1 <- rnorm(n, mean = 0, sd = 2)
X <- cbind(1, X1)
beta_true <- c(beta_0, beta_1)
eta <- X %*% beta_true
p <- 1 / (1 + exp(-eta))
y <- rbinom(n, size = 1, prob = p)
df <- data.frame(
y = y,
X1 = X1
)
  
return(df)
}
data <- generate_data(50, 0.5, -0.3)

head(data)
