#add a manual intercept, which takes a value of 0 for existing data
df$intercept = 1

#we will fit the glm using the weight method, assign a weight of 1 to existing data
#(or as appropriate based on existing grouping)

df$weight = 1
#select our penalty parameter m
m=3

#add_row function is from tibble/tidyverse
df <- add_row(df,
y=0.5, X1 = 1,
X2 = 0,
intercept = 0,
weight = m)

#y represents the proportion of successes out of weight trials
#fit a glm function with -1, which removes the automatic intercept fit by R

model <- glm(y ~ -1 + intercept + X1 + X2,
family=binomial,
data=df,
weight = weight)
