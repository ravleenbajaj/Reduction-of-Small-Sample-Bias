library(brglm2)
set.seed(2025)


fit_MLE <- glm(y ~ X1, family=binomial(link="logit"), data=data)
summary(fit_MLE)

fit_Firth <- glm(y ~ X1, family=binomial(link="logit"), data=data,
method=brglmFit, type="AS_mean")
summary(fit_Firth)





