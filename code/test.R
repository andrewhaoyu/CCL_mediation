n = 1000
C = rnorm(n)
X = rnorm(n, mean = 0, sd = 1)
M = 0.3*X + 0.5* C  + rnorm(n)
score = -3+M*0.2+0.1*X + 0.3* C 
P = exp(score)/(1+exp(score))
Y = rbinom(P,size = 1,prob = P)
sum(Y)
data = data.frame(Y, X, M, C)
out_model = glm(Y ~ X + M + C, data = data, family = binomial(link = 'logit'))
med_model = lm(M~X + C, data = data)
total_model = glm(Y ~ X + C, data = data, family = binomial(link = 'logit'))
summary(out_model)
summary(med_model)
summary(total_model) 
library(mediation)
fit_model = mediate(med_model, out_model, treat="X", mediator= "M", 
                    #outcome = c("censor_days_cancer_ignore", "case_control_cancer_ignore"),
                    boot=T, boot.ci.type = "bca")
summary(fit_model)

install.packages("lavaan")
