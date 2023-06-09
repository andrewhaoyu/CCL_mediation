n = 1000
A = rnorm(n)
C = rnorm(n)
alpha = 0.3
alpha_c = 0.5
M = A*alpha+C*alpha_c 
logit_fun = function(x){
  return(exp(x)/(1+exp(x)))
}
beta_0 = 0.1
beta = 0.5
beta_c = 0.3
P_Y = logit_fun(beta_0+beta*M+beta_c*C)
Y = rbinom(n,1,P_Y)

total_model = glm(Y~A+C,family = binomial())
med_model = lm(M~A+C)
out_model = glm(Y~A+M+C, family = binomial())

Mediation(out_model,med_model,total_model)
