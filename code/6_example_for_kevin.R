library(Matrix)
set.seed(12334)
n = 1000

C = rnorm(n)
A = rnorm(n)
beta_1 = 0.3
#beta_2 is beta_c
beta_2 = 0.25
theta_1 = 0.5
theta_2 = 0.3
theta_c = 0.21


#in theory
#theta_1 is the direct effect (0.5)
#theta_2*beta_1 is the indirect effect (0.09)
#theta_1+theta_2*beta_1 is the total effect (0.59)

library(regmedint)
probit = function(x){exp(x)/(1+exp(x))}
# result_list = list()
# result_list2 = list()
# result_list3 = list()
#for(i in 1:1000){
  
  #P_M = probit(beta_1 *A  +beta_2*C )
  #M = rbinom(n,1,P_M)
  M = beta_1 *A  +beta_2*C + rnorm(n)
  P = probit(theta_1*A + theta_2*M + theta_c*C )
  Y = rbinom(n,1,P)
  data = data.frame(Y = Y, M = M, A = A, C = C)
  #standardized the PRS using scale()
  regmedint_obj1 <- regmedint(data = data,
                              ## Variables
                              yvar = "Y",
                              avar = "A",
                              mvar = "M",
                              #age, PC1-PC10
                              cvar = c("C"),
                              eventvar = NULL,
                              ## Values at which effects are evaluated
                              a0 = 0,
                              a1 = 1,
                              m_cde = 0,
                              #c_cond means the baseline covariates value
                              #if the covariate is continous, use mean(covariate) as baseline
                              #if the covariate is binary, use the reference group
                              c_cond = mean(C),
                              ## Model types
                              #note if the mediator is binary use "logistic"
                              #if the mediator is continuous, use "linear"
                              mreg = "linear",
                              yreg = "logistic",
                              ## Additional specification
                              interaction = FALSE,
                              casecontrol = FALSE)
  
  summary(regmedint_obj1)
  
#}
#In the result
#row 2 pnde is the direct effect
#row 3 tnie is the indirect effect
#row 6 te is the total effect
#row 7 pm is the proportion of indirect effect over total effect
#covariates to adjust for Devika: age, PC1-10
#covariates to adjust for Kevin: age, sex, PC1-10
pfun = function(z){
  2*pnorm(-abs(z),lower.tail = T)
}
pfun(5.164887)
