library(Matrix)
df=iris
set.seed(12334)
n = 5000
C = rnorm(n)
A = rnorm(n)
NDE = 0.5
beta_1 = 0.3
beta_2 = 0.25
theta_1 = 0.5
theta_2 = 0.3
theta_c = 0.21
theta_0 = 0
beta_0 = 0
#CC test
result_list = list()
for(i in 1:1000){
  
  M = beta_1 *A  +beta_2*C + rnorm(n)
  Y = theta_1*A + theta_2*M + theta_c*C + rnorm(n)
  med_model = lm( M ~ A + C)
  out_model = lm(Y ~ A + M + C)
  results = MediationCC(out_model,med_model, A0 = 0, A1 = 1,C = c(0), M = NULL, Interaction = NULL)  
  result_list[[i]] = results[[2]]
}

library(data.table)
final_result = rbindlist(result_list)
NDE = theta_1
NIE = theta_2 * beta_1
TE = NDE + NIE

mean(final_result$NDE) - NDE
sum(final_result$NDE_low <= NDE & final_result$NDE_high >= NDE )/1000
mean(final_result$NIE) - NIE
sum(final_result$NIE_low <= NIE & final_result$NIE_high >= NIE )/1000
mean(final_result$TE) - TE
sum(final_result$TE_low <= TE & final_result$TE_high >= TE )/1000



#CB test
probit = function(x){exp(x)/(1+exp(x))}
result_list = list()
for(i in 1:1000){
  
  M = beta_1 *A  +beta_2*C + rnorm(n)
  P = probit(theta_1*A + theta_2*M + theta_c*C)
  Y = rbinom(n,1,P)
  med_model = lm( M ~ A + C)
  out_model = glm(Y ~ A + M + C, family = binomial())
  results = MediationCB(out_model,med_model, A0 = 0, A1 = 1,C = c(0), M = NULL, Interaction = NULL)  
  result_list[[i]] = results[[2]]
}
OR_NDE = exp(theta_2)
OR_NIE = exp(theta_1 * beta_1)
OR_TE = OR_NDE*OR_NIE
library(data.table)
final_result = rbindlist(result_list)

mean(final_result$OR_NDE) - OR_NDE
sum(final_result$OR_NDE_low <= OR_NDE & final_result$OR_NDE_high >= OR_NDE )/1000
mean(final_result$OR_NIE) - OR_NIE
sum(final_result$OR_NIE_low <= OR_NIE & final_result$OR_NIE_high >= OR_NIE )/1000
mean(final_result$OR_TE) - OR_TE
sum(final_result$OR_TE_low <= OR_TE & final_result$OR_TE_high >= OR_TE )/1000



#BC test
library(data.table)
probit = function(x){exp(x)/(1+exp(x))}
result_list = list()
result_list2 = list()
for(i in 1:1000){
  
  P = probit(beta_1 *A  +beta_2*C )
  M = rbinom(n,1,P)
  Y = theta_1*A + theta_2*M + theta_c*C + rnorm(n)
  med_model = glm( M ~ A + C,family = binomial())
  out_model = lm(Y ~ A + M + C)
  results = MediationBC(out_model,med_model, A0 = 0, A1 = 1,C = c(0), M = NULL, Interaction = NULL)  
  result_list[[i]] = results[[2]]
  result_list2[[i]] = results[[1]]
}
final_result = rbindlist(result_list)
test_result =rbindlist(result_list2)
A1 = 1
A0 = 0
NDE = theta_1
NIE = theta_2* (probit(beta_1*A1+beta_2*0)-probit(beta_1*A0+beta_2*0))
TE = NDE + NIE

mean(final_result$NDE) - NDE
sum(final_result$NDE_low <= NDE & final_result$NDE_high >= NDE )/1000
mean(final_result$NIE) - NIE
sum(final_result$NIE_low <= NIE & final_result$NIE_high >= NIE )/1000
mean(final_result$TE) - TE
sum(final_result$TE_low <= TE & final_result$TE_high >= TE )/1000


mean(test_result$SE_NIE)
sd(final_result$NIE)



#BB test

#CB test
Probit = function(x){exp(x)/(1+exp(x))}
result_list = list()
result_list2 = list()
result_list3 = list()
for(i in 1:1000){
  
  P_M = probit(beta_1 *A  +beta_2*C )
  M = rbinom(n,1,P_M)
  P = probit(theta_1*A + theta_2*M + theta_c*C )
  Y = rbinom(n,1,P)


  data = data.frame(Y = Y, M = M, A = A, C = C)
  regmedint_obj1 <- regmedint(data = data,
                              ## Variables
                              yvar = "Y",
                              avar = "A",
                              mvar = "M",
                              cvar = c("C"),
                              eventvar = NULL,
                              ## Values at which effects are evaluated
                              a0 = 0,
                              a1 = 1,
                              m_cde = 1,
                              c_cond = 3,
                              ## Model types
                              mreg = "logistic",
                              yreg = "logistic",
                              ## Additional specification
                              interaction = FALSE,
                              casecontrol = FALSE)
  
  summary(regmedint_obj1)
  med_model = glm( M ~ A + C, family = binomial())
  out_model = glm(Y ~ A + M + C, family = binomial())
  results = MediationBB(out_model,med_model, A0 = 0, A1 = 1,C = c(0), M = NULL, Interaction = NULL)  
  new_results = MediationBB_true(theta_0,theta_1,theta_2, beta_0, beta_1, beta_2,out_model, med_model, C = c(0))
  result_list[[i]] = results[[2]]
  result_list2[[i]] = results[[1]]
  result_list3[[i]] = new_results
}
OR_NDE = as.numeric((exp(theta_1*A1)*(1+exp(theta_2+theta_3*A1+beta_0+beta_1*A0+ crossprod(beta_2,C))))/
  (exp(theta_1*A0)*(1+exp(theta_2+theta_3*A0+beta_0+beta_1*A0+crossprod(beta_2,C)))))
OR_NIE =   as.numeric(((1+exp(beta_0+beta_1*A0+crossprod(beta_2,C)))*(1+exp(theta_2+theta_3*A1+beta_0+beta_1*A1+crossprod(beta_2,C))))/
                        ((1+exp(beta_0+beta_1*A1+crossprod(beta_2,C)))*(1+exp(theta_2+theta_3*A1+beta_0+beta_1*A0+crossprod(beta_2,C)))
                        ))
OR_TE = OR_NDE*OR_NIE
library(data.table)
final_result = rbindlist(result_list)
test_result =rbindlist(result_list2)
test_se = mean((rbindlist(result_list3)$V1))
mean(final_result$OR_NDE) - OR_NDE
sum(final_result$OR_NDE_low <= OR_NDE & final_result$OR_NDE_high >= OR_NDE )/1000
mean(final_result$OR_NIE) - OR_NIE
sum(final_result$OR_NIE_low <= OR_NIE & final_result$OR_NIE_high >= OR_NIE )/1000
mean(final_result$OR_TE) - OR_TE
sum(final_result$OR_TE_low <= OR_TE & final_result$OR_TE_high >= OR_TE )/1000


mean(test_result$SE_logOR_NIE)
sd(log(final_result$OR_NIE))



