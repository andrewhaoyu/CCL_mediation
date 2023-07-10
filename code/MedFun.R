#Mediation function without interaction
library(Matrix)
Probit = function(x){
  as.numeric(exp(x)/(1+exp(x)))
}
#A is the exposure
#M is the mediator
#Y is the outcome
#C are the covariates
#out_model is Y ~ A + M + C
#med_model is Y ~ M + C
#if Y/M is continous, then the model is linear model
#if Y/M is binary, then the model is logit model
#Binary Mediator Binary Outcome
#A0 is the fixing value of exposure in the control group 
#A1 is the fixing value of exposure in the case group
#The averaged natural direct effect is defined as E(Y|M_{A0},A1) - E(Y|M_{A0},A0)
#The averaged natural indirect effect is defined as E(Y|M_{A1},A1) - E(Y|M_{A0},A1)
#The averaged total effect is defined as E(Y|M_{A1},A1) - E(Y|M_{A0},A0)
#C here is the controlled value for other covariates.

MediationBB_true = function(theta_0,theta_1,theta_2, beta_0, beta_1, beta_2,out_model, med_model , C ){
 
  theta_3 = 0
  cov_theta <- vcov(out_model)
  
  
  cov_beta <- vcov(med_model)
  Sigma = bdiag(list(cov_beta,cov_theta))
  A = Probit(theta_2+theta_3*A1+beta_0+beta_1*A1+crossprod(beta_2,C))
  B = Probit(theta_2+theta_3*A1+beta_0+beta_1*A0+crossprod(beta_2,C))
  K = Probit(beta_0+beta_1*A1+crossprod(beta_2,C))
  D = Probit(beta_0+beta_1*A0+crossprod(beta_2,C))
  b1 = (D+A)-(K+B)
  b2 = A0*(D-B)+A1*(A-K)
  b3 = (D+A-K-B)*C
  b4 = 0
  b5 = 0
  b6 = as.numeric(A-B)
  #b7 = A1*(as.numeric(A-B))
  b8 = 0*C
  #Gamma = c(b1,b2,b3,b4,b5,b6,b7,b8)
  Gamma = c(b1,b2,b3,b4,b5,b6,b8)
  SE_logOR_NIE = as.numeric(sqrt(t(Gamma)%*%Sigma%*%Gamma))
  
  return(list(SE_logOR_NIE))
  
  
}




MediationBB = function(out_model,med_model, A0, A1,C, M = NULL, Interaction = NULL){
  out_coef = coef(out_model)
  theta_0 = out_coef[1]
  theta_1 = out_coef[2]
  theta_2 = out_coef[3]
  
  #since no interaction is assumed in the model
  #theta_3 = 0
  theta_3 = 0
  cov_theta <- vcov(out_model)
  
  med_coef = coef(med_model)
  beta_0 = med_coef[1]
  beta_1 = med_coef[2]
  beta_2 = med_coef[3:length(med_coef)]
  
  cov_beta <- vcov(med_model)
  
  OR_NDE = (exp(theta_1*A1)*(1+exp(theta_2+theta_3*A1+beta_0+beta_1*A0+ crossprod(beta_2,C))))/
              (exp(theta_1*A0)*(1+exp(theta_2+theta_3*A0+beta_0+beta_1*A0+crossprod(beta_2,C))))
  logOR_NDE = log(OR_NDE)
  OR_NIE =  as.numeric(((1+exp(beta_0+beta_1*A0+crossprod(beta_2,C)))*(1+exp(theta_2+theta_3*A1+beta_0+beta_1*A1+crossprod(beta_2,C))))/
                         ((1+exp(beta_0+beta_1*A1+crossprod(beta_2,C)))*(1+exp(theta_2+theta_3*A1+beta_0+beta_1*A0+crossprod(beta_2,C)))
                         ))
  logOR_NIE = log(OR_NIE)          
            
  OR_TE = OR_NDE*OR_NIE
  logOR_TE = log(OR_TE)
  #proportion of indirect effect over total effect
  PM = logOR_NIE/logOR_TE
  
  #covariance matrix of the regression coefficients from two models
  Sigma = bdiag(list(cov_beta,cov_theta))
  
  A = Probit(theta_2+theta_3*A1+beta_0+beta_1*A0+crossprod(beta_2,C))
  B = Probit(theta_2+theta_3*A0+beta_0+beta_1*A0+crossprod(beta_2,C))
  d1 = as.numeric(A-B)
  d2 = A0*(as.numeric(A-B))
  d3 = as.numeric(A-B)*C
  d4 = 0
  d5 = A1-A0
  d6 = as.numeric(A-B)
  #d7 = A1*A-A0*B
  d8 = 0*C
  #Gamma = c(d1,d2,d3,d4,d5,d6,d7,d8)
  #we have no interaction in this model
  #so we don't include d7 in the model
  Gamma = c(d1,d2,d3,d4,d5,d6,d8)
  SE_logOR_NDE = as.numeric(sqrt(t(Gamma)%*%Sigma%*%Gamma))
  OR_NDE_low = exp(logOR_NDE-1.96*SE_logOR_NDE)
  OR_NDE_high = exp(logOR_NDE+1.96*SE_logOR_NDE)
  NDE_p = 2*pnorm(-abs(logOR_NDE/SE_logOR_NDE), lower.tail = T)
  
  A = Probit(theta_2+theta_3*A1+beta_0+beta_1*A1+crossprod(beta_2,C))
  B = Probit(theta_2+theta_3*A1+beta_0+beta_1*A0+crossprod(beta_2,C))
  K = Probit(beta_0+beta_1*A1+crossprod(beta_2,C))
  D = Probit(beta_0+beta_1*A0+crossprod(beta_2,C))
  b1 = (D+A)-(K+B)
  b2 = A0*(D-B)+A1*(A-K)
  b3 = (D+A-K-B)*C
  b4 = 0
  b5 = 0
  b6 = as.numeric(A-B)
  #b7 = A1*(as.numeric(A-B))
  b8 = 0*C
  #Gamma = c(b1,b2,b3,b4,b5,b6,b7,b8)
  Gamma = c(b1,b2,b3,b4,b5,b6,b8)
  SE_logOR_NIE = as.numeric(sqrt(t(Gamma)%*%Sigma%*%Gamma))
  OR_NIE_low = exp(logOR_NIE-1.96*SE_logOR_NIE)
  OR_NIE_high = exp(logOR_NIE+1.96*SE_logOR_NIE)
  NIE_p = 2*pnorm(-abs(logOR_NIE/SE_logOR_NIE), lower.tail = T)
  
  l1 = d1 + b1
  l2 = d2 + b2
  l3 = d3 + b3
  l4 = d4 + b4
  l5 = d5 + b5
  l6 = d6 + b6
  #l7 = d7 + b7
  l8 = d8 + b8
  Gamma = c(l1,l2,l3,l4,l5,l6,l8)
  SE_logOR_TE = as.numeric(sqrt(t(Gamma)%*%Sigma%*%Gamma))
  OR_TE_low = exp(logOR_TE-1.96*SE_logOR_TE)
  OR_TE_high = exp(logOR_TE+1.96*SE_logOR_TE)
  TE_p = 2*pnorm(-abs(logOR_TE/SE_logOR_TE), lower.tail = T)
  result_core = data.frame(logOR_NIE,SE_logOR_NIE,NIE_p,
                           logOR_NDE,SE_logOR_NDE,NDE_p,
                           logOR_TE,SE_logOR_TE,TE_p)
    
    
    result_summary = data.frame(OR_NDE,OR_NDE_low,OR_NDE_high,NDE_p,
                      OR_NIE,OR_NIE_low,OR_NIE_high,NIE_p,
                      OR_TE,OR_TE_low,OR_TE_high,TE_p,PM)
    
    result = list(result_core,
                  result_summary)
    
  
  
}









#Continuous Mediator Binary Outcome
MediationCB = function(out_model,med_model, A0, A1,C, M = NULL, Interaction = NULL){
  out_coef = coef(out_model)
  theta_0 = out_coef[1]
  theta_1 = out_coef[2]
  theta_2 = out_coef[3]
  theta_4 = out_coef[4:length(out_coef)]
  #since no interaction is assumed in the model
  #theta_3 = 0
  theta_3 = 0
  cov_theta <- vcov(out_model)
  
  med_coef = coef(med_model)
  beta_0 = med_coef[1]
  beta_1 = med_coef[2]
  beta_2 = med_coef[3:length(med_coef)]

  cov_beta <- vcov(med_model)
  #estimate error variance
  rss <- sum(resid(med_model)^2)
  df = med_model$df.residual
  var_error = rss/df
  
  OR_NDE = exp(
    (theta_1+theta_3*(beta_0+beta_1*A0+crossprod(C,beta_2)+theta_2*var_error))*(A1-A0)+
      0.5*theta_3^2*(A1^2-A0^2))
  logOR_NDE = log(OR_NDE)
  
  OR_NIE = exp((theta_2*beta_1+theta_3*beta_1*A1)*(A1-A0))
  logOR_NIE = log(OR_NIE)          
  
  OR_TE = OR_NDE*OR_NIE
  logOR_TE = log(OR_TE)
  #proportion of indirect effect over total effect
  PM = logOR_NIE/log(OR_TE)
  #covariance matrix of the regression coefficients from two models
  Sigma = bdiag(list(cov_beta,cov_theta,var_error))
  
  d1 = theta_3
  d2 = theta_3*A0
  d3 = (theta_3)*C
  d4 = 0
  d5 = 1
  d6 = theta_3*var_error
  #d7 = beta_0 + beta_1 * A0 + beta_2*C + theta_2 * var_error + theta_3 * var_error * (A0 + A1)
  d8 = 0* C
  d9 = theta_2 * theta_3 + 0.5 * theta_3^2 *( A0 + A1)
  #Gamma = c(d1,d2,d3,d4,d5,d6,d7,d8)
  #we have no interaction in this model
  #so we don't include d7 in the model
  Gamma = c(d1,d2,d3,d4,d5,d6,d8,d9)
  SE_logOR_NDE = as.numeric(sqrt(t(Gamma)%*%Sigma%*%Gamma*(A1-A0)^2))
  OR_NDE_low = exp(logOR_NDE-1.96*SE_logOR_NDE)
  OR_NDE_high = exp(logOR_NDE+1.96*SE_logOR_NDE)
  NDE_p = 2*pnorm(-abs(logOR_NDE/SE_logOR_NDE), lower.tail = T)
  
  b1 = 0
  b2 = theta_2 + theta_3 * A1
  b3 = 0*C
  b4 = 0
  b5 = 0
  b6 = beta_1
  #b7 = beta_1 * A1
  b8 = 0*C
  b9 = 0
  #Gamma = c(b1,b2,b3,b4,b5,b6,b7,b8)
  Gamma = c(b1,b2,b3,b4,b5,b6,b8,b9)
  SE_logOR_NIE = as.numeric(sqrt(t(Gamma)%*%Sigma%*%Gamma*(A1-A0)^2))
  OR_NIE_low = exp(logOR_NIE-1.96*SE_logOR_NIE)
  OR_NIE_high = exp(logOR_NIE+1.96*SE_logOR_NIE)
  NIE_p = 2*pnorm(-abs(logOR_NIE/SE_logOR_NIE), lower.tail = T)
  
  l1 = d1 + b1
  l2 = d2 + b2
  l3 = d3 + b3
  l4 = d4 + b4
  l5 = d5 + b5
  l6 = d6 + b6
  #l7 = d7 + b7
  l8 = d8 + b8
  l9 = d9 + b9
  Gamma = c(l1,l2,l3,l4,l5,l6,l8,l9)
  SE_logOR_TE = as.numeric(sqrt(t(Gamma)%*%Sigma%*%Gamma*(A1-A0)^2))
  OR_TE_low = exp(logOR_TE-1.96*SE_logOR_TE)
  OR_TE_high = exp(logOR_TE+1.96*SE_logOR_TE)
  TE_p = 2*pnorm(-abs(logOR_TE/SE_logOR_TE), lower.tail = T)
  result_core = data.frame(logOR_NIE,SE_logOR_NIE,NIE_p,
                           logOR_NDE,SE_logOR_NDE,NDE_p,
                           logOR_TE,SE_logOR_TE,TE_p)
  
  
  result_summary = data.frame(OR_NDE,OR_NDE_low,OR_NDE_high,NDE_p,
                              OR_NIE,OR_NIE_low,OR_NIE_high,NIE_p,
                              OR_TE,OR_TE_low,OR_TE_high,TE_p,PM)
  
  result = list(result_core,
                result_summary)
  
  
  
}








#Continuous Mediator Continuous Outcome
MediationCC = function(out_model,med_model, A0, A1,C, M = NULL, Interaction = NULL){
  out_coef = coef(out_model)
  theta_0 = out_coef[1]
  theta_1 = out_coef[2]
  theta_2 = out_coef[3]
  
  #since no interaction is assumed in the model
  #theta_3 = 0
  theta_3 = 0
  cov_theta <- vcov(out_model)
  
  med_coef = coef(med_model)
  beta_0 = med_coef[1]
  beta_1 = med_coef[2]
  beta_2 = med_coef[3:length(med_coef)]
  
  cov_beta <- vcov(med_model)
  #estimate error variance
  rss <- sum(resid(med_model)^2)
  df = med_model$df.residual
  var_error = rss/df
  
  NDE = (theta_1 +theta_3*beta_0+theta_3*beta_1*A0+theta_3*crossprod(beta_2, C))*(A1 - A0)
    
  
  NIE = (theta_2*beta_1 + theta_3*beta_1*A1)*(A1-A0)
  
  
  TE = NDE + NIE
  
  #proportion of indirect effect over total effect
  PM = NIE/TE
  #covariance matrix of the regression coefficients from two models
  Sigma = bdiag(list(cov_beta,cov_theta))
  
  d1 = theta_3
  d2 = theta_3*A0
  d3 = (theta_3)*C
  d4 = 0
  d5 = 1
  d6 = theta_3*var_error
  #d7 = beta_0 + beta_1 * A0 + beta_2*C + theta_2 * var_error + theta_3 * var_error * (A0 + A1)
  d8 = 0* C
  #d9 = theta_2 * theta_3 + 0.5 * theta_3^2 *( A0 + A1)
  #Gamma = c(d1,d2,d3,d4,d5,d6,d7,d8)
  #we have no interaction in this model
  #so we don't include d7 in the model
  Gamma = c(d1,d2,d3,d4,d5,d6,d8)
  SE_NDE = as.numeric(sqrt(t(Gamma)%*%Sigma%*%Gamma*(A1-A0)^2))
  NDE_low = NDE-1.96*SE_NDE
  NDE_high = NDE+1.96*SE_NDE
  NDE_p = 2*pnorm(-abs(NDE/SE_NDE), lower.tail = T)
  
  b1 = 0
  b2 = theta_2 + theta_3 * A1
  b3 = 0*C
  b4 = 0
  b5 = 0
  b6 = beta_1
  #b7 = beta_1 * A1
  b8 = 0*C
  #b9 = 0
  #Gamma = c(b1,b2,b3,b4,b5,b6,b7,b8)
  Gamma = c(b1,b2,b3,b4,b5,b6,b8)
  SE_NIE = as.numeric(sqrt(t(Gamma)%*%Sigma%*%Gamma*(A1-A0)^2))
  NIE_low = NIE-1.96*SE_NIE
  NIE_high = NIE+1.96*SE_NIE
  NIE_p = 2*pnorm(-abs(NIE/SE_NIE), lower.tail = T)
  
  l1 = d1 + b1
  l2 = d2 + b2
  l3 = d3 + b3
  l4 = d4 + b4
  l5 = d5 + b5
  l6 = d6 + b6
  #l7 = d7 + b7
  l8 = d8 + b8
  Gamma = c(l1,l2,l3,l4,l5,l6,l8)
  SE_TE = as.numeric(sqrt(t(Gamma)%*%Sigma%*%Gamma*(A1-A0)^2))
  TE_low = TE-1.96*SE_TE
  TE_high = TE+1.96*SE_TE
  TE_p = 2*pnorm(-abs(TE/SE_TE), lower.tail = T)
  result_core = data.frame(NIE,SE_NIE,NIE_p,
                           NDE,SE_NDE,NDE_p,
                           TE,SE_TE,TE_p)
  
  
  result_summary = data.frame(NDE,NDE_low,NDE_high,NDE_p,
                              NIE,NIE_low,NIE_high,NIE_p,
                              TE,TE_low,TE_high,TE_p,PM)
  
  result = list(result_core,
                result_summary)
  
  
  
}






#Binary Mediator Continuous Outcome
MediationBC = function(out_model,med_model, A0, A1,C, M = NULL, Interaction = NULL){
  out_coef = coef(out_model)
  theta_0 = out_coef[1]
  theta_1 = out_coef[2]
  theta_2 = out_coef[3]
  
  #since no interaction is assumed in the model
  #theta_3 = 0
  theta_3 = 0
  cov_theta <- vcov(out_model)
  
  med_coef = coef(med_model)
  beta_0 = med_coef[1]
  beta_1 = med_coef[2]
  beta_2 = med_coef[3:length(med_coef)]
  
  cov_beta <- vcov(med_model)
  
  
  NDE = theta_1*(A1-A0) + theta_3*(A1-A0)*Probit(beta_0 + beta_1*A0+crossprod(beta_2,C))
  
  
  NIE = (theta_2 + theta_3*A1)*(Probit(beta_0+beta_1*A1+crossprod(beta_2,C))-Probit(beta_0+beta_1*A0+crossprod(beta_2,C)))
  
  
  TE = NDE + NIE
  
  #proportion of indirect effect over total effect
  PM = NIE/TE
  #covariance matrix of the regression coefficients from two models
  Sigma = bdiag(list(cov_beta,cov_theta))
  
  d1 = theta_3 *Probit(beta_0+beta_1*A0+crossprod(beta_2,C))-theta_3*Probit(beta_0+beta_1*A0+crossprod(beta_2,C))^2
  d2 = theta_3*(A0*Probit(beta_0+beta_1*A0+crossprod(beta_2,C))-Probit(beta_0+beta_1*A0+crossprod(beta_2,C))^2)
  d3 = as.numeric(theta_3*((Probit(beta_0+beta_1*A0+crossprod(beta_2,C)))*C-Probit(beta_0+beta_1*A0+crossprod(beta_2,C))^2))
  d4 = 0
  d5 = 1
  d6 = 0
  #d7 = Probit(beta_0+beta_1*A0+crossprod(beta_2,C))
  d8 = 0* C
  #d9 = theta_2 * theta_3 + 0.5 * theta_3^2 *( A0 + A1)
  #Gamma = c(d1,d2,d3,d4,d5,d6,d7,d8)
  #we have no interaction in this model
  #so we don't include d7 in the model
  Gamma = c(d1,d2,d3,d4,d5,d6,d8)
  SE_NDE = as.numeric(sqrt(t(Gamma)%*%Sigma%*%Gamma*(A1-A0)^2))
  NDE_low = NDE-1.96*SE_NDE
  NDE_high = NDE+1.96*SE_NDE
  NDE_p = 2*pnorm(-abs(NDE/SE_NDE), lower.tail = T)
  
  A = Probit(beta_0+beta_1*A1+crossprod(beta_2,C))-Probit(beta_0+beta_1*A1+crossprod(beta_2,C))^2
  B = Probit(beta_0+beta_1*A0+crossprod(beta_2,C))-Probit(beta_0+beta_1*A0+crossprod(beta_2,C))^2
  K = Probit(beta_0+beta_1*A1+crossprod(beta_2,C))
  D = Probit(beta_0+beta_1*A0+crossprod(beta_2,C))
  
  b1 = (theta_2 + theta_3*A1)*as.numeric(A-B)
  b2 = (theta_2 + theta_3 * A1)*(A1*A-A0*B)
  b3 = (theta_2 + theta_3 * A1)*as.numeric(A-B)*C
  b4 = 0
  b5 = 0
  b6 = K-D
  #b7 = a*(K-D)
  b8 = 0*C
  #b9 = 0
  #Gamma = c(b1,b2,b3,b4,b5,b6,b7,b8)
  Gamma = c(b1,b2,b3,b4,b5,b6,b8)
  SE_NIE = as.numeric(sqrt(t(Gamma)%*%Sigma%*%Gamma))
  NIE_low = NIE-1.96*SE_NIE
  NIE_high = NIE+1.96*SE_NIE
  NIE_p = 2*pnorm(-abs(NIE/SE_NIE), lower.tail = T)
  
  l1 = d1 + b1
  l2 = d2 + b2
  l3 = d3 + b3
  l4 = d4 + b4
  l5 = d5 + b5
  l6 = d6 + b6
  #l7 = d7 + b7
  l8 = d8 + b8
  Gamma = c(l1,l2,l3,l4,l5,l6,l8)
  SE_TE = as.numeric(sqrt(t(Gamma)%*%Sigma%*%Gamma*(A1-A0)^2))
  TE_low = TE-1.96*SE_TE
  TE_high = TE+1.96*SE_TE
  TE_p = 2*pnorm(-abs(TE/SE_TE), lower.tail = T)
  result_core = data.frame(NIE,SE_NIE,NIE_p,
                           NDE,SE_NDE,NDE_p,
                           TE,SE_TE,TE_p)
  
  
  result_summary = data.frame(NDE,NDE_low,NDE_high,NDE_p,
                              NIE,NIE_low,NIE_high,NIE_p,
                              TE,TE_low,TE_high,TE_p,PM)
  
  result = list(result_core,
                result_summary)
  
  
  
}


