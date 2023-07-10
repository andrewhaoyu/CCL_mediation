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
beta_a = 0.5
beta = 2.9
beta_c = 0.3
P_Y = logit_fun(beta_0+beta*M+beta_c*C)
Y = rbinom(n,1,P_Y)

total_model = glm(Y~A+C,family = binomial())
med_model = lm(M~A+C)
out_model = glm(Y~A+M+C, family = binomial())

Mediation(out_model,med_model,total_model)

library(data.table)
library(bit64)
library(devtools)
install_github("RajLabMSSM/echotabix")
library(echotabix)
new_file <- gzfile('/data/BB_Bioinformatics/HZ_RQ/data/21001_irnt.gwas.imputed_v3.female.tsv.bgz',
                   'rt') 
path = '/data/BB_Bioinformatics/HZ_RQ/data/21001_irnt.gwas.imputed_v3.female.tsv.bgz'
data <- echotabix::read_bgz(path=path, method="utils")

data = read.table(file = new_file, header =T,
                sep = '\t')
