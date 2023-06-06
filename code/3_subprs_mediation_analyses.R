args = commandArgs(trailingOnly = T)
#i1 for different variable
i1 = as.numeric(args[[1]])
#i2 for different prs
i2 = as.numeric(args[[2]])
#goal: mediation analyses for cll projects
setwd("/data/zhangh24/CLL_mediation/")
library(survival)
library(mediation)
library(data.table)
library(dplyr)
#load data with 436784 subjects
data = readRDS("./data/subPRS_updated.rds")
#436361 controls, 423 cases
#remove 2241 subjects (2236 controls, 5 cases) with missing smoking status
data = data %>% filter(smoke_NFC!=9)
#create data matrix for the smoking status
smoke_bin = model.matrix(~as.factor(smoke_NFC), data = data)[,-1]
colnames(smoke_bin) = c("Former", "Current")
data = cbind(data,smoke_bin)
#load prs file
#finalized data: 434125 controls, 418 cases
data_com = data
if(i2 == 1){
  data_com$SCORESUM = data_com$prs1
}else{
  data_com$SCORESUM = data_com$prs2
}
data_com_control = data_com[data_com$case_control_cancer_control==0,]
mean_prs = mean(data_com_control$SCORESUM,na.rm = T)
se_prs = sd(data_com_control$SCORESUM,na.rm = T)
data_com$SCORESUM_sd = (data_com$SCORESUM-mean_prs)/se_prs
#variable_list
var_list = c("YRI_scale","ASN_scale","Former","Current","white_blood_cell_count",
             "monocyte_percentage","neutrophil_percentage",
             "autosome_mosaic")
#binary variable uses logistic regression
bin_var = c("Former","Current","autosome_mosaic")
#continuous variable uses linear regression
med_var_name = var_list[i1]
med_var = data_com[,med_var_name,drop=F]
colnames(med_var) = c("med_var")
data_clean = cbind(data_com,med_var)
#fit the total_effect model
if(med_var_name=="Former"){
  total_model = glm(case_control_cancer_ignore~ SCORESUM_sd  + age + age2 + YRI_scale + ASN_scale + Current + sex_new + white_blood_cell_count, 
                    data = data_clean,
                    family = "binomial")
  
}else if(med_var_name=="Current"){
  total_model = glm(case_control_cancer_ignore~SCORESUM_sd + age + age2 + YRI_scale + ASN_scale + Former + sex_new + white_blood_cell_count, 
                    data = data_clean,
                    family = "binomial")
}else if(med_var_name=="YRI_scale"){
  total_model = glm(case_control_cancer_ignore~SCORESUM_sd  + age + age2 + ASN_scale + Former + Current + sex_new + white_blood_cell_count, 
                    data = data_clean,
                    family = "binomial")
}else if(med_var_name=="ASN_scale"){
  total_model = glm(case_control_cancer_ignore~SCORESUM_sd  + age + age2 + YRI_scale + Former + Current + sex_new + white_blood_cell_count, 
                    data = data_clean,
                    family = "binomial")
}else if(med_var_name=="white_blood_cell_count"){
  total_model = glm(case_control_cancer_ignore~SCORESUM_sd  + age + age2 + YRI_scale + ASN_scale + Former + Current + sex_new, 
                    data = data_clean,
                    family = "binomial")
}else{
  total_model = glm(case_control_cancer_ignore~SCORESUM_sd  + age + age2 + YRI_scale + ASN_scale + Former + Current + sex_new + white_blood_cell_count, 
                    data = data_clean,
                    family = "binomial")
  
}

#fit the mediator model
if(med_var_name%in%bin_var){
  #separate different exiting covariates
  if(med_var_name=="Former"){
    med_model = glm(med_var ~ SCORESUM_sd + age + age2 + YRI_scale + ASN_scale  + Current + white_blood_cell_count, 
                    data = data_clean,
                    family = "binomial")
  }else if(med_var_name=="Current"){
    med_model = glm(med_var ~ SCORESUM_sd + age + age2 + YRI_scale + ASN_scale  + Former + white_blood_cell_count, 
                    data = data_clean,
                    family = "binomial")
  }else{
    med_model = glm(med_var ~ SCORESUM_sd + age + age2 + YRI_scale + ASN_scale  + Former + Current + white_blood_cell_count, 
                    data = data_clean,
                    family = "binomial")
  }
  
}else{
  #separate different exiting covariates
  if(med_var_name=="YRI_scale"){
    med_model = lm(med_var ~ SCORESUM_sd + age + age2  + ASN_scale  + Former + Current + white_blood_cell_count, 
                   data = data_clean)
  }else if(med_var_name=="ASN_scale"){
    med_model = lm(med_var ~ SCORESUM_sd + age + age2 + YRI_scale + ASN_scale  + Former + Current + white_blood_cell_count, 
                   data = data_clean)
  }else if(med_var_name=="white_blood_cell_count"){
    med_model = lm(med_var ~ SCORESUM_sd + age + age2 + YRI_scale + ASN_scale  + Former + Current , 
                   data = data_clean)
  }else{
    med_model = lm(med_var ~ SCORESUM_sd + age + age2 + YRI_scale + ASN_scale  + Former + Current + white_blood_cell_count, 
                   data = data_clean)
  }
}
#fit the output model
if(med_var_name=="Former"){
  out_model = glm(case_control_cancer_ignore~SCORESUM_sd + med_var + age + age2 + YRI_scale + ASN_scale + Current + sex_new + white_blood_cell_count, 
                  data = data_clean,
                  family = "binomial")
  
}else if(med_var_name=="Current"){
  out_model = glm(case_control_cancer_ignore~SCORESUM_sd + med_var + age + age2 + YRI_scale + ASN_scale + Former + sex_new + white_blood_cell_count, 
                  data = data_clean,
                  family = "binomial")
}else if(med_var_name=="YRI_scale"){
  out_model = glm(case_control_cancer_ignore~SCORESUM_sd + med_var + age + age2 + ASN_scale + Former + Current + sex_new + white_blood_cell_count, 
                  data = data_clean,
                  family = "binomial")
}else if(med_var_name=="ASN_scale"){
  out_model = glm(case_control_cancer_ignore~SCORESUM_sd + med_var + age + age2 + YRI_scale + Former + Current + sex_new + white_blood_cell_count, 
                  data = data_clean,
                  family = "binomial")
}else if(med_var_name=="white_blood_cell_count"){
  out_model = glm(case_control_cancer_ignore~SCORESUM_sd + med_var + age + age2 + YRI_scale + ASN_scale + Former + Current + sex_new, 
                  data = data_clean,
                  family = "binomial")
  
}else{
  out_model = glm(case_control_cancer_ignore~SCORESUM_sd + med_var + age + age2 + YRI_scale + ASN_scale + Former + Current + sex_new + white_blood_cell_count, 
                  data = data_clean,
                  family = "binomial")
  
}

Mediation = function(out_model, med_model, total_model){
  out_coef = coefficients(summary(out_model))
  log_NDE = out_coef[2,1]
  log_NDE_se = out_coef[2,2]
  NDE_p = out_coef[2,4]
  OR_NDE = exp(log_NDE)
  OR_NDE_low = exp(log_NDE-1.96*log_NDE_se)
  OR_NDE_high = exp(log_NDE+1.96*log_NDE_se)
  med_coef = coefficients(summary(med_model))
  log_NIE = out_coef[3,1]*med_coef[2,1]
  log_NIE_se = sqrt(out_coef[3,1]^2*med_coef[2,2]^2+
                      med_coef[2,1]^2*out_coef[3,2]^2)
  NIE_p = 2*pnorm(-abs(log_NIE/log_NIE_se), lower.tail = T)
  OR_NIE = exp(log_NIE)
  OR_NIE_low = exp(log_NIE-1.96*log_NIE_se)
  OR_NIE_high = exp(log_NIE+1.96*log_NIE_se)
  log_TDE = log_NDE + log_NIE
  log_TDE_se = sqrt(log_NIE_se^2+log_NDE_se^2)
  TDE_p = 2*pnorm(-abs(log_TDE/log_TDE_se), lower.tail = T)
  OR_TDE = exp(log_TDE)
  OR_TDE_low = exp(log_TDE-1.96*log_TDE_se)
  OR_TDE_high = exp(log_TDE+1.96*log_TDE_se)
  OR_TDE = exp(log_TDE)
  proportion = log_NIE/log_TDE
  
  result = data.frame(OR_NDE,OR_NDE_low,OR_NDE_high,NDE_p,
                      OR_NIE,OR_NIE_low,OR_NIE_high,NIE_p,
                      OR_TDE,OR_TDE_low,OR_TDE_high,TDE_p,proportion)
  
  return(result)
  
  
}
result = Mediation(out_model,med_model,
                   total_model)


# fit_model = mediate(med_model, out_model, treat='SCORESUM', mediator= "med_var", 
#                     #outcome = c("censor_days_cancer_ignore", "case_control_cancer_ignore"),
#                     boot=T, boot.ci.type = "bca")
#result_list = list(med_model, out_model, fit_model,total_model)
#save(result_list, file = paste0("./result/mediation_result_",i1,".rdata"))
save(result, file = paste0("./result/mediation_result_sub_",i1,"_",i2,".rdata"))
