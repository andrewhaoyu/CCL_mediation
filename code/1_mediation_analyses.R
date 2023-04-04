args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
#goal: mediation analyses for cll projects
setwd("/data/zhangh24/CLL_mediation/")
library(survival)
library(mediation)
library(data.table)
library(dplyr)
#load data with 436784 subjects
data = readRDS("./data/mediation1.rds")
#436361 controls, 423 cases
#remove 2241 subjects (2236 controls, 5 cases) with missing smoking status
data = data %>% filter(smoke_NFC!=9)
#create data matrix for the smoking status
smoke_bin = model.matrix(~as.factor(smoke_NFC), data = data)[,-1]
colnames(smoke_bin) = c("Former", "Current")
data = cbind(data,smoke_bin)
#load prs file
prs = fread("./data/CLL_PRS_info/CLL_score.profile")
#finalized data: 434125 controls, 418 cases
data_com = left_join(data,prs, by = c("f.eid"="IID"))
#variable_list
var_list = c("YRI_scale","ASN_scale","Former","Current","white_blood_cell_count",
             "monocyte_percentage","neutrophil_percentage",
             "autosome_mosaic")
#binary variable uses logistic regression
bin_var = c("Former","Current","autosome_mosaic")
#continuous variable uses linear regression
treat_var_name = var_list[i1]
treat_var = data_com[,treat_var_name,drop=F]
colnames(treat_var) = c("treat_var")
data_clean = cbind(data_com,treat_var)
#fit the mediator model
if(treat_var_name%in%bin_var){
  #separate different exiting covariates
  if(treat_var_name=="Former"){
    med_model = glm(treat_var ~ SCORESUM + age + age2 + YRI_scale + ASN_scale  + Current, 
                    data = data_clean,
                    family = "binomial")
  }else if(treat_var_name=="Current"){
    med_model = glm(treat_var ~ SCORESUM + age + age2 + YRI_scale + ASN_scale  + Former, 
                    data = data_clean,
                    family = "binomial")
  }else{
    med_model = glm(treat_var ~ SCORESUM + age + age2 + YRI_scale + ASN_scale  + Former + Current, 
                    data = data_clean,
                    family = "binomial")
  }
  
}else{
  #separate different exiting covariates
  if(treat_var_name=="YRI_scale"){
    med_model = lm(treat_var ~ SCORESUM + age + age2  + ASN_scale  + Former + Current, 
                    data = data_clean)
  }else if(treat_var_name=="ASN_scale"){
    med_model = lm(treat_var ~ SCORESUM + age + age2 + YRI_scale + ASN_scale  + Former + Current, 
                    data = data_clean)
  }else{
    med_model = lm(treat_var ~ SCORESUM + age + age2 + YRI_scale + ASN_scale  + Former + Current, 
                    data = data_clean)
  }
}
#fit the output model
if(treat_var_name=="Former"){
  out_model = glm(case_control_cancer_ignore~SCORESUM + treat_var + age + age2 + YRI_scale + ASN_scale + Current + sex_new, 
                  data = data_clean,
                  family = "binomial")
  
}else if(treat_var_name=="Current"){
  out_model = glm(case_control_cancer_ignore~SCORESUM + treat_var + age + age2 + YRI_scale + ASN_scale + Former + sex_new, 
                  data = data_clean,
                  family = "binomial")
}else if(treat_var_name=="YRI_scale"){
  out_model = glm(case_control_cancer_ignore~SCORESUM + treat_var + age + age2 + ASN_scale + Former + Current + sex_new, 
                  data = data_clean,
                  family = "binomial")
}else if(treat_var_name=="ASN_scale"){
  out_model = glm(case_control_cancer_ignore~SCORESUM + treat_var + age + age2 + YRI_scale + Former + Current + sex_new, 
                  data = data_clean,
                  family = "binomial")
}else{
  out_model = glm(case_control_cancer_ignore~SCORESUM + treat_var + age + age2 + YRI_scale + ASN_scale + Former + Current + sex_new, 
                  data = data_clean,
                  family = "binomial")
  
}

fit_model = mediate(med_model, out_model, treat='SCORESUM', mediator='sex_new', 
                    #outcome = c("censor_days_cancer_ignore", "case_control_cancer_ignore"),
                    boot=T, boot.ci.type = "bca")
result_list = list(med_model, out_model, fit_model)
save(result_list, file = paste0("./result/mediation_result_",i1,".rdata"))
