args = commandArgs(trailingOnly = T)
#i1 for different variables
i1 = as.numeric(args[[1]])
#i2 for different prs
i2 = as.numeric(args[[2]])
#goal: mediation analyses for cll projects
setwd("/data/zhangh24/CLL_mediation/")
source("./code/MedFun.R")
library(survival)
library(regmedint)
#library(mediation)
library(data.table)
library(dplyr)
#load data with 436784 subjects
#data = readRDS("./data/mediation1.rds")
data = readRDS("./data/mediation_prscomp_chip.rds")
data$y = as.numeric(gsub(" days","",data$censor_days_cancer_ignore))
#436361 controls, 423 cases
#remove 2241 subjects (2236 controls, 5 cases) with missing smoking status
data = data %>% filter(smoke_NFC!=9)
#create data matrix for the smoking status
smoke_bin = model.matrix(~as.factor(smoke_NFC), data = data)[,-1]
colnames(smoke_bin) = c("Former", "Current")
data = cbind(data,smoke_bin)
#load prs file
prs = fread("./data/CLL_PRS_info/CLL_score.profile")
#combined data: 433257 controls, 418 cases
data_com = inner_join(data,prs, by = c("f.eid"="IID"))
#assign different PRS
if(i2 == 1){
  data_com$SCORESUM = data_com$SCORESUM
}else if(i2 ==2){
  data_com$SCORESUM = data_com$prs1
}else{
  data_com$SCORESUM = data_com$prs2
}
data_com_control = data_com[data_com$case_control_cancer_control==0,]
mean_prs = mean(data_com_control$SCORESUM,na.rm = T)
se_prs = sd(data_com_control$SCORESUM,na.rm = T)
data_com$SCORESUM_sd = (data_com$SCORESUM-mean_prs)/se_prs
#remove anyone people with missing PRS
#finalized dataset: 433362 controls, 313 cases
data_com_new = data_com %>% filter(!is.na(data_com$SCORESUM))

#variable_list
var_list = c("YRI_scale","ASN_scale","Former","Current","white_blood_cell_count",
             "monocyte_percentage","neutrophil_percentage",
             "autosome_mosaic","ch_chip","lymphoid","myeloid","lymphoid_chip",
             "myeloid_chip","both","lymphoid_myeloid")
#binary variable uses logistic regression
bin_var = c("Former","Current","autosome_mosaic","ch_chip","lymphoid","myeloid","lymphoid_chip",
            "myeloid_chip","both","lymphoid_myeloid")
#continuous variable uses linear regression
med_var_name = var_list[i1]
med_var = data_com_new[,med_var_name,drop=F]
colnames(med_var) = c("med_var")
data_clean = cbind(data_com_new,med_var)




#remove 13046 with missing white_blood_cell_count
data_clean = data_clean %>% filter(!is.na(white_blood_cell_count))

if(med_var_name%in%bin_var){
  if(med_var_name=="Former"){
    #set up the baseline value for other covariates
    C = c(mean(data_clean$age), mean(data_clean$age2), mean(data_clean$YRI_scale),  
          mean(data_clean$ASN_scale), 0, 0, mean(data_clean$white_blood_cell_count,na.rm = T))
    
    result <- regmedint(data = data_clean,
                                ## Variables
                                yvar = "y",
                                avar = "SCORESUM_sd",
                                mvar = "med_var",
                                cvar = c("age", "age2","YRI_scale", "ASN_scale","Current",
                                         "sex_new", "white_blood_cell_count"),
                                eventvar = "case_control_cancer_ignore",
                                ## Values at which effects are evaluated
                                a0 = 0,
                                a1 = 1,
                                m_cde = 0,
                                c_cond = C,
                                ## Model types
                                mreg = "logistic",
                                yreg = "survCox",
                                ## Additional specification
                                interaction = TRUE,
                                casecontrol = FALSE)
    
    
  }else if(med_var_name=="Current"){
    C = c(mean(data_clean$age), mean(data_clean$age2), mean(data_clean$YRI_scale),  
          mean(data_clean$ASN_scale), 0, 0, mean(data_clean$white_blood_cell_count,na.rm = T))
    
    result <- regmedint(data = data_clean,
                                ## Variables
                                yvar = "y",
                                avar = "SCORESUM_sd",
                                mvar = "med_var",
                                cvar = c("age", "age2","YRI_scale", "ASN_scale","Former",
                                         "sex_new", "white_blood_cell_count"),
                                eventvar = "case_control_cancer_ignore",
                                ## Values at which effects are evaluated
                                a0 = 0,
                                a1 = 1,
                                m_cde = 0,
                                c_cond = C,
                                ## Model types
                                mreg = "logistic",
                                yreg = "survCox",
                                ## Additional specification
                                interaction = TRUE,
                                casecontrol = FALSE)
  }else{
    C = c(mean(data_clean$age), mean(data_clean$age2), mean(data_clean$YRI_scale),  
          mean(data_clean$ASN_scale), 0, 0, 0, mean(data_clean$white_blood_cell_count,na.rm = T))
    result <- regmedint(data = data_clean,
                                ## Variables
                                yvar = "y",
                                avar = "SCORESUM_sd",
                                mvar = "med_var",
                                cvar = c("age", "age2","YRI_scale", "ASN_scale","Former","Current",
                                         "sex_new", "white_blood_cell_count"),
                                eventvar = "case_control_cancer_ignore",
                                ## Values at which effects are evaluated
                                a0 = 0,
                                a1 = 1,
                                m_cde = 0,
                                c_cond = C,
                                ## Model types
                                mreg = "logistic",
                                yreg = "survCox",
                                ## Additional specification
                                interaction = TRUE,
                                casecontrol = FALSE)
  }
}else{
  #continous mediator
  if(med_var_name=="YRI_scale"){
    C = c(mean(data_clean$age), mean(data_clean$age2), 
          mean(data_clean$ASN_scale), 0, 0, 0, mean(data_clean$white_blood_cell_count,na.rm = T))
    
    result <- regmedint(data = data_clean,
                                ## Variables
                                yvar = "y",
                                avar = "SCORESUM_sd",
                                mvar = "med_var",
                                cvar = c("age", "age2", "ASN_scale","Former","Current",
                                         "sex_new", "white_blood_cell_count"),
                                eventvar = "case_control_cancer_ignore",
                                ## Values at which effects are evaluated
                                a0 = 0,
                                a1 = 1,
                                m_cde = 0,
                                c_cond = C,
                                ## Model types
                                mreg = "linear",
                                yreg = "survCox",
                                ## Additional specification
                                interaction = TRUE,
                                casecontrol = FALSE)
    
  }else if(med_var_name=="ASN_scale"){
    C = c(mean(data_clean$age), mean(data_clean$age2), mean(data_clean$YRI_scale),  
          0, 0, 0, mean(data_clean$white_blood_cell_count,na.rm = T))
    
    result <- regmedint(data = data_clean,
                                ## Variables
                                yvar = "y",
                                avar = "SCORESUM_sd",
                                mvar = "med_var",
                                cvar = c("age", "age2","YRI_scale","Former","Current",
                                         "sex_new", "white_blood_cell_count"),
                                eventvar = "case_control_cancer_ignore",
                                ## Values at which effects are evaluated
                                a0 = 0,
                                a1 = 1,
                                m_cde = 0,
                                c_cond = C,
                                ## Model types
                                mreg = "linear",
                                yreg = "survCox",
                                ## Additional specification
                                interaction = TRUE,
                                casecontrol = FALSE)
    
  }else if(med_var_name=="white_blood_cell_count"){
    C = c(mean(data_clean$age), mean(data_clean$age2), mean(data_clean$YRI_scale),  
          mean(data_clean$ASN_scale), 0, 0, 0)
    
    result <- regmedint(data = data_clean,
                                ## Variables
                                yvar = "y",
                                avar = "SCORESUM_sd",
                                mvar = "med_var",
                                cvar = c("age", "age2","YRI_scale", "ASN_scale","Former","Current",
                                         "sex_new"),
                                eventvar = "case_control_cancer_ignore",
                                ## Values at which effects are evaluated
                                a0 = 0,
                                a1 = 1,
                                m_cde = 0,
                                c_cond = C,
                                ## Model types
                                mreg = "linear",
                                yreg = "survCox",
                                ## Additional specification
                                interaction = TRUE,
                                casecontrol = FALSE)
    
  }else{
    C = c(mean(data_clean$age), mean(data_clean$age2), mean(data_clean$YRI_scale),  
          mean(data_clean$ASN_scale), 0, 0, 0, mean(data_clean$white_blood_cell_count,na.rm = T))
    
    result <- regmedint(data = data_clean,
                                ## Variables
                                yvar = "y",
                                avar = "SCORESUM_sd",
                                mvar = "med_var",
                                cvar = c("age", "age2","YRI_scale", "ASN_scale","Former","Current",
                                         "sex_new", "white_blood_cell_count"),
                                eventvar = "case_control_cancer_ignore",
                                ## Values at which effects are evaluated
                                a0 = 0,
                                a1 = 1,
                                m_cde = 0,
                                c_cond = C,
                                ## Model types
                                mreg = "linear",
                                yreg = "survCox",
                                ## Additional specification
                                interaction = TRUE,
                                casecontrol = FALSE)
    
  }
  
}

save(result, file = paste0("./result/mediation_result_delta_",i1,"_",i2,".rdata"))
#data_new = readRDS("./data/mediation_prscomp_chip.rds")


#fit the mediator model
# if(med_var_name%in%bin_var){
#   #separate different exiting covariates
#   if(med_var_name=="Former"){
#     med_model = glm(med_var ~ SCORESUM_sd + age + age2 + YRI_scale + ASN_scale  + Current + sex_new + white_blood_cell_count, 
#                     data = data_clean,
#                     family = "binomial")
#   }else if(med_var_name=="Current"){
#     med_model = glm(med_var ~ SCORESUM_sd + age + age2 + YRI_scale + ASN_scale  + Former + sex_new + white_blood_cell_count, 
#                     data = data_clean,
#                     family = "binomial")
#   }else{
#     med_model = glm(med_var ~ SCORESUM_sd + age + age2 + YRI_scale + ASN_scale  + Former + Current + sex_new + white_blood_cell_count, 
#                     data = data_clean,
#                     family = "binomial")
#   }
#   
# }else{
#   #separate different exiting covariates
#   if(med_var_name=="YRI_scale"){
#     med_model = lm(med_var ~ SCORESUM_sd + age + age2  + ASN_scale  + Former + Current + sex_new + white_blood_cell_count, 
#                    data = data_clean)
#   }else if(med_var_name=="ASN_scale"){
#     med_model = lm(med_var ~ SCORESUM_sd + age + age2 + YRI_scale + ASN_scale  + Former + Current + sex_new + white_blood_cell_count, 
#                    data = data_clean)
#   }else if(med_var_name=="white_blood_cell_count"){
#     med_model = lm(med_var ~ SCORESUM_sd + age + age2 + YRI_scale + ASN_scale  + Former + Current + sex_new, 
#                    data = data_clean)
#   }else{
#     med_model = lm(med_var ~ SCORESUM_sd + age + age2 + YRI_scale + ASN_scale  + Former + Current + sex_new + white_blood_cell_count, 
#                    data = data_clean)
#   }
# }
# #fit the output model
# if(med_var_name=="Former"){
#   out_model = glm(case_control_cancer_ignore~SCORESUM_sd + med_var + age + age2 + YRI_scale + ASN_scale + Current + sex_new + white_blood_cell_count, 
#                   data = data_clean,
#                   family = "binomial")
#   
# }else if(med_var_name=="Current"){
#   out_model = glm(case_control_cancer_ignore~SCORESUM_sd + med_var + age + age2 + YRI_scale + ASN_scale + Former + sex_new + white_blood_cell_count, 
#                   data = data_clean,
#                   family = "binomial")
# }else if(med_var_name=="YRI_scale"){
#   out_model = glm(case_control_cancer_ignore~SCORESUM_sd + med_var + age + age2 + ASN_scale + Former + Current + sex_new + white_blood_cell_count, 
#                   data = data_clean,
#                   family = "binomial")
# }else if(med_var_name=="ASN_scale"){
#   out_model = glm(case_control_cancer_ignore~SCORESUM_sd + med_var + age + age2 + YRI_scale + Former + Current + sex_new + white_blood_cell_count, 
#                   data = data_clean,
#                   family = "binomial")
# }else if(med_var_name=="white_blood_cell_count"){
#   out_model = glm(case_control_cancer_ignore~SCORESUM_sd + med_var + age + age2 + YRI_scale + ASN_scale + Former + Current + sex_new, 
#                   data = data_clean,
#                   family = "binomial")
#   
# }else{
#   out_model = glm(case_control_cancer_ignore~SCORESUM_sd + med_var + age + age2 + YRI_scale + ASN_scale + Former + Current + sex_new + white_blood_cell_count, 
#                   data = data_clean,
#                   family = "binomial")
#   
# }
# 
# 
# if(med_var_name%in%bin_var){
#   #binary mediator binary outcome
#   result = MediationBB(out_model,med_model, A0 = 0, A1 = 1,C = C, M = NULL, Interaction = NULL)  
# }else{
#   #continous mediator binary outcome
#   result = MediationCB(out_model,med_model, A0 = 0, A1 = 1,C = C, M = NULL, Interaction = NULL)  
# }


# fit_model = mediate(med_model, out_model, treat='SCORESUM', mediator= "med_var", 
#                     #outcome = c("censor_days_cancer_ignore", "case_control_cancer_ignore"),
#                     boot=T, boot.ci.type = "bca")
#result_list = list(med_model, out_model, fit_model,total_model)
#save(result_list, file = paste0("./result/mediation_result_",i1,".rdata"))

# 
# 
# #fit the total_effect model
# if(med_var_name=="Former"){
#   total_model = glm(case_control_cancer_ignore~ SCORESUM_sd  + age + age2 + YRI_scale + ASN_scale + Current + sex_new + white_blood_cell_count, 
#                     data = data_clean,
#                     family = "binomial")
#   
# }else if(med_var_name=="Current"){
#   total_model = glm(case_control_cancer_ignore~SCORESUM_sd + age + age2 + YRI_scale + ASN_scale + Former + sex_new + white_blood_cell_count, 
#                     data = data_clean,
#                     family = "binomial")
# }else if(med_var_name=="YRI_scale"){
#   total_model = glm(case_control_cancer_ignore~SCORESUM_sd  + age + age2 + ASN_scale + Former + Current + sex_new + white_blood_cell_count, 
#                     data = data_clean,
#                     family = "binomial")
# }else if(med_var_name=="ASN_scale"){
#   total_model = glm(case_control_cancer_ignore~SCORESUM_sd  + age + age2 + YRI_scale + Former + Current + sex_new + white_blood_cell_count, 
#                     data = data_clean,
#                     family = "binomial")
# }else if(med_var_name=="white_blood_cell_count"){
#   total_model = glm(case_control_cancer_ignore~SCORESUM_sd  + age + age2 + YRI_scale + ASN_scale + Former + Current + sex_new, 
#                     data = data_clean,
#                     family = "binomial")
# }else{
#   total_model = glm(case_control_cancer_ignore~SCORESUM_sd  + age + age2 + YRI_scale + ASN_scale + Former + Current + sex_new + white_blood_cell_count, 
#                     data = data_clean,
#                     family = "binomial")
#   
# }