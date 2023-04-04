#goal: mediation analyses for cll projects
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
data_com = left_join(data,prs, by = c("f.eid"="IID"))
#finalized data: 434125 controls, 418 cases
#fit the output model
out_model = coxph(Surv(censor_days_cancer_ignore, case_control_cancer_ignore) ~ 
                    SCORESUM + age + age2 + YRI_scale + ASN_scale + Former + Current + sex_new, data = data_com)
med_model = glm(sex_new ~ SCORESUM + age + age2 + YRI_scale + ASN_scale + Former + Current, data = data_com,
                family = "binomial")
fit_model = mediate(med_model, out_model, treat='SCORESUM', mediator='sex_new', 
                    outcome = c("censor_days_cancer_ignore", "case_control_cancer_ignore"),
                    boot=F, data = data_com)
