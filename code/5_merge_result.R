var_list = c("YRI_scale","ASN_scale","Former","Current","white_blood_cell_count",
             "monocyte_percentage","neutrophil_percentage",
             "autosome_mosaic","ch_chip","lymphoid","myeloid","lymphoid_chip",
             "myeloid_chip","both","lymphoid_myeloid")
setwd("/data/zhangh24/CLL_mediation/")
library(writexl)
library(data.table)
result_list = list()

SummarizeResult = function(result){
  sum_result = summary(result)
  coef = sum_result$summary_myreg
  logNDE = coef[2,1]
  NDE = exp(logNDE)
  NDE_low = exp(coef[2,5])
  NDE_high = exp(coef[2,6])
  NDE_p = 2*pnorm(-abs(coef[2,3]),lower.tail = T)
  logNIE = coef[3,1]
  NIE = exp(logNIE)
  NIE_low = exp(coef[3,5])
  NIE_high = exp(coef[3,6])
  NIE_p = 2*pnorm(-abs(coef[3,3]),lower.tail = T)
  logTE = coef[6,1]
  TE = exp(logTE)
  TE_low = exp(coef[6,5])
  TE_high = exp(coef[6,6])
  TE_p = 2*pnorm(-abs(coef[6,3]),lower.tail = T)
  pm = coef[7,1]
  pm_low = coef[7,5]
  pm_high = coef[7,6]
  pm_p = 2*pnorm(-abs(coef[7,3]),lower.tail = T)
  result = list(NDE,NDE_low,NDE_high,NDE_p,NIE,NIE_low,NIE_high,NIE_p,TE,TE_low,TE_high,TE_p,pm,pm_low,
             pm_high,pm_p)
  return(result)
}
result_list = list()
i2 = 1
for(i1 in 1:length(var_list)){
  load(paste0("./result/mediation_result_delta_",i1,"_",i2,".rdata"))
  result_list[[i1]] = SummarizeResult(result)
  
}
final_result1 = cbind(var_list,rbindlist(result_list))
i2 = 2
for(i1 in 1:length(var_list)){
  load(paste0("./result/mediation_result_delta_",i1,"_",i2,".rdata"))
  result_list[[i1]] = SummarizeResult(result)
  
}
final_result2 = cbind(var_list,rbindlist(result_list))

i2 = 3
for(i1 in 1:length(var_list)){
  load(paste0("./result/mediation_result_delta_",i1,"_",i2,".rdata"))
  result_list[[i1]] = SummarizeResult(result)
  
}
final_result3 = cbind(var_list,rbindlist(result_list))

result_names = c("Mediator",
                 "OR for Natural Director Effect (NDE)",
                 "95% CI low for OR of NDE",
                 "95% CI high for OR of NDE",
                 "P-value for OR of NDE",
                 "OR for Natural Indirector Effect (NIE)",
                 "95% CI low for OR of NIE",
                 "95% CI high for OR of NIE",
                 "P-value for OR of NIE",
                 "OR for Total Effect (TE)",
                 "95% CI low for OR of TE",
                 "95% CI high for OR of TE",
                 "P-value for OR of TE",
                 "Proportion of effect explained by indrect effect",
                 "95% CI low for proportion of effect explained by indrect effect",
                 "95% CI high for proportion of effect explained by indrect effect",
                 "P-value for proportion of effect explained by indrect effect")
colnames(final_result1) = colnames(final_result2) = 
  colnames(final_result3) = result_names

# Assuming df1 and df2 are your data frames
write_xlsx(list("PRS" = final_result1, "PRS1" = final_result2,
                "PRS2" = final_result3), path = "./result/mediation_result_071223.xlsx")


# write.csv(final_result, file = "./result/mediation_result_060523.csv")
# library(mediation)
# i1 = 6
# load(paste0("./result/mediation_result_",i1,".rdata"))
# fit_model = result_list[[3]]
# out_model = result_list[[2]]
# med_model = result_list[[1]]
# total_model = result_list[[4]]
# 
# summary(fit_model)
# #out model Y ~ M + X + C
# summary(out_model)
# #med model M ~ X + C
# summary(med_model)
# #total model Y ~  X + C
# summary(total_model)


pnorm(sqrt(2.7/2))
