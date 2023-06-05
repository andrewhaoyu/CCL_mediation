var_list = c("YRI_scale","ASN_scale","Former","Current","white_blood_cell_count",
             "monocyte_percentage","neutrophil_percentage",
             "autosome_mosaic")
setwd("/data/zhangh24/CLL_mediation/")
result_list = list()
for(i1 in 1:length(var_list)){
  load(paste0("./result/mediation_result_delta_",i1,".rdata"))
  result_list[[i1]] = result
}
final_result = cbind(var_list,rbindlist(result_list))
colnames(final_result) = c("Mediator",
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
                           "Proportion of effect explained by indrect effect")
write.csv(final_result, file = "./result/mediation_result_060523.csv")
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
