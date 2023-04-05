var_list = c("YRI_scale","ASN_scale","Former","Current","white_blood_cell_count",
             "monocyte_percentage","neutrophil_percentage",
             "autosome_mosaic")
setwd("/data/zhangh24/CLL_mediation/")
library(mediation)
i1 = 6
load(paste0("./result/mediation_result_",i1,".rdata"))
fit_model = result_list[[3]]
out_model = result_list[[2]]
med_model = result_list[[1]]
summary(fit_model)
summary(out_model)
summary(med_model)
