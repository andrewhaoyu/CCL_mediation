setwd("/Users/zhangh24/Desktop")
data = read.csv("IVW_Final_Results.csv")
data_index = read.csv("UKBB_Data_Index.csv")
data$Pheno = data_index$Phenotype.Description
idx_include = which(data$P.Value!="N/A")
data_clean = data[idx_include,]
data_clean$P.Value = as.numeric(data_clean$P.Value)
range(data_clean$P.Value)
data_clean_new = data_clean[order(data_clean$P.Value),]
dim(data_clean_new)
0.05/1440
#Bonferroni correction
idx_sig = which(data_clean_new$P.Value<=0.05/1440)
data_clean_sig = data_clean_new[idx_sig,]
#write.csv(data_clean_sig, file = "UKBB_Data_IVW_sig.csv")




#p_clean_BF = p.adjust(data_clean_new$P.value, method = "bonferroni")
p_clean = p.adjust(data_clean_new$P.Value, method = "BH")
idx_sig2 = which(p_clean<=0.05)
data_clean_sig2 = data_clean_sig[idx_sig2,]
# idx_remove = which(data$P.Value=="N/A")
# data_remove = data[idx_remove,]
# table(data_remove$N)
