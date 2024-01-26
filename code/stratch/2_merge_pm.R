library(data.table)
data_list = list()
for(i in 1:1000){
  load(paste0("/data/BB_Bioinformatics/HZ/example_mediation/pm_matrix_",i1,".rdata"))
  data_list[[i]]  = as.data.frame(pm_matrix)
}
pm_result = rbindlist(data_list)

#95% CI and p-value for the results with functions
