args = commandArgs(trailingOnly = T)
#i1 for different
i1 = as.numeric(args[[1]])

library(data.table)
library(readxl)
library(vroom)
library(dplyr)
library(regmedint)
library(tidyverse)
library(boot)
library(doParallel)

probit = function(x){exp(x)/(1+exp(x))}

icogs <- fread(file = "/data/BB_Bioinformatics/DG/725/zhang_icogs_pheno_v15_02_age_used.txt")
onco <- fread("/data/BB_Bioinformatics/DG/725/zhang_onco_pheno_v15_02_corrected_age_used.txt")
prs_icogs <- fread("/data/BB_Bioinformatics/DG/725/icogsPRS.sscore")
prs_onco  <- fread("/data/BB_Bioinformatics/DG/725/oncoPRS.sscore")

prs_icogs$'#IID' <- sapply(strsplit(prs_icogs$`#IID`, "_"), `[`, 1)
prs_onco$'#IID' <- sapply(strsplit(prs_onco$`#IID`, "_"), `[`, 1)
names(prs_icogs)[names(prs_icogs) == "#IID"] <- "id"
names(prs_onco)[names(prs_onco) == "#IID"] <- "id"
names(onco)[names(onco) == "Onc_ID"] <- "id"
names(icogs)[names(icogs) == "SG_ID"] <- "id"

prs_icogs$id <- as.integer(prs_icogs$id)
prs_onco$id <- as.integer(prs_onco$id)
icogs_with_prs <- left_join(icogs, prs_icogs, by = "id")
onco_with_prs <- left_join(onco, prs_onco, by = "id")

icogs_with_prs <- icogs_with_prs %>%
  rename(
    PC_1 = pc1,
    PC_2 = pc2,
    PC_3 = pc3,
    PC_4 = pc4,
    PC_5 = pc5,
    PC_6 = pc6,
    PC_7 = pc7,
    PC_8 = pc8,
    PC_9 = pc9,
    PC_10 = pc10,
    PC_11 = pc11,
    PC_12 = pc12,
    PC_13 = pc13,
    PC_14 = pc14,
    PC_15 = pc15
  ) %>%
  select(-c(PC_LMBC, Expr1057))

onco_with_prs <- onco_with_prs %>% select(-c(Expr1056))
identical(names(icogs_with_prs), names(onco_with_prs))
combined_data <- rbind(icogs_with_prs, onco_with_prs)

##cleaning data based on status variable
combined_data$SCORE1_SUM <- scale(combined_data$SCORE1_SUM)
combined_data <- as.data.frame(combined_data)
combined_data <- combined_data[combined_data[, "status"] %in% c(0, 1), ]
pm_fun <- function(nde, nie){
  result = (exp(nde)*(exp(nie)-1))/(exp(nde)*exp(nie)-1)
  return(result)
}
mediators <- c("ageMenarche", "mensAgeLast")
mediate <- function(data) {
    covariates <- c("ageInt","PC_1", "PC_2", "PC_3", "PC_4", "PC_5", "PC_6", "PC_7", "PC_8", "PC_9", "PC_10")
    #mediators <- c("ageMenarche", "mensAgeLast", "breastFed", "breastMos", "BMI", "OCEver", "HRTEver", "HRTCurrent", "EPEver", "smokingEver", "smokingCurrent", "physActBefore30", "physAct30_50", "physAct50plus")
    mediators <- c("ageMenarche", "mensAgeLast")
    covariate_columns <- data[, covariates]
  means <- colMeans(covariate_columns, na.rm = TRUE)
  mediation_results <- list()
  
  for (mediator in mediators){
    columns_to_subset <- c("status", "SCORE1_SUM", mediator, covariates)
    sub <- data %>% select(all_of(columns_to_subset)) %>% filter(!.data[[mediator]] %in% c(NA, 777, 888)) %>% filter(complete.cases(.))
    
    if(mediator %in% c("ageMenarche", "mensAgeLast", "breastMos", "BMI", "physActBefore30", "physAct30_50", "physAct50plus")){
      regmedint_obj1 <- regmedint(data = sub,
                                  yvar = "status",
                                  avar = "SCORE1_SUM",
                                  mvar = mediator,
                                  cvar = covariates,
                                  a0 = 0,
                                  a1 = 1,
                                  m_cde = 0,
                                  c_cond = means,
                                  mreg = "linear",
                                  yreg = "logistic",
                                  interaction = FALSE,
                                  casecontrol = FALSE)
    }
    else{
      sub <- sub[sub[, mediator] %in% c(0, 1), ]
      regmedint_obj1 <- regmedint(data = sub,
                                  yvar = "status",
                                  avar = "SCORE1_SUM",
                                  mvar = mediator,
                                  cvar = covariates,
                                  a0 = 0,
                                  a1 = 1,
                                  m_cde = 0,
                                  c_cond = means,
                                  mreg = "logistic",
                                  yreg = "logistic",
                                  interaction = FALSE,
                                  casecontrol = FALSE)
    }
    s <- summary(regmedint_obj1)
    mediation_results <- append(mediation_results, mediator)
    mediation_results <- append(mediation_results, list(s$summary_myreg[c(2,3,6),]))
  }
  
  df_list <- lapply(seq(1, length(mediation_results), by = 2), function(i) {
    mediator <- mediation_results[[i]]
    coefficients <- as.data.frame(mediation_results[[i+1]])
    coefficients$x <- rownames(coefficients)
    rownames(coefficients) <- NULL
    cbind(mediator, coefficients)
  })
  combined_df <- bind_rows(df_list)
  rownames(combined_df) <- NULL
  wider_df <- pivot_wider(combined_df, names_from = x, values_from = est:upper)
  return(wider_df)
}

LogoddsMetaAnalysis <- function(logodds1,sigma1,logodds2,sigma2){
  sigma1.inv <- solve(sigma1)
  sigma2.inv <- solve(sigma2)
  sigma.meta <- solve(sigma1.inv+sigma2.inv)
  logodds.meta <- sigma.meta%*%(sigma1.inv%*%logodds1+sigma2.inv%*%logodds2)
  
  return(list(logodds.meta = logodds.meta,
              sigma.meta = sigma.meta))
}

n1 <- nrow(onco_with_prs)
n2 <- nrow(icogs_with_prs)

bootstrap_meta <- function(data) {
  indices <- sample(1:(n1 + n2), size = (n1 + n2), replace = TRUE)
  sample1 <- data[indices[1:n1], ]
  sample2 <- data[indices[(n1+1):(n1+n2)], ]
  
  fit1 <- mediate(sample1)
  fit2 <- mediate(sample2)
  
  results_list <- list()
  
  values <- c("pnde", "tnie", "te")
  for (var in values) {
    sigmai <- diag(fit2[[paste0("se_", var)]]^2)
    sigmao <- diag(fit1[[paste0("se_", var)]]^2)
    logoddsi <- matrix(fit2[[paste0("est_", var)]], ncol = 1)
    logoddso <- matrix(fit1[[paste0("est_", var)]], ncol = 1)
    result <- LogoddsMetaAnalysis(logoddsi, sigmai, logoddso, sigmao)
    results_list[[var]] <- data.frame(
      Mediator = var,
      Est = result$logodds.meta,
      SE = diag(result$sigma.meta),
      lower = result$logodds.meta - 1.96 * diag(result$sigma.meta),
      upper = result$logodds.meta + 1.96 * diag(result$sigma.meta)
    )
  }
  
  combined_results <- do.call(rbind, results_list)
  combined_results$var <- rep(fit1$mediator)
  
  wide_combined_results <- combined_results %>%
    pivot_wider(
      names_from = Mediator,
      values_from = c(Est, SE, lower, upper)
    )
  
  pm <- pm_fun(wide_combined_results$Est_pnde,wide_combined_results$Est_tnie)
  return(pm)
}

set.seed(i1)
num_bootstrap <- 2
pm_matrix <- matrix(numeric(0), nrow = num_bootstrap, ncol = length(mediators))

for (i in 1:num_bootstrap) {
  print(i)
  pm_matrix[,i] <- t(bootstrap_meta(combined_data))
}

# save(pm_matrix,
#      file = paste0("/data/BB_Bioinformatics/DG/725/pm_matrix_",i1,".rdata"))
save(pm_matrix, file = 
       paste0("/data/BB_Bioinformatics/HZ/example_mediation/pm_matrix_",i1,".rdata")
       )


# cl <- makeCluster(parallelly::availableCores())
# clusterSetRNGStream(cl, 456)
# registerDoParallel(cl)
# 
# pm_matrix <- foreach(1:num_bootstrap, .combine=cbind, .packages=c("dplyr", "regmedint", "tidyr")) %dopar% {
#   bootstrap_meta(combined_data)
# }
# 
# stopCluster(cl)



# # Confience interval
# conf_int <- quantile(results$t, c(0.025, 0.975))
# 
# # Two-tailed p-value for combined effect size being different from 0
# p_val <- mean(abs(results$t) >= abs(mean(results$t)))
# 
# list(confidence_interval = conf_int, p_value = p_val)
