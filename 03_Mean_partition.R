#############################################################
library(dplyr)
library(foreach)
library(tibble)
library(tidyr)
library(MM4LMM)
#############################################################
#A) Paths
############################################################# 

submission      = "submission1"

path            = "/home/tasha/breedfuture/"
path_input      = paste0(path, "prepared_data/v2/")
path_output     = paste0(path, "results/v2/model_output/mean_partition/")
path_kinships   = paste0(path, "results/v2/kinships/")

#############################################################
#B) Inputs
############################################################# 

quantiles       = c(paste0(seq(30, 90, by = 20), "%"), "95%", "99%")

GRM_G           = readRDS(paste0(path_kinships, "random_effects/GRM_baseline_v2.rds"))
GRM_G_upd       = as.matrix(Matrix::nearPD(GRM_G)$mat)

E_kernels       = readRDS(paste0(path_input, "environment_kernels.rds"))

data_list       = readRDS(paste0(path_input, "phenotypes.rds")) %>% 
                  dplyr::filter(analysis == "GP") %>% 
                  dplyr::mutate_at(vars("id", "year", "region", "country"), as.character) %>% split(.$trait)

loads_unpermuted= readRDS(paste0(path_kinships, "fixed_effects/loads_unpermuted_v2.rds"))
loads_permuted  = readRDS(paste0(path_kinships, "fixed_effects/loads_permuted_v2.rds"))

#############################################################
#C) Mean partition (unpermuted)
#############################################################    
 
mean_partition  = function(analysis, loads) {
  
  output        = data.frame()
  
  for(current_quantile in quantiles){
  print(current_quantile)
  
  loads_sub       = loads %>% dplyr::filter(quantile == current_quantile)
  
    for(current_trait in names(data_list)) {
    
    info          = list("analysis" = analysis, "trait" = current_trait, "quantile" = current_quantile, "permute" = "no", "permutation" = "none")
      
    data_subset   = data_list[[current_trait]] %>% 
                    left_join(loads_sub, by = "id") %>% 
                    dplyr::mutate(country = factor(country, levels = c("DK", "DE"))) 
                  
    KC            = E_kernels[["KC"]][as.character(data_subset$country), as.character(data_subset$country)]
    KR            = E_kernels[["KR"]][data_subset$region, data_subset$region]
    KY            = E_kernels[["KY"]][data_subset$year, data_subset$year] 
  
    GRM_G_ext     = GRM_G_upd[data_subset$id, data_subset$id]
                                   
                if(analysis == "marginal") {            
                                  
                X_input  = as.formula('effect ~ 1 + G + B + D + country + G:country') 
  
         } else if(analysis == "conditional") {
  
                X_input  = as.formula('effect ~ 1 + G + B + D + country + G:country + B:country + D:country')  } 
            
                V_list   = list("G" = GRM_G_ext, "GxC" = GRM_G_ext * KC, "CxY" = KC * KY, "error" = diag(1, nrow(data_subset)))

                if (length(unique(data_subset$region)) > 2) { V_list[["CxR"]] = KC * KR }


   tryCatch({
  
   SK_model       = MMEst(Y= data_subset$effect, Cofactor = data_subset, formula = X_input, Method = "Reml", VarList = V_list, NbCores=30) 
      
   Sigma          = SK_model$NullModel$VarBeta
   colnames(Sigma)= rownames(Sigma)
   Beta           = SK_model$NullModel$Beta
       
   if(analysis == "marginal") {  
   
   terms_list     = list(BD = c("B", "D"), B = "B", D = "D")

   wald_pvals     = sapply(terms_list, function(terms) {
                    aod::wald.test(Sigma = Sigma, b = Beta, Terms = which(names(Beta) %in% terms))$result$chi2["P"] })
   
   temp_output    = data.frame("B" = Beta["B"], "D" = Beta["D"], "B_SE" = sqrt(Sigma["B", "B"]), "D_SE" = sqrt(Sigma["D", "D"]),
                               "Joint_sig" = wald_pvals["BD.P"], "B_sig" = wald_pvals["B.P"], "D_sig" = wald_pvals["D.P"]) %>% 
                    dplyr::mutate(!!!info)
   
   output         = rbind(output, temp_output) 
                                    
    } else if(analysis == "conditional") {
  
   terms_list     = list(Joint = c("B", "D", "B:countryDE", "D:countryDE"), B = "B", D = "D", BxE = "B:countryDE", DxE = "D:countryDE")
 
   wald_pvals     = sapply(terms_list, function(terms) {
                    aod::wald.test(Sigma = Sigma, b = Beta, Terms = which(names(Beta) %in% terms))$result$chi2["P"] })

   temp_output    = data.frame("B" = Beta["B"], "D" = Beta["D"],
                               "BxE" = Beta["B:countryDE"], "DxE" = Beta["D:countryDE"],
                               "B_SE" = sqrt(Sigma["B", "B"]),"D_SE" = sqrt(Sigma["D", "D"]),
                               "BxE_SE" = sqrt(Sigma["B:countryDE", "B:countryDE"]), "DxE_SE" = sqrt(Sigma["D:countryDE", "D:countryDE"]),
                               "Joint_sig" = wald_pvals["Joint.P"], "B_sig" = wald_pvals["B.P"], "D_sig" = wald_pvals["D.P"], 
                               "BxE_sig" = wald_pvals["BxE.P"], "DxE_sig" = wald_pvals["DxE.P"]) %>%
                    dplyr::mutate(!!!info)

   output         = rbind(output, temp_output) } }, error = function(e) { cat("Error: ", e$message, "\n") }) }} 
  
   saveRDS(output, paste0(path_output,"mean_partition_", analysis, "_unpermuted.rds")) }
   
mean_partition(analysis = "marginal", loads = loads_unpermuted) 
mean_partition(analysis = "conditional", loads = loads_unpermuted) 

#############################################################
#D) Mean partition (permuted)
#############################################################    

mean_partition_perm = function(analysis, loads) {
  
  output        = data.frame()
  
  for(current_quantile in quantiles){
  print(current_quantile)
  
  loads_sub     = loads %>% dplyr::filter(quantile == current_quantile) %>% split(.$permutation)
  
    for(current_trait in names(data_list)) {
    
    info        = list("analysis" = analysis, "trait" = current_trait, "quantile" = current_quantile, "permute" = "yes")
      
    data_subset = data_list[[current_trait]] %>% 
                  dplyr::mutate(country = factor(country, levels = c("DK", "DE"))) 
                  
    KC          = E_kernels[["KC"]][as.character(data_subset$country), as.character(data_subset$country)]
    KR          = E_kernels[["KR"]][data_subset$region, data_subset$region]
    KY          = E_kernels[["KY"]][data_subset$year, data_subset$year] 
  
    GRM_G_ext   = GRM_G_upd[data_subset$id, data_subset$id]
    
    perm_list     = lapply(loads_sub, function(df) {
    rownames(df)  = NULL    
    out           = df %>% select(c(G, B, D, id)) %>% column_to_rownames("id") %>% as.matrix()                             
    out           = out[data_subset$id, ] })
                                                 
                if(analysis == "marginal") {            
                                  
                X_input  = as.formula('effect ~ 1 + G + B + D + country + G:country') 
  
         } else if(analysis == "conditional") {
  
                X_input  = as.formula('effect ~ 1 + G + B + D + country + G:country + B:country + D:country')  } 
  
          
                V_list          = list("G" = GRM_G_ext, "GxC"  = GRM_G_ext * KC, "CxY" = KC * KY, "error" = diag(1, nrow(data_subset)))

                if (length(unique(data_subset$region)) > 2) { V_list[["CxR"]] = KC * KR }

   tryCatch({
  
   SK_model     = MMEst(Y= data_subset$effect, 
                        formula = X_input, X = perm_list, Cofactor = data_subset, Method = "Reml", VarList = V_list, NbCores=cores)
   
   for(permutation in names(SK_model)){
 
   Beta           = SK_model[[permutation]]$Beta
   Sigma          = SK_model[[permutation]]$VarBeta
   colnames(Sigma)= rownames(Sigma)
       
          if(analysis == "marginal") {  
   
   terms_list     = list(BD = c("B", "D"), B = "B", D = "D")

   wald_pvals     = sapply(terms_list, function(terms) {
                    aod::wald.test(Sigma = Sigma, b = Beta, Terms = which(names(Beta) %in% terms))$result$chi2["P"] })
   
   temp_output    = data.frame("B" = Beta["B"], "D" = Beta["D"], "B_SE" = sqrt(Sigma["B", "B"]), "D_SE" = sqrt(Sigma["D", "D"]),
                             "Joint_sig" = wald_pvals["BD.P"], "B_sig" = wald_pvals["B.P"], "D_sig" = wald_pvals["D.P"]) %>% 
                    dplyr::mutate(!!!info) %>% mutate("permutation" = permutation)
   
   output         = rbind(output, temp_output) 
                                    
    } else if(analysis == "conditional") {
  
   terms_list     = list(Joint = c("B", "D", "B:countryDE", "D:countryDE"), B = "B", D = "D", BxE = "B:countryDE", DxE = "D:countryDE")
 
   wald_pvals     = sapply(terms_list, function(terms) {
                    aod::wald.test(Sigma = Sigma, b = Beta, Terms = which(names(Beta) %in% terms))$result$chi2["P"]})

   temp_output    = data.frame("B" = Beta["B"], "D" = Beta["D"],
                               "BxE" = Beta["B:countryDE"], "DxE" = Beta["D:countryDE"],
                               "B_SE" = sqrt(Sigma["B", "B"]),"D_SE" = sqrt(Sigma["D", "D"]),
                               "BxE_SE" = sqrt(Sigma["B:countryDE", "B:countryDE"]), "DxE_SE" = sqrt(Sigma["D:countryDE", "D:countryDE"]),
                               "Joint_sig" = wald_pvals["Joint.P"], "B_sig" = wald_pvals["B.P"], "D_sig" = wald_pvals["D.P"], 
                               "BxE_sig" = wald_pvals["BxE.P"], "DxE_sig" = wald_pvals["DxE.P"]) %>%
                    dplyr::mutate(!!!info) %>% mutate("permutation" = permutation)

   output        = rbind(output, temp_output) }} }, error = function(e) { cat("Error:", e$message, "\n") }) }} 
  
   saveRDS(output, paste0(path_output,"mean_partition_", analysis, "_permuted.rds")) } 

mean_partition_perm(analysis = "marginal", loads = loads_permuted) 
mean_partition_perm(analysis = "conditional", loads = loads_permuted)

#############################################################