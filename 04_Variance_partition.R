#############################################################
library(dplyr)
library(foreach)
library(tibble)
library(tidyr)
library(MM4LMM)
library(qgg)
#############################################################
#A) Paths
############################################################# 

submission      = "submission1"

path            = "/home/tasha/breedfuture/"
path_input      = paste0(path, "prepared_data/v2/")
path_output     = paste0(path, "results/v2/model_output/variance_partition/")
path_kinships   = paste0(path, "results/v2/kinships/")

#############################################################
#B) Inputs
############################################################# 

quantiles       = c(paste0(seq(30, 90, by = 20), "%"), "95%", "99%")

X_input         = as.formula('effect ~ 1 + G + country + G:country')

############################################################# 

GRM_G           = readRDS(paste0(path_kinships, "random_effects/GRM_baseline_v2.rds"))

E_kernels       = readRDS(paste0(path_input, "environment_kernels.rds"))

data_list       = readRDS(paste0(path_input, "phenotypes.rds")) %>% 
                  dplyr::filter(analysis == "GP") %>% 
                  dplyr::mutate_at(vars("id", "year", "region", "country"), as.character) %>% split(.$trait)

loads           = readRDS(paste0(path_kinships, "fixed_effects/loads_unpermuted_v2.rds")) %>% distinct(id, G)

############################################################# 
#C) Variance partition- single GRM (unpermuted)
############################################################# 
   
   SK_var       = data.frame()
   SK_llr       = data.frame()
   
   for(current_trait in names(data_list)) {
   
   meta         = list("trait" = current_trait, random_model = "M1", "quantile" = "none", "kernels" = "SK", "permuted" = "no", "permutation" = 0) 
          
   data_subset  = data_list[[current_trait]] %>% left_join(loads, by = "id") 
   
   y            = setNames(data_subset$effect, rownames(data_subset))
   
   X_form       = model.matrix(as.formula(X_input), data = data_subset)
                                    
   KC           = E_kernels[["KC"]][as.character(data_subset$country), as.character(data_subset$country)]
   KR           = E_kernels[["KR"]][data_subset$region, data_subset$region]
   KY           = E_kernels[["KY"]][data_subset$year, data_subset$year]   
   
   GRM_G_ext    = GRM_G[data_subset$id, data_subset$id]
   
   V_M1_list    = list("G" = GRM_G_ext, "GxC" = GRM_G_ext * KC, "CxY" = KC * KY)

                  if (length(unique(data_subset$region)) > 2) { V_M1_list[["CxR"]] = KC * KR }
                         
   tryCatch({
  
   SK_model     = greml(y = y, X = X_form, GRM = V_M1_list, ncores = 25)
                   
   SK_var_temp  = SK_model$theta %>% as.data.frame() %>% 
                  rownames_to_column("Component") %>% 
                  dplyr::rename("Variance" = ".") %>% 
                  dplyr::mutate(!!!meta) 
                   
   SK_llr_temp  = SK_model$llik %>% as.data.frame() %>%
                  dplyr::rename("llik" = ".") %>% 
                  dplyr::mutate(!!!meta) %>% mutate("iter" = SK_model$niter)
                                                          
   SK_var       = rbind(SK_var, SK_var_temp) 
   SK_llr       = rbind(SK_llr, SK_llr_temp) }, 
   
   error = function(e) { cat("Error ", e$message, "\n") }) }  
   
   saveRDS(SK_var, paste0(path_output, "SK_variance_components_unpermuted.rds")) 
   saveRDS(SK_llr, paste0(path_output, "SK_likelihoods_unpermuted.rds")) 
   
############################################################# 
#D) Variance partition-multiple GRM (unpermuted)
############################################################# 
   
   MK_var        = data.frame()
   MK_llr        = data.frame()
   
   for(current_quantile in quantiles){
   print(paste0("quantile is ", current_quantile))

   G_kernels     = readRDS(paste0(path_kinships, "random_effects/GRMs_unpermuted_", current_quantile,"_v2.rds"))[[1]][["kinships"]]
    GRM_S        = G_kernels[["Sel"]]
    GRM_N        = G_kernels[["Neu"]]
   
   for(current_trait in names(data_list)) {
   
   meta          = list("trait" = current_trait, "quantile" = current_quantile, "kernels" = "MK", "permuted" = "no", "permutation" = 0) 
          
   data_subset   = data_list[[current_trait]] %>% left_join(loads, by = "id") 
   
   y             = setNames(data_subset$effect, rownames(data_subset))
   
   X_form        = model.matrix(as.formula(X_input), data = data_subset)
                     
   KC            = E_kernels[["KC"]][as.character(data_subset$country), as.character(data_subset$country)]
   KR            = E_kernels[["KR"]][data_subset$region, data_subset$region]
   KY            = E_kernels[["KY"]][data_subset$year, data_subset$year] 
  
   GRM_G_ext     = GRM_G[data_subset$id, data_subset$id]
   GRM_S_ext     = GRM_S[data_subset$id, data_subset$id]
   GRM_N_ext     = GRM_N[data_subset$id, data_subset$id]

   V_M1_list     = list("G"  = GRM_G_ext, "GxC"= GRM_G_ext * KC, "CxY" = KC * KY)
   V_M2_list     = list("GN" = GRM_N_ext, "GS" = GRM_S_ext, "GxC" = GRM_G_ext * KC, "CxY" = KC * KY)
   V_M3_list     = list("GN" = GRM_N_ext, "GS" = GRM_S_ext, "GNxC"= GRM_N_ext * KC, "GSxC" = GRM_S_ext * KC, "CxY" = KC * KY)

   if(length(unique(data_subset$region)) > 2) {
  
   V_M2_list[["CxR"]] = KC * KR
   V_M3_list[["CxR"]] = KC * KR }

   tryCatch({
  
   models_list   = list("M2" = V_M2_list, "M3" = V_M3_list)

   results       = lapply(names(models_list), function(model_name) {
  
   MK_model      = greml(y = y, X = X_form, GRM = models_list[[model_name]], ncores = 25)
                   
   MK_var_temp   = MK_model$theta %>% as.data.frame() %>% 
                   rownames_to_column("Component") %>% 
                   dplyr::rename("Variance" = ".") %>%  
                   dplyr::mutate("random_model" = model_name) %>% 
                   dplyr::mutate(!!!meta) 
                   
   MK_llr_temp   = MK_model$llik %>% as.data.frame() %>%
                   dplyr::rename("llik" = ".") %>%  
                   dplyr::mutate("random_model" = model_name) %>% 
                   dplyr::mutate(!!!meta) %>% mutate("iter" = MK_model$niter)
  
   list("llrt" = MK_llr_temp, "var" = MK_var_temp) })

   MK_llr        = rbind(MK_llr, do.call(rbind, lapply(results, `[[`, "llrt")))
      
   MK_var        = rbind(MK_var, do.call(rbind, lapply(results, `[[`, "var"))) }, 
                   
   error = function(e) { cat("Error for trait", current_trait, ":", e$message, "\n") }) }}

   saveRDS(MK_llr, paste0(path_output, "MK_likelihoods_unpermuted.rds"))
   saveRDS(MK_var, paste0(path_output, "MK_variance_components_unpermuted.rds")) 
 
############################################################# 
#E) Variance partition (permuted)
############################################################# 
   
   MK_var        = data.frame()
   MK_llr        = data.frame()
   
   for(current_quantile in quantiles){
   print(paste("quantile is ", current_quantile))

   G_kernels     = readRDS(paste0(path_kinships, "random_effects/GRMs_permuted_", current_quantile,"_v2.rds"))
   
   for(current_trait in names(data_list)) {
          
   data_subset   = data_list[[current_trait]] %>% left_join(loads, by = "id") 
   
   y             = setNames(data_subset$effect, rownames(data_subset))
   
   X_form        = model.matrix(as.formula(X_input), data = data_subset)
                  
   KC            = E_kernels[["KC"]][as.character(data_subset$country), as.character(data_subset$country)]
   KR            = E_kernels[["KR"]][data_subset$region, data_subset$region]
   KY            = E_kernels[["KY"]][data_subset$year, data_subset$year] 
   
     for(current_permute in names(G_kernels)) {
     
   meta          = list("trait" = current_trait, "quantile" = current_quantile, "kernels" = "MK", "permuted" = "yes", "permutation" = current_permute) 
   
   GRM_S         = G_kernels[[current_permute]][["kinships"]][["Sel"]]
   GRM_N         = G_kernels[[current_permute]][["kinships"]][["Neu"]]
  
   GRM_G_ext     = GRM_G[data_subset$id, data_subset$id]
   GRM_S_ext     = GRM_S[data_subset$id, data_subset$id]
   GRM_N_ext     = GRM_N[data_subset$id, data_subset$id]

   V_M2_list     = list("GN" = GRM_N_ext, "GS" = GRM_S_ext, "GxC" = GRM_G_ext * KC, "CxY" = KC * KY)
   V_M3_list     = list("GN" = GRM_N_ext, "GS" = GRM_S_ext, "GNxC"= GRM_N_ext * KC, "GSxC" = GRM_S_ext * KC, "CxY" = KC * KY)

   if(length(unique(data_subset$region)) > 2) {
  
   V_M2_list[["CxR"]] = KC * KR
   V_M3_list[["CxR"]] = KC * KR }

   tryCatch({
  
   models_list   = list("M2" = V_M2_list, "M3" = V_M3_list)

   results       = lapply(names(models_list), function(model_name) {
  
   MK_model      = greml(y = y, X = X_form, GRM = models_list[[model_name]], ncores = 25)
                   
   MK_var_temp   = MK_model$theta %>% as.data.frame() %>% 
                   rownames_to_column("Component") %>% 
                   dplyr::rename("Variance" = ".") %>%  
                   dplyr::mutate("random_model" = model_name) 
                   
   MK_llr_temp   = MK_model$llik %>% as.data.frame() %>%
                   dplyr::rename("llik" = ".") %>%  
                   dplyr::mutate("random_model" = model_name) 
  
   list("llrt" = MK_llr_temp, "var" = MK_var_temp) })

   MK_llr        = rbind(MK_llr, do.call(rbind, lapply(results, `[[`, "llrt")) %>% mutate(!!!meta))
      
   MK_var        = rbind(MK_var, do.call(rbind, lapply(results, `[[`, "var")) %>% mutate(!!!meta)) }, 
                   
   error = function(e) { cat("Error for trait", current_trait, ":", e$message, "\n") }) }}}

   saveRDS(MK_llr, paste0(path_output, "MK_likelihoods_permuted.rds"))
   saveRDS(MK_var, paste0(path_output, "MK_variance_components_permuted.rds"))
    
#############################################################
#############################################################