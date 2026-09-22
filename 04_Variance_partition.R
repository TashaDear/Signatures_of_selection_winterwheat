#############################################################
library(dplyr)
library(foreach)
library(tibble)
library(tidyr)
library(MM4LMM)
library(qgg)
library(sommer)
#############################################################
#A) Paths
############################################################# 

path            = "/home/tasha/breedfuture/"
path_input      = paste0(path, "prepared_data/submission2/")
path_output     = paste0(path, "output/submission2/")

#############################################################
#B) Inputs
############################################################# 

cores           = 30

quantiles       = c(paste0(seq(30, 90, by = 20), "%"), "95%", "99%")

X_input         = as.formula('effect ~ 1 + G + country + G:country')

GRM_G           = readRDS(paste0(path_input, "kinships/GRM_baseline.rds"))

E_kernels       = readRDS(paste0(path_input, "kinships/ENV_kernels.rds"))

data_list       = readRDS(paste0(path_input, "phenotypes/phenotypes_df.rds")) %>% dplyr::filter(analysis == "GP") %>% 
                  dplyr::mutate_at(vars("id", "year", "region", "country"), as.character) %>% split(.$trait)

loads           = readRDS(paste0(path_input, "loads/loads_unpermuted.rds")) %>% distinct(id, G)

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
#D) Variance partition-multi GRM (unpermuted)
############################################################# 
   
   MK_var        = data.frame()
   MK_llr        = data.frame()
   
   for(current_quantile in quantiles){
   print(paste0("quantile is ", current_quantile))

   G_kernels     = readRDS(paste0(path_input, "kinships/GRMs_unpermuted_", current_quantile,".rds"))[[1]][["kinships"]]
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
  
   MK_model      = greml(y = y, X = X_form, GRM = models_list[[model_name]], ncores = cores)
                   
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
#E) Variance partition- multi GRM (permuted)
############################################################# 
   
   MK_var        = data.frame()
   MK_llr        = data.frame()
   
   for(current_quantile in quantiles){
   print(paste("quantile is ", current_quantile))

   G_kernels     = readRDS(paste0(path_input, "kinships/GRMs_permuted_", current_quantile,".rds"))
   
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
  
   MK_model      = greml(y = y, X = X_form, GRM = models_list[[model_name]], ncores = cores)
                   
   MK_var_temp   = MK_model$theta %>% as.data.frame() %>% 
                   rownames_to_column("Component") %>% 
                   dplyr::rename("Variance" = ".") %>%  
                   dplyr::mutate("random_model" = model_name) 
                   
   MK_llr_temp   = MK_model$llik %>% as.data.frame() %>%
                   dplyr::rename("llik" = ".") %>%  
                   dplyr::mutate("random_model" = model_name) 
  
     list("llrt" = MK_llr_temp, "var" = MK_var_temp) })

   MK_llr        = rbind(MK_llr, do.call(rbind, lapply(results, `[[`, "llrt")) %>% dplyr::mutate(!!!meta))
      
   MK_var        = rbind(MK_var, do.call(rbind, lapply(results, `[[`, "var")) %>% dplyr::mutate(!!!meta)) }, 
                   
   error = function(e) { cat("Error for trait", current_trait, ":", e$message, "\n") }) }}}

   saveRDS(MK_llr, paste0(path_output, "MK_likelihoods_permuted.rds"))
   saveRDS(MK_var, paste0(path_output, "MK_variance_components_permuted.rds"))
    
#############################################################
#F)  Reliability of predicted genetic effects
#############################################################

reliability_df  = data.frame()

for(current_trait in names(data_list)) {
   
   data_subset  = data_list[[current_trait]] %>% left_join(loads, by = "id") %>% 
                  dplyr::mutate(id_country     = paste(.$id, .$country, sep = ":"),
                                country_year   = paste(.$country, .$year, sep = ":"),
                                country_region = paste(.$country, .$region, sep = ":")) 

   GRM_G_ext    = GRM_G_upd[as.character(data_subset$id), as.character(data_subset$id)]
  
   GxKC         = kronecker(GRM_G, E_kernels[["KC"]], make.dimnames = T)[data_subset$id_country, data_subset$id_country]
   CxY          = kronecker(E_kernels[["KC"]], E_kernels[["KY"]], make.dimnames = T)[data_subset$country_year, data_subset$country_year]
   CxR          = kronecker(E_kernels[["KC"]], E_kernels[["KR"]], make.dimnames = T)[data_subset$country_region, data_subset$country_region]

     if (length(unique(data_subset$region)) > 2) {
   
    fit         = sommer::mmer(fixed  = effect ~ 1 + G + country + G:country,
                               random = ~ vsr(id, Gu = GRM_G_ext) + vsr(id_country, Gu = GxKC) + 
                                          vsr(country_year, Gu = CxY) +vsr(country_region, Gu = CxR),
                               rcov   = ~ units,
                               data   = data_subset,
                               getPEV = TRUE,
                               dateWarning = FALSE) 
                                
                                                   } else {
    
  fit           = sommer::mmer(fixed  = effect ~ 1 + G + country + G:country,
                               random = ~ vsr(id, Gu = GRM_G_ext) + vsr(id_country, Gu = GxKC) + vsr(country_year, Gu = CxY),
                               rcov   = ~ units,
                               data   = data_subset,
                               getPEV = TRUE,
                               dateWarning = FALSE)  }
                                
  PEV_u_diag     = diag(fit$PevU$`u:id`$effect)

  temp_output    = data.frame("trait" = current_trait, "R" = 1 - mean(PEV_u_diag, na.rm = TRUE) / fit$sigma$`u:id`[1, 1])
  
  reliability_df = rbind(temp_output, reliability_df) }

   saveRDS(reliability_df, paste0(path_output, "SK_reliability_unpermuted.rds"))

#############################################################
#############################################################
