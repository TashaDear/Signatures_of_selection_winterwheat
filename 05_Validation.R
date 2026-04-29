#############################################################
library("dplyr")
library("foreach") 
library("tidyr")
library("tibble")
library("stringr")
library("qgg")
library("purrr")
#############################################################
#Paths 
#############################################################

submission      = "submission1"

path            = "/home/tasha/breedfuture/"
path_input      = paste0(path, "prepared_data/v2/")
path_output     = paste0(path, "results/v2/model_output/validation/")
path_kinships   = paste0(path, "results/v2/kinships/")

#############################################################
#B) Inputs
############################################################# 

quantiles          = c(paste0(seq(30, 90, by = 20), "%"), "95%", "99%")

GRM_G              = readRDS(paste0(path_kinships, "random_effects/GRM_baseline_v2.rds"))
GRM_G_upd          = as.matrix(Matrix::nearPD(GRM_G)$mat)

E_kernels          = readRDS(paste0(path_input, "environment_kernels.rds"))

data_list          = readRDS(paste0(path_input, "phenotypes.rds")) %>% 
                     dplyr::filter(analysis == "GP") %>% 
                     dplyr::mutate_at(vars("id", "year", "region", "country"), as.character) %>% split(.$trait)

loads_unpermuted   = readRDS(paste0(path_kinships, "fixed_effects/loads_unpermuted_v2.rds")) 
loads_permuted     = readRDS(paste0(path_kinships, "fixed_effects/loads_permuted_v2.rds")) 

#############################################################
#B) Helper function 
#############################################################

#Prepare model combinations
fixed_models       = c("X1", "X2", "X3")
random_models      = c("M1", "M2", "M3")
model_grid         = expand.grid(fixed_model = fixed_models, random_model = random_models, stringsAsFactors = FALSE)

#function
run_prediction     = function(data_input, fixed_model, random_model, V_input, validate_input) {

y                  = setNames(data_input$effect, rownames(data_input))
  
X_formula          = switch(fixed_model,
                            X1 = "~ 1 + G + country + G:country",
                            X2 = "~ 1 + G + B + D + country + G:country", 
                            X3 = "~ 1 + G + B + D + country + G:country + B:country + D:country",
                            stop("Unsupported fixed_model"))
                                   
X_input            = model.matrix(as.formula(X_formula), data = data_input)

GRM_list           = V_input[[random_model]]
 
model              = tryCatch({ greml(y = y, X = X_input, GRM = GRM_list, validate = validate_input, ncores = 25) }, 
                     error = function(e) { cat("Error ", fixed_model, "/", random_model, ":", conditionMessage(e), "\n")
                     return(NULL) })
  
if (is.null(model)) return(NULL)
  
temp_output        = tryCatch({ 
                     model$accuracy %>% as.data.frame() %>% 
                     rownames_to_column("fold") %>%
                     dplyr::mutate("fixed_model" = fixed_model, "random_model" = random_model)  }, 

error = function(e) { cat("Error for", fixed_model, "/", random_model, ":", conditionMessage(e), "\n")

return(NULL) })  

return(temp_output) }

#############################################################
#C) Validation (unpermuted)
#############################################################      

validation_unperm   = function() {                  
  
  final_output      = data.frame()

for(current_quantile in quantiles) {
print(paste("quantile is ", current_quantile))
    
  G_kernels         = readRDS(paste0(path_kinships, "random_effects/GRMs_unpermuted_", current_quantile,"_v2.rds"))[["notperm_0"]][["kinships"]] 
    GRM_S           = as.matrix(Matrix::nearPD(G_kernels[["Sel"]])$mat) 
    GRM_N           = as.matrix(Matrix::nearPD(G_kernels[["Neu"]])$mat) 
  
  loads_sub         = loads_unpermuted %>% dplyr::filter(quantile == current_quantile)
  
     for(current_trait in names(data_list)) {
     print(current_trait)
  
  data_subset       = data_list[[current_trait]] %>% 
                      left_join(loads_sub, by = "id") %>% 
                      dplyr::mutate(genetic_cluster = as.factor(paste0("cluster_", cluster)))
                  
  validate_list     = split(seq_len(nrow(data_subset)), data_subset$genetic_cluster) 
   
    KC              = E_kernels[["KC"]][data_subset$country, data_subset$country]
    KR              = E_kernels[["KR"]][data_subset$region, data_subset$region]
    KY              = E_kernels[["KY"]][data_subset$year, data_subset$year] 
  
    GRM_G_ext       = GRM_G_upd[data_subset$id, data_subset$id]
    GRM_S_ext       = GRM_S[data_subset$id, data_subset$id]
    GRM_N_ext       = GRM_N[data_subset$id, data_subset$id] 
  
    V_M1_list       = list("G"  = GRM_G_ext, "GxC"= GRM_G_ext * KC, "CxY" = KC * KY)                             
    V_M2_list       = list("GN" = GRM_N_ext, "GS" = GRM_S_ext, "GxC" = GRM_G_ext * KC, "CxY" = KC * KY)
    V_M3_list       = list("GN" = GRM_N_ext, "GS" = GRM_S_ext, "GNxC"= GRM_N_ext * KC, "GSxC" = GRM_S_ext * KC, "CxY" = KC * KY)
    
    if(length(unique(data_subset$region)) > 2) {
   
      V_M1_list[["CxR"]]= KC * KR
      V_M2_list[["CxR"]]= KC * KR
      V_M3_list[["CxR"]]= KC * KR }
  
      V_list        = list("M1" = V_M1_list, "M2" = V_M2_list, "M3" = V_M3_list)
                           
      out           = pmap_dfr(model_grid, function(fixed_model, random_model) {
                      run_prediction("data_input" = data_subset, "fixed_model" = fixed_model, "random_model" = random_model, 
                                     "V_input" = V_list, "validate_input" = validate_list)})

      temp_output   = out %>% dplyr::mutate("trait" = current_trait, "quantile" = current_quantile, "permuted" = "no", "permutation" = 0)
      
      final_output  = rbind(final_output, temp_output) }}
      
saveRDS(final_output, paste0(path_output, "validation_genetic_cluster_unpermuted.rds")) }

#############################################################
#Run
validation_unperm()

#############################################################
#C) Validation (permuted)
#############################################################  

validation_perm    = function() {

 final_output      = data.frame()

for(current_quantile in quantiles) {
print(paste("quantile is ", current_quantile))
    
  G_kernels         = readRDS(paste0(path_kinships, "random_effects/GRMs_permuted_", current_quantile,"_v2.rds")) 
  
  loads_sub         = loads_permuted %>% dplyr::filter(quantile == current_quantile)
  
     for(current_trait in names(data_list)) {
     print(current_trait)
  
  data_subset       = data_list[[current_trait]] %>% 
                      dplyr::mutate(genetic_cluster = as.factor(paste0("cluster_", cluster)))
                  
  validate_list     = split(seq_len(nrow(data_subset)), data_subset$genetic_cluster) 
      
    KC              = E_kernels[["KC"]][data_subset$country, data_subset$country]
    KR              = E_kernels[["KR"]][data_subset$region, data_subset$region]
    KY              = E_kernels[["KY"]][data_subset$year, data_subset$year] 
    
    GRM_G_ext       = GRM_G_upd[data_subset$id, data_subset$id]
    
    V_M1_list       = list("G"  = GRM_G_ext, "GxC"= GRM_G_ext * KC, "CxY" = KC * KY)    
    
    for(current_permute in unique(loads_perm$permutation)) {
    
    loads_final     = loads_sub %>% dplyr::filter(permutation == current_permute)
    
    GRM_S_ext       = G_kernels[[current_permute]][["kinships"]][["Sel"]][data_subset$id, data_subset$id]
    GRM_N_ext       = G_kernels[[current_permute]][["kinships"]][["Neu"]][data_subset$id, data_subset$id]
                               
    V_M2_list       = list("GN" = GRM_N_ext, "GS" = GRM_S_ext, "GxC" = GRM_G_ext * KC, "CxY" = KC * KY)
    V_M3_list       = list("GN" = GRM_N_ext, "GS" = GRM_S_ext, "GNxC"= GRM_N_ext * KC, "GSxC" = GRM_S_ext * KC, "CxY" = KC * KY)
    
    if(length(unique(data_subset$region)) > 2) {
   
      V_M1_list[["CxR"]]= KC * KR
      V_M2_list[["CxR"]]= KC * KR
      V_M3_list[["CxR"]]= KC * KR }
  
      V_list        = list("M1" = V_M1_list, "M2" = V_M2_list, "M3" = V_M3_list)
      
      data_final    = data_subset %>% 
                      dplyr::mutate(.row_order = row_number()) %>% left_join(loads_final, by = "id") %>% arrange(.row_order) 
                                               
      out           = pmap_dfr(model_grid, function(fixed_model, random_model) {
                      run_prediction("data_input" = data_final, "fixed_model" = fixed_model, "random_model" = random_model, 
                                     "V_input" = V_list, "validate_input" = validate_list)})

      temp_output   = out %>% 
                      dplyr::mutate("trait" = current_trait, "quantile" = current_quantile, "permuted" = "yes", "permutation" = current_permute)
      
      final_output  = rbind(final_output, temp_output) }}
      
saveRDS(final_output, paste0(path_output, "validation_genetic_cluster_permuted.rds")) }

#############################################################
#Run
validation_perm()

#############################################################
#############################################################