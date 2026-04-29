############################################################# 
library("data.table")
library("dplyr")
library("tidyverse")
library("tidyr")
library("tibble")

library(parallel)    
library(doParallel)  

n_workers = 4
cl = makeCluster(n_workers, type = "FORK")
registerDoParallel(cl)
#############################################################
#A) Paths 
#############################################################

path              = "/faststorage/project/breedfuture/results/v2/"
path_input        = paste0(path, "fitness_estimates/")
path_output       = paste0(path, "kinships/")

#############################################################
#B) Load input
#############################################################

geno_input        = readRDS(paste0(path_input, "geno.rds")) 

SSA_output        = readRDS(paste0(path_input, "SSA_output.rds"))

data_list         = readRDS(paste0(path_pheno, "phenotypes.rds")) 

permutations      = 1:250

quantiles         = c(seq(0.30, 0.90, by = 0.20), 0.95, 0.99)

thresholds        = SSA_output %>% 
                    dplyr::mutate(S_abs = abs(S)) %>% pull(S_abs) %>% 
                    quantile(probs = quantiles) %>% as.data.frame() %>% 
                    setNames("threshold") %>% rownames_to_column("quantile")

#############################################################
#C) Function (Mean partition/Fixed effects)
#############################################################

Prepare_loads     = function(input, geno, permute, seeds, thresholds_vec){ 
  print(permute)
              
  loads_df        = data.frame()
             
  load_MA         = rowSums(geno) %>% as.data.frame() %>% rownames_to_column("id") %>% rename("G" = ".") 
                                                                   
  for(current_quantile in unique(thresholds_vec$quantile)){
  print(paste("threshold", current_quantile))
  
  current_threshold     = thresholds_vec %>% dplyr::filter(quantile == current_quantile) %>% pull(threshold)
  
     for(current_seed in seeds){  
     set.seed(current_seed)
     
     current_seed       = paste0("perm_", current_seed)
  
     input$S_coef       = if (permute == "permuted") { sample(input$S, replace = FALSE) } else { input$S }
    
     Selected_positive  = input %>% dplyr::filter(S_coef >=  current_threshold)
     Selected_negative  = input %>% dplyr::filter(S_coef <= -current_threshold) 
      
     B_load             = rowSums(geno[, Selected_positive$SNP]) %>% as.data.frame() %>% rownames_to_column("id") %>% rename("B" = ".") 
     D_load             = rowSums(geno[, Selected_negative$SNP]) %>% as.data.frame() %>% rownames_to_column("id") %>% rename("D" = ".")

     temp_df            = load_MA %>% 
                          left_join(B_load, by = "id") %>% 
                          left_join(D_load, by = "id") %>% 
                          dplyr::mutate("permute" = permute, "permutation" = current_seed, 
                                        "quantile" = current_quantile, "threshold" = round(current_threshold, 2))
                          
     loads_df           = rbind(temp_df, loads_df) }
  
     loads_df           = loads_df %>% dplyr::mutate(across(c(G, B, D), as.numeric)) 

     saveRDS(loads_df, paste0(path_output, "fixed_effects/loads_", permute,"_v2.rds")) }}
  
#############################################################

Prepare_loads(SSA_output, permute = "unpermuted", seeds = c(0), geno = geno_input, thresholds_vec = thresholds)  
Prepare_loads(SSA_output, permute = "permuted", seeds = permutations, geno = geno_input, thresholds_vec = thresholds)  

#############################################################
#D) Function (Variance partition/Random effects)
#############################################################

Prepare_GRMs      = function(input, geno, permute, seeds, thresholds_vec) { 
   print(permute)

   geno_full        = geno %>% as.data.frame() %>% dplyr::mutate(across(everything(), as.numeric)) %>% as.matrix()
                           
   for (current_quantile in unique(thresholds_vec$quantile)) {
   print(current_quantile)

   current_threshold= thresholds_vec %>% dplyr::filter(quantile == current_quantile) %>% dplyr::pull(threshold)

   seed_results     = foreach(current_seed = seeds, .packages = "dplyr", .combine = 'c') %dopar% {
   
   set.seed(current_seed)
   
   seed_name        = if (permute == "permuted") { paste0("perm_", current_seed)} else { paste0("notperm_", current_seed)}

   input_upd        = input
   input_upd$S_coef = if (permute == "permuted") { sample(input_upd$S, replace = FALSE) } else { input_upd$S }

   Selected         = input_upd %>% dplyr::filter(abs(S_coef) >= current_threshold)
   sel_idx          = colnames(geno_full) %in% Selected$SNP

   seed_kinships    = list()
   seed_scaling     = list()

      for (geno_name in c("Sel", "Neu")) {

        G        = if (geno_name == "Sel") { geno_full[, sel_idx, drop = FALSE] } else { geno_full[, !sel_idx, drop = FALSE] }  

        X_center = scale(G, center = TRUE, scale = FALSE)

        p        = colMeans(G) / 2
        
        sum_2pq  = sum(2 * p * (1 - p))
        
        GRM      = tcrossprod(X_center) / sum_2pq

        seed_kinships[[geno_name]] = GRM
        
        if (permute == "unpermuted") {
        
        seed_scaling[[length(seed_scaling)+1]]  = data.frame("GRM" = geno_name, "scaling" = sum_2pq, "threshold" = current_threshold,
                                                             "quantile" = current_quantile, "permute" = permute, "permutation" = seed_name,
                                                             "no_snps" = ncol(G)) }}
      
       list_result = list()
       list_result[[seed_name]] <- list(kinships = seed_kinships, scaling = seed_scaling)
       list_result }
        
       saveRDS(seed_results, paste0(path_output, "random_effects/GRMs_", permute, "_", current_quantile, "_v2.rds"))
    
       gc() }

        if(permute == "unpermuted") {
        X_center = scale(geno_full, center = TRUE, scale = FALSE)
        p        = colMeans(geno_full) / 2
        sum_2pq  = sum(2 * p * (1 - p))
        GRM_full = tcrossprod(X_center) / sum_2pq

        saveRDS(GRM_full, paste0(path_output, "random_effects/GRM_baseline_v2.rds")) }}
  
#############################################################

Prepare_GRMs(SSA_output, permute = "unpermuted", seeds = c(0), geno = geno_input, thresholds_vec = thresholds)  
Prepare_GRMs(SSA_output, permute = "permuted", seeds = permutations, geno = geno_input, thresholds_vec = thresholds)  

#############################################################
#E) Function prepare kernels for environments
#############################################################

#Function
create_kernels    = function(input_data, column_name) {

labels            = unique(input_data[[column_name]])
  
K                 = diag(1, nrow = length(labels), ncol = length(labels))

rownames(K) <- colnames(K) <- labels
  
return(K)}

#Prepare kernels
kernel_cols       = c("country", "year", "region")

kernels           = lapply(kernel_cols, function(col) create_kernels(data_list, col))

names(kernels)    = c("KC", "KY", "KR") 

saveRDS(kernels, paste0(path_pheno, "environment_kernels.rds")) #maybe save at kinships next time
                 
#############################################################
stopCluster(cl)
#############################################################