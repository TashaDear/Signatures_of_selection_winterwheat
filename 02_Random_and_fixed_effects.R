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

path              = "/faststorage/project/breedfuture/prepared_data/submission2/"
path_raw          = "/faststorage/project/breedfuture/raw_data/WW/geno/plink/all_data/"

#Output
path_random       = paste0(path, "kinships/")
path_fixed        = paste0(path, "loads/")

#############################################################
#B) Load input
#############################################################

SNPs_exclude      = readLines(paste0(path_raw, "missing_SNPs_final.txt"))

SSA_output        = readRDS(paste0(path, "SSA/SSA_output_GMMAT_threshold_50_final.rds")) %>% 
                    dplyr::filter(converged == TRUE, !SNP %in% SNPs_exclude) 
                    
phenotypes_df     = readRDS(paste0(path, "phenotypes/prepared_data_list.rds"))[["GP"]] 

samples_GP        = phenotypes_df %>% pull(id) %>% unique()                    
                    
geno_input        = readRDS(paste0(path, "genotypes/geno_with_IDs_validation.rds"))[samples_GP, SSA_output$SNP] %>% as.matrix()

print(paste("SNPs in SSA output", nrow(SSA_output), "and SNPs in geno input", ncol(geno_input)))
stopifnot(all(SSA_output$SNP %in% colnames(geno_input)))
print(paste("Number of samples in geno input", nrow(geno_input)))

perms             = 1:500

quantiles         = c(seq(0.30, 0.90, by = 0.20), 0.95, 0.99)

thresholds        = SSA_output %>% 
                    dplyr::mutate(S_abs = abs(S)) %>% pull(S_abs) %>% 
                    quantile(probs = quantiles) %>% as.data.frame() %>% 
                    setNames("threshold") %>% rownames_to_column("quantile")
                    
#############################################################
#Check
#############################################################

SNPs_in_SSA       = readRDS(paste0(path, "SSA/SSA_output_GMMAT_threshold_50_final.rds")) %>%  
                    dplyr::filter(!SNP %in% SNPs_exclude) %>% pull(SNP) %>% length()
                    
SNPs_no_converge  = readRDS(paste0(path, "SSA/SSA_output_GMMAT_threshold_50_final.rds")) %>%  
                    dplyr::filter(converged != "TRUE", !SNP %in% SNPs_exclude) %>% pull(SNP) %>% length()
                    
final_SNPs        = readRDS(paste0(path, "SSA/SSA_output_GMMAT_threshold_50_final.rds")) %>%  
                    dplyr::filter(converged == TRUE, !SNP %in% SNPs_exclude) %>% pull(SNP) %>% length()
                    
print(paste("SNPs in SSA", SNPs_in_SSA, ", SNPs which do not converge", SNPs_no_converge, ", SNPs final", final_SNPs))

#############################################################
#C) Function (Mean partition/Fixed effects)
#############################################################

Prepare_loads     = function(input, geno, permute, seeds, thresholds_vec){ 
                    print(permute)
              
  loads_df        = data.frame()
             
  load_G          = rowSums(geno) %>% as.data.frame() %>% rownames_to_column("id") %>% rename("G" = ".") 
                                                                   
  for(current_quantile in unique(thresholds_vec$quantile)){
  print(paste("--- threshold", current_quantile))
  
  current_threshold     = thresholds_vec %>% dplyr::filter(quantile == current_quantile) %>% pull(threshold)
  
     for(current_seed in seeds){  
     set.seed(current_seed)
     
     current_seed       = paste0("perm_", current_seed)
  
     input$S_coef       = if (permute == "permuted") { sample(input$S, replace = FALSE) } else { input$S }
    
     Selected_positive  = input %>% dplyr::filter(S_coef >=  current_threshold)
     Selected_negative  = input %>% dplyr::filter(S_coef <= -current_threshold) 
      
    if(nrow(Selected_positive) > 0){
   
    B_load = rowSums(geno[, Selected_positive$SNP, drop = FALSE]) } else {
    B_load = rep(0, nrow(geno)) }

    if(nrow(Selected_negative) > 0){
   
    D_load = rowSums(geno[, Selected_negative$SNP, drop = FALSE]) } else {
    D_load = rep(0, nrow(geno)) }
    
    B_load = B_load %>% as.data.frame() %>% rownames_to_column("id") %>% rename("B" = ".")
    D_load = D_load %>% as.data.frame() %>% rownames_to_column("id") %>% rename("D" = ".")

     temp_df            = load_G %>% 
                          left_join(B_load, by = "id") %>% 
                          left_join(D_load, by = "id") %>% 
                          dplyr::mutate("permute" = permute, "permutation" = current_seed, 
                                        "quantile" = current_quantile, "threshold" = round(current_threshold, 2)) %>% 
                          dplyr::mutate(across(c(G, B, D), as.numeric)) 
                          
     loads_df           = rbind(temp_df, loads_df) }}
     
     saveRDS(loads_df, paste0(path_fixed, "loads_", permute,".rds"))  }
  
#############################################################

Prepare_loads(SSA_output, permute = "unpermuted", seeds = c(0),  geno = geno_input, thresholds_vec = thresholds)  
Prepare_loads(SSA_output, permute = "permuted",   seeds = perms, geno = geno_input, thresholds_vec = thresholds)  

#############################################################
#D) Function (Variance partition/Random effects)
#############################################################

Prepare_GRMs      = function(input, geno, permute, seeds, thresholds_vec) { 
   print(permute)
              
   for (current_quantile in unique(thresholds_vec$quantile)) {
   print(paste("--- threshold", current_quantile))
   
   current_threshold= thresholds_vec %>% dplyr::filter(quantile == current_quantile) %>% dplyr::pull(threshold)

   seed_results     = foreach(current_seed = seeds, .packages = "dplyr", .combine = 'c') %dopar% {
   
   set.seed(current_seed)
   
   seed_name        = if (permute == "permuted") { paste0("perm_", current_seed)} else { paste0("notperm_", current_seed) }
   
   input_upd        = input

   input_upd$S_coef = if (permute == "permuted") { sample(input_upd$S, replace = FALSE) } else { input_upd$S }

   Selected         = input_upd %>% dplyr::filter(abs(S_coef) >= current_threshold)
   sel_idx          = colnames(geno) %in% Selected$SNP

   seed_kin         = list()
   scaling          = list()

      for (geno_name in c("Sel", "Neu")) {

        X_input               = if (geno_name == "Sel") { geno[, sel_idx, drop = FALSE] } else { geno[, !sel_idx, drop = FALSE] }  

        X_center              = scale(X_input, center = TRUE, scale = FALSE)

        p                     = colMeans(X_input) / 2
        
        sum_2pq               = sum(2 * p * (1 - p))

        seed_kin[[geno_name]] = tcrossprod(X_center) / sum_2pq
        
        if (permute == "unpermuted") {
        
        scaling[[length(scaling)+1]] = data.frame("GRM" = geno_name, "scaling" = sum_2pq, 
                                                  "threshold" = current_threshold, "quantile" = current_quantile, 
                                                  "permute" = permute, "permutation" = seed_name, "no_snps" = ncol(X_input)) }}
      
       list_result                   = list()
       list_result[[seed_name]]      = list(kinships = seed_kin, scaling = scaling)
       list_result }
        
       saveRDS(seed_results, paste0(path_random, "GRMs_", permute, "_", current_quantile, ".rds"))
    
       gc() }

        if(permute == "unpermuted") {
        X_center = scale(geno, center = TRUE, scale = FALSE)
        p        = colMeans(geno) / 2
        sum_2pq  = sum(2 * p * (1 - p))
        GRM_full = tcrossprod(X_center) / sum_2pq

        saveRDS(GRM_full, paste0(path_random, "GRM_baseline.rds")) }}
  
#############################################################

Prepare_GRMs(SSA_output, permute = "unpermuted", seeds = c(0),  geno = geno_input, thresholds_vec = thresholds)  
Prepare_GRMs(SSA_output, permute = "permuted",   seeds = perms, geno = geno_input, thresholds_vec = thresholds)

#############################################################
#E) Function prepare kernels for environments
#############################################################

#Function
create_kernels    = function(input_data, column_name) {

labels            = unique(input_data[[column_name]])
  
K                 = diag(1, nrow = length(labels), ncol = length(labels))

rownames(K) <- colnames(K) <- labels
  
return(K) }

#Prepare kernels
kernel_cols       = c("country", "year", "region")

kernels           = lapply(kernel_cols, function(col) create_kernels(phenotypes_df, col))

names(kernels)    = c("KC", "KY", "KR") 

saveRDS(kernels, paste0(path_random, "ENV_kernels.rds"))
                 
#############################################################
stopCluster(cl)
#############################################################
