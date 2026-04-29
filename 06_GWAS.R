#############################################################
library("data.table")
library("dplyr")
library("Matrix")
library("MM4LMM") 
library("tidyverse")
#############################################################
#Part A) Paths
#############################################################

submission      = "submission1"

path            = "/home/tasha/breedfuture/"
path_input      = paste0(path, "prepared_data/v2/")
path_output     = paste0(path, "results/v2/model_output/GWAS/")
path_geno       = paste0(path, "results/v2/fitness_estimates/")
path_kinships   = paste0(path, "results/v2/kinships/")

#############################################################
#Part B) Load input 
#############################################################

data_list       = readRDS(paste0(path_input, "phenotypes.rds")) %>% 
                  dplyr::filter(analysis == "GWAS") %>% 
                  dplyr::mutate_at(vars("id", "year", "region", "country"), as.character) %>% split(.$trait)

E_kernels       = readRDS(paste0(path_input, "environment_kernels.rds"))
                  
loads           = readRDS(paste0(path_kinships, "fixed_effects/loads_unpermuted_v2.rds")) %>% distinct(id, G)
                  
geno            = readRDS(paste0(path_geno, "geno.rds")) 

GRM_G           = readRDS(paste0(path_kinships, "random_effects/GRM_baseline_v2.rds"))
GRM_G_upd       = as.matrix(Matrix::nearPD(GRM_G)$mat)

#############################################################
#Part C) GWAS
#############################################################

run_GWAS        = function(ver) {

   output       = list()
  
   for(current_trait in names(data_list)){
   print(current_trait)

   data_subset      = data_list[[current_trait]] %>% left_join(loads, by = "id") %>% 
                      dplyr::filter(country == "DK", !is.na(effect))
   
   samples          = data_subset$id %>% unique()
   
   geno_sub         = geno[samples, ]

   maf              = colMeans(geno_sub, na.rm = TRUE) / 2

   geno_final       = geno_sub[, maf >= 0.05 & maf <= 0.95, drop = FALSE]
                
   X_input          = as.formula('~ 1 + G + Xeffect + PC1 + PC2 + PC3')  
   
   KR               = E_kernels[["KR"]][data_subset$region, data_subset$region]
   KY               = E_kernels[["KY"]][data_subset$year, data_subset$year] 
   GRM_G_ext        = GRM_G_upd[data_subset$id, data_subset$id]           
   
          if(ver == "A") {         
   
   V_list           = list("G" = GRM_G_ext, "R" = KR, "Y" = KY, "error" = diag(1, nrow(data_subset))) 
   
   } else if(ver == "B") {
   
   V_list           = list("G" = GRM_G_ext, "R" = KR, "Y" = KY, "GxY" = GRM_G_ext * KY, "error" = diag(1, nrow(data_subset))) 
   
   } else if (ver == "C") {
   
   V_list           = list("G" = GRM_G_ext, "R" = KR, "Y" = KY, "GxY" = GRM_G_ext * KY, "RxY" = KR * KY, "error" = diag(1, nrow(data_subset))) }  
                                                                                             
   tryCatch({ 
   
   output[[current_trait]] = MMEst(Y = data_subset$effect, Cofactor= data_subset, X = geno_final, 
                             formula = X_input, VarList = V_list, CritVar = 10e-5, NbCores = 25) }, error = function(e) { print(paste("Error:", e)) }) }
   
   saveRDS(output, paste0(path_output,"GWAS_marginal_past_", ver,"_nearPD.rds")) }
   
#############################################################
#Part D) Run function
#############################################################

run_GWAS("A")   

#############################################################
#############################################################
