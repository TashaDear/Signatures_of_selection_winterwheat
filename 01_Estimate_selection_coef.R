#############################################################
library(dplyr)
library(gridExtra)
library(tidyverse)
library(tidyr)
#############################################################
#A) Path(s)
#############################################################

path      = "/faststorage/project/breedfuture/results/v2/fitness_estimates/"

#############################################################
#B) Inputs
#############################################################

SSA_input = readRDS(paste0(path, "SSA_input.rds"))                                        

#############################################################
#C) Create function for inference of selection coefficients
#############################################################

Estimate_fitness_effects   = function(input, ploidy){

         if (ploidy == 1) {
  
  input_upd                = input %>% mutate(across(everything(), ~ replace(., . == 1, NA)),  across(everything(), ~ replace(., . == 2, 1)))                   
  } else if (ploidy == 2) {
  
  input_upd                = input %>% mutate(across(everything(), ~ replace(., . == 1, 0.5)), across(everything(), ~ replace(., . == 2, 1)))  
  
  }
  
  All_SNPs                 = setdiff(colnames(input_upd), c("id", "crossing_year", "PC1", "PC2", "PC3"))
  
  final_output             = data.frame()
  
  for(SNP in All_SNPs) {
  
  cat("SNP number:", match(SNP, colnames(input_upd)), "\n")  
    
  subset_data             = input_upd %>% rename("SNP" = SNP) %>% 
                            dplyr::mutate(SNP = as.integer(SNP), crossing_year = as.numeric(crossing_year)) 
                                   
  weight_vec              = rep(ploidy, nrow(subset_data))                          
    
  GLM_PC                  = tryCatch({glm(SNP ~ crossing_year + PC1 + PC2 + PC3, family = binomial(link = "logit"), 
                            weights = weight_vec, data = subset_data) }, error = function(e)
     
  { cat("Error occurred for SNP:", SNP, "\n") })
    
  GLM_PC_term             = which(names(GLM_PC$coefficients) %in% c("crossing_year"))
    
  GLM_PC_pval             = aod::wald.test(Sigma = vcov(GLM_PC), b = coef(GLM_PC), Terms = GLM_PC_term)$result$chi2["P"]
    
  temp_output             = data.frame("SNP"     = SNP,
                                       "SE"      = sqrt(vcov(GLM_PC)["crossing_year", "crossing_year"]),
                                       "Effect"  = coef(GLM_PC)["crossing_year"],
                                       "Pval"    = GLM_PC_pval) %>%
                           dplyr::mutate("S" = Effect/SE)
        
  final_output            = rbind(temp_output, final_output) }

  saveRDS(final_output, paste0(path, "SSA_output.rds")) }   
   
#############################################################
#D) Run function
#############################################################

Estimate_fitness_effects(SSA_input, ploidy = 1)   
 
############################################################# 