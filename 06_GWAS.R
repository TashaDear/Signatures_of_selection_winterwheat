#############################################################
library("data.table")
library("dplyr")
library("Matrix")
library("MM4LMM") 
library("tidyverse")
#############################################################
#A) Paths
#############################################################

path            = "/home/tasha/breedfuture/"
path_input      = paste0(path, "prepared_data/submission2/")
path_output     = paste0(path, "output/submission2/")

path_raw        = "/faststorage/project/breedfuture/raw_data/WW/geno/plink/all_data/"

#############################################################
#B) Input 
#############################################################

SNPs_exclude    = readLines(paste0(path_raw, "missing_SNPs_final.txt"))

data_list       = readRDS(paste0(path_input, "phenotypes/phenotypes_df.rds")) %>% dplyr::filter(analysis == "GWAS") %>% 
                  dplyr::mutate_at(vars("id", "year", "region", "country"), as.character) %>% split(.$trait)

SSA_output      = readRDS(paste0(path_input, "SSA/SSA_output_GMMAT_threshold_50_final.rds")) %>% 
                  dplyr::filter(converged == TRUE, !SNP %in% SNPs_exclude)

#############################################################
#C) Prepare kernels
#############################################################
                
phenotypes_df    = readRDS(paste0(path_input, "phenotypes/prepared_data_list.rds"))[["GWAS"]] 

samples_GWAS     = phenotypes_df %>% pull(id) %>% unique() 
                                
geno             = readRDS(paste0(path_input, "genotypes/geno_with_IDs_validation.rds"))[samples_GWAS, SSA_output$SNP] %>% as.matrix()

G_load           = rowSums(geno) %>% as.data.frame() %>% rownames_to_column("id") %>% rename("G" = ".") 

#Create GRM for GWAS samples
X_center         = scale(geno, center = TRUE, scale = FALSE)
p                = colMeans(geno) / 2
sum_2pq          = sum(2 * p * (1 - p))
GRM_G            = tcrossprod(X_center) / sum_2pq
GRM_G_upd        = as.matrix(Matrix::nearPD(GRM_G)$mat)

#Kernels ENG
create_kernels    = function(input_data, column_name) {

labels            = unique(input_data[[column_name]])
  
K                 = diag(1, nrow = length(labels), ncol = length(labels))

rownames(K) <- colnames(K) <- labels
  
return(K) }

#Prepare kernels
kernel_cols       = c("country", "year", "region")

E_kernels         = lapply(kernel_cols, function(col) create_kernels(phenotypes_df, col))

names(E_kernels)  = c("KC", "KY", "KR") 

#############################################################
#D) Run GWAS
#############################################################

   output       = list()
  
   for(current_trait in names(data_list)){
   print(current_trait)

   data_subset      = data_list[[current_trait]] %>% left_join(G_load, by = "id") %>% 
                      dplyr::filter(country == "DK", !is.na(effect))
      
   geno_upd         = geno[data_subset$id, ]
                
   X_input          = as.formula('~ 1 + G + Xeffect + PC1 + PC2 + PC3') 
   
   KR               = E_kernels[["KR"]][data_subset$region, data_subset$region]
   KY               = E_kernels[["KY"]][data_subset$year, data_subset$year] 
   GRM_G_ext        = GRM_G_upd[data_subset$id, data_subset$id]           
         
   V_list           = list("G" = GRM_G_ext, "R" = KR, "Y" = KY, "error" = diag(1, nrow(data_subset))) 
                                                                                   
   tryCatch({ 
   
   output[[current_trait]] = MMEst(Y = data_subset$effect, Cofactor= data_subset, X = geno_upd, Method = "Reml",
                             formula = X_input, VarList = V_list, CritVar = 10e-5, NbCores = 25) }, error = function(e) { print(paste("Error:", e)) }) }
   
   saveRDS(output, paste0(path_output, "GWAS_marginal_final_submission2.rds")) 

#############################################################
#############################################################
