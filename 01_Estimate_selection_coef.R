#############################################################
library(gaston)
library(GMMAT)
library(dplyr)
library(gridExtra)
library(tidyverse)
library(tidyr)
library(future.apply)

options(future.globals.maxSize = 100 * 1024^3)
#############################################################
#A) Path(s)
#############################################################

set.seed(0)

path              = "/home/tasha/breedfuture/"
path_input        = paste0(path, "prepared_data/submission2/genotypes/")
path_output       = paste0(path, "prepared_data/submission2/SSA/")

#############################################################
#B) Inputs
#############################################################

SSA_input          = readRDS(paste0(path_out, "SSA/SSA_input_final.rds"))

#Function
prepare_GRM       = function(input) {

   X              = input %>% dplyr::select(-c("crossing_year", "PC1", "PC2", "PC3")) %>% column_to_rownames("id") %>% as.matrix()  
   X[X == 1]      = 2
   X[X == 0.5]    = 1    
   X_center       = scale(X, center = TRUE, scale = FALSE)
   p              = colMeans(X) / 2
   sum_2pq        = sum(2 * p * (1 - p))
   G_GRM_upd      = tcrossprod(X_center) / sum_2pq

   return(G_GRM_upd) }
   
SSA_GRM           = SSA_input %>% prepare_GRM()

#############################################################
#C) Create function for inference of selection coefficients
#############################################################

SSA_analysis       = function(input, model, GRM, block_size = 50, workers = 20, threshold) {

  pipeline_start   = Sys.time()

  All_SNPs         = setdiff(colnames(input), c("id","crossing_year","PC1","PC2","PC3"))
  
  input_matrix     = input %>% column_to_rownames("id")

  X                = model.matrix(~ crossing_year + PC1 + PC2 + PC3, data = input_matrix)
  coef_idx         = match("crossing_year", colnames(X))
  
  dat              = input %>% dplyr::select(c("id","crossing_year","PC1","PC2","PC3"))
  
  input            = input[match(dat$id, input$id), ]
  
  weight_vec       = rep(2, nrow(dat))

  GRM_upd          = GRM[dat$id, dat$id]
  
  GRM_upd          = as.matrix(Matrix::nearPD(GRM_upd)$mat)

  snp_blocks       = split(All_SNPs, ceiling(seq_along(All_SNPs) / block_size))

  plan(multicore, workers = workers)

  results          = future_lapply(seq_along(snp_blocks), function(b) {

    block          = snp_blocks[[b]]

    cat("Block", b, "of", length(snp_blocks), "\n")

    block_results  = vector("list", length(block))

    for(j in seq_along(block)) {

      SNP        = block[j]

      y          = input[[SNP]]
      
      dat$y      = y
      
      fit     = tryCatch({

               if(model == "GLM") {

      fit        = glm.fit(x = X, y = y, weights = weight_vec, family = binomial())

        } else if(model == "GMMAT") {
        
     fit         = GMMAT::glmmkin(fixed = y ~ crossing_year + PC1 + PC2 + PC3, data = dat, kins = GRM_upd, id = "id", 
                                  family = binomial(), weights = weight_vec, maxiter = 1000)
                               
        } else if(model == "GASTON") {
        
     y[y == 0.5] = NA

     keep        = !is.na(y)

     fit         = gaston::logistic.mm.aireml(Y = y[keep], X = X[keep, ], K = GRM_upd[keep, keep], max.iter = 500,
                   verbose = FALSE) } }, error = function(e) NULL)

      if(is.null(fit))
        next

      if(model == "GLM") {

        beta = fit$coefficients[coef_idx]
        vcov = summary.glm(fit)$cov.scaled
        SE   = sqrt(vcov[coef_idx, coef_idx])
        con  = NA

      } else if(model == "GMMAT") {

        beta = fit$coefficients["crossing_year"]
        SE   = sqrt(diag(fit$cov))["crossing_year"]
        con  = fit$converged

      } else if(model == "GASTON") {

        beta = fit$BLUP_beta[coef_idx]
        SE   = sqrt(diag(fit$varbeta))[coef_idx]
        con  = fit$niter

      }

      z      = beta / SE

    block_results[[j]] = data.frame(model = model, converged = con, SNP = SNP, Effect = beta, SE = SE, S = z, Pval = 2 * pnorm(-abs(z))) }
    
    cat("Finished block", b, "of", length(snp_blocks), "\n")

    data.table::rbindlist(block_results, fill = TRUE) }, future.seed = TRUE)

  plan(sequential)

  final_output         = data.table::rbindlist(results)
  
  saveRDS(final_output, paste0(path_output, "SSA_output_", model, "_threshold_50_final.rds")) }
  
#############################################################
#D) Run function
#############################################################

SSA_analysis(SSA_input, "GMMAT", SSA_GRM, block_size = 200, workers = 41, threshold = 50) 

#############################################################
#############################################################
