#############################################################
library(dplyr)
library(foreach)
library(doParallel)
library(tibble)
library(tidyr)
#############################################################
#A) Paths 
#############################################################

path               = "/home/tasha/breedfuture/"
path_input         = paste0(path, "prepared_data/submission2/")
path_output        = paste0(path, "output/submission2/")
path_map           = paste0(path, "raw_data/WW/geno/plink/all_data/")
path_raw           = "/faststorage/project/breedfuture/raw_data/WW/geno/plink/all_data/"


#############################################################
#B) Input
#############################################################
SNPs_exclude       = readLines(paste0(path_raw, "missing_SNPs_final.txt"))

samples            = readRDS(paste0(path_input, "phenotypes/phenotypes_df.rds")) %>% pull(id) %>% unique()

SSA_output         = readRDS(paste0(path_input, "SSA/SSA_output_GMMAT_threshold_50_final.rds")) %>% 
                     dplyr::filter(converged == TRUE,  !SNP %in% SNPs_exclude)
                     
geno_input         = readRDS(paste0(path_input, "genotypes/geno_with_IDs_validation.rds"))[samples, SSA_output$SNP] %>% as.matrix()

map                = read.table(paste0(path_map, "ww_updated_final.map"), header= F,col.names=c("CHROM","SNP", "DIST", "POS")) %>% as.data.frame() %>% 
                     dplyr::filter(SNP %in% colnames(geno_input))
                                
chromosomes        = unique(map$CHROM) %>% na.omit() %>% as.vector()

#############################################################
#C) Create GRM
#############################################################
 
X_center           = scale(geno_input, center = TRUE, scale = FALSE)
p                  = colMeans(geno_input) / 2
sum_2pq            = sum(2 * p * (1 - p))
GRM_G              = tcrossprod(X_center) / sum_2pq
GRM_G_upd          = as.matrix(Matrix::nearPD(GRM_G)$mat)
                                          
#############################################################        
#D) Functions
#############################################################

invSqrt            = function(mat) {
stopifnot(isSymmetric(mat))
  
ED                 = eigen(mat) 
L                  = ED$vectors %*% diag(1/sqrt(ED$values)) %*% t(ED$vectors)
  
return(L) }              

#############################################################        
#E) Prepare matrices
#############################################################

Q                  = prcomp(geno_input, center = T, scale = F)$x[,1:3]
Q                  = cbind(Intercept = 1, Q)

genotypes          = rownames(Q)
n                  = length(genotypes)
GRM_G_upd          = GRM_G_upd[genotypes, genotypes]
diag(GRM_G_upd)    = diag(GRM_G_upd) + 1e-5

# Projection matrices
P                  = list()
P[["Int"]]         = diag(n) - matrix(1, nrow=n, ncol=n)/n
P[["Q"]]           = diag(n) - Q %*% solve(crossprod(Q)) %*% t(Q)

Ginv               = solve(GRM_G_upd)
H                  = Q %*% solve(t(Q) %*% Ginv %*% Q) %*% t(Q) %*% Ginv
P[["Q+K"]]         = invSqrt(GRM_G_upd) %*% (diag(n) - H)

max_dist           = 10000000

############################################################
#C) LDdecay
#############################################################

output             = list()

for (chrom in chromosomes) {
print(chrom)

current_SNPs        = map %>% dplyr::filter(CHROM == chrom) %>% pull(SNP)

X_raw               = geno_input[, current_SNPs] %>% as.matrix() 

temp_results        = foreach(i = current_SNPs, .combine = 'rbind', .packages = "SNPRelate") %dopar% {
    
pos_i               = map %>% dplyr::filter(SNP == i) %>% pull(POS)
    
snps                = map %>% dplyr::filter(CHROM == chrom, POS >= pos_i - max_dist, POS <= pos_i + max_dist, SNP != i)

ld_list             = lapply(names(P), function(adj) {

X_adj               = P[[adj]] %*% X_raw

X_adj               = scale(X_adj, center = FALSE, scale = sqrt(apply(X_adj, 2, crossprod)))

temp                = lapply(snps$SNP, function(j) {

pos_j               = snps %>% dplyr::filter(SNP == j) %>% pull(POS)

r                   = sum(X_adj[, i] * X_adj[, j])

data.frame("adjustment" = adj, "CHROM" = chrom, "snp" = i, "D" = abs(pos_i - pos_j), "r" = min(max(r, -1), 1), "r2" = r^2)})

bind_rows(temp)})

bind_rows(ld_list) }

temp_split          = split(temp_results, temp_results$adjustment)
  
for (adj in names(temp_split)) {

output[[adj]]       = bind_rows(output[[adj]], temp_split[[adj]]) } }

saveRDS(output, paste0(path_output, "LDdecay.rds"))

#############################################################
