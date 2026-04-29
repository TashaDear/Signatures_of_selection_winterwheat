#############################################################
library(dplyr)
library(foreach)
library(doParallel)
library(tibble)
library(tidyr)
#############################################################
#A) Paths 
#############################################################

submission         = "submission1"

path               = "/home/tasha/breedfuture/"

path_pheno         = paste0(path, "prepared_data/v2/")

path_geno          = paste0(path, "results/v2/fitness_estimates/")

path_map           = paste0(path, "raw_data/WW/geno/plink/all_data/")

path_kinships      = paste0(path, "results/v2/kinships/random_effects/")

path_output        = paste0(path, "results/v2/model_output/LDdecay/")

#############################################################
#B) Input
#############################################################

all_samples        = readRDS(paste0(path_pheno, "phenotypes.rds")) %>% pull(id) %>% unique()
                     
geno               = readRDS(paste0(path_geno, "geno.rds")) 
geno_subset        = geno[all_samples, ]

map                = read.table(paste0(path_map, "ww_updated_final.map"), header=TRUE,col.names=c("CHROM","SNP","POS")) %>% as.data.frame() %>% 
                     dplyr::filter(SNP %in% colnames(geno))
                     
G                  = readRDS(paste0(path_kinships, "GRM_baseline_v2.rds"))[all_samples, all_samples]
                     
chromosomes        = unique(map$CHROM) %>% na.omit() %>% as.vector()
                                            
#############################################################        
#C) Functions
#############################################################

invSqrt            = function(mat) {
stopifnot(isSymmetric(mat))
  
ED                 = eigen(mat) 
L                  = ED$vectors %*% diag(1/sqrt(ED$values)) %*% t(ED$vectors)
  
return(L) }              

#############################################################        
#D) Prepare matrices
#############################################################

Q                  = prcomp(geno_subset, center = T, scale = F)$x[,1:3]
Q                  = cbind(Intercept = 1, Q)

G_upd              = as.matrix(Matrix::nearPD(G)$mat)
genotypes          = rownames(Q)
n                  = length(genotypes)
G_upd              = G_upd[genotypes, genotypes]
diag(G_upd)        = diag(G_upd) + 1e-5

# Projection matrices
P                  = list()
P[["Int"]]         = diag(n) - matrix(1, nrow=n, ncol=n)/n
P[["Q"]]           = diag(n) - Q %*% solve(crossprod(Q)) %*% t(Q)

Ginv               = solve(G_upd)
H                  = Q %*% solve(t(Q) %*% Ginv %*% Q) %*% t(Q) %*% Ginv
P[["Q+K"]]         = invSqrt(G_upd) %*% (diag(n) - H)

max_dist           = 10000000

############################################################
#C) LDdecay
#############################################################

output             = list()

for (chrom in chromosomes) {
print(chrom)

current_SNPs        = map %>% dplyr::filter(CHROM == chrom) %>% pull(SNP)

X_raw               = geno[all_samples, current_SNPs] %>% as.matrix() 

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

saveRDS(output, paste0(path_output, "LDdecay_",submission,".rds"))

#############################################################
