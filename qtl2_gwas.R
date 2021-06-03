################################################################################
# Test whether HS/Npt SNP database works with qtl2.
# Daniel Gatti
# Mar. 1, 2021
# dmgatti@coa.edu
################################################################################

#library(qtl2convert)
library(qtl2)

base_dir  = '/media/dmgatti/hdb/projects/ColoState'
# Original file from Elijah.
#data_file = file.path(base_dir, 'data', 'HZEproject.Rdata')
# New file in qtl2 format.
data_file = file.path(base_dir, 'data', 'HZEproject_qtl2.Rdata')
linkage_dir = file.path(base_dir, 'results', 'linkage')
gwas_dir   = file.path(base_dir, 'results', 'gwas')
binary_dir = file.path(gwas_dir, 'binary_model')
linear_dir = file.path(gwas_dir, 'linear_model')
snp_file   = file.path('/media/dmgatti/hda/data/MUGA', 'hsnpt_variants.sqlite')

# Create SNP DB function.
snp_func = create_variant_query_func(snp_file)

# Read in data.
load(data_file)

pheno_bin = as.matrix(pheno[,8:ncol(pheno)])

# # Convert from DOQTL format to qtl2 format.
# colnames(markers)[1:3] = c('marker', 'chr', 'pos')
# probs = qtl2convert::probs_doqtl_to_qtl2(probs, markers)
# map   = qtl2convert::map_df_to_list(markers, pos_column = 'pos')
# 
# # Subset phenotypes.
# pheno_bin = as.matrix(pheno[,8:ncol(pheno)])
# rownames(pheno_bin) = pheno$rownames
# 
# addcovar = matrix(addcovar, ncol = 1)
# rownames(addcovar) = pheno$rownames
# 
# # Check for sample alignment.
# samples = intersect(rownames(pheno_bin), rownames(probs[[1]]))
# 
# pheno     = pheno[samples,, drop = FALSE]
# pheno_bin = pheno_bin[samples,, drop = FALSE]
# addcovar  = addcovar[samples,, drop = FALSE]
# for(i in 1:length(probs)) {
# 
#   probs[[i]] = probs[[i]][samples,,]
# 
# } # for(i)
# 
# K = calc_kinship(probs, type = 'loco', cores = 4)
# 
# stopifnot(all(rownames(addcovar)   == rownames(pheno_bin)))
# stopifnot(all(rownames(K[[1]])     == rownames(pheno_bin)))
# stopifnot(all(rownames(probs[[1]]) == rownames(pheno_bin)))
# 
# # Save new Rdata file.
# save(pheno, addcovar, probs, K, map, file = file.path(base_dir, 'data', 'HZEproject_qtl2.Rdata'))

######
# Heritability estimates for each phenotype.
Kall = calc_kinship(probs, type = 'overall', cores = 4)
herit = qtl2::est_herit(pheno = pheno_bin, kinship = Kall, addcovar = addcovar)
write.csv(herit, file.path(binary_dir, 'herit_sex_only.csv'))

######
# Map phenotypes using GWAS with sex as only covariate.
results_dir = file.path(binary_dir, 'sex_covar_only')

for(i in 1:ncol(pheno_bin)) {

  pheno_name = colnames(pheno_bin)[i]
  
  print(pheno_name)
  
  gwas = scan1snps(genoprobs = probs, map = map, pheno = pheno_bin[,i,drop = FALSE],
                   addcovar = addcovar, query_func = snp_func, model = 'binary', 
                   keep_all_snps = FALSE, cores = 2)
  saveRDS(gwas, file = file.path(results_dir, paste0(pheno_name, '_gwas.rds')))
  
  png(file.path(results_dir, paste0(pheno_name, '_gwas.png')), width = 2000, height = 1000, res = 128)
  plot_snpasso(gwas$lod, gwas$snpinfo, altcol = 'blue', gap = 0, main = pheno_name)
  dev.off()
  
  rm(gwas)

} # for(i)

######
# Permutations using GWAS with sex as only covariate.
# perms = matrix(0, nrow = nperm, ncol = ncol(pheno_bin),
#                dimnames = list(paste0('p', 1:nperm), colnames(pheno_bin)))
# for(p in 1:nperm) {
#   
#     print(p)
#     print(Sys.time())
# 
#     # Resample phenotypes.
#     pheno_perm = pheno_bin
#     rownames(pheno_perm) = sample(rownames(pheno_perm))
#     rownames(addcovar)  = rownames(pheno_perm)w
#     
#     gwas = scan1snps(genoprobs = probs, map = map, pheno = pheno_perm,
#                      addcovar = addcovar, query_func = snp_func, model = 'binary', 
#                      keep_all_snps = FALSE, cores = 2)
#     
#     perms[p,] = apply(gwas$lod, 2, max, na.rm = T)
#     
#     saveRDS(perms, file = file.path(results_dir, 'perms.rds'))
# 
#     rm(gwas)
#     
# } # for(p)


######
# Map phenotypes using GWAS with sex and radiation as additive covariates.
results_dir = file.path(binary_dir, 'radiation_additive')
pheno$group = factor(pheno$group, levels = c('Unirradiated', 'Gamma', 'HZE'))
addcovar = model.matrix(~sex + group, data = pheno)[,-1]

for(i in 1:ncol(pheno_bin)) {
  
  pheno_name = colnames(pheno_bin)[i]
  
  print(pheno_name)
  
  gwas = scan1snps(genoprobs = probs, map = map, pheno = pheno_bin[,i,drop = FALSE],
                   addcovar = addcovar, query_func = snp_func, model = 'binary', 
                   keep_all_snps = FALSE, cores = 2)
  saveRDS(gwas, file = file.path(results_dir, paste0(pheno_name, '_gwas.rds')))
  
  png(file.path(results_dir, paste0(pheno_name, '_gwas.png')), width = 2000, height = 1000, res = 128)
  plot_snpasso(gwas$lod, gwas$snpinfo, altcol = 'blue', gap = 0, main = pheno_name)
  dev.off()
  
  rm(gwas)
  
} # for(i)


######
# Map phenotypes using GWAS with sex as additive covariate and radiation as interactive covariates.
results_dir = file.path(binary_dir, 'radiation_interactive')
pheno$group = factor(pheno$group, levels = c('Unirradiated', 'Gamma', 'HZE'))
addcovar = model.matrix(~sex + group, data = pheno)[,-1]

for(i in 1:ncol(pheno_bin)) {
  
  pheno_name = colnames(pheno_bin)[i]
  
  print(pheno_name)
  
  
  gwas = scan1snps(genoprobs = probs, map = map, pheno = pheno_bin[,i,drop = FALSE],
                   addcovar = addcovar, intcovar = addcovar[,-1], query_func = snp_func, model = 'binary', 
                   keep_all_snps = FALSE, cores = 2)
  saveRDS(gwas, file = file.path(results_dir, paste0(pheno_name, '_gwas.rds')))
  
  png(file.path(results_dir, paste0(pheno_name, '_gwas.png')), width = 2000, height = 1000, res = 128)
  plot_snpasso(gwas$lod, gwas$snpinfo, altcol = 'blue', gap = 0, main = pheno_name)
  dev.off()
  
  rm(gwas)
  
} # for(i)

#####
# Once all of the mapping is done, plot the additive LOD, the interactive LOD and the LOD difference
# for each phenotype.
add_dir = file.path(binary_dir, 'radiation_additive')
int_dir = file.path(binary_dir, 'radiation_interactive')

add_files = dir(path = add_dir, pattern = 'rds$', full.names = FALSE)
int_files = file.path(int_dir, add_files)
add_files = file.path(add_dir, add_files)

for(i in 1:length(add_files)) {
  
  a_gwas = readRDS(add_files[i])
  i_gwas = readRDS(int_files[i])
  d_gwas = i_gwas
  d_gwas$lod = d_gwas$lod - a_gwas$lod
  
  pheno_name = sub(add_dir,       '', add_files[i])
  pheno_name = gsub('^/|_gwas\\.rds$', '', pheno_name)
  
  max_lod = max(i_gwas$lod)

  png(file.path(binary_dir, paste0(pheno_name, '_gwas_all.png')), width = 1600, height = 1600, res = 200)  
  par(mfrow = c(3, 1))
  plot_snpasso(a_gwas$lod, a_gwas$snpinfo, ylim = c(0, max_lod), main = paste(pheno_name, 'Additive'))
  plot_snpasso(i_gwas$lod, i_gwas$snpinfo, ylim = c(0, max_lod), main = paste(pheno_name, 'Interactive'))
  plot_snpasso(d_gwas$lod, d_gwas$snpinfo, ylim = c(0, max_lod), main = paste(pheno_name, 'Interactive - Additive'))
  dev.off()
  
} # for(i)

###### Harvest qtl peaks.
# sex-sdditive perms.
perms = readRDS(file.path(binary_dir, 'perms.rds'))

# perms = readRDS(file.path(binary_dir, 'perms1.rds'))
# for(i in 2:4) {
#   perms = rbind(readRDS(file.path(binary_dir, paste0('perms', i, '.rds'))))
# }
# perms = subset(perms, perms[,1] > 0)

thr = quantile(perms, probs = 0.95)

gwas_files = dir(file.path(binary_dir, 'sex_covar_only'), pattern = '_gwas\\.rds$')

result = NULL

for(f in gwas_files) {
  
  pheno_name = sub('_gwas\\.rds$', '', f)
  
  print(pheno_name)
  
  gwas = readRDS(file.path(binary_dir, 'sex_covar_only', f))

  local_result = NULL
  
  # Get peaks from each chromosome.
  unique_chr = unique(gwas$snpinfo$chr)
  for(chr in unique_chr) {
  
    ts = qtl2::top_snps(gwas$lod, gwas$snpinfo, chr = chr)
    local_result = rbind(local_result, subset(ts, lod >= thr))
    
  } # for(chr)
  
  # Compile results
  if(nrow(local_result) > 0) {
    local_result = cbind(phenotype = pheno_name, local_result)
    result = rbind(result, local_result)
  } # if(nrow(local_result) > 0)
  
} # for(f)

write.table(result, file = file.path(binary_dir, 'sex_covar_only', 'gwas_top_snps.tsv'),
          sep = '/t', row.names = F)




