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
gwas_dir  = file.path(base_dir, 'results', 'gwas')
snp_file  = file.path(base_dir, 'data', 'hsnpt_variants.sqlite')
# Number of permutations.
nperm     = 1000

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
#   K[[i]]     = K[[i]][samples, samples]
#   
# } # for(i)
# 
# stopifnot(all(rownames(addcovar)   == rownames(pheno_bin)))
# stopifnot(all(rownames(K[[1]])     == rownames(pheno_bin)))
# stopifnot(all(rownames(probs[[1]]) == rownames(pheno_bin)))
# 
# # Save new Rdata file.
# save(pheno, addcovar, probs, K, map, file = file.path(base_dir, 'data', 'HZEproject_qtl2.Rdata'))

######
# Map phenotypes using GWAS with sex as only covariate.
results_dir = file.path(gwas_dir, 'sex_covar_only')

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

} # for(i)

