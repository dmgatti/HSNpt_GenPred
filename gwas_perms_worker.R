################################################################################
# GWAS permutations for Colorado State radiation phenotypes.
# Mar 5, 2021
# Daniel Gatti
# dmgatti@coa.edu
################################################################################

args = commandArgs(trailingOnly = TRUE)

# The permuation file to write to. Integer.
perm_num = as.numeric(args[1])
# The number of permutations to run. Integer.
nperm    = as.numeric(args[2])

print(paste(perm_num, ',', nperm))

library(qtl2)

base_dir  = '/media/dmgatti/hdb/projects/ColoState'
# New file in qtl2 format.
data_file   = file.path(base_dir, 'data', 'HZEproject_qtl2.Rdata')
results_dir = file.path(base_dir, 'results', 'gwas')
snp_file  = '/media/dmgatti/hdb/projects/ColoState/data/hsnpt_variants_old.sqlite'

# Create SNP DB function.
snp_func = create_variant_query_func(snp_file)

# Read in data.
load(data_file)

pheno_bin = as.matrix(pheno[,8:ncol(pheno)])

perms = matrix(0, nrow = nperm, ncol = ncol(pheno_bin),
               dimnames = list(paste0('p', 1:nperm), colnames(pheno_bin)))

for(p in 1:nperm) {

  print(p)
  print(Sys.time())

  # Resample phenotypes.
  pheno_perm = pheno_bin
  rownames(pheno_perm) = sample(rownames(pheno_perm))
  rownames(addcovar)  = rownames(pheno_perm)

  gwas = scan1snps(genoprobs = probs, map = map, pheno = pheno_perm,
                   addcovar = addcovar, query_func = snp_func, model = 'binary',
                   keep_all_snps = FALSE, cores = 1)

  perms[p,] = apply(gwas$lod, 2, max, na.rm = T)

  saveRDS(perms, file = file.path(results_dir, paste0('perms', perm_num, '.rds')))

  rm(gwas)

} # for(p)


