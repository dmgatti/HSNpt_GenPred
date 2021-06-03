# Test how many cores to use for perms.
# Vary number and check time and memory.

####################################################################################
# Mar 5, 2021
# I found that the cores argument has no affect on the timing or number of 
# cores uses. I looked at the system monitor and never saw more than one
# core being used. This took 40 minutes per mapping run, regardless of the 
# number of cores.
####################################################################################


library(qtl2)

base_dir  = '/media/dmgatti/hdb/projects/ColoState'
# Original file from Elijah.
#data_file = file.path(base_dir, 'data', 'HZEproject.Rdata')
# New file in qtl2 format.
data_file = file.path(base_dir, 'data', 'HZEproject_qtl2.Rdata')
linkage_dir = file.path(base_dir, 'results', 'linkage')
gwas_dir  = file.path(base_dir, 'results', 'gwas')
snp_file  = '/media/dmgatti/hda/data/MUGA/hsnpt_variants.sqlite'
# Number of permutations.
nperm     = 1000

# Create SNP DB function.
snp_func = create_variant_query_func(snp_file)

# Read in data.
load(data_file)

pheno_bin = as.matrix(pheno[,8:ncol(pheno)])

addcovar = model.matrix(~sex + group, data = pheno)[,-1]

num_cores = data.frame(num_cores = 1:10,
                       time      = rep(0, 10),
                       mem       = rep(0, 10))

for(nc in num_cores$num_cores) {
  
  print(nc)
  start_time = proc.time()[3]
  
  gwas = scan1snps(genoprobs = probs, map = map, pheno = pheno_bin,
                   addcovar = addcovar, intcovar = addcovar[,-1], query_func = snp_func, 
                   model = 'binary', keep_all_snps = FALSE, cores = nc,
                   chr = '1', start = 0, end = 200)
  
  num_cores$time[nc] = proc.time()[3] - start_time
  num_cores$mem[nc]  = sum(gc()[,2])
  
  print(num_cores$time[nc])
  
  rm(gwas)
  gc()
  
} # for(nc)
