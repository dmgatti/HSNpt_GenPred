################################################################################
# Make a new qtl2-style Rdata file.
# Daniel Gatti.
# Jan. 25, 2021
################################################################################

# HS founders: A/J, AKR/J, BALB/cJ, C3H/HeJ, C57BL/6J, DBA/2J, CBA/2J and LP/J.

options(stringsAsFactors = FALSE)

# LIBRARIES
library(readxl)
library(qtl2convert)
library(qtl2)
library(tidyverse)

# CONSTANTS
base_dir = '/media/dmgatti/hdb/projects/ColoState'
data_dir = file.path(base_dir, 'data')
pheno_file = file.path(data_dir, 'GRSD.pheno.powercalc.xlsx')
rdata_file = file.path(data_dir, 'HZEproject.Rdata')
new_rdata_file = file.path(data_dir, 'HZEproject_qtl2.Rdata')

###########################################################
# Run this first to make a new qtl2-style probs, map and K.
# Loads in pheno, addcovar, K, markers and probs in DOQTL format.
load(rdata_file)

# Read in a new phenotype file.
pheno = read_xlsx(path = pheno_file, sheet = 'GRSD.pheno')
pheno = as.data.frame(pheno)
rownames(pheno) = paste0('X', pheno$`Animal number`)

# Intersect and synch samples.
samples = sort(intersect(rownames(pheno), rownames(probs)))

# Samples NOT in both pheno & probs.
setdiff(rownames(pheno), rownames(probs))

pheno = pheno[samples,]
probs = probs[samples,,]

# Convert markers to qtl2 format.
colnames(markers)[colnames(markers) == 'SNP_ID']    = 'marker'
colnames(markers)[colnames(markers) == 'Chr']       = 'chr'
colnames(markers)[colnames(markers) == 'Mb_NCBI38'] = 'pos'
markers$pos = markers$pos * 1e-6
map = qtl2convert::map_df_to_list(markers, pos_column = 'pos')

# Convert probs to qtl2 format.
# Sex is coded as F = 0 and M = 1 in addcovar.
new_probs = qtl2convert::probs_doqtl_to_qtl2(probs = probs, map = markers, is_female = addcovar[,1] == 0)

# Calculate new kinship matrices.
K = calc_kinship(probs = new_probs, type = 'loco', cores = 4)

# Clean up.
probs = new_probs
rm(new_probs)
gc()

# Save a new Rdata file.
save(pheno, addcovar, map, K, probs, file = new_rdata_file)
