################################################################################
# Map solid tumors in the HS/Npt mice using different models and covariates.
# Daniel Gatti
# dmgatti@coa.edu
# 2021-05-14
################################################################################
options(stringsAsFactors = FALSE)

library(tidyverse)
library(ROCR)       # Precision/Recall plot.
library(glmnet)     # Elastic net & support fxns.
library(readxl)
library(doMC)       # Parallel execution.
library(qtl2)

base_dir   = '/media/dmgatti/hdb/projects/ColoState'
data_dir   = file.path(base_dir,   'data')
pheno_file = file.path(data_dir,   'GRSD.pheno.powercalc.xlsx')
qtl2_file  = file.path(data_dir,   'HZEproject_qtl2.Rdata')
result_dir = file.path(base_dir,   'results')
gwas_dir   = file.path(result_dir, 'gwas')
binary_dir = file.path(gwas_dir,   'binary_model')
linear_dir = file.path(gwas_dir,   'linear_model')
prs_dir    = file.path(result_dir, 'polygenic_risk_score')
snp_file = '/media/dmgatti/hda/data/MUGA/hsnpt_variants.sqlite'

# SNP lookup function for qtl2.
snp_func = create_variant_query_func(dbfile = snp_file)

# Set up for parallel execution.
registerDoMC(cores = 2)

################################################################################
# Data input and setup.

# Read in qtl2 data.
load(qtl2_file)

# Read in phenotypes.
pheno = read_xlsx(pheno_file, sheet = 'GRSD.pheno')
pheno = pheno[,c('Animal number', 'sex', 'family', 'group', 'Solid Tumor')]
pheno = as.data.frame(pheno)
colnames(pheno) = c('mouse', 'sex', 'family', 'group', 'solid_tumor')
pheno$mouse = paste0('X', pheno$mouse)
rownames(pheno) = pheno$mouse
pheno$sex = factor(pheno$sex, levels = c('F', 'M'))
pheno$group = sub('-(Fe|Si)', '', pheno$group)
pheno$group = sub('^sham ', 'un', pheno$group)
pheno$group = factor(pheno$group, levels = c('unirradiated', 'gamma', 'HZE'))

# Synch samples.
samples = intersect(pheno$mouse, rownames(probs[[1]]))
pheno = pheno[samples,]
for(i in seq_along(probs)) {
  
  probs[[i]] = probs[[i]][samples,,]
  K[[i]] = K[[i]][samples, samples]
  
} # for(i)

covar = model.matrix(~ sex + group, data = pheno)[,-1,drop = FALSE]

################################################################################
# GWAS using qtl2.

# GWAS using linear model.
lod_linear = scan1snps(genoprobs = probs, map = map, pheno = pheno[,'solid_tumor', drop = FALSE], 
                       kinship = K, addcovar = covar, model = 'normal', query_func = snp_func, 
                       keep_all_snps = FALSE, cores = 4)
saveRDS(lod_linear, file = file.path(linear_dir, 'solid_tumor_gwas.Rdata'))

png(file.path(linear_dir, 'solid_tumor_gwas.png'), width = 2000, height = 1000, res = 128)
plot_snpasso(lod_linear$lod, lod_linear$snpinfo, main = 'HS/Npt, Solid Tumor, Linear Model, sex + group covar')
dev.off()

# GWAS using binary model.
lod_binary = scan1snps(genoprobs = probs, map = map, pheno = pheno[,'solid_tumor', drop = FALSE], 
                       addcovar = covar, model = 'binary', query_func = snp_func, 
                       keep_all_snps = FALSE,  cores = 4)
saveRDS(lod_binary, file = file.path(binary_dir, 'solid_tumor_gwas.Rdata'))
lod_binary = readRDS(file = file.path(binary_dir, 'solid_tumor_gwas.Rdata'))

png(file.path(binary_dir, 'solid_tumor_gwas.png'), width = 2000, height = 1000, res = 128)
plot_snpasso(lod_binary$lod, lod_binary$snpinfo, main = 'HS/Npt, Solid Tumor, Binary Model, sex + group covar')
dev.off()

# NOTE: The binary model produces higher LOD scores. This may be the correct model to use.

rm(lod_linear)

################################################################################
# ANOVA for sex and radiation group.
mod = glm(solid_tumor ~ sex + group + sex:group, data = pheno, family = binomial(link = 'logit'))
summary(mod)

# Males have more tumors. The irradiated groups have more tumors.
# There may be some interaction between sex and irradiation as well.

# Fitting simplified model with only 'unirradiated' & 'irradiated' groups.
tmp = pheno[,c('sex', 'group', 'solid_tumor')]
tmp$group = sub('gamma|HZE', 'irradiated', tmp$group)
tmp$group = factor(tmp$group, levels = c('unirradiated', 'irradiated'))
  
mod2 = glm(solid_tumor ~ sex + group + sex:group, data = tmp, family = binomial(link = 'logit'))
summary(mod2)

# It may be worth thinking abotut a model like this. Holding off for now.

rm(tmp)

################################################################################
# Use glmnet to build models.

# family = binomial(link = "logit")
# alpha = 1 --> LASSO, alpha = 0 --> ridge.

# Reformat the qtl2 SNPs into a matrix for glmnet.
snps = NULL
for(chr in names(probs)) {
  
  # Get SNPs from HS/Npt database.
  s = snp_func(chr = chr, start = 0, end = 1e9)
  s = s[,c('snp_id', 'chr', 'pos', 'sdp')]
  si = index_snps(map, s)
  s = qtl2::genoprob_to_snpprob(probs, si)[[1]][,1,]
  
  snps = cbind(snps, s)

} # for(chr)

# Remove probs to save memory.
rm(probs)
gc()

# Remove SNPs with MAF < 5%.
maf = colMeans(snps)
maf[maf > 0.5] = 1.0 - maf[maf > 0.5]
dim(snps)
snps = snps[,maf > 0.05]
dim(snps)

# Add sex and treatment to the front of the snps.
covar = model.matrix(~sex + group, data = pheno)[,-1]
stopifnot(all(rownames(covar) == rownames(snps)))
snps = cbind(covar, snps)

saveRDS(snps, file  = file.path(prs_dir, 'snps.rds'))
snps = readRDS(file = file.path(prs_dir, 'snps.rds'))

# We have 1784 samples in three roughly 600 sample groups.
table(pheno$group)

# About 47% have solid tumors.
mean(pheno$solid_tumor)

# Use elastic net with logistic regression model.

mod = glmnet(snps, pheno[,'solid_tumor'], alpha = 0.5, family = 'binomial')
plot(mod)

cvmod = cv.glmnet(snps, pheno[,'solid_tumor'], alpha = 0.5, family = 'binomial',
                  type.measure = 'auc', parallel = TRUE)
plot(cvmod)

cvmod$lambda.1se

mod = glmnet(snps, pheno[,'solid_tumor'], alpha = 0.5, family = 'binomial', lambda = cvmod$lambda.1se)

pred = predict(mod, newx = snps, s = cvmod$lambda.1se, type = 'class')

tmp = coef(mod, lambda = cvmod$lambda.1se)
tmp = tmp[tmp[,1] != 0,]

lrmod = glm(pheno$solid_tumor ~ snps[,names(tmp)[-1]], 
            family = binomial(link = 'logit'))
lrpred = predict(lrmod, type = 'response')
plot(lrpred ~ pheno$solid_tumor)
lrpred = lrpred >= 0.5

table(lrpred, pheno$solid_tumor)

#              Actual Values
#                T   |   F
#            --+-----+------
# Predicted  T | TP  |  FP
#            --+-----+------
#   Values   F|  FN  |  TN
#            --+-----+------
#
# Sensitivity = TP / (TP + FN)   aka Recall
# Specificity = TN / (TN + FP) Prop. of TN out of actual F values.
# Precision = TP / (TP + FP)   Prop. of TP out of predicted T values.
# Recall = TP / (TP + FN)      Prop. of TP out of actual T values.

cvmod = cv.glmnet(snps, pheno[,'solid_tumor'], alpha = 0.5, family = 'binomial',
                  parallel = TRUE)
preds = predict(cvmod, newx = snps, type = 'response')
perf = performance(prediction(preds, pheno$solid_tumor), measure = 'tpr', x.measure = 'fpr')
plot(perf)
abline(0, 1, lty = 2)

# Testing glmnet ROC code. It's slow.
mod = glmnet(snps, pheno$solid_tumor, alpha = 0.5, family = 'binomial', lambda = cvmod$lambda.1se)
# This produces AUC, MSE, etc.
summry = assess.glmnet(mod, newy = pheno$solid_tumor, newx = snps, family = 'binomial')

confusion.glmnet(mod, newx = snps, newy = pheno$solid_tumor, family = 'binomial')

roc = roc.glmnet(mod, newx = snps, newy = pheno$solid_tumor)
plot(roc, type = 'l')
abline(0, 1, lty = 2)


###### Split data into training and test sets to estimate AUC.
set.seed(2021051613)
nsim = 100
prop_train = 0.8
# Sample factor is the grouping for each sex & radiation group.
sample_fctr = factor(pheno$sex:pheno$group)
samples = rownames(pheno)
samples = split(samples, sample_fctr)
sample_n = sapply(samples, length)
# Get genome positions for GWAS plots.
chrlen = sapply(split(lod_binary$snpinfo$pos, lod_binary$snpinfo$chr), max)
# Warning is OK here. It's chr X not converting to a number.
chrlen = chrlen[order(as.numeric(names(chrlen)))]
# qtl2 add a gap between each chr that is 1% of genome length.
# See https://github.com/kbroman/qtl2/blob/master/R/plot_scan1.R
# But this doesn't work... Adjusting by -5.
gap = sum(chrlen) / 100 - 5
chrend = cumsum(c(0, chrlen[-length(chrlen)] + gap))
chrend[1] = gap / 2
names(chrend) = names(chrlen)

auc_results = data.frame(sim = 1:nsim, 
                         auc = rep(0, nsim), 
                         num_snps = rep(0, nsim),
                         tp = rep(0, nsim),
                         fp = rep(0, nsim),
                         tn = rep(0, nsim),
                         fn = rep(0, nsim))

for(i in 1:nsim) {

  print(paste('Sim', i))

  # Sample evenly from each sex and radiation group.
  training_samples = unlist(mapply(x = samples, y = sample_n, 
                                   FUN = function(x, y) { sample(x, round(prop_train * y)) }))
  test_samples     = setdiff(unlist(samples), training_samples)
  
  save(training_samples, test_samples, 
       file = file.path(prs_dir, paste0('sim_', i, '_train_test_samples.Rdata')))
  
  # Training.
  cvmod = cv.glmnet(x = snps[training_samples,], 
                    y = pheno[training_samples, 'solid_tumor'], 
                    alpha = 0.5, family = 'binomial', parallel = TRUE)
  
  ### Test set
  # Fit model using all training set data with lambda form CV.
  mod = glmnet(x = snps[training_samples,], 
               y = pheno[training_samples, 'solid_tumor'],
               lambda = cvmod$lambda.1se, alpha = 0.5, 
               family = 'binomial', parallel = TRUE)

  # Save this model.
  saveRDS(mod, file = file.path(prs_dir, paste0('sim_', i, '_lasso_model.rds')))
  
  # Get AUC.
  auc = assess.glmnet(mod, newy = pheno[test_samples, 'solid_tumor'], 
                      newx = snps[test_samples,], family = 'binomial')
  auc_results$auc[i] = auc$auc
  
  # Get confusion matrix.
  conf_tbl = confusion.glmnet(mod, newy = pheno[test_samples, 'solid_tumor'], 
                              newx = snps[test_samples,], family = 'binomial')
  auc_results$tp[i] = conf_tbl[2, 2]
  auc_results$fp[i] = conf_tbl[2, 1]
  auc_results$tn[i] = conf_tbl[1, 1]
  auc_results$fn[i] = conf_tbl[1, 2]

  # Get test set coefs.
  lasso_coef = coef(mod)
  lasso_coef = lasso_coef[which(lasso_coef != 0),]
  write.csv(lasso_coef, file = file.path(prs_dir, paste0('sim_', i, '_lasso_coef.csv')),
            quote = F)
  num_snps = length(grep('^([0-9]|X)', names(lasso_coef)))
  auc_results$num_snps[i] = num_snps
  
  # Genome plot of coef locations overlaid on binary model GWAS.
  #####################################################
  # NOTE: The SNPS may not be in the correct location!
  #####################################################
  # png(file.path(prs_dir, 'lasso_coef_genome_plots', paste0('solid_tumor_gwas_sim_', i, '.png')), 
  #     width = 2000, height = 1000, res = 128)
  # plot_snpasso(lod_binary$lod, lod_binary$snpinfo, main = 'HS/Npt, Solid Tumor, Binary Model, sex + group covar')
  # snploc = subset(lasso_coef, grepl('^([0-9]|X)', names(lasso_coef)))
  # snploc = sub('_[ACGT]/[ACGT]$', '', names(snploc))
  # snploc = strsplit(snploc, split = ':')
  # snploc = data.frame(chr = sapply(snploc, '[', 1),
  #                  pos = as.numeric(sapply(snploc, '[', 2)) * 1e-6)
  # m = match(snploc$chr, names(chrend))
  # snploc$pos = snploc$pos + chrend[m]
  # usr = par('usr')
  # points(snploc$pos, rep(usr[4] - 0.1 * diff(usr[3:4]), nrow(snploc)), pch = 3, col = 2)
  # rm(snploc)
  # dev.off()

  # Sensitivity / Specificity plot.
  png(file = file.path(prs_dir, 'sens_spec_plots', paste0('prs_sim_', i,'_auc.png')),
      width = 1600, height = 800, res = 128)
  layout(matrix(1:2, 1, 2))
  roc = roc.glmnet(mod, newy = pheno[test_samples, 'solid_tumor'], newx = snps[test_samples,])
  plot(roc, main = paste('Sim', i, ', Num SNPs =', num_snps,
       ', AUC =', round(auc_results$auc[i], digits = 3)),
       xlim = c(0, 1), ylim = c(0, 1), type = 'l')
  abline(0, 1, lty = 2)
  
  preds = predict.glmnet(mod, newx = snps[test_samples,], type = 'class')
  perf = performance(prediction(preds, pheno[test_samples, 'solid_tumor']), 
                     measure = 'prec', x.measure = 'rec')
  plot(perf, main = paste('Sim', i, 'Precision/Recall'),
       xlim = c(0, 1), ylim = c(0, 1))
  abline(h = mean(pheno$solid_tumor), lty = 2)
  dev.off()
  
  rm(cvmod, mod, mod_all, preds, perf)
  gc()

} # for(i)

write.csv(auc_results, file = file.path(prs_dir, 'prs_sims_test_auc.csv'),
          row.names = F, quote = F)


################################################################################
# Use the SNPs in each model to estimate the proportion of phenotypic
# variance explained by the SNPs selected by the model.
# TBD: How do we do this? Variance of binary trait? np(1-p)
# Not sure if this even makes sense with binary trait.
coef_files = dir(prs_dir, pattern = '_lasso_coef.csv$')

for(f in coef_files) {
  
  coefs = read.csv(file.path(prs_dir, f))
  
} # for(f)


# Look at the correlation between SNPs in each model.




################################################################################
# Simulate random data and see how prediction is.
set.seed(2021051613)
nperm = 100
prop_train = 0.8
# Sample factor is the grouping for each sex & radiation group.
sample_fctr = factor(pheno$sex:pheno$group)
samples = rownames(pheno)
samples = split(samples, sample_fctr)
sample_n = sapply(samples, length)

auc_results = data.frame(perm = 1:nperm, auc = rep(0, nperm), num_snps = rep(0, nperm))

pheno_perm = pheno[,'solid_tumor', drop = FALSE]

for(i in 1:nperm) {
  
  print(paste('Perm', i))
  
  # Sample evenly from each sex and radiation group.
  training_samples = unlist(mapply(x = samples, y = sample_n, 
                                   FUN = function(x, y) { sample(x, round(prop_train * y)) }))
  test_samples     = setdiff(unlist(samples), training_samples)
  
  # Permutation of samples.
  pheno_perm$solid_tumor = sample(pheno$solid_tumor)
  
  # Training.
  cvmod = cv.glmnet(x = snps[training_samples,], 
                    y = pheno_perm[training_samples, 'solid_tumor'], 
                    alpha = 0.5, family = 'binomial', parallel = TRUE)
  
  ### Test set
  # fit test set model.
  mod = glmnet(x = snps[test_samples,], 
               y = pheno_perm[test_samples, 'solid_tumor'],
               lambda = cvmod$lambda.1se, alpha = 0.5, 
               family = 'binomial', parallel = TRUE)
  # Get test set coefs.
  lasso_coef = coef(mod)
  lasso_coef = lasso_coef[which(lasso_coef != 0),]
  num_snps = length(grep('^([0-9]|X)', names(lasso_coef)))
  auc_results$num_snps[i] = num_snps
  
  # Get predictions for test set.
  preds = predict(mod, newx = snps[test_samples,], type = 'response')
  
  rocr_pred = prediction(preds, pheno[test_samples, 'solid_tumor'])
  
  # Get model performance from ROCR.
  perf = performance(rocr_pred, measure = 'tpr', x.measure = 'fpr')
  perf_auc = performance(rocr_pred, measure = 'auc', fpr.stop = 1.0)
  auc_results$auc[i] = perf_auc@y.values[[1]]

} # for(i)

write.csv(auc_results, file = file.path(prs_dir, 'prs_random_perms_auc.csv'),
          row.names = F, quote = F)



