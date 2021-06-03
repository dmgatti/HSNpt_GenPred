#!/bin/bash

BASE_DIR=/media/dmgatti/hdb/projects/ColoState

R CMD BATCH --no-save --no-restore '--args 1 250' gwas_perms_worker.R perms1.Rout &
R CMD BATCH --no-save --no-restore '--args 2 250' gwas_perms_worker.R perms2.Rout &
R CMD BATCH --no-save --no-restore '--args 3 250' gwas_perms_worker.R perms3.Rout & 
R CMD BATCH --no-save --no-restore '--args 4 250' gwas_perms_worker.R perms4.Rout &

