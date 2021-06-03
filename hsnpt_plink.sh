#!/bin/bash

# Uses PLINK PLINK v1.90b6.21 64-bit (19 Oct 2020)
BASE_DIR=/media/dmgatti/hdb/projects/ColoState
DATA_DIR=${BASE_DIR}/data/plink
INPUT_FILE=${DATA_DIR}/hsnpt_hze
PHENO_FILE=${DATA_DIR}/hsnpt_hze_pheno.txt
RESULT_DIR=${BASE_DIR}/results/plink

# Convert PLINK text files to binary. Run once, them comment out.
# Filter out SNPs with MAF < 0.05
# plink --tfile ${INPUT_FILE} \
#        --maf 0.05 \
#        --make-bed \
#        --out ${INPUT_FILE}

# Make genetic relatedness matrix.
plink --tfile ${INPUT_FILE} \
      --maf 0.05 \
      --make-rel square\
      --out ${RESULT_DIR}/hsnpt_hze

# NOTE: PLINK doesn't have a heritability estimate right now....

# Association mapping for all phenotypes.
plink --bfile ${INPUT_FILE} \
      --pheno ${PHENO_FILE} \
      --all-pheno \
      --assoc \
      --out ${RESULT_DIR}/hsnpt_hze

# NOTE: the files now end with '.assoc' in them.

# Add a column called BETA that is the log(OR).
Rscript --vanilla hsnpt_plink_add_beta.R `ls ${RESULT_DIR}/*.assoc`

# NOTE: the files now end with '.trans' in them.

# Clumping
for F in `ls ${RESULT_DIR}/*.trans`
do
   plink \
       --bfile ${INPUT_FILE} \
       --clump-p1 1 \
       --clump-r2 0.1 \
       --clump-kb 250 \
       --clump ${F} \
       --clump-snp-field SNP \
       --clump-field P \
       --out ${F}
done

# NOTE: the files now end with '.clumped' in them.

# Extract the clumped SNPs.
for F in `ls ${RESULT_DIR}/*.clumped`
do
   awk 'NR!=1{print $3}' ${F} > ${F/clumped/valid_snp}
done

# Create a file called 'range_list.txt' with different PRS p-value ranges
# to use.
echo $'0.001 0 0.001\n0.005 0 0.005\n0.01 0 0.01\n0.05 0 0.05\n0.01 0 0.01\n' > ${RESULT_DIR}/range_list.txt

# Calculate PRS
for F in `ls ${RESULT_DIR}/*.clumped`
do
   TRANS_FILE=${F/.clumped/}
   awk '{print $2,$9}' ${TRANS_FILE} > ${TRANS_FILE}.snppvalue
   plink \
      --bfile ${INPUT_FILE} \
      --score ${TRANS_FILE} 2 4 11 header \
      --q-score-range ${RESULT_DIR}/range_list.txt ${TRANS_FILE}.snppvalue \
      --extract ${TRANS_FILE}.valid_snp \
      --out ${TRANS_FILE}
done
