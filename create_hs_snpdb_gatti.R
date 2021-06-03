################################################################################
# Build a qtl2-style SQLLite data base of Northport/HS SNPs.
#
# NOTE: I'm only including SNPs where all strains have homozygous calls.
#       TBD: Handling heterozygotes, no-calls and low quality SNPs.
# 
# Daniel Gatti
# Feb. 22, 2021
# dmgatti@coa.edu
################################################################################

# HS/Npt founder strains: A/J, AKR/J, BALBc/J, CBA/J, C3H/HeJ, C57BL/6J, DBA/2J, LP/J

library(qtl2)
library(RSQLite)
library(VariantAnnotation)

base_dir    = '/media/dmgatti/hdb/projects/ColoState'
sanger_dir  = '/media/dmgatti/hda/data/Sanger/REL-1807-SNPs_Indels'
sanger_file = file.path(sanger_dir, 'mgp.v6.merged.norm.snp.indels.sfiltered.vcf.gz')
founder_names = c('A/J', 'AKR/J', 'BALB/cJ', 'CBA/J', 'C3H/HeJ', 'C57BL/6J', 'DBA/2J', 'LP/J')
cc_var_file = '/media/dmgatti/hda/data/MUGA/cc_variants.sqlite'
hs_var_file = '/media/dmgatti/hda/data/MUGA/hsnpt_variants.sqlite'

##### FUNCTIONS #####

# Extract consequence information and return unique Ensembl genes
# and ensmebl/consequence combinations.
# NOTE: This may not be the fastest way.
process_csq = function(csq) {

  # Split the data by '|'.  
  csq = sapply(csq, strsplit, split = '\\|')
  # Keep the first 5 columns of each CSQ.
  # Column 2 has the ENSMUSG and column 5 has the consequence.
  csq = lapply(csq, lapply, '[', c(2, 5))
  # Unlist. There are often multiplt CSQ per SNP.
  csq = lapply(csq, unlist)
  # Create matrices for each CSQ to make it easier to parse.
  csq = lapply(csq, matrix, ncol = 2, byrow = TRUE)
  # Get unique ENSMUSG ids and paste together using ','.
  eg = lapply(csq, function(z) { unique(z[,1]) })
  eg = sapply(eg, paste0, collapse = ',')

  # Paste ENSMUSG and CSQ with ':'.
  csq = lapply(csq, function(z) { paste(z[,1], z[,2], sep = ':') })
  # Keep unique ENSMUSG:CSQ combinations.
  csq = lapply(csq, unique)
  # Paste the ENSMUSG:CSQ values together with ','.
  csq = sapply(csq, paste, collapse = ',')

  return(data.frame(ensembl_gene = eg, consequence = csq))
  
} # process_csq()

# Extract the genotypes and convert them to qtl2 output format.
# We expect only homozygotes right now.
# TBD: This will need to be overhauled if we handle hets or no-calls.
process_gt = function(gt) {
  
  # Add C57BL/6J.
  gt = cbind(gt[,1:4], C57BL_6J = rep('0/0', nrow(gt)), gt[,5:7])
  rn = rownames(gt)
  cn = colnames(gt)

  # Get first character of each call.
  gt = substring(gt, 1, 1)
  gt = matrix(as.numeric(gt), nrow = nrow(gt), dimnames = list(rn, cn))
  # Add one because qtl2 uses 1 as the ref call.
  gt = gt + 1

  return(gt)

} # process_gt()


# Calculate the SDP. Karl's qtl2::calc_sdp() requires 1/3 calls.
calc_sdp = function(gt) {
  
  # Convert to 1/3 format.
  gt[gt != 1] = 3
  return(qtl2::calc_sdp(gt))
    
} # calc_sdp()


##### ANALYSIS #####

# Get the header information from the Sanger file.
header = scanVcfHeader(file = sanger_file)

# Are all of the HS/Npt founders in the Sanger SNP file?
sum(sub('/', '_', founder_names) %in% samples(header))

founder_names[which(!sub('/', '_', founder_names) %in% samples(header))]

# We expect C57BL/6J to be missing since it's the reference.
# So we have all eight founder strains.
founder_names = sort(sub('/', '_', founder_names))

# Get the chromosome names and keep the standard ones.
chr_info = keepStandardChromosomes(seqinfo(header))

# Build the description table.
description = data.frame(description = c('SNPs & Indels in HS/Npt founders', 
                                         'Large deletions in HS/Npt founders',
                                         'Large insertions in HS/Npt founders'), 
                         source      = c('Sanger Mouse Genomes Project', 'Sanger Mouse Genomes Project', 'Sanger Mouse Genomes Project'),
                         url         = c('ftp://ftp-mouse.sanger.ac.uk/REL-2004-v7-SNPs_Indels/mgp_REL2005_snps_indels.vcf.gz', 
                                         'ftp://ftp-mouse.sanger.ac.uk/REL-1606-SV/mgpv5.SV_deletions.bed.gz',
                                         'ftp://ftp-mouse.sanger.ac.uk/REL-1606-SV/mgpv5.SV_insertions.bed.gz'),
                         date_created = c('2021-02-22', '2021-02-22', '2021-02-22'), 
                         date_source  = c('2021-02-22', '2021-02-22', '2021-02-22'), 
                         genom_build  = c('GRCm38', 'GRCm38', 'GRCm38'))

# Create HS/Npt SNP database.
hsdb = dbConnect(RSQLite::SQLite(), hs_var_file)

# Write out 'description' table.
dbWriteTable(hsdb, 'description', description)

# Build the 'variants' table.
# Size of data chunks to pull from Sanger DB.
chunk_size = 1e7

# For each chromosome.
for(chr in seqnames(chr_info)) {
  
  print(paste('CHR', chr, '****************************'))
  
  # Get data in chunks of 10 Mb.
  curr_start = 0
  curr_end   = 0
  
  while(curr_start < seqlengths(chr_info)[chr]) {

    print(paste('   ', curr_start / 1e6))
    curr_end = curr_start + chunk_size
    
    # Read SNPs from Sanger VCF.
    gr = GRanges(seqnames = chr, ranges = IRanges(start = curr_start, end = curr_end))
    param = ScanVcfParam(fixed = c('ALT'), info = c('INDEL', 'CSQ', 'MQ'), geno = c('GT', 'FT', 'GQ'), 
                         samples = founder_names, which = gr)
    snps = readVcf(file = sanger_file, param = param)
    
    # Filter to keep SNPS for v1. (i.e. !INDEL is a SNP)
    snps = snps[!info(snps)$INDEL,]
    
    # Filter to keep polymorphic SNPs.
    # For now, I define this as all strains != '0/0'.
    # This might not be correct all of the time.
    snps = snps[rowMeans(geno(snps)$GT == '0/0') < 1.0,]

    # Filter to keep only homozygotes.
    # TBD: How to handle heterozygotes, missing calls, and low quality SNPs.
    snps = snps[rowMeans(matrix(geno(snps)$GT %in% c('0/0', '1/1', '2/2', '3/3', '4/4', '5/5', '6/6', '7/7'), 
                                nrow = nrow(geno(snps)$GT))) == 1,]

    if(length(snps) != 0) {
    
      # Get the alleles for each variant.
      ref = as.character(fixed(snps)$REF)
      alt = CharacterList(fixed(snps)$ALT)
      alt = unstrsplit(alt, sep = '/')
      alleles = paste(ref, alt, sep = '|')
      
      # Convert consequences to ENSMUSG and consequences.
      csq_result = process_csq(as.list(info(snps)$CSQ))
      
      # Convert genotypes to qtl2 format.
      gt = process_gt(geno(snps)$GT)
      # Get SDP for each SNP.
      sdps = calc_sdp(gt)
      
      # Build type string.
      type = rep('snp', nrow(snps))
      type[info(snps)$INDEL] = 'indel'
      
      # Build table to write to DB.
      output = data.frame(snp_id  = names(rowRanges(snps)),
                          chr     = seqnames(rowRanges(snps)),
                          pos     = start(rowRanges(snps)),
                          alleles = alleles,
                          sdp     = sdps,
                          ensembl_gene = csq_result$ensembl_gene,
                          consequence  = csq_result$consequence,
                          gt,
                          type    = type)
  
      # Write to DB file.
      dbWriteTable(conn = hsdb, name = 'variants', value = output, append = TRUE)
      
      # Clean up a bit.
      rm(snps, gt, output)

    } else {
      rm(snps)
    } # else
    
    # Increment genome chunk location.
    curr_start = curr_start + chunk_size

  } # while(curr_start < seqlengths(chr_info)[chr])
  
} # for(chr)

# Disconnect from HS/Npt database.
dbDisconnect(hsdb)


################################################
### CC Variants file for example data.
ccdb = dbConnect(RSQLite::SQLite(), cc_var_file)
dbListTables(ccdb)
dbGetQuery(ccdb, 'SELECT * FROM description')

cc_var_example = dbGetQuery(ccdb, "SELECT * FROM variants WHERE chr='1' And pos>17000000 AND pos<172000000")

dbDisconnect(ccdb)



