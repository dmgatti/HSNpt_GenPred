# Create sqlite database with variants (SNPs/indels/SVs)
# in the 8 founder strains of the Collaborative Cross (CC)
# Karl Broman
# from https://github.com/rqtl/qtl2/blob/master/inst/scripts/create_ccvariants.R
#
######################################################################
# sources:
#   SNPs:
#     ftp://ftp-mouse.sanger.ac.uk/REL-2004-v7-SNPs_Indels/mgp_REL2005_snps_indels.vcf.gz
#     ftp://ftp-mouse.sanger.ac.uk/REL-2004-v7-SNPs_Indels/mgp_REL2005_snps_indels.vcf.gz.tbi
#   Indels:
#     ftp://ftp-mouse.sanger.ac.uk/REL-2004-v7-SNPs_Indels/mgp_REL2005_snps_indels.vcf.gz
#     ftp://ftp-mouse.sanger.ac.uk/REL-2004-v7-SNPs_Indels/mgp_REL2005_snps_indels.vcf.gz.tbi
#   SVs:
#     ftp://ftp-mouse.sanger.ac.uk/REL-1606-SV/mgpv5.SV_deletions.bed.gz
#     ftp://ftp-mouse.sanger.ac.uk/REL-1606-SV/mgpv5.SV_deletions.bed.gz.tbi
#     ftp://ftp-mouse.sanger.ac.uk/REL-1606-SV/mgpv5.SV_insertions.bed.gz
#     ftp://ftp-mouse.sanger.ac.uk/REL-1606-SV/mgpv5.SV_insertions.bed.gz.tbi
#
######################################################################

##############################
### download files
##############################
# site <- c("ftp://ftp-mouse.sanger.ac.uk",
#           "ftp://ftp-mouse.sanger.ac.uk",
#           "ftp://ftp-mouse.sanger.ac.uk",
#           "ftp://ftp-mouse.sanger.ac.uk")
# subdir <- c("REL-2004-v7-SNPs_Indels",
#             "REL-2004-v7-SNPs_Indels",
#             "REL-1606-SV",
#             "REL-1606-SV")
# files <- c("mgp_REL2005_snps_indels.vcf.gz",
#            "mgp_REL2005_snps_indels.vcf.gz",
#            "mgpv5.SV_deletions.bed.gz",
#            "mgpv5.SV_insertions.bed.gz")
# date_source <- c("2020-05-20",
#                  "2020-05-20",
#                  "2016-06-20",
#                  "2016-06-20")
# genome_build <- rep("GRCm38/mm10", 4)

# for(i in seq_along(files)) {
#   file <- files[i]
#   url <- paste0(site[i], "/", subdir[i], "/", file)
#   tbi_file <- paste0(file, ".tbi")
#   tbi_url <- paste0(site[i], "/", subdir[i], "/", tbi_file)
#   
#   if(!file.exists(file)) {
#     cat(" -Downloading", file, "\n")
#     download.file(url, file)
#   }
#   if(!file.exists(tbi_file)) {
#     cat(" -Downloading", tbi_file, "\n")
#     download.file(tbi_url, tbi_file)
#   }
# }

library(VariantAnnotation)
library(RSQLite)

# function for re-formating snp/indel "consequence"
# Consequences: vectors of |-separated valued
#    2nd value is gene ID
#    3rd is transcript ID (ignored)
#    5th is consequence (potentially with multiple separated by &'s)
format_consequence <-
  function(csq_record)
  {
    x <- strsplit(csq_record, "|", fixed=TRUE)
    genes <- sapply(x, "[", 2)
    csq <- sapply(x, "[", 5)
    csq <- unlist(lapply(seq_along(csq),
                         function(i) {
                           if(genes[i] == "") {
                             return(csq[i])
                           } else {
                             tmp <- unlist(strsplit(csq[i], "&", fixed=TRUE))
                             return(paste(genes[i], tmp, sep=":"))
                           }}))
    c(paste(unique(genes), collapse=","),
      paste(unique(csq), collapse=","))
  }





##############################
### SNPs
##############################
files = '/media/dmgatti/hda/data/Sanger/REL-2005-SNPs_Indels/mgp_REL2005_snps_indels.vcf.gz'

header = scanVcfHeader(files[1])
chr <- seqnames(seqinfo(header))
chr = chr[order(as.numeric(chr))]
hs_founders <- c('A/J', 'AKR/J', 'BALB/cJ', 'C3H/HeJ', 'C57BL/6J', 'CBA/J', 'DBA/2J', 'LP/J')
strains <- sub("/", "_", hs_founders[hs_founders != 'C57BL/6J'])
n_strains <- length(strains)

db_file <- "/media/dmgatti/hda/data/MUGA/hsnpt_variants.sqlite"
db <- dbConnect(SQLite(), dbname=db_file)

dbExecute(db, paste0("ATTACH '", db_file, "' AS NEW"))

cat(" -SNPs\n")
tabfile <- TabixFile(files[1], paste0(files[1], ".tbi"))

db_started <- FALSE

for(thechr in chr) {
  for(left in seq(0, 190, by=10)) {
    cat(thechr, left, "\n")
    
    # 10 Mbp range
    gr <- GRanges(seqnames=thechr, ranges=IRanges(start=left*1e6, end=(left+10)*1e6-1))
    
    # grab data
    param <- ScanVcfParam(geno = c("GT", "FI"), samples = strains,
                          which = gr)
    snps <- readVcf(tabfile, genome = "mm10", param = param)
    if(nrow(snps)==0) next
    
    # drop snps with any quality < 1
    fi <- geno(snps)$FI
    snps <- snps[rowSums(!is.na(fi) & fi==1) == n_strains]
    
    # drop snps that are indels (for now).
    indel = info(snps)$INDEL
    snps = snps[!indel,]
    if(nrow(snps)==0) next
    
    # drop snps that are all 0/0
    g <- geno(snps)$GT
    snps <- snps[rowSums(is.na(g)) == 0 & rowSums(g=="0/0") < n_strains]
    g <- geno(snps)$GT
    if(nrow(snps)==0) next
    
    # grab genotypes
    g <- geno(snps)$GT
    # add B6 genotypes (reference) and change column names
    g <- cbind(g[,1:4], C57BL_6J="0/0", g[,5:7])
    colnames(g) <- hs_founders
    
    # alleles
    major <- as.character(ref(snps))
    minor <- CharacterList(alt(snps))
    alleles <- matrix("", nrow=nrow(snps), ncol=4)
    alleles[,1] <- major
    for(i in 2:4) {
      alleles[,i] <- sapply(minor, "[", i-1)
    }
    
    # numbers of each genotype
    rs <- sapply(c("0/0", "1/1", "2/2", "3/3"),
                 function(a) rowSums(g==a))
    
    # rows should all sum to 8
    stopifnot(all(rowSums(rs)==8))
    
    # find the most common allele and swap with major
    for(i in ncol(rs):2) {
      # this allele most common?
      wh <- (rs[,i] > rs[,1] & rowSums(rs <= rs[,i]) == ncol(rs))
      
      if(any(wh)) {
        # swap alleles
        tmp <- alleles[wh,i]
        alleles[wh,i] <- alleles[wh,1]
        alleles[wh,1] <- tmp
        
        # swap genotypes
        pat <- paste0(i-1, "/", i-1)
        gg <- g[wh,,drop=FALSE]
        gg[gg==pat] <- "x/x"
        gg[gg=="0/0"] <- pat
        gg[gg=="x/x"] <- "0/0"
        g[wh,] <- gg
      }
    }
    
    # version with A/C/G/T
    glet <- gnum <- matrix(nrow=nrow(g), ncol=ncol(g))
    dimnames(glet) <- dimnames(gnum) <- dimnames(g)
    for(i in 1:4) {
      pat <- paste0(i-1, "/", i-1)
      for(j in 1:ncol(g)) {
        wh <- (g[,j] == pat)
        glet[wh,j] <- alleles[wh,i]
        gnum[wh,j] <- i
      }
    }
    
    # NAs in alleles when not seen
    for(i in 1:4) {
      wh <- which(rowSums(gnum==i)==0)
      alleles[wh,i] <- NA
    }
    
    # alleles as a pattern A|C/G/T
    alleles_char <- paste(alleles[,1],
                          apply(alleles[,-1,drop=FALSE], 1,
                                function(a) paste(a[!is.na(a)], collapse="/")),
                          sep="|")
    
    # gnum: turn it into consecutive numbers
    #   (if "2" is missing make 3 = 2)
    for(i in 3:2) {
      wh <- is.na(alleles[,i])
      if(any(wh)) {
        tmp <- gnum[wh,]
        tmp[tmp >= i] <- tmp[tmp >= i] - 1
        gnum[wh,] <- tmp
      }
    }
    
    # convert to 1/3
    gbin <- gnum
    gbin[gbin > 1] <- 3
    
    # first row = gene; 2nd row = consequence
    csq <- sapply(info(snps)$CSQ, format_consequence)
    
    # create full table of info
    snps <- data.frame(snp_id=rownames(g),
                       chr=as.vector(seqnames(snps)),
                       pos=start(snps),
                       alleles=alleles_char,
                       sdp=qtl2::calc_sdp(gbin),
                       ensembl_gene=csq[1,],
                       consequence=csq[2,],
                       gnum,
                       type="snp",
                       stringsAsFactors=FALSE)
    # make sure column names are what we want
    #colnames(snps)[8:15] <- c(strains[1], "C57BL_6J", strains[-1])
    colnames(snps) = gsub('\\.', '_', colnames(snps))
    
    dbWriteTable(db, "variants", snps, row.names = FALSE, overwrite = !db_started,
                 append = db_started, field.types = NULL)
    db_started <- TRUE
    
  } # for(left)
} # for(thechr)

##############################
### indels
##############################
cat(" -InDels\n")
chr <- c(1:19, "X", "Y") # drop MT because not present in indel file
tabfile <- TabixFile(files[2], paste0(files[2], ".tbi"))

for(thechr in chr) {
  for(left in seq(0, 190, by=10)) {
    cat(thechr, left, "\n")
    
    # 10 Mbp range
    gr <- GRanges(seqnames=thechr, ranges=IRanges(start=left*1e6, end=(left+10)*1e6-1))
    
    # grab data
    param <- ScanVcfParam(geno = c("GT", "FI"), samples = strains,
                          which = gr)
    indels <- readVcf(tabfile, genome = "mm10", param = param)
    if(nrow(indels)==0) next
    
    # drop indels with any quality < 1
    fi <- geno(indels)$FI
    indels <- indels[rowSums(!is.na(fi) & fi==1) == n_strains]
    if(nrow(indels)==0) next
    
    # drop indels that are all 0/0
    g <- geno(indels)$GT
    indels <- indels[rowSums(is.na(g)) == 0 & rowSums(g=="0/0") < n_strains]
    g <- geno(indels)$GT
    if(nrow(indels)==0) next
    
    # add B6 genotypes (reference) and change column names
    g <- cbind(g[,1,drop=FALSE], C57BL_6J="0/0", g[,-1,drop=FALSE])
    colnames(g) <- hs_founders
    
    # alleles
    major <- as.character(ref(indels))
    minor <- CharacterList(alt(indels))
    max_minor <- max(sapply(minor, length))
    alleles <- matrix("", nrow=nrow(g), ncol=max_minor+1)
    alleles[,1] <- major
    for(i in 2:(max_minor+1))
      alleles[,i] <- sapply(minor, "[", i-1)
    
    # numbers of each genotype
    pat <- paste0(0:max_minor, "/", 0:max_minor)
    rs <- sapply(pat, function(a) rowSums(g==a))
    if(!is.matrix(rs)) {
      rs <- matrix(rs, nrow=1)
      dimnames(rs) <- list(rownames(g), pat)
    }
    
    # rows should all sum to 8
    stopifnot(all(rowSums(rs)==8))
    
    # find the most common allele and swap with major
    for(i in ncol(rs):2) {
      # this allele most common?
      wh <- (rs[,i] > rs[,1] & rowSums(rs <= rs[,i]) == ncol(rs))
      
      if(any(wh)) {
        # swap alleles
        tmp <- alleles[wh,i]
        alleles[wh,i] <- alleles[wh,1]
        alleles[wh,1] <- tmp
        
        # swap genotypes
        pat <- colnames(rs)[i]
        gg <- g[wh,,drop=FALSE]
        gg[gg==pat] <- "x/x"
        gg[gg=="0/0"] <- pat
        gg[gg=="x/x"] <- "0/0"
        g[wh,] <- gg
      }
    }
    
    # version with alleles
    glet <- gnum <- matrix(nrow=nrow(g), ncol=ncol(g))
    dimnames(glet) <- dimnames(gnum) <- dimnames(g)
    for(i in 1:ncol(rs)) {
      pat <- paste0(i-1, "/", i-1)
      for(j in 1:ncol(g)) {
        wh <- (g[,j] == pat)
        glet[wh,j] <- alleles[wh,i]
        gnum[wh,j] <- i
      }
    }
    
    # NAs in alleles when not seen
    for(i in 1:ncol(alleles)) {
      wh <- which(rowSums(gnum==i)==0)
      alleles[wh,i] <- NA
    }
    
    # alleles as a pattern A|C/G/T
    alleles_char <- paste(alleles[,1],
                          apply(alleles[,-1,drop=FALSE], 1,
                                function(a) paste(a[!is.na(a)], collapse="/")),
                          sep="|")
    
    # gnum: turn it into consecutive numbers
    #   (if "2" is missing make 3 = 2)
    if(ncol(alleles) >= 3) {
      for(i in 3:2) {
        wh <- is.na(alleles[,i])
        if(any(wh)) {
          tmp <- gnum[wh,]
          tmp[tmp >= i] <- tmp[tmp >= i] - 1
          gnum[wh,] <- tmp
        }
      }
    }
    
    # convert to 1/3
    gbin <- gnum
    gbin[gbin > 1] <- 3
    
    # first row = gene; 2nd row = consequence
    csq <- sapply(info(indels)$CSQ, format_consequence)
    
    # create full table of info
    indels <- data.frame(snp_id=rownames(g),
                         chr=thechr,
                         pos=start(indels),
                         alleles=alleles_char,
                         sdp=qtl2::calc_sdp(gbin),
                         ensembl_gene=csq[1,],
                         consequence=csq[2,],
                         gnum,
                         type="indel",
                         stringsAsFactors=FALSE)
    # make sure column names are what we want
    colnames(indels)[8:15] <- c(strains[1], "C57BL_6J", strains[-1])
    
    dbWriteTable(db, "variants", indels, row.names=FALSE, overwrite=FALSE,
                 append=TRUE, field.types=NULL)
  }
}


##############################
# structural variants (SVs)
##############################
cat(" -Stuctural variants\n")
tmpfile <- tempfile()
system(paste0("gunzip -c ", files[3], " > ", tmpfile))
svs <- data.table::fread(tmpfile, data.table=FALSE)
unlink(tmpfile)

# pull out genotypes
g <- svs[,colnames(svs) %in% strains]
g <- g[,strains]

# add B6 ref genotype
g <- cbind(g[,1,drop=FALSE], C57BL_6J="0", g[,-1,drop=FALSE], stringsAsFactors=FALSE)

# drop monomorphic ones
n_allele <- apply(g, 1, function(a) length(unique(a)))
svs <- svs[n_allele > 1,]
g <- g[n_allele > 1,]
g[g=="0"] <- "-"

# alleles
alleles <- apply(g, 1, function(a) {
  tab <- table(a)
  result <- names(sort(tab, decreasing=TRUE))
  if(length(result) < 8)
    result <- c(result, rep(NA, 8-length(result)))
  result })
alleles <- t(alleles)

# numeric version
gnum <- matrix(nrow=nrow(g), ncol=ncol(g))
dimnames(gnum) <- dimnames(g)
for(i in 1:ncol(alleles)) {
  for(j in 1:ncol(g)) {
    wh <- (!is.na(alleles[,i]) & g[,j] == alleles[,i])
    gnum[wh,j] <- i
  }
}

# alleles as a pattern A|C/G/T
alleles_char <- paste(alleles[,1],
                      apply(alleles[,-1,drop=FALSE], 1,
                            function(a) paste(a[!is.na(a)], collapse="/")),
                      sep="|")

# convert to 1/3
gbin <- gnum
gbin[gbin > 1] <- 3

# create full table of info
svs <- data.frame(snp_id=paste0("SV_", svs[,"#CHROM"], "_", svs[,"START"], "_", svs[,"END"]),
                  chr=svs[,"#CHROM"],
                  pos=round((svs[,"START"] + svs[,"END"])/2), # average of start and end
                  alleles=alleles_char,
                  sdp=qtl2::calc_sdp(gbin),
                  ensembl_gene=NA,
                  consequence=NA,
                  gnum,
                  type="SV",
                  stringsAsFactors=FALSE)
# make sure column names are what we want
colnames(svs)[8:15] <- c(strains[1], "C57BL_6J", strains[-1])

dbWriteTable(db, "variants", svs, row.names=FALSE, overwrite=FALSE,
             append=TRUE, field.types=NULL)


##############################
### Add description table
##############################

# description <- data.frame(description=c("SNPs in Collaborative Cross founders",
#                                         "Indels in Collaborative Cross founders",
#                                         "SVs in Collaborative Cross founders"),
#                           source=c("Mouse Genome Informatics (MGI), Jackson Lab",
#                                    "Mouse Genome Informatics (MGI), Jackson Lab",
#                                    "Sanger"),
#                           url=paste0(site, "/", subdir, "/", files),
#                           date_created=rep(as.character(Sys.Date()), 3),
#                           date_source=date_source,
#                           genome_build=genome_build,
#                           stringsAsFactors=FALSE)
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
dbWriteTable(db, "description", description, append=TRUE)


##############################
### add index
##############################
dbExecute(db, "CREATE INDEX chr_pos ON variants(chr, pos)")
dbDisconnect(db)

