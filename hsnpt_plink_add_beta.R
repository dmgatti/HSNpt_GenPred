################################################################################
# Add an effect size column to the given PLINK assoc file.
# April 11, 2021
# Daniel Gatti
# dmgatti@coa.edu
################################################################################
options(stringsAsFactors = FALSE)
args = commandArgs(trailingOnly = TRUE)

library(data.table)

if(length(args) == 0) {
  stop('Need at least one file')
}

for(i in seq_along(args)) {
   
  input_file = args[i]
  output_file = paste0(input_file, '.trans')
  
  print(paste('Adding BETA to', input_file))

  data = fread(input_file, header = TRUE)
  data$BETA = log(data$OR)
  fwrite(data, output_file, sep = ' ', quote = FALSE, row.names = FALSE)
  
} # for(i)
