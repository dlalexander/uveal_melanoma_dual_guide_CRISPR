suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character"),
  make_option(c("-c", "--counts"), type = "character",
              help = "full path to count matrix", metavar = "character")
);

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

############################################################
# GENERAL                                                  #
############################################################

# Set top level directory
repo_path <- opt$dir

# Check top level directory exists
if (!dir.exists(repo_path)) {
  stop(sprintf("Repository directory not exist: %s", repo_path))
}

# Add helper functions
source(file.path(repo_path, 'SCRIPTS', 'preprocessing', 'subscripts', 'helper.R'))

############################################################
# Count matrix                                             #
############################################################

check_file_exists(opt$counts)

print(paste("Reading count matrix:", opt$counts))

# Read in count matrix
count_matrix <- read.delim(opt$counts, sep = "\t", header = T, check.names = F)

############################################################
# Normalise using BAGEL method                             #
############################################################

print("Adding pseudocount...")

# Add a pseudocount (5)
count_matrix.pseudo <- count_matrix %>% mutate(across(c(16:ncol(count_matrix)), ~ . + 5))

print("Total normalisation with scaling factor...")

# Total normalisation with scaling factor (10 million)
count_matrix.norm <- count_matrix.pseudo %>%
  mutate(across(16:ncol(count_matrix.pseudo), ~ (. / sum(.)) * 10000000))

# Check same number of guides as we started with
if (nrow(count_matrix) != nrow(count_matrix.norm)) {
  print(paste('Number of rows pre-normalisation:', nrow(count_matrix)))
  print(paste('Number of rows post-normalisation:', nrow(count_matrix.norm)))
  stop('Number of rows is not consistent.')
}

############################################################
# Outputs                                                  #
############################################################

# Write normalised counts to file
count_matrix.norm.path <- file.path(repo_path, 'DATA', 'preprocessing', 'count_matrix.norm.tsv')
write.table(count_matrix.norm, count_matrix.norm.path, sep = "\t", quote = F, row.names = F)

print(paste("Normalised count matrix written to:", count_matrix.norm.path))

# Save as RDS
count_matrix.norm.rds.path <- file.path(repo_path, 'DATA', 'RDS', 'preprocessing', 'count_matrix.norm.rds')
saveRDS(count_matrix.norm, file = count_matrix.norm.rds.path )

print(paste("Normalised count matrix RDS written to:", count_matrix.norm.rds.path))

print("DONE.")