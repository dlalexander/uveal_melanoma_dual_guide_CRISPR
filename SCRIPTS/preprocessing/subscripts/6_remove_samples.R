suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character"),
  make_option(c("-s", "--samples"), type = "character",
              help = "comma-separated list of samples to remove", metavar = "character"),
  make_option(c("-f", "--fc"), type = "character",
              help = "full path to LFC matrix", metavar = "character"),
  make_option(c("-b", "--batch"), type = "character",
              help = "batch", metavar = "character")
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

# Process sample names
if (is.null(opt$samples)) {
  stop("Please provide sample names.")
}
samples_to_remove <- as.vector(str_split(opt$samples, ",", simplify = T))
if (0 == length(samples_to_remove) || grepl(',', samples_to_remove[1])) {
  stop(paste("Problem splitting sample_names:", opt$samples))
}

############################################################
# LFC matrix                                               #
############################################################

check_file_exists(opt$fc)

print(paste("Reading LFC matrix:", opt$fc))

# Read in LFC matrix
lfc_matrix <- read.delim(opt$fc, sep = "\t", header = T, check.names = F)


############################################################
# Specify batch.                                           #
############################################################

print(paste("Batch processing:", opt$batch))

# Create batch variable
batch <- opt$batch


############################################################
# Remove samples                                           #
############################################################

print('Validating samples...')

# Check samples exist in column names 

if (samples_to_remove != "") {
  missing_samples <- setdiff(samples_to_remove, colnames(lfc_matrix))
  if (length(missing_samples) > 0) {
    stop(paste('Samples not found in LFC matrix:', paste(missing_samples, collapse = ", ")))
  }
  print('Removing samples...')
  # Remove samples
  lfc_matrix.filt <- lfc_matrix %>% select(-all_of(samples_to_remove))
  # Check column names have been removed
  if (ncol(lfc_matrix) - ncol(lfc_matrix.filt) != length(samples_to_remove)) {
    print(paste('Pre-filter column number:', ncol(lfc_matrix)))
    print(paste('Post-filter column number:', ncol(lfc_matrix.filt)))
    stop('Incorrect number of columns post-filtering.')
  }
} else {
  print('No samples to remove.')
  lfc_matrix.filt <- lfc_matrix
}


############################################################
# Outputs                                                  #
############################################################

# Write filtered LFCs to file
lfc_matrix.filt.path <- file.path(repo_path, 'DATA', 'preprocessing', batch, 'lfc_matrix.unscaled.samples_removed.tsv')
write.table(lfc_matrix.filt,lfc_matrix.filt.path, sep = "\t", quote = F, row.names = F)

print(paste("Filtered LFC matrix written to:", lfc_matrix.filt.path))

# Save as RDS
lfc_matrix.filt.rds.path <- file.path(repo_path, 'DATA', 'RDS', 'preprocessing',batch, 'lfc_matrix.unscaled.samples_removed.rds')
saveRDS(lfc_matrix.filt, file = lfc_matrix.filt.rds.path )

print(paste("Filtered LFC matrix RDS written to:", lfc_matrix.filt.rds.path))

print("DONE.")