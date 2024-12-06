suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character"),
  make_option(c("-c", "--counts"), type = "character",
              help = "full path to count matrix", metavar = "character"),
  make_option(c("-m", "--mapping"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
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

############################################################
# Count matrix                                             #
############################################################

check_file_exists(opt$counts)

print(paste("Reading count matrix:", opt$counts))

# Read in count matrix
count_matrix <- read.delim(opt$counts, sep = "\t", header = T, check.names = F)

############################################################
# Sample mapping                                           #
############################################################

check_file_exists(opt$mapping)

print(paste("Reading sample annotations from:", opt$mapping))

# Read in sample mapping
sample_mapping <- read.delim(opt$mapping, header = T, sep = "\t")

############################################################
# Specify batch                                            #
############################################################


print(paste("Batch processing:", opt$batch))

# Create batch variable
batch <- opt$batch


############################################################
# Identify control samples                                 #
############################################################

print("Identifying control samples...")

# Identify indices of control samples
control_samples <-  sample_mapping %>%
  filter(grepl('Control', sample_label))

control_sample_labels <- control_samples %>%
  pull(sample_label)

print('Control sample labels:')
print(control_sample_labels)

control_sample_indices <- which(colnames(count_matrix) %in% control_sample_labels)

############################################################
# Filter low count guides (in controls)                    #
############################################################

print('Identifying low count guides...')

# Identify low count guides
# Edited this to specify the correct count columns
filt_guides <- get_guides_failing_filter(count_matrix,
                                         id_column = 1,
                                         count_column = 16:ncol(count_matrix),
                                         filter_indices = control_sample_indices,
                                         filter_method = 'mean',
                                         min_reads = 20)

print(paste('Number of low count guides ids:', length(filt_guides)))

print('Removing low count guides...')

# Remove low count guides from normalised counts
count_matrix.filt <- count_matrix %>%
  filter(!id %in% filt_guides)

# Check that the correct number of guides have been removed
num_guide_to_remove <- length(filt_guides)
print(paste('Count matrix rows:', nrow(count_matrix)))
print(paste('Filtered count matrix rows:', nrow(count_matrix.filt)))
if (nrow(count_matrix) - nrow(count_matrix.filt) != num_guide_to_remove) {
  print('Unexpected number of guides removed')
}

############################################################
# Outputs                                                  #
############################################################

# Write filtered guides to file
filt_guides.path <- file.path(repo_path, 'DATA', 'preprocessing',batch,'filtered_guides.txt')
write.table(filt_guides, filt_guides.path, sep = "\t", quote = F, row.names = F, col.names = F)

print(paste("Filtered guide ids written to:", filt_guides.path))

# Write filtered counts to file
count_matrix.filt.path <- file.path(repo_path, 'DATA', 'preprocessing',batch,'count_matrix.norm.filt.tsv')
write.table(count_matrix.filt, count_matrix.filt.path, sep = "\t", quote = F, row.names = F)

print(paste("Filtered count matrix written to:", count_matrix.filt.path))

# Save as RDS
count_matrix.filt.rds.path <- file.path(repo_path, 'DATA', 'RDS', 'preprocessing',batch,'count_matrix.norm.filt.rds')
saveRDS(count_matrix.filt, file = count_matrix.filt.rds.path )

print(paste("Filtered count matrix RDS written to:", count_matrix.filt.rds.path))

print("DONE.")
