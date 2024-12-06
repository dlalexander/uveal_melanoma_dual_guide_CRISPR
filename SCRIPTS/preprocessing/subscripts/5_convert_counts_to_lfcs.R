suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))


############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character"),
  make_option(c("-m", "--mapping"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
  make_option(c("-c", "--counts"), type = "character",
              help = "full path to count matrix", metavar = "character"),
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
# Specify batch.                                           #
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
# Calculate control mean                                   #
############################################################

print('Calculating control mean...')

# Get control mean
count_matrix.mean_control <- count_matrix %>%
  mutate('control_mean' = rowMeans(count_matrix[,control_sample_indices])) %>%
  select(everything(), -contains('control'), control_mean)

# Write filtered counts with mean control to file
count_matrix.mean_control.path <- file.path(repo_path, 'DATA', 'preprocessing',batch,'count_matrix.norm.filt.mean_control.tsv')
write.table(count_matrix.mean_control, count_matrix.mean_control.path, sep = "\t", quote = F, row.names = F)

print(paste("Count matrix with control mean written to:", count_matrix.mean_control.path))

############################################################
# Calculate log fold changes with respect to control mean  #
############################################################

print("Calculating unscaled LFCs...")

# Set annotation column names
annotation_colnames <- c(
  'id', 'sgrna_ids', 'sgrna_seqs', 'gene_pair_id',
  'sorted_gene_pair', 'targetA', 'targetB', 
  'sgrna_symbols', 'sgrna_symbol_A', 'sgrna_symbol_B', 
  'sgrna_libraries', 'sgrna_group', 'guide_type', 
  'guide_orientation', 'singles_target_gene')

# Calculate log fold changes
lfc_matrix <- suppressMessages(calculate_lfc(count_matrix.mean_control,
                            control_indices = which(colnames(count_matrix.mean_control) == 'control_mean'),
                            treatment_indices = which(! colnames(count_matrix.mean_control) %in% c(annotation_colnames, 'control_mean')),
                            pseudocount =  0.5))

# Add annotations back into matrix
lfc_matrix.ann <- lfc_matrix %>%
 left_join(count_matrix.mean_control[,1:15], by = c('id', 'sgrna_ids')) %>%
 relocate(annotation_colnames[3:15], .after = 'sgrna_ids')

# Check number of rows in count and LFC matrix are the same
if (nrow(count_matrix) != nrow(lfc_matrix.ann)) {
  print(paste('Number of count matrix rows:', nrow(count_matrix)))
  print(paste('Number of LFC matrix rows:', nrow(lfc_matrix.ann)))
  stop('Number of rows is not consistent.')
}

############################################################
# Outputs                                                  #
############################################################

# Write normalised/filtered log2 fold changes to file
lfc_matrix.path <- file.path(repo_path, 'DATA', 'preprocessing',batch, 'lfc_matrix.unscaled.all.tsv')
write.table(lfc_matrix.ann, lfc_matrix.path, sep = "\t", quote = F, row.names = F)

print(paste("Unscaled LFC matrix written to:", lfc_matrix.path))

# Save as RDS
lfc_matrix.rds.path <- file.path(repo_path, 'DATA', 'RDS', 'preprocessing',batch, 'lfc_matrix.unscaled.all.rds')
saveRDS(lfc_matrix.ann, file = lfc_matrix.rds.path )

print(paste("Unscaled LFC matrix RDS written to:", lfc_matrix.rds.path ))

print("DONE.")
