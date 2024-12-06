suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character"),
  make_option(c("-f", "--fc"), type = "character",
              help = "full path to lfc matrix", metavar = "character"),
  make_option(c("-m", "--mapping"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
  make_option(c("-e", "--essentials"), type = "character",
              help = "full path to essentials gene list for scaling", metavar = "character"),
  make_option(c("-s", "--suffix"), type = "character",
              help = "file name suffix", metavar = "character"),
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


# Determine file extensions with suffix
if (is.null(opt$suffix)) {
  tsv_ext <- 'tsv'
  rds_ext <- 'rds'
} else {
  tsv_ext <- paste0(opt$suffix, '.tsv')
  rds_ext <- paste0(opt$suffix, '.rds')
}

# Set annotation column names
annotation_colnames <- c(
  'id', 'sgrna_ids', 'sgrna_seqs', 'gene_pair_id',
  'sorted_gene_pair', 'targetA', 'targetB', 
  'sgrna_symbols', 'sgrna_symbol_A', 'sgrna_symbol_B', 
  'sgrna_libraries', 'sgrna_group', 'guide_type', 
  'guide_orientation', 'singles_target_gene')

############################################################
# LFC matrix                                               #
############################################################

check_file_exists(opt$fc)

print(paste("Reading LFC matrix:", opt$fc))

# Read in LFC matrix
lfc_matrix <- read.delim(opt$fc, sep = "\t", header = T, check.names = F)

############################################################
# Sample mapping                                           #
############################################################

check_file_exists(opt$mapping)

print(paste("Reading sample annotations from:", opt$mapping))

# Read in sample mapping
sample_mapping <- read.delim(opt$mapping, header = T, sep = "\t")


############################################################
# Essentials for scaling                                   #
############################################################

check_file_exists(opt$essentials)

print(paste("Reading essential genes from:", opt$essentials))

# Read in sample mapping
essentials <- read.delim(opt$essentials, header = T, sep = "\t")


############################################################
# Specify batch.                                           #
############################################################


print(paste("Batch processing:", opt$batch))

# Create batch variable
batch <- opt$batch


############################################################
# Prepare LFC matrix                                       #
############################################################

print('Preparing LFC matrix...')

# Gather LFC matrix
lfc.narrow <- lfc_matrix %>%
  gather(sample_label, LFC, -all_of(annotation_colnames)) 

# Annotate with sample information
lfc.narrow <- lfc.narrow %>%
  left_join(sample_mapping %>% select(sample_label, cell_line_label, replicate, cancer_type), by = 'sample_label') %>%
  relocate('LFC', .after = 'cancer_type')

############################################################
# Scale LFC matrix                                         #
############################################################

# Add additional annotation column using bagel essentials list
print('Annotating bagel essentials')
lfc.narrow <- lfc.narrow %>%
  mutate(essential_status = case_when(
    sgrna_group == "safe_safe|safe_safe" ~ 'safe_safe',
    sgrna_group == "Non_essential|Non_essential"  ~ 'non_essential',
    singles_target_gene %in% essentials$GENE ~ 'essential',
    TRUE ~ 'single_unknown' #For when other conditions are not met
  ))

# Get median of all replicates in a cell line by library type (e.g. safe_safe)
print('Calculating replicate medians per cell line by essentiality status..')
lfc.narrow.medians <- lfc.narrow %>%
  group_by(cell_line_label, essential_status) %>%
  mutate(median_lfc_cell_line_essentiality = median(LFC)) %>%
  ungroup()

# Calculate safe and essential medians
print('Calculating safe and essential medians...')
ess_group_medians <- lfc.narrow.medians %>%
  filter(essential_status %in% c('essential', 'safe_safe')) %>%
  group_by(cell_line_label, essential_status) %>%
  summarise(median = median(LFC), .groups = 'keep') %>%
  spread(essential_status, median) %>%
  rename( safe_safe_median = `safe_safe`, essential_median = `essential`)

# Add a new column with the median of the safe_safes for that cell line
print('Combining median data frames...')
lfc.narrow.medians <- lfc.narrow.medians %>%
  left_join(ess_group_medians, by = 'cell_line_label')

# Calculated scaled LFC
# There will be a new essential median after minusing safe_safe median, so need to account for that by adding safe_safe_median to denominator
print('Scaling LFCs...')
lfc.narrow.scaled <- lfc.narrow.medians %>%
  mutate(scaled_LFC = ((LFC - safe_safe_median) / (safe_safe_median - essential_median)))

# Spread scaled LFCs
lfc.scaled.wide <- lfc.narrow.scaled %>%
  select(all_of(annotation_colnames), sample_label, scaled_LFC) %>%
  spread(sample_label, scaled_LFC)

# Check that the new matrix is the correct length
if (nrow(lfc_matrix) != nrow(lfc.scaled.wide)) {
  print(paste('Number of rows pre-scaling:', nrow(lfc_matrix)))
  print(paste('Number of rows post-scaling:', nrow(lfc.scaled.wide)))
  stop("Incorrect number of rows after scaling.")
}
  
# Set all positive LFCs to 0
print('Setting positive LFCs to 0...')
lfc.narrow.no_positives <- lfc.narrow.scaled %>%
  mutate(scaled_LFC = ifelse(scaled_LFC > 0, 0, scaled_LFC))
print(paste('Number of positive ids pre-transform:', lfc.narrow.scaled %>% filter(scaled_LFC > 0) %>% nrow()))
print(paste('Number of positive ids post-transform:', lfc.narrow.no_positives %>% filter(scaled_LFC > 0) %>% nrow()))

# Spread transformed LFCs
lfc.no_positives.wide <- lfc.narrow.no_positives %>%
  select(all_of(annotation_colnames), sample_label, scaled_LFC) %>%
  spread(sample_label, scaled_LFC)

# Check that the new matrix is the correct length
if (nrow(lfc_matrix) != nrow(lfc.no_positives.wide)) {
  print(paste('Number of rows before setting positive LFCs to zero:', nrow(lfc_matrix)))
  print(paste('Number of rows after setting positive LFCs to zero:', nrow(lfc.no_positives.wide)))
  stop("Incorrect number of rows  after setting positive LFCs to zero.")
}

############################################################
# RDS outputs                                              #
############################################################

# Write LFC medians to RDS
lfc_medians.rds.path <- file.path(repo_path, 'DATA', 'RDS', 'preprocessing', batch, paste0('lfc_medians.unscaled.', rds_ext))
saveRDS(lfc.narrow.medians, file = lfc_medians.rds.path)
print(paste("Unscaled LFC median RDS written to:", lfc_medians.rds.path))

# Write scaled LFCs to RDS
scaled_lfc.rds.path <- file.path(repo_path, 'DATA', 'RDS', 'preprocessing', batch, paste0('lfc_matrix.scaled.', rds_ext))
saveRDS(lfc.narrow.scaled, file = scaled_lfc.rds.path)
print(paste("Scaled LFC RDS written to:", scaled_lfc.rds.path))

# Write scaled and positives to zero LFCs to RDS
scaled_pos_zero_lfc.rds.path <- file.path(repo_path, 'DATA', 'RDS', 'preprocessing', batch, paste0('lfc_matrix.scaled_pos_zero.', rds_ext))
saveRDS(lfc.narrow.no_positives, file = scaled_pos_zero_lfc.rds.path)
print(paste("Scaled LFC with positives to zero RDS written to:", scaled_pos_zero_lfc.rds.path))

############################################################
# TSV outputs                                              #
############################################################

# Write scaled LFCs to TSV
scaled_lfc.path <- file.path(repo_path, 'DATA', 'preprocessing', batch, paste0('lfc_matrix.scaled.', tsv_ext))
write.table(lfc.scaled.wide, scaled_lfc.path, sep = "\t", quote = F, row.names = F)
print(paste("Scaled LFC TSV written to:", scaled_lfc.path))

# Write scaled and positives to zero LFCs to TSV
scaled_pos_zero_lfc.path <- file.path(repo_path, 'DATA', 'preprocessing', batch, paste0('lfc_matrix.scaled_pos_zero.', tsv_ext))
write.table(lfc.no_positives.wide, scaled_pos_zero_lfc.path, sep = "\t", quote = F, row.names = F)
print(paste("Scaled LFC with positives to zero TSV written to:", scaled_pos_zero_lfc.path))

print('Done.')
