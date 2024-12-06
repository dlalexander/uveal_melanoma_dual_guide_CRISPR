suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(stringi)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character"),
  make_option(c("-m", "--mapping"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
  make_option(c("-a", "--annotation"), type = "character",
              help = "full path to annotation", metavar = "character")
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
# Sample mapping                                           #
############################################################

check_file_exists(opt$mapping)

print(paste("Reading sample annotations from:", opt$mapping))

# Read in sample mapping
sample_mapping <- read.delim(opt$mapping, header = T, sep = "\t")

############################################################
# pyCROQUET library                                        #
############################################################

check_file_exists(opt$annotation)

print(paste("Reading pyCROQUET library:", opt$annotation))

# Read in library
# Want to skip all 4 header lines - including the newly specified dual-orientation
pycroquet_library <- read.delim(opt$annotation, sep = "\t", skip = 4, header = F,
                                col.names = c('id', 'sgrna_ids', 'sgrna_seqs', 'gene_pair_id', 'sgrna_symbols', 'sgrna_libraries', 'sgrna_group'))

############################################################
# Expand library                                           #
############################################################

print('Expanding library...')

expanded_library <- pycroquet_library %>% 
  separate(gene_pair_id, sep = "\\|", into = c('targetA', 'targetB'), remove = F) %>% # Split the genes targeted in each pair into separate columns
  rowwise() %>%
  mutate('sorted_gene_pair' = str_c(str_sort(str_split(gene_pair_id, "\\|", simplify = T), numeric = T), collapse = "|"), .after = 'gene_pair_id') %>% # Sort gene pair alphabetically
  separate(sgrna_symbols, sep = "\\|", into = c('sgrna_symbol_A', 'sgrna_symbol_B'), remove = F) %>% # Split sgrna targets
  mutate('guide_type' = case_when(grepl("^F[0-9]+$", sgrna_symbol_A) & grepl("^F[0-9]+$", sgrna_symbol_B) ~ "safe_control|safe_control", # Add annotation for library type
                                !grepl("^F[0-9]+$", sgrna_symbol_A) & grepl("^F[0-9]+$", sgrna_symbol_B) ~ "gene|safe_control",
                                grepl("^F[0-9]+$", sgrna_symbol_A) & !grepl("^F[0-9]+$", sgrna_symbol_B) ~ "safe_control|gene",
                                !grepl("^F[0-9]+$", sgrna_symbol_A) & !grepl("^F[0-9]+$", sgrna_symbol_B) ~ "gene|gene",
                                TRUE ~ "UNKNOWN")) %>%
  mutate('singles_target_gene' = case_when(guide_type == 'gene|safe_control' ~ sgrna_symbol_A, # Add annotation for gene targeted in single pairs
                                         guide_type == 'safe_control|gene' ~ sgrna_symbol_B)) %>%
  mutate('guide_orientation' = case_when(guide_type == 'gene|gene' & sgrna_symbol_A == targetA & sgrna_symbol_B == targetB ~ 'geneA|geneB', # Add annotation for orientation
                                         guide_type == 'gene|gene' & sgrna_symbol_A == targetB & sgrna_symbol_B == targetA ~ 'geneB|geneA',
                                         guide_type == 'gene|safe_control' & sgrna_symbol_A == targetA ~ 'geneA|safe_control',
                                         guide_type == 'gene|safe_control' & sgrna_symbol_A == targetB ~ 'geneB|safe_control',
                                         guide_type == 'safe_control|gene' & sgrna_symbol_B == targetB ~ 'safe_control|geneB',
                                         guide_type == 'safe_control|gene' & sgrna_symbol_B == targetA ~ 'safe_control|geneA',
                                         guide_type == 'safe_control|safe_control' ~ 'safe_control|safe_control'), .after = 'guide_type')

############################################################
# Build dual guide matrix                                  #
############################################################

# For each dual (gene|gene) guide we need to know the corresponding single guides (e.g. safe|gene and gene|safe)
# Function is in helper.R
dual_guide_matrix <- get_guide_matrix(pycroquet_library)

############################################################
# Expanded library output                                  #
############################################################

# Write expanded library to output file
expanded_library.path <- file.path(repo_path, 'METADATA', 'libraries', 'expanded_library.tsv')
write.table(expanded_library, expanded_library.path, sep = "\t", quote = F, row.names = F)

print(paste("Expanded library written to:", expanded_library.path))

# Save as RDS
expanded_library.rds.path <- file.path(repo_path, 'DATA', 'RDS', 'preprocessing', 'expanded_library.rds')
saveRDS(expanded_library, file = expanded_library.rds.path )

print(paste("Expanded library RDS written to:", expanded_library.rds.path))

############################################################
# Guide matrix output                                      #
############################################################

# Write guide matrix to output file
dual_guide_matrix.path <- file.path(repo_path, 'DATA', 'dual_guide', 'dual_guide_matrix.tsv')
write.table(dual_guide_matrix, dual_guide_matrix.path, sep = "\t", quote = F, row.names = F)

print(paste("Guide matrix written to:", dual_guide_matrix.path))

# Save as RDS
dual_guide_matrix.rds.path <- file.path(repo_path, 'DATA', 'RDS', 'preprocessing', 'dual_guide_matrix.rds')
saveRDS(dual_guide_matrix, file = dual_guide_matrix.rds.path )

print(paste("Guide matrix RDS written to:", dual_guide_matrix.rds.path))

print('Done.')