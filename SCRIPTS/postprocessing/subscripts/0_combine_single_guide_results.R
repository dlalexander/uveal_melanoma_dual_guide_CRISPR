suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "repository directory path", metavar = "character"),
  make_option(c("-m", "--mapping"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
  make_option(c("-a", "--annotation"), type = "character",
              help = "full path to annotation", metavar = "character")
);

opt_parser <- OptionParser( option_list = option_list );
opt <- parse_args( opt_parser );

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

# Set datasets
datasets <- c('unscaled', 'scaled', 'scaled_pos_zero')

############################################################
# Sample mapping                                           #
############################################################

check_file_exists(opt$mapping)

print(paste("Reading sample annotations from:", opt$mapping))

# Read in sample mapping
sample_mapping <- read.delim(opt$mapping, header = T, sep = "\t")

############################################################
# Expanded library                                         #
############################################################

check_file_exists(opt$annotation)

print(paste("Reading expanded library:", opt$annotation))

# Read in library
expanded_library <- read.delim(opt$annotation, sep = "\t", header = T)

############################################################
# MAGeCK results                                           #
############################################################

# Create empty results list
mageck_results_list <- list()

# Loop over datasets
for (dataset in datasets) {
  # Get a list of MAGeCK gene level results
  mageck_result_files <- list.files(path = file.path(opt$dir, 'DATA', 'single_guide', 'combined', dataset), pattern = 'MAGeCK.treatment_vs_control.gene_summary.txt', full.names = F, recursive = T)
  # Create empty data frame and loop over datasets
  mageck_results_list[[dataset]] <- data.frame()
  for (result_file in mageck_result_files) {
    # Get cell line label
    cell_line <- str_split(result_file, pattern = "/", simplify = T)[2]
    cell_line_label <- sample_mapping %>% filter(stripped_cell_line_name == cell_line) %>% pull(cell_line_label) %>% unique()
    # Read in MAGeCK gene result file
    mageck_tmp <- read.delim(file.path(opt$dir, 'DATA', 'single_guide', 'combined', dataset, result_file))
    # Add cell line label and rename id
    mageck_tmp <- mageck_tmp %>%
      mutate('dataset' = dataset,
             'cell_line_label' = cell_line_label,
             .before = 'id') %>%
      rename('gene' = 'id')
    # Add binary for significantly depleted and enriched genes
    mageck_tmp <- mageck_tmp %>%
      mutate('is_depleted_mageck' = ifelse(neg.fdr < 0.05, 1, 0),
             'is_enriched_mageck' = ifelse(pos.fdr < 0.05, 1, 0))
    # Add to existing dataset results
    if (nrow(mageck_results_list[[dataset]]) == 0) {
      mageck_results_list[[dataset]] <- mageck_tmp
    } else {
      mageck_results_list[[dataset]] <- bind_rows(mageck_results_list[[dataset]], mageck_tmp)
    }
  }
}

############################################################
# BAGEL results                                            #
############################################################

# Create empty results list
bagel_results_list <- list()

# Loop over datasets
for (dataset in datasets) {
  # Get a list of BAGEL gene level results
  bagel_binary_result_files <- list.files(path = file.path(opt$dir, 'DATA', 'single_guide', 'combined', dataset), pattern = 'BF.scaled_depletions_matrix.gene.BF.treatment_vs_control.tsv', full.names = F, recursive = T)
  bagel_bf_result_files <- list.files(path = file.path(opt$dir, 'DATA', 'single_guide', 'combined', dataset), pattern = 'BF.scaled.gene.BF.treatment_vs_control.tsv', full.names = F, recursive = T)
  # Create empty data frame and loop over datasets
  bagel_results_list[[dataset]] <- data.frame()
  for (i in 1:length(bagel_binary_result_files)) {
    # Get file paths
    bf_file <- bagel_bf_result_files[i]
    binary_file <- bagel_binary_result_files[i]
    # Get cell line label
    cell_line <- str_split(bf_file, pattern = "/", simplify = T)[2]
    cell_line_label <- sample_mapping %>% filter(stripped_cell_line_name == cell_line) %>% pull(cell_line_label) %>% unique()
    # Read in BAGEL gene result files
    bagel_bf_tmp <- read.delim(file.path(opt$dir, 'DATA', 'single_guide', 'combined', dataset, bf_file))
    bagel_binary_tmp <- read.delim(file.path(opt$dir, 'DATA', 'single_guide', 'combined', dataset, binary_file))
    # Combine results files
    bagel_tmp <- bagel_bf_tmp %>%
      left_join(bagel_binary_tmp, by = 'gene')
    # Add cell line and dataset
    bagel_tmp <- bagel_tmp %>%
      mutate('dataset' = dataset,
             'cell_line_label' = cell_line_label,
             .before = 'gene') 
    # Rename is_depleted
    bagel_tmp <- bagel_tmp %>%
      rename('is_depleted_bagel' = 'is_depleted')
    # Add to existing dataset results
    if (nrow(bagel_results_list[[dataset]]) == 0) {
      bagel_results_list[[dataset]] <- bagel_tmp
    } else {
      bagel_results_list[[dataset]] <- bind_rows(bagel_results_list[[dataset]], bagel_tmp)
    }
  }
}

############################################################
# Combine MAGeCK and BAGEL results                         #
############################################################

# Create empty results list
singles_results_list <- list()

# Loop over datasets and combine MAGeCK and BAGEL results
for (dataset in datasets) {
  singles_results_list[[dataset]] <- mageck_results_list[[dataset]] %>%
    left_join(bagel_results_list[[dataset]], by = c('cell_line_label', 'dataset', 'gene'))
}

############################################################
# Expand dataset by library metadata and cell line         #
############################################################

# Create empty results list
binary_results_list <- list()

# Loop over dataset and expand results across library
for (dataset in datasets) {
  print(paste("Expanding dataset:", dataset))

  print('Expanding library by cell line...')

  binary_results_list[[dataset]] <- expand.grid(expanded_library$sorted_gene_pair, sample_mapping$cell_line_label) %>%
    unique() %>%
    select('sorted_gene_pair' = 'Var1', 'cell_line_label' = 'Var2' ) %>%
    filter(cell_line_label %in% singles_results_list[[dataset]]$cell_line_label)

  print('Adding additional library metadata...')

  binary_results_list[[dataset]] <- binary_results_list[[dataset]] %>%
    left_join(expanded_library %>% select(gene_pair_id, sorted_gene_pair, targetA, targetB), by = 'sorted_gene_pair') %>%
    unique()

  print('Adding single essentiality for targetA...')

  binary_results_list[[dataset]] <- binary_results_list[[dataset]] %>%
    left_join(singles_results_list[[dataset]] %>% 
                select(dataset, cell_line_label, gene, is_depleted_mageck, is_enriched_mageck, is_depleted_bagel) %>%
                rename_at(vars(-gene, -cell_line_label), ~ paste('targetA', ., sep = '__')), by = c('targetA' = 'gene', 'cell_line_label')) %>%
    unique()

  print('Adding single essentiality for targetB...')

  binary_results_list[[dataset]] <- binary_results_list[[dataset]] %>%
    left_join(singles_results_list[[dataset]] %>% 
                select(dataset, cell_line_label, gene, is_depleted_mageck, is_enriched_mageck, is_depleted_bagel) %>%
                rename_at(vars(-gene, -cell_line_label), ~ paste('targetB', ., sep = '__')), by = c('targetB' = 'gene', 'cell_line_label')) %>%
    unique()
  
  print('Adding sample mapping...')

  binary_results_list[[dataset]] <- binary_results_list[[dataset]] %>%
    left_join(sample_mapping %>% select(cell_line_label, cancer_type), by = 'cell_line_label') %>%
    unique()

  print('Postprocessing...')

  binary_results_list[[dataset]] <- binary_results_list[[dataset]] %>%
    select(cell_line_label, cancer_type, 
          gene_pair_id, sorted_gene_pair, 
          targetA, targetB, 
          contains('is_depleted_mageck'), contains('is_enriched_mageck'), contains('is_depleted_bagel'))
}

############################################################
# Save MAGeCK and BAGEL combined results                   #
############################################################

# Loop over datasets and write combined MAGeCK and BAGEL results
for (dataset in datasets) {
  # Save singles results as TSV
  singles_results.path <- file.path(repo_path, 'DATA', 'postprocessing', 'intermediate_tables', paste0('combined_singles_results.', dataset, '.tsv'))
  write.table(singles_results_list[[dataset]], singles_results.path, sep = "\t", quote = F, row.names = F)
  print(paste("Combined MAGeCK and BAGEL results TSV written to:", singles_results.path))

  # Save binary results as TSV
  binary_results.path <- file.path(repo_path, 'DATA', 'postprocessing', 'intermediate_tables', paste0('combined_singles_results.binary.', dataset, '.tsv'))
  write.table(binary_results_list[[dataset]], binary_results.path, sep = "\t", quote = F, row.names = F)
  print(paste("Combined MAGeCK and BAGEL binary results TSV written to:", binary_results.path))
  
  # Save as RDS
  singles_results.rds.path <- file.path(repo_path, 'DATA', 'RDS', 'postprocessing', paste0('combined_singles_results.full.', dataset, '.rds'))
  saveRDS(singles_results_list[[dataset]], file = singles_results.rds.path)
  print(paste("Combined MAGeCK and BAGEL results RDS written to:", singles_results.rds.path))

  # Save binary results as RDS
  binary_results.rds.path <- file.path(repo_path, 'DATA', 'RDS', 'postprocessing', paste0('combined_singles_results.binary.', dataset, '.rds'))
  saveRDS(binary_results_list[[dataset]], file = binary_results.rds.path)
  print(paste("Combined MAGeCK and BAGEL binary results RDSwritten to:", binary_results.path))
} 

print('Done.')