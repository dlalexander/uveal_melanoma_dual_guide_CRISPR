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
# Download DepMap relative copy number                     #
############################################################

print('Downloading DepMap relative copy number...')

# Read in DepMap 22Q2 gene copy number per cell line
# Gene level copy number data, log2 transformed with a pseudo count of 1. 
tmpfile <- tempfile()
options(timeout = 500)
download.file(url = 'https://ndownloader.figshare.com/files/34989937', destfile = tmpfile, quiet = T)

print('Reading DepMap relative copy number...')

depmap_cn <- read.delim(tmpfile, sep = ",", header = T)
colnames(depmap_cn)[1] <- 'depMapID'

############################################################
# Prepare DepMap CN                                        #
############################################################

print('Filtering DepMap copy number...')

# Filter to only include screened cell lines 
# That are present in DepMap
depmap_cn <- depmap_cn %>% 
  filter(depMapID %in% sample_mapping$depMapID)

# Gather and split gene information
depmap_cn <- depmap_cn %>% 
  gather(gene_label, relative_copy_number_log2ps1, -depMapID) %>%
  separate(gene_label, sep = "\\.\\.", into = c('gene', 'entrez_id')) %>%
  mutate(entrez_id = str_replace_all(entrez_id, "\\.", ""))

############################################################
# Expand CN by library metadata and cell line              #
############################################################

print('Expanding library by cell line...')

depmap_cn.ann <- expand.grid(expanded_library$sorted_gene_pair, sample_mapping$depMapID) %>%
  unique() %>%
  select('sorted_gene_pair' = 'Var1', 'depMapID' = 'Var2' ) %>%
  filter(!is.na(depMapID))

print('Adding additional library metadata...')

depmap_cn.ann <- depmap_cn.ann %>%
  left_join(expanded_library %>% select(gene_pair_id, sorted_gene_pair, targetA, targetB), by = 'sorted_gene_pair') %>%
  unique()
  
print('Adding copy number for targetA...')

depmap_cn.ann <- depmap_cn.ann %>%
  left_join(depmap_cn %>% 
              select(-entrez_id) %>%
              rename_at(vars(-gene, -depMapID), ~ paste('targetA', ., sep = '__')), by = c('targetA' = 'gene', 'depMapID')) %>%
  unique()

print('Adding copy number for targetB...')

depmap_cn.ann <- depmap_cn.ann %>%
  left_join(depmap_cn %>% 
              select(-entrez_id) %>%
              rename_at(vars(-gene, -depMapID), ~ paste('targetB', ., sep = '__')), by = c('targetB' = 'gene', 'depMapID')) %>%
  unique()

print('Add cell line metadata to DepMap copy number...')

depmap_cn.ann <- depmap_cn.ann %>%
  left_join(sample_mapping %>% select(depMapID, cell_line_label, cancer_type), by = 'depMapID') %>%
  unique() %>%
  select(depMapID, cell_line_label, cancer_type, 
         gene_pair_id, sorted_gene_pair, 
         targetA, targetB, 
         targetA__relative_copy_number_log2ps1, 
         targetB__relative_copy_number_log2ps1)

############################################################
# Write outputs to file                                    #
############################################################

print('Writing to file...')

# Save as a TSV
cn.filepath <- file.path(opt$dir, 'DATA', 'postprocessing', 'intermediate_tables', 'depmap_relative_copy_number.tsv')
write.table(depmap_cn.ann, file = cn.filepath, sep = "\t", row.names = F, quote = F)
print(paste0('DepMap relative copy number TSV written to:', cn.filepath))

# Save as an Rdata object
cn.rds.filepath <- file.path(opt$dir, 'DATA', 'RDS', 'postprocessing', 'depmap_relative_copy_number.tsv')
saveRDS(depmap_cn.ann, file = cn.rds.filepath)
print(paste0('DepMap relative copy number RDS written to:', cn.rds.filepath))

print('Done.')