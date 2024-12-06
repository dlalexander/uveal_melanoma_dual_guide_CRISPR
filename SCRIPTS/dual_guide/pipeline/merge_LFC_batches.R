suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(optparse)))

############################################################
# OPTIONS                                                  # 
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "repository directory path", metavar = "character"),
  make_option(c("-a", "--batchA"), type = "character",
              help = "full path to batchA LFC matrix", metavar = "character"),
  make_option(c("-b", "--batchB"), type = "character",
              help = "full path to batchB LFC matrix", metavar = "character"),
  make_option(c("-s", "--suffix"), type = "character",
              help = "suffix", metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character",
              help = "full path to output directory", metavar = "character")
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

# Check output directory exists
check_dir_exists(opt$outdir)


############################################################
# LFC matrix                                               #
############################################################

check_file_exists(opt$batchA)
check_file_exists(opt$batchB)

print(paste("Reading LFC matrix:", opt$batchA))
print(paste("Reading LFC matrix:", opt$batchB))

# Read in LFC matrix
lfc_matrixA <- read.delim(opt$batchA, sep = "\t", header = T, check.names = F)
lfc_matrixB <- read.delim(opt$batchB, sep = "\t", header = T, check.names = F)


# Set annotation column names
annotation_colnames <- c(
  'id', 'sgrna_ids', 'sgrna_seqs', 'gene_pair_id',
  'sorted_gene_pair', 'targetA', 'targetB', 
  'sgrna_symbols', 'sgrna_symbol_A', 'sgrna_symbol_B', 
  'sgrna_libraries', 'sgrna_group', 'guide_type', 
  'guide_orientation', 'singles_target_gene')


# Join LFC matrices 
# Include a check of the dimensions of the original matrices and final matrix
# Therefore requires that the two batches had the same library components

# NB because of the filtering for low count guides there will be some differences in the LFC matrices in their rows
# Keeping only the guide pairs that did not have low counts 

print('Joining LFC matrices from batchA and batchB...')
print('Retaining guide pairs without low counts in both batches ...')

lfc_matrix <- lfc_matrixA  %>%
  inner_join(lfc_matrixB, by = annotation_colnames)


############################################################
# Saving LFC matrix                                    #
############################################################

print('Writing LFC matrix...')

filepath <- file.path(opt$outdir, paste('lfc_matrix.', opt$suffix, '.tsv',sep=""))
write.table(lfc_matrix, filepath, sep = "\t", quote = F, row.names = F)

print('Done.')