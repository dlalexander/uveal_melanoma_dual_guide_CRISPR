suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(optparse)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "repository directory path", metavar = "character"),
  make_option(c("-g", "--guide"), type = "character",
              help = "full path to file with list of guide ids to remove", metavar = "character"),
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
# User-defined guides                                      #
############################################################

# Check file with guide ids exists
check_file_exists(opt$guide)

# Read in guide ids to remove
guides_to_remove <- scan(opt$guide, what = character(), strip.white = T, quiet = T)

# Check we have one or more ids
if (is.null(guides_to_remove) || length(guides_to_remove) == 0) {
  print(paste('Error reading guides to remove:', opt$guide))
  stop()
}

############################################################
# Expanded library                                         #
############################################################

check_file_exists(opt$annotation)

print(paste("Reading expanded library:", opt$annotation))

# Read in library
expanded_library <- read.delim(opt$annotation, sep = "\t", header = T)

############################################################
# Remove user-defined guides from library                  #
############################################################

print("Removing user defined guides from library...")

lib.guides_removed <- expanded_library %>%
  filter(!id %in% guides_to_remove)

print(paste('Number of guides removed:', length(guides_to_remove)))

# Check that the new matrix is the correct length
if (nrow(expanded_library) - nrow(lib.guides_removed) != length(guides_to_remove)) {
  print(paste('Number of rows pre-removal:', nrow(expanded_library)))
  print(paste('Number of rows post-removal:', nrow(lib.guides_removed)))
  stop("Incorrect number of rows after removing guides.")
}

############################################################
# Extract single guides from library                       #
############################################################

print("Extracting single guides...")

# Get guide pairs containing safe-targeting guides by guide type
singles_guide_type <- c('gene|safe_control', 'safe_control|gene', 'safe_control|safe_control')
single_guides <- list()
single_guides[['combined']] <- lib.guides_removed %>%
  filter(guide_type %in% singles_guide_type) 

print("Adding gene annotation for singles...")

# Add a gene annotation for safe controls
single_guides[['combined']] <- single_guides[['combined']] %>%
  mutate('gene' = ifelse(guide_type == "safe_control|safe_control", 'safe_control', singles_target_gene)) 

print(paste('Number of single guides:', nrow(single_guides[['combined']])))

############################################################
# Extract guides which are gene|safe_control               #
############################################################

single_guides[['A']] <- single_guides[['combined']] %>%
  filter(!guide_type == 'safe_control|gene') 

############################################################
# Extract guides which are safe_control|gene               #
############################################################

single_guides[['B']] <- single_guides[['combined']] %>%
  filter(guide_type == 'gene|safe_control') 

############################################################
# Convert library for MAGeCK RRA                           #
############################################################

print("Converting singles library for use with MAGeCK...")

# NB: Could switch to just doing combined 
datasets <- c('combined', 'A', 'B')

# Select required columns
single_guides_mageck <- list()
for (d in datasets) {
  single_guides_mageck[[d]] <- single_guides[[d]] %>% 
    select(sgrna_ids, sgrna_seqs, gene) %>% 
    unique()
}

############################################################
# Singles library output                                   #
############################################################

# Write singles libraries to output file
for (d in datasets) {
  singles_mageck.path <- file.path(repo_path, 'DATA', 'single_guide', d, paste0('singles_library.', d, '.tsv'))
  write.table(single_guides_mageck[[d]], singles_mageck.path, sep = "\t", quote = F, row.names = F)
  print(paste("Singles library (", d,  ") written to:", singles_mageck.path))
}

# Write singles libraries to RDS
for (d in datasets) {
  singles_mageck.rds.path <- file.path(repo_path, 'DATA', 'RDS', 'single_guide', paste0('singles_library.', d, '.rds'))
  saveRDS(single_guides_mageck[[d]], file = singles_mageck.rds.path)
  print(paste("Singles library (", d,  ") RDS written to:", singles_mageck.rds.path))
}

print('Done.')