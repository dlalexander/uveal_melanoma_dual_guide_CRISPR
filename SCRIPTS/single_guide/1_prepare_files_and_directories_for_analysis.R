suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(optparse)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "repository directory path", metavar = "character"),
  make_option(c("-a", "--annotation"), type = "character",
              help = "relative path to annotation", metavar = "character"),
  make_option(c("-m", "--mapping"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
  make_option(c("-c", "--counts"), type = "character",
              help = "full path to count matrix", metavar = "character"),
  make_option(c("-e", "--exclude"), type = "character",
              help = "cell line labels to exclude", metavar = "character"),
  make_option(c("-l", "--log"), type = "character",
              help = "relative path for log files", metavar = "character"),
  make_option(c("-s", "--suffix"), type = "character",
              help = "suffix", metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character",
              help = "relative path to output directory", metavar = "character"),
  make_option(c("-b", "--batch"), type = "character",
              help = "batch", metavar = "character")
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
# Singles library                                          #
############################################################

check_file_exists(opt$annotation)

print(paste("Reading singles library:", opt$annotation))

# Read in library
singles_library <- read.delim(opt$annotation, sep = "\t", header = T)

############################################################
# Count matrix                                             #
############################################################

check_file_exists(opt$counts)

print(paste("Reading count matrix:", opt$counts))

# Read in count matrix
count_matrix <- read.delim(opt$counts, sep = "\t", header = T, check.names = F)

# Set annotation column names
annotation_colnames <- c(
  'id', 'sgrna_ids', 'sgrna_seqs', 'gene_pair_id',
  'sorted_gene_pair', 'targetA', 'targetB', 
  'sgrna_symbols', 'sgrna_symbol_A', 'sgrna_symbol_B', 
  'sgrna_libraries', 'sgrna_group', 'guide_type', 
  'guide_orientation', 'singles_target_gene')

############################################################
# Sample mapping                                           #
############################################################

check_file_exists(opt$mapping)

print(paste("Reading sample annotations from:", opt$mapping))

# Read in sample mapping
sample_mapping <- read.delim(opt$mapping, header = T, sep = "\t")

############################################################
# Filter sample mapping                                    #
############################################################

print('Filtering sample mapping...')

# Only keep samples which are in the count matrix
sample_mapping.filt <- sample_mapping %>%
  filter(sample_label %in% colnames(count_matrix))

# Remove sample labels to exclude
if (!is.null(opt$exclude)) {
  cell_lines_to_exclude <- as.vector(str_split(opt$exclude, ",", simplify = T))
  if (0 == length(cell_lines_to_exclude) || grepl(',', cell_lines_to_exclude[1])) {
    stop(paste("Problem splitting cell line labels:", opt$exclude))
  }
  sample_mapping.filt <- sample_mapping.filt %>%
    filter(!qc_group %in% cell_lines_to_exclude)
}

############################################################
# Prepare count matrix                                     #
############################################################

print('Filtering count matrix...')

# Only keep guides which are in the singles library
count_matrix.filt <- count_matrix %>%
  filter(sgrna_ids %in% singles_library$sgrna_ids)

print('Adding gene to count matrix...')

# Add gene annotation 
count_matrix.filt <- count_matrix.filt  %>%
  left_join(singles_library, by = c('sgrna_ids', 'sgrna_seqs'))

print('Getting unique counts...')

# Get unique singles
count_matrix.filt <- count_matrix.filt  %>%
  select(sgrna_ids, gene, control_mean, all_of(sample_mapping.filt$sample_label)) %>%
  unique() 

print('Replacing sample labels...')

# Replace sample labels with stripped sample labels in sample mapping
count_matrix.filt <- count_matrix.filt  %>%
  gather(sample_label, count, -sgrna_ids, -gene, -control_mean) %>%
  left_join(sample_mapping.filt %>% select(sample_label, stripped_sample_label), by = 'sample_label') %>%
  select(-sample_label) %>%
  spread(stripped_sample_label, count)

############################################################
# Filter library                                           #
############################################################

print('Filtering library...')

singles_library.filt <- singles_library %>%
  filter(sgrna_ids %in% count_matrix.filt$sgrna_ids)

############################################################
# Get cell line directory names                            #
############################################################

print('Getting directory names...')

# Get stripped cell line directory names
cell_lines <- unique(sample_mapping.filt$stripped_cell_line_name)

############################################################
# Set up analysis list                                     #
############################################################

print('Preparing analysis list...')

analysis_list <- list()
for (cl in cell_lines) {
  analysis_list[[cl]][['directory']] <- file.path(opt$outdir, cl)
  analysis_list[[cl]][['manifest']]  <- file.path(opt$outdir, cl, 'sample_manifest.tsv')
  analysis_list[[cl]][['counts']]    <- file.path(opt$outdir, cl, 'sample_counts.tsv')
  analysis_list[[cl]][['library']]   <- str_replace(opt$annotation, ".tsv", ".filt.tsv")
  analysis_list[[cl]][['lsf']][['job_name']] <- paste0("csar_", cl)
  analysis_list[[cl]][['lsf']][['queue']]    <- 'normal'
  analysis_list[[cl]][['lsf']][['memory']]   <- '-R "select[mem>4000] rusage[mem=4000] span[hosts=1]" -M 4000'
  analysis_list[[cl]][['lsf']][['cpu']]      <- '-n 6'
  analysis_list[[cl]][['lsf']][['outfile']]  <- file.path("${REPO_PATH}", opt$log, paste0('csar_', cl, '.o'))
  analysis_list[[cl]][['lsf']][['errfile']]  <- file.path("${REPO_PATH}", opt$log, paste0('csar_', cl, '.e'))
}

############################################################
# Create directories                                       #
############################################################

print('Creating cell line directories...')

for (cl in names(analysis_list)) {
  create_directory(analysis_list[[cl]][['directory']])
  check_dir_exists(analysis_list[[cl]][['directory']])
}

############################################################
# Generating libraries                                     #
############################################################

print('Writing filtered singles library...')

write.table(singles_library.filt, analysis_list[[cl]][['library']], row.names = F, sep = "\t", quote = F)
check_file_exists(analysis_list[[cl]][['library']])

############################################################
# Generating sample manifests                              #
############################################################

print('Creating cell line sample manifests...')

#Suppressing warnings because funs is deprecated (but still works)
suppressWarnings({
for (cl in names(analysis_list)) {
  # Get subset of manifest in C-SAR format and add control mean
  sample_manifest.subset <- sample_mapping.filt %>%
    filter(stripped_cell_line_name == cl) %>%
    mutate('filename' = basename(analysis_list[[cl]][['counts']]), 'plasmid' = 0, 'control' = 0, 'treatment' = 1) %>%
    select(filename, stripped_sample_label, plasmid, control, treatment, qc_group, replicate, sample_label) %>%
    rbind(c(basename(analysis_list[[cl]][['counts']]), 'control_mean', 0, 1, 0, 'control', 'R1', 'control_mean'))
    if(cl == "921"){
      sample_manifest.subset <- sample_manifest.subset %>%
      mutate_all(funs(str_replace(., "921", "X921")))
      # Replaced 921 with X921 in the sample manifest to avoid error due to cell line name starting with a number
    }
  # Write subset manifest to file
  write.table(sample_manifest.subset, file.path(opt$dir, analysis_list[[cl]][['manifest']]), sep = "\t", quote = F, row.names = F)
  check_file_exists(analysis_list[[cl]][['manifest']])
}
})
  
############################################################
# Generating count matrices                                #
############################################################

print('Creating cell line count matrices...')

for (cl in names(analysis_list)) {
  # Get sample labels
  samples_in_cell_line <- sample_mapping.filt %>%
    filter(stripped_cell_line_name == cl) %>%
    pull(stripped_sample_label)
  # Get subset of count matrix
  count_matrix.filt.subset <- count_matrix.filt %>%
    select(sgrna_ids, gene, control_mean, all_of(samples_in_cell_line))
  # Set sample count indices
  analysis_list[[cl]][['count_count_column_index']] <- paste(3:ncol(count_matrix.filt.subset), collapse = ',')
  # Write subset count matrix to file
  write.table(count_matrix.filt.subset, file.path(opt$dir, analysis_list[[cl]][['counts']]), sep = "\t", quote = F, row.names = F)
  check_file_exists(analysis_list[[cl]][['counts']])
}

############################################################
# Generating bsub commands                                 #
############################################################

print('Generating bsub commands...')

bsub_command_list <- vector()
for (cl in names(analysis_list)) {
  csar_command <- paste0( 'c-sar /opt/wsi-t113/c-sar/main.nf -c "${REPO_PATH}/METADATA/single_guide/csar.config" ',
                          '--counts "${REPO_PATH}/', analysis_list[[cl]][['counts']], '" ',
                          '--count_count_column_index "', analysis_list[[cl]][['count_count_column_index']], '" ',
                          '--library "${REPO_PATH}/', analysis_list[[cl]][['library']], '" ',
                          '--info "${REPO_PATH}/', analysis_list[[cl]][['manifest']], '" ',
                          '--outdir "${REPO_PATH}/', analysis_list[[cl]][['directory']], '"')

  bsub_command <- paste(paste0('cd "${REPO_PATH}/', analysis_list[[cl]][['directory']], '"; '),
                        'bsub',
                        '-q', analysis_list[[cl]][['lsf']][['queue']],
                        '-J', analysis_list[[cl]][['lsf']][['job_name']],
                        '-o', paste0('"', analysis_list[[cl]][['lsf']][['outfile']], '"'),
                        '-e', paste0('"', analysis_list[[cl]][['lsf']][['errfile']], '"'),
                        analysis_list[[cl]][['lsf']][['memory']],
                        analysis_list[[cl]][['lsf']][['cpu']])
  bsub_command <- paste(bsub_command, csar_command)
  bsub_command_list <- c(bsub_command_list, bsub_command)
}

############################################################
# Write bsub commands                                      #
############################################################

print('Writing bsub commands...')

filepath <- file.path(opt$dir, 'SCRIPTS', 'single_guide', paste0('bsub_commands.', opt$suffix, '.sh'))
fileConn <- file(filepath)
write_lines(bsub_command_list, fileConn, append = F)

print('Done.')