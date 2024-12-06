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
# Expanded library                                         #
############################################################

check_file_exists(opt$annotation)

print(paste("Reading expanded library:", opt$annotation))

# Read in library
expanded_library <- read.delim(opt$annotation, sep = "\t", header = T)

############################################################
# Sample mapping                                           #
############################################################

check_file_exists(opt$mapping)

print(paste("Reading sample annotations from:", opt$mapping))

# Read in sample mapping
sample_mapping <- read.delim(opt$mapping, header = T, sep = "\t")

############################################################
# pyCROQUET counts per lane                                #
############################################################

# Get list of count file paths
pycroquet_dir <- file.path(repo_path, 'DATA', 'pyCROQUET', 'COUNTS')
pycroquet_count_files <- list.files(path = pycroquet_dir, pattern = '.counts.tsv', full.names = T)

print(paste("Reading pyCROQUET count files from:", pycroquet_dir))
print(paste("Number of count files:", length(pycroquet_count_files)))

# Read in count files. Organise by sample and lane. Takes a while to run.
pycroquet_counts_per_lane <- data.frame()
for (pcf in pycroquet_count_files) {
  rlt <- basename(pcf) %>% str_replace('.counts.tsv', '')
  sn <- read_lines(pcf, n_max=1, skip = 2) %>% substring(58) # Identify the sample name
  tmp_counts <- read.delim(pcf, header = F, sep = "\t", skip = 3) %>% select(id = V1, counts = V6) 
  tmp_counts <- tmp_counts %>% mutate('sample_name' = sn, 'run_info' = rlt)
  if (nrow(pycroquet_counts_per_lane) == 0) {
    pycroquet_counts_per_lane <- tmp_counts
  } else {
    pycroquet_counts_per_lane <- bind_rows(pycroquet_counts_per_lane, tmp_counts)
  }
}

print(paste("Number of samples found:", length(unique(pycroquet_counts_per_lane$sample_name))))

############################################################
# pyCROQUET counts per sample                              #
############################################################

print("Calculating total counts per sample...")

# Get total counts per sample
pycroquet_counts <- pycroquet_counts_per_lane %>%
  group_by(sample_name, id) %>%
  summarise('sample_counts' = sum(counts), .groups = 'keep') %>%
  ungroup()

#Replacing accession numbers with sample names - in some cases accession numbers were present in the counts files
accession_mapping_path <- file.path(repo_path, 'METADATA','accession_sample_mapping.tsv')

check_file_exists(accession_mapping_path)

print(paste("Reading accession to sample mapping from:",accession_mapping_path))

# Read in sample mapping
accession_mapping <- read.delim(accession_mapping_path, header = T, sep = "\t")

#Create a table with the current and mapping
accession_mapping <- accession_mapping %>% filter(accession_number != "NULL")
row_count <- nrow(accession_mapping)

sample_name_or_accession <- pycroquet_counts[,'sample_name']

#Replacing accession with sample name
for (i in 1:row_count) {
  search_keyword <- accession_mapping[i,"accession_number"]
  replace_keyword <- accession_mapping[i,"supplier_name"]
  sample_name_or_accession <- sapply(sample_name_or_accession, function(x){
    x[x==search_keyword] <- replace_keyword
    return(x)
  })
}

pycroquet_counts[,"sample_name"]<-sample_name_or_accession



print("Building count matrix...")

# Spread into sample count matrix
pycroquet_counts.wide <- sample_mapping %>% select(sanger_sample_name, sample_label) %>%
  left_join(pycroquet_counts, by = c('sanger_sample_name' = 'sample_name')) %>%
  select(-sanger_sample_name) %>%
  spread(sample_label, sample_counts, fill = 0)

print("Adding annotations to counts...")

# Add annotations into count matrix
pycroquet_counts.wide.ann <- expanded_library %>%
  right_join(pycroquet_counts.wide, by = c('id'))

############################################################
# Outputs                                                  #
############################################################

# Write raw counts per sample to output file
pycroquet_counts.path <- file.path(repo_path, 'DATA', 'preprocessing', 'count_matrix.tsv')
write.table(pycroquet_counts.wide.ann, pycroquet_counts.path, sep = "\t", quote = F, row.names = F)

print(paste("Unprocessed count matrix written to:", pycroquet_counts.path))

# Save as RDS
pycroquet_counts.rds.path <- file.path(repo_path, 'DATA', 'RDS', 'preprocessing', 'count_matrix.rds')
saveRDS(pycroquet_counts.wide.ann, file = pycroquet_counts.rds.path )

print(paste("Unprocessed count matrix RDS written to:", pycroquet_counts.rds.path))

print('Done.')
