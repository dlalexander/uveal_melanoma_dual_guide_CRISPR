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
              help = "full path to sample mapping", metavar = "character")
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
# Batch Splitting                                          #
############################################################


#Splitting sample annotations file so can run batches separately

print(paste("Splitting sample annotations file by batch:", opt$mapping))

sample_mapping_batchA <- sample_mapping %>% filter(batch=="A")
sample_mapping_batchB <- sample_mapping %>% filter(batch=="B")


#Split filtered counts files into batches, within batch specific folders so these can be picked up by downstream steps
print("Splitting normalised counts matrix by batch")

#Keeping all labelling columns, and sample columns that are in the different batches
#Simplest to drop the columns that are in the other batch
count_matrix_batchA <- count_matrix %>% select(-c(sample_mapping_batchB$sample_label))
count_matrix_batchB <- count_matrix %>% select(-c(sample_mapping_batchA$sample_label))

############################################################
# Outputs                                                  #
############################################################

#Create directories to store split results so can run pre-processing on different directories
batchA_path <- paste(repo_path,"/DATA/preprocessing","/batchA",sep="")
batchB_path <- paste(repo_path,"/DATA/preprocessing","/batchB",sep="")

if (!dir.exists(batchA_path)) {
  dir.create(batchA_path)
}

if (!dir.exists(batchB_path)) {
  dir.create(batchB_path)
}

batchA_RDS_path <- paste(repo_path,"/DATA/RDS/preprocessing","/batchA",sep="")
batchB_RDS_path <- paste(repo_path,"/DATA/RDS/preprocessing","/batchB",sep="")

if (!dir.exists(batchA_RDS_path)) {
  dir.create(batchA_RDS_path)
}

if (!dir.exists(batchB_RDS_path)) {
  dir.create(batchB_RDS_path)
}


# Write split annotations to files
sample_mapping_batchA.path <- file.path(repo_path, 'METADATA', '5429_sample_annotations_batchA.tsv')
write.table(sample_mapping_batchA, sample_mapping_batchA.path, sep = "\t", quote = F, row.names = F)

sample_mapping_batchB.path <- file.path(repo_path, 'METADATA', '5429_sample_annotations_batchB.tsv')
write.table(sample_mapping_batchB, sample_mapping_batchB.path, sep = "\t", quote = F, row.names = F)

print(paste("Split annotation files written to:", sample_mapping_batchA.path, "and:",sample_mapping_batchB.path))

# Write split counts matrices to files
count_matrix_batchA.path <- file.path(batchA_path,'count_matrix.norm.batchA.tsv')
write.table(count_matrix_batchA, count_matrix_batchA.path, sep = "\t", quote = F, row.names = F)

count_matrix_batchB.path <- file.path(batchB_path,'count_matrix.norm.batchB.tsv')
write.table(count_matrix_batchB, count_matrix_batchB.path, sep = "\t", quote = F, row.names = F)

print(paste("Split count matrices written to:", count_matrix_batchA.path,"and ", count_matrix_batchB.path))

# Save as RDS
count_matrix.batchA.rds.path <- file.path(batchA_RDS_path,'count_matrix.norm.batchA.rds')
saveRDS(count_matrix_batchA, file = count_matrix.batchA.rds.path )

count_matrix.batchB.rds.path <- file.path(batchB_RDS_path,'count_matrix.norm.batchB.rds')
saveRDS(count_matrix_batchB, file = count_matrix.batchB.rds.path )

print(paste("Split count matrices RDS files written to:", count_matrix.batchA.rds.path,"and ", count_matrix.batchB.rds.path))

print("DONE.")

###