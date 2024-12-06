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
  make_option(c("-s", "--samples"), type = "character",
              help = "full path to samples", metavar = "character")
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

print(paste("Reading DepMap sample annotations from:", opt$mapping))

# Read in sample mapping
sample_mapping <- read.delim(opt$mapping, header = T, sep = ",")

############################################################
# Sample IDs                                               #
############################################################

check_file_exists(opt$samples)

print(paste("Reading sample IDs from:", opt$samples))

# Read in sample mapping
sample_information <- read.delim(opt$samples, header = T, sep = "\t")

############################################################
# Annotating cell lines.                                   #
############################################################

print(paste("Annotating cell lines"))

#The start of sample_supplier_name gives the stripped cell line names from DepMap
#But NB: - OMM2.3, MP38, MP41 are not in DepMap so do not have full sample annotations.

sample_annotations <- sample_information %>% 
  rename(sanger_sample_id=sample,sanger_sample_name=sample_supplier_name) %>%
  select(sanger_sample_id,sanger_sample_name) %>%
  mutate(sample_label = str_replace_all(sanger_sample_name, "_", " "))%>%
  mutate(stripped_sample_label = str_replace_all(sample_label, " ", "."))%>%
  mutate(replicate = str_sub(stripped_sample_label,-1,-1))%>%
  rename(replicate_letter=replicate)%>%
  mutate(replicate=case_when(replicate_letter=="A"~"R1",
                             replicate_letter=="B"~"R2",
                             replicate_letter=="C"~"R3"))%>%
  #QC group - the sample label without replicate, if it is a control want to specify this
  mutate(qc_group=str_sub(sample_label,end=-2))%>%
  mutate(time_point_days=as.character(lapply(strsplit(as.character(qc_group), split="D"),tail, n=1)))%>%
  mutate(cell_line_label=as.character(lapply(strsplit(as.character(sample_label), split=" "),head, n=1)))%>%
  mutate(cell_line_name=cell_line_label)%>%
  mutate(stripped_cell_line_name=toupper(cell_line_label))%>%
  #Since not all in DepMap req. manually annotating cancer type for this study
  mutate(batch=case_when(stripped_cell_line_name=="921"|stripped_cell_line_name=="MP41"|stripped_cell_line_name=="MP38"|stripped_cell_line_name=="OMM1"|stripped_cell_line_name=="OMM23"~"B",
                         stripped_cell_line_name=="MEL202"|stripped_cell_line_name=="MEL270"|stripped_cell_line_name=="MEL285"|stripped_cell_line_name=="OMM25"|stripped_cell_line_name=="MP46"~"A"))%>%
  mutate(batch_control=case_when(batch=="A"& grepl("Control",sample_label)~"control_A",
                                 batch=="B"& grepl("Control",sample_label)~"control_B"))%>%
  mutate(cancer_type="Uveal melanoma")

sample_annotations <- merge(sample_annotations,sample_mapping[ , c("stripped_cell_line_name","CCLE_Name","DepMap_ID","COSMICID", "RRID","cell_line_name")],by.x="stripped_cell_line_name",by.y = "stripped_cell_line_name",all.x = TRUE)

sample_annotations <- sample_annotations %>%
  rename(cell_line_name_depmap=cell_line_name.y)%>%
  rename(cell_line_name=cell_line_name.x)%>%
  rename(depMapID=DepMap_ID)

#Removing duplicate rows that are there as a result of the number of counts files
sample_annotations <- distinct(sample_annotations)

sample_annotations.path <- file.path(repo_path, 'METADATA','5429_sample_annotations.tsv')
write.table(sample_annotations, sample_annotations.path, sep = "\t", quote = F, row.names = F)

print(paste("Sample annotations written to:", sample_annotations.path))



