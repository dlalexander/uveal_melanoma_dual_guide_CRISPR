suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(GGally)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character"),
  make_option(c("-m", "--mapping"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
  make_option(c("-d", "--data"), type = "character",
              help = "data type", metavar = "character")
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

QC_path <- file.path(repo_path, 'DATA', 'QC')

# Creating all output directories
# Create QC output directory
if (!dir.exists(QC_path)) {
  dir.create(QC_path)
}

# Create QC output directory
if (!dir.exists(file.path(QC_path,opt$data,"unscaled","density"))) {
  dir.create(file.path(QC_path,opt$data,"unscaled","density"),recursive=T)
}


############################################################
# Sample mapping                                           #
############################################################

check_file_exists(opt$mapping)

print(paste("Reading sample annotations from:", opt$mapping))

# Read in sample mapping
sample_mapping <- read.delim(opt$mapping, header = T, sep = "\t")

############################################################
# Data matrix                                             #
############################################################

# Read in count matrix or LFC matrix depending on data type 
# Create merged LFC matrix 

if (opt$data == "normalised_counts") {
  count_matrix_path <- file.path(repo_path,"DATA/preprocessing/count_matrix.norm.tsv")
  print(paste("Reading normalised count matrix:",count_matrix_path))
  # Read in count matrix
  data_matrix <- read.delim(count_matrix_path, sep = "\t", header = T, check.names = F)
}


if (opt$data == "lfcs") {
  print(paste("Reading batch A and batch B unscaled LFC matrices"))
  # Read in count matrix
  data_matrix_A <- read.delim(file.path(repo_path,"DATA/preprocessing/batchA/lfc_matrix.unscaled.all.tsv"), sep = "\t", header = T, check.names = F)
  data_matrix_B <- read.delim(file.path(repo_path,"DATA/preprocessing/batchB/lfc_matrix.unscaled.all.tsv"), sep = "\t", header = T, check.names = F)
  print(paste("Merging batch LFC matrices, retaining guides common to both")) 
  data_matrix <- inner_join(data_matrix_A,data_matrix_B)
  # NB slightly different guides in batch A and batch B LFCs because of different guides being filtered out due to low counts
  # Only retaining those common to both LFC matrices 
}

############################################################
# Unscaled LFC density for the different replicates, tiled by cell line#
############################################################

data_matrix_longer <- data_matrix %>%
  pivot_longer(cols=c(16:72), names_to="sample_label")%>%
  left_join(sample_mapping%>%select(c("sample_label","qc_group")),by="sample_label")


for(sample_i in unique(data_matrix_longer$sample_label)){
   unscaled_sko_safe <- data_matrix_longer %>% 
    filter(sample_label == sample_i) %>% 
    filter(guide_type %in% c('gene|safe_control', 'safe_control|gene', 'safe_control|safe_control')) %>% 
    mutate(essential_status = case_when(
      sgrna_group == "safe_safe|safe_safe" ~ 'safe_safe',
      sgrna_group == "Essential|Essential" ~ 'essential',
      sgrna_group == "Non_essential|Non_essential"  ~ 'non_essential',
      TRUE ~ 'single_unknown' #For when other conditions are not met, but still will be singles that get paired in duals
    ))
  unscaled_sko_safe_plot <- ggplot(unscaled_sko_safe,aes(x = value, fill = essential_status)) +
      geom_density(alpha = 0.7) +
      theme_bw() +
      labs(title = paste('Unscaled LFCs for SKO and safe_safe guides in sample:', sample_i))
  pdf(file.path(QC_path,opt$data,"unscaled","density",paste(sample_i,'single_safe_unscaled_lfc_density_plots.pdf')), height = 5, width = 8)
  print(unscaled_sko_safe_plot)
  dev.off()
}

# Use different classifications for the essentials 
# From DepMap Public 22Q4 CRISPRInferredCommonEssentials.csv
# Create comparison data directory
if (!dir.exists(file.path(QC_path,"comparison_data"))) {
  dir.create(file.path(QC_path,"comparison_data"),recursive=T)
}

depmap_essentials_path <- file.path(QC_path,"comparison_data/CRISPRInferredCommonEssentials.csv")

check_file_exists(depmap_essentials_path)

print(paste("Reading DepMap post Chronos essentials from:", depmap_essentials_path))

depmap_essentials <- read_csv(depmap_essentials_path)

depmap_essentials <- depmap_essentials %>% separate(Essentials,into=c("gene_symbol","gene_ID"),sep=" ")

for(sample_i in unique(data_matrix_longer$sample_label)){
  unscaled_sko_safe_depmap <- data_matrix_longer %>% 
    filter(sample_label == sample_i) %>% 
    filter(guide_type %in% c('gene|safe_control', 'safe_control|gene', 'safe_control|safe_control')) %>% 
    mutate(essential_status = case_when(
      sgrna_group == "safe_safe|safe_safe" ~ 'safe_safe',
      sgrna_group == "Essential|Essential" & !(singles_target_gene %in% depmap_essentials$gene_symbol) ~ 'library_essential_not_depmap',
      sgrna_group == "Non_essential|Non_essential"  ~ 'non_essential',
      singles_target_gene %in% depmap_essentials$gene_symbol ~ 'depmap_common_essential',
      TRUE ~ 'single_unknown' #For when other conditions are not met, but still will be singles that get paired in duals
    ))
  unscaled_sko_safe_plot_depmap <- ggplot(unscaled_sko_safe_depmap,aes(x = value, fill = essential_status)) +
    geom_density(alpha = 0.7) +
    theme_bw() +
    labs(title = paste('Unscaled LFCs for SKO and safe_safe guides in sample:', sample_i))
  pdf(file.path(QC_path,opt$data,"unscaled","density",paste(sample_i,'single_safe_unscaled_lfc_density_depmap_essentials.pdf')), height = 5, width = 8)
  print(unscaled_sko_safe_plot_depmap)
  dev.off()
}


print("Comparing to published uveal melanoma datasets")

mel202_fitness_path <- file.path(QC_path,"comparison_data/depmap-fitness-table-mel-202.csv")
check_file_exists(mel202_fitness_path)
print("Reading available Uveal Mel202 SKO fitness scores")
mel202_fitness <- read_csv(mel202_fitness_path)

pancan_essentials_mel202 <- mel202_fitness %>% filter(isPanCancer == TRUE)

# Split into Pan cancer essentials and uveal essentials
# Scaled BF < 0 are classed as essential
mel202_essentials_uveal <- mel202_fitness %>% filter(bf_scaled < 0 & is.na(isPanCancer))

# Plotting the uveal essentials and the pan cancer essentials 
for(sample_i in unique(data_matrix_longer$sample_label)){
  unscaled_sko_safe_mel202 <- data_matrix_longer %>% 
    filter(sample_label == sample_i) %>% 
    filter(guide_type %in% c('gene|safe_control', 'safe_control|gene', 'safe_control|safe_control')) %>% 
    mutate(essential_status = case_when(
      sgrna_group == "safe_safe|safe_safe" ~ 'safe_safe',
      sgrna_group == "Non_essential|Non_essential"  ~ 'non_essential',
      singles_target_gene %in% pancan_essentials_mel202$geneSymbol ~ 'pancan_essential',
      singles_target_gene %in% mel202_essentials_uveal$geneSymbol ~ 'mel202_essential',
      TRUE ~ 'single_unknown' #For when other conditions are not met, but still will be singles that get paired in duals
    ))
  unscaled_sko_safe_plot_mel202 <- ggplot(unscaled_sko_safe_mel202,aes(x = value, fill = essential_status)) +
    geom_density(alpha = 0.7) +
    theme_bw() +
    labs(title = paste('Unscaled LFCs for SKO and safe_safe guides in sample:', sample_i))
  pdf(file.path(QC_path,opt$data,"unscaled","density",paste(sample_i,'single_safe_unscaled_lfc_density_mel202_essentials.pdf')), height = 5, width = 8)
  print(unscaled_sko_safe_plot_mel202)
  dev.off()
}


# Read in the other uveal melanoma cell line dataset 
# Check what the overlap is in the essentials in the two uveal cell lines 

mel285_fitness_path <- file.path(QC_path,"comparison_data/depmap-fitness-table-mel-285.csv")
check_file_exists(mel285_fitness_path)
print("Reading available Uveal Mel285 SKO fitness scores")
mel285_fitness <- read_csv(mel285_fitness_path)

pancan_essentials_mel285 <- mel285_fitness %>% filter(isPanCancer == TRUE)

# Split into Pan cancer essentials and uveal essentials
# Scaled BF < 0 are classed as essential
mel285_essentials_uveal <- mel285_fitness %>% filter(bf_scaled < 0 & is.na(isPanCancer))

# Check what the fitness scores are of these genes
print(paste("Mean fitness scores for mel285 non-pancancer essentials:",round(mean(mel285_essentials_uveal$bf_scaled),digits=3)))
print(paste("Mean fitness scores for pancancer essentials in Mel285:",round(mean(pancan_essentials_mel285$bf_scaled),digits=3)))

# Check overlap between uveal line essentials 
print(paste(nrow(mel202_essentials_uveal[mel202_essentials_uveal$geneSymbol %in% mel285_essentials_uveal$geneSymbol,]), "essentials in both",nrow(mel202_essentials_uveal),"Mel202 essentials and",nrow(mel285_essentials_uveal),"Mel285 essentials"))

for(sample_i in unique(data_matrix_longer$sample_label)){
  unscaled_sko_safe_mel285 <- data_matrix_longer %>% 
    filter(sample_label == sample_i) %>% 
    filter(guide_type %in% c('gene|safe_control', 'safe_control|gene', 'safe_control|safe_control')) %>% 
    mutate(essential_status = case_when(
      sgrna_group == "safe_safe|safe_safe" ~ 'safe_safe',
      sgrna_group == "Non_essential|Non_essential"  ~ 'non_essential',
      singles_target_gene %in% pancan_essentials_mel285$geneSymbol ~ 'pancan_essential',
      singles_target_gene %in% mel285_essentials_uveal$geneSymbol ~ 'mel285_essential',
      TRUE ~ 'single_unknown' #For when other conditions are not met, but still will be singles that get paired in duals
    ))
  unscaled_sko_safe_plot_mel285 <- ggplot(unscaled_sko_safe_mel285,aes(x = value, fill = essential_status)) +
    geom_density(alpha = 0.7) +
    theme_bw() +
    labs(title = paste('Unscaled LFCs for SKO and safe_safe guides in sample:', sample_i))
  pdf(file.path(QC_path,opt$data,"unscaled","density",paste(sample_i,'single_safe_unscaled_lfc_density_mel285_essentials.pdf')), height = 5, width = 8)
  print(unscaled_sko_safe_plot_mel285)
  dev.off()
}

