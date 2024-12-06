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
if (opt$data == "normalised_counts" & !dir.exists(file.path(QC_path,opt$data,"correlation","controls"))) {
  dir.create(file.path(QC_path,opt$data,"correlation","controls"),recursive = T)
}
if (!dir.exists(file.path(QC_path,opt$data,"correlation","test_lines"))) {
  dir.create(file.path(QC_path,opt$data,"correlation","test_lines"),recursive=T)
}


############################################################
# Sample mapping                                           #
############################################################

# For normalised counts 
opt$mapping <- file.path(repo_path,paste("METADATA/5429_sample_annotations.tsv",sep=""))

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
  ###
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
# Initial QC comparing the counts of the controls          #
############################################################

# Only want to do this for the normalised counts since don't have LFCs for the controls 
if (opt$data == "normalised_counts") {
  print("Plotting correlation of the pyCROQUET normalised counts for all control cell lines ")
  
  # Identify indices of control samples
  control_samples <-  sample_mapping %>%
    filter(grepl('Control', sample_label))
  
  control_sample_labels <- control_samples %>%
    pull(sample_label)
  
  print('Control sample labels:')
  print(control_sample_labels)
  
  control_sample_indices <- which(colnames(data_matrix) %in% control_sample_labels)
  
  normalised_counts_correlation_plot_all <- data_matrix %>%
    ggpairs(
      columns = control_sample_indices,
      xlab =("normalised count"),
      ylab =("normalised count"),
      upper = list(continuous = wrap("cor", size = 7))
    ) + 
    ggtitle("Correlation of normalised counts between control replicates")+
    theme(strip.text = element_text(size = 15),axis.text = element_text(size=11),plot.title = element_text(size = 20, face = "bold"))
  
  ggsave(filename = file.path(QC_path,"normalised_counts","correlation","controls","correlation_normalised_count_controls.png"),plot = normalised_counts_correlation_plot_all, height = 15, width = 15)
  
  
}


#####################
# Correlation between samples QC
#####################


# Function to get MDS plot from a correlation matrix
get_mds_plot <- function(correlation_matrix) {
  
  # Using correlation matrix, measure distance and plot in 2 dimensions
  distance_matrix <- sqrt(2 * (1 - correlation_matrix))
  distance_matrix <- as.dist(distance_matrix)
  mds_model <- cmdscale(distance_matrix)
  mds_model <- as.data.frame(mds_model)
  colnames(mds_model) <- c('x', 'y')
  mds_model <- mds_model %>% rownames_to_column('sample_id')
  
  # Draw MDS plot
  mds_plot <- mds_model %>% 
    ggplot(aes(x= x, y = y, label = sample_id)) +
    geom_point() +
    geom_text(size=5)+
    theme_bw()
  
  return(mds_plot)
}


print(paste("Generating correlation matrix of ",opt$data," for all cell lines"))

# Heatmap of correlation of all samples normalised counts 
# Generate correlation matrix between all samples (using Pearson) and draw heatmap
correlation_matrix_all_samples <- data_matrix %>% 
  select(-c(1:15)) %>% 
  cor(method = 'pearson')

write.table(correlation_matrix_all_samples,file=file.path(QC_path,opt$data,"correlation",'all_samples_correlation_matrix.txt'))
write.csv(as.data.frame(correlation_matrix_all_samples),file=file.path(QC_path,opt$data,"correlation",'all_samples_correlation_matrix.csv'))

#Finding the average R value for each line, and their range 
corr_matrix_no_controls <- as.data.frame(correlation_matrix_all_samples)  %>%
  select(-contains("Control"))%>%
  filter(!str_detect(rownames(correlation_matrix_all_samples), "Control"))

control_samples <-  sample_mapping %>%
  filter(grepl('Control', sample_label))

control_sample_labels <- control_samples %>%
  pull(sample_label)

print('Control sample labels:')
print(control_sample_labels)

control_sample_indices <- which(colnames(data_matrix) %in% control_sample_labels)


# MDS plot
mds_plot_all_samples <- get_mds_plot(correlation_matrix_all_samples)

# Define heatmap colour palette
colour_palette <- colorRampPalette(c("#F0F0F0", "#DB5634"))(10)

print("Generating Heatmap of correlation between cell lines and replicates")

# Write out all plots
pdf(file.path(QC_path,opt$data,"correlation","test_lines",paste("correlation_heatmap_",opt$data,".pdf",sep="")), height = 19, width = 19)
heatmap(correlation_matrix_all_samples, col = colour_palette, scale = "none", sym = TRUE, margins = c(7, 7))
dev.off()

print("Generating MDS plot to show similarity of cell lines and replicates")

pdf(file.path(QC_path,opt$data,"correlation","test_lines", paste("correlation_mds_plot_",opt$data,".pdf",sep="")), height = 19, width = 19)
mds_plot_all_samples
dev.off()

#####################
# Correlation of control and day 28 samples, excluding the later time point
#####################

# Filter to exclude later time points 
data_matrix_day28_control <- data_matrix %>%
  select(-c(1:15)) %>% # Remove guide pair details 
  select(contains("D28") | contains("Control"))

# Heatmap of correlation of all samples normalised counts 
# Generate correlation matrix between all samples (using Pearson) and draw heatmap
correlation_matrix_day28_control <- data_matrix_day28_control %>% 
  cor(method = 'pearson')

# MDS plot
mds_plot_day28_control <- get_mds_plot(correlation_matrix_day28_control)

# Define heatmap colour palette
colour_palette <- colorRampPalette(c("#F0F0F0", "#DB5634"))(10)

print("Generating Heatmap of correlation between cell lines and replicates")

# Write out all plots
pdf(file.path(QC_path,opt$data,"correlation","test_lines",paste("correlation_heatmap_day28_controls",opt$data,".pdf",sep="")), height = 19, width = 19)
heatmap(correlation_matrix_day28_control, col = colour_palette, scale = "none", sym = TRUE, margins = c(20, 20),cexRow=2,cexCol=2)
dev.off()

#Alternative way of plotting 
jpeg(filename=file.path(QC_path,opt$data,"correlation","test_lines",paste("correlation_heatmap_day28_controls_with_legend_",opt$data,".jpg",sep="")),width = 220, height = 220,units='mm', res = 300)
heatmap.2(correlation_matrix_day28_control,
          trace = "none",
          dendrogram = "both",
          xlab="Cell line and replicate",
          ylab="Cell line and replicate",
          col=colour_palette,
          key=TRUE,
          density.info="none",
          keysize=0.75,
          margins=c(8.5,8.5),
          key.title = "",
          key.xlab = "Pearson correlation\n coefficient",
          key.ylab="",
          key.par = list(cex=0.7,font=1),
          family="Arial",
          offsetRow = 0,
          offsetCol = 0,
          cexRow = 0.95,
          cexCol = 0.95,
          #If wanted to include grid lines
          #sepwidth=c(0.000001,0.000001),
          #sepcolor="white",
          #colsep=1:ncol(correlation_matrix_day28_control),
          #rowsep=1:nrow(correlation_matrix_day28_control),
          scale="none")
dev.off()


print("Generating MDS plot to show similarity of cell lines and replicates")

pdf(file.path(QC_path,opt$data,"correlation","test_lines", paste("correlation_mds_plot_day28_controls",opt$data,".pdf",sep="")), height = 19, width = 19)
mds_plot_day28_control
dev.off()


#####################
# Loop through to plot correlation matrix and MDS plot by cell line, with control for that batch 
#####################

print("Generating Heatmap of correlation between samples by cell line")

# For loop to plot by cell line 

for(cell_line in unique(sample_mapping$cell_line_label)){
  #Filter to focus on the cell line and control, or just cell line for LFCs
  
  plot_batch <- unique(sample_mapping$batch[sample_mapping$cell_line_label == cell_line])
  
  if (opt$data=="normalised_counts") {
    plot_samples <-  sample_mapping %>%
      filter(cell_line_label==cell_line | batch==plot_batch&grepl('control', batch_control))%>%
      pull(sample_label)
  } else if (opt$data=="lfcs") {
    plot_samples <-  sample_mapping %>%
      filter(cell_line_label==cell_line & is.na(sample_mapping$batch_control))%>%
      pull(sample_label)
  }
  
  plot_correlation_matrix <- correlation_matrix_all_samples[plot_samples,plot_samples]

  pdf(file.path(QC_path,opt$data,"correlation","test_lines",paste(cell_line,"correlation_heatmap_",opt$data,".pdf",sep="")))
  heatmap(plot_correlation_matrix, col = colour_palette, scale = "none", sym = TRUE, margins = c(10, 10))
  dev.off()
}


#####################
# Excluding outlier samples to see separation of remaining samples better
#####################

if (opt$data == "normalised_counts") {
  
  print("Generating MDS plots with outliers excluded to show sample separation in better detail")
  
  #Finding indicies of samples that want to exclude 
  exclude_samples <-  sample_mapping %>%
    filter(grepl('Mel285', sample_label)) %>%
    pull(sample_label)
  
  exclude_sample_indices <- which(colnames(data_matrix) %in% exclude_samples)
  
  # Generate correlation matrix between all samples (using Pearson) and draw heatmap
  correlation_matrix_excl_mel285 <- data_matrix %>% 
    select(-c(1:15,exclude_sample_indices)) %>% 
    cor(method = 'pearson')
  
  # MDS plot
  mds_plot_excl_mel285 <- get_mds_plot(correlation_matrix_excl_mel285)
  
  pdf(file.path(QC_path,"normalised_counts","correlation","test_lines", 'correlation_mds_plot_normalised_counts_excl_mel285.pdf'), height = 19, width = 19)
  mds_plot_excl_mel285
  dev.off()
  
  pdf(file.path(QC_path,"normalised_counts","correlation","test_lines",'correlation_heatmap_normalised_counts_excl_mel285.pdf'), height = 19, width = 19)
  heatmap(correlation_matrix_excl_mel285 , col = colour_palette, scale = "none", sym = TRUE, margins = c(7, 7))
  dev.off()
  
  # Excluding additional samples for MDS plot to see better separation 
  
  # Finding indicies of samples that want to exclude 
  exclude_samples2 <-  sample_mapping %>%
    filter(grepl('Mel285|MP46 D64|MP38 D65A|Control|921', sample_label)) %>%
    pull(sample_label)
  
  exclude_sample_indices2 <- which(colnames(data_matrix) %in% exclude_samples2)
  
  # Generate correlation matrix between all samples (using Pearson) and draw heatmap
  correlation_matrix_excl <- data_matrix %>% 
    select(-c(1:15,exclude_sample_indices2)) %>% 
    cor(method = 'pearson')
  
  # MDS plot
  mds_plot_excl <- get_mds_plot(correlation_matrix_excl)
  
  pdf(file.path(QC_path,"normalised_counts","correlation","test_lines", 'correlation_mds_plot_normalised_counts_excl.pdf'), height = 19, width = 19)
  mds_plot_excl
  dev.off()
  
  pdf(file.path(QC_path,"normalised_counts","correlation","test_lines",'correlation_heatmap_normalised_counts_excl.pdf'), height = 19, width = 19)
  heatmap(correlation_matrix_excl , col = colour_palette, scale = "none", sym = TRUE, margins = c(7, 7))
  dev.off()
  
}


print("DONE.")



