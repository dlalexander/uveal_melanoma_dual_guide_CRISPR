
suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(scales)))
suppressPackageStartupMessages(suppressWarnings(library(ggpubr)))
suppressPackageStartupMessages(suppressWarnings(library(ggsci)))
suppressPackageStartupMessages(suppressWarnings(library(GGally)))
suppressPackageStartupMessages(suppressWarnings(library(ggrepel)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character"),
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

QC_path <- file.path(repo_path, 'DATA', 'QC')

# Creating all output directories
# Create QC output directory
if (!dir.exists(QC_path)) {
  dir.create(QC_path)
}


nnmd_plot_path <- file.path(QC_path,'nnmd_ssmd')

if (!dir.exists(nnmd_plot_path)) {
  dir.create(nnmd_plot_path)
}


############################################################
# Sample mapping                                           #
############################################################

print(paste("Reading sample annotations from:", opt$mapping))

# Read in sample mapping
sample_mapping <- read.delim(opt$mapping, header = T, sep = "\t")

############################################################
# LFC matrix.                                              #
############################################################

message(paste("Reading fold change matrix from:", opt$lfc))

# Read in lfc matrix
lfc_matrix <- read.delim(file.path(repo_path,"DATA/preprocessing/lfc_matrix.unscaled.all.tsv"), sep = "\t", header = T, check.names = F)

message("Gathering LFC matrix...")

annotation_colnames <- c(
  'id', 'sgrna_ids', 'sgrna_seqs', 'gene_pair_id',
  'sorted_gene_pair', 'targetA', 'targetB', 
  'sgrna_symbols', 'sgrna_symbol_A', 'sgrna_symbol_B', 
  'sgrna_libraries', 'sgrna_group', 'guide_type', 
  'guide_orientation', 'singles_target_gene')

lfc_matrix.narrow <- lfc_matrix %>%
   gather(sample_label, lfc, -all_of(annotation_colnames)) %>%
   left_join(sample_mapping, by = 'sample_label') 

############################################################
# Essential/nonessential separation assessment             #
############################################################

# Loading NNMD calculation functions 

calculate_nnmd <- function(data = NULL) {
  results <- data %>%
    group_by(sample_label, sgrna_group) %>%
    summarise(.groups = 'keep',
              'mean' = mean(lfc),
              'median' = median(lfc),
              'sd' = sd(lfc)) %>%
    gather(key, value, -sample_label, -sgrna_group) %>%
    pivot_wider(names_from = c(key,sgrna_group), values_from = c(value)) %>%
    mutate('NNMD' = (`mean_Essential|Essential` - `mean_Non_essential|Non_essential`) / `sd_Non_essential|Non_essential`) %>%
    arrange(-NNMD)
  return(results)
}

plot_nnmd <- function(data = NULL) {
  p <-
    ggplot(data, aes(x = sample_label, y = NNMD,fill = stripped_cell_line_name)) +
    geom_col(color = 'gray5', linewidth = 0.2, alpha = 0.8) +
    scale_y_reverse()+
    theme_classic()+
    #scale_y_reverse(breaks = pretty_breaks(12)) +
    facet_grid(. ~ cell_line_label, scales = 'free_x', space = 'free') +
    #geom_hline(yintercept = -2, color = 'gray30', linetype = 'dashed', linewidth = 0.5) +
    #theme_pubr(base_size = 14) +
    scale_fill_brewer(palette = "Spectral")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          panel.border = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()) +
    labs(color = "Category",  x = '', y = 'NNMD',fill="Uveal Cell Line")
  return(p)
}

message("Calculating NNMD...")

# Calculate NNMD (separation of Essential and Non-essential)
nnmd_results <- calculate_nnmd(lfc_matrix.narrow) %>%
  left_join(sample_mapping, by = 'sample_label')%>%
  arrange(sanger_sample_name)

# Ordering the labels 
nnmd_results$sanger_sample_name <- factor(nnmd_results$sample_label, levels = sample_mapping$sample_label)

# Plot NMMD
message("Plotting NNMD...")
nnmd_ordered_path <- file.path(nnmd_plot_path, 'normalised_LFC_NNMD_ordered.png')
nnmd_plot <- plot_nnmd(nnmd_results)
ggsave(filename = nnmd_ordered_path, plot = nnmd_plot, device = 'png', dpi = 300 , width = 4000, height = 2000, units = 'px')
message(paste("NNMD plot written to:", nnmd_ordered_path))


############################################################
# SSMDs                                                    #
############################################################

# Initially calculating this at a guide level 

calculate_ssmd <- function(data=NULL){
  # Selecting the positive and negative controls 
  results <- data %>%
    group_by(sample_label, sgrna_group) %>%
    summarise(.groups = 'keep',
              'mean' = mean(lfc),
              # Identifying the number of lfcs that are positive and negative controls
              'count'= length(lfc),
              # Find the standard deviation of each positive and negative control 
              'sd' = sd(lfc))%>%
    # Group by different vaule type 
    gather(key, value, -sample_label, -sgrna_group)%>%
    pivot_wider(names_from = c(key,sgrna_group), values_from = c(value))%>%
    mutate('SSMD' = (`mean_Essential|Essential` - `mean_Non_essential|Non_essential`) / sqrt(`sd_Essential|Essential`**2 + `sd_Non_essential|Non_essential`**2)) %>%
    arrange(-SSMD)
  return(results)
}

plot_ssmd <- function(data = NULL) {
  p <-
    ggplot(data, aes(x = sample_label, y = SSMD,fill = stripped_cell_line_name)) +
    geom_col(color = 'gray5', linewidth = 0.2, alpha = 0.8) +
    scale_y_reverse()+
    theme_classic()+
    #scale_y_reverse(breaks = pretty_breaks(12)) +
    facet_grid(. ~ cell_line_label, scales = 'free_x', space = 'free') +
    #geom_hline(yintercept = -2, color = 'gray30', linetype = 'dashed', linewidth = 0.5) +
    #theme_pubr(base_size = 14) +
    scale_fill_brewer(palette = "Spectral")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          panel.border = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()) +
    labs(color = "Category",  x = '', y = 'SSMD',fill="Uveal Cell Line")
  return(p)
}

message("Calculating SSMD...")

# Calculate SSMD
ssmd_results <- calculate_ssmd(lfc_matrix.narrow) %>%
  left_join(sample_mapping, by = 'sample_label')%>%
  arrange(sanger_sample_name)

# Ordering the labels 
ssmd_results$sanger_sample_name <- factor(ssmd_results$sample_label, levels = sample_mapping$sample_label)

# Plot SSMD
message("Plotting SSMD...")
ssmd_plot_path <- file.path(nnmd_plot_path, 'normalised_LFC_SSMD_ordered.png')
ssmd_plot <- plot_ssmd(ssmd_results)
ggsave(filename = ssmd_plot_path, plot = ssmd_plot, device = 'png', dpi = 300 , width = 4000, height = 2000, units = 'px')
message(paste("SSMD plot written to:", ssmd_plot_path))

############################################################
# Comparing NNMDs to those from Project Score              #
############################################################

# Project Score LFCs
score_lfcs <- read.delim(file.path(repo_path, 'DATA', 'QC','comparison_data','Project_Score_corrected_LFC.tsv', sep = "\t", header = T, check.names = F)

lfc_matrix_ess_only <- lfc_matrix %>%
  filter(sgrna_group == "Essential|Essential")
lfc_matrix_noness_only <- lfc_matrix %>%
  filter(sgrna_group == "Non_essential|Non_essential")

# Classify the genes as they have been classified in the library 
# Label as whether essential or non-essential 
score_lfcs_controls <- score_lfcs %>%
  filter(Gene %in% lfc_matrix.narrow$singles_target_gene)

score_lfcs.narrow <- score_lfcs_controls %>%
  pivot_longer(cols=!"Gene",names_to = "sample_label",values_to = "lfc")%>%
  mutate(sgrna_group= case_when(
    Gene %in% lfc_matrix_ess_only$singles_target_gene  ~ "Essential|Essential",
    Gene %in% lfc_matrix_noness_only$singles_target_gene ~ "Non_essential|Non_essential"))

# Calculate NNMDs per cell line
score_nnmd_results <- calculate_nnmd(score_lfcs.narrow)
score_ssmd_results <- calculate_ssmd(score_lfcs.narrow)

# Re-compute the uveal nnmds after having found the mean gene LFCs across reps
# So this is comparable to the Score NNMDs 
lfc_matrix.control_narrow_gene_means <- lfc_matrix.narrow %>%
  filter(sgrna_group %in% c("Essential|Essential","Non_essential|Non_essential"))%>%
  group_by(singles_target_gene,cell_line_name,sgrna_group)%>%
  summarise('mean_line_gene_lfc'=mean(lfc))

lfc_matrix.control_narrow_gene_means_lfcs <- lfc_matrix.control_narrow_gene_means %>%
  rename("lfc"="mean_line_gene_lfc")%>%
  rename("sample_label"="cell_line_name")

# Mean-level rather than guide-level lfcs influences NNMD calculated
gene_mean_lfc_nnmd_results <- calculate_nnmd(lfc_matrix.control_narrow_gene_means_lfcs)
gene_mean_lfc_ssmd_results <- calculate_ssmd(lfc_matrix.control_narrow_gene_means_lfcs)

# Checking by time point in addition to cell line
lfc_matrix.control_narrow_gene_means_by_time_point <- lfc_matrix.narrow %>%
  filter(sgrna_group %in% c("Essential|Essential","Non_essential|Non_essential"))%>%
  group_by(singles_target_gene,qc_group,sgrna_group)%>%
  summarise('mean_time_point_gene_lfc'=mean(lfc))

lfc_matrix.control_narrow_gene_means_lfcs_by_time_point <- lfc_matrix.control_narrow_gene_means_by_time_point %>%
  rename("lfc"="mean_time_point_gene_lfc")%>%
  rename("sample_label"="qc_group")

gene_mean_lfc_nnmd_results_time_point <- calculate_nnmd(lfc_matrix.control_narrow_gene_means_lfcs_by_time_point)
gene_mean_lfc_ssmd_results_time_point <- calculate_ssmd(lfc_matrix.control_narrow_gene_means_lfcs_by_time_point)

# Join cell line label 
gene_mean_lfc_nnmd_results_time_point <- gene_mean_lfc_nnmd_results_time_point %>%
  left_join(unique((sample_mapping %>% select(cell_line_name,qc_group))),by=c("sample_label"="qc_group"))%>%
  separate(col=sample_label,into=c("cell_line","time_point"),sep=" ")
# Split to separate time points
gene_mean_lfc_ssmd_results_time_point <- gene_mean_lfc_ssmd_results_time_point %>%
  left_join(unique((sample_mapping %>% select(cell_line_name,qc_group))),by=c("sample_label"="qc_group"))%>%
  separate(col=sample_label,into=c("cell_line","time_point"),sep=" ")

# Violin plot showing the range of NNMDs and SSMDs from Project Score 
# Overlay with dotplot showing the NNMDs from the cell lines in this study 
ggplot(score_nnmd_results,aes(x= "",y = NNMD,fill="NNMD")) +
  geom_violin(trim=FALSE)+
  #geom_point(gene_mean_lfc_nnmd_results,aes(x="",y=NNMD),position = position_jitter(0.2))+
  geom_jitter(data=gene_mean_lfc_nnmd_results, aes(x="", y=NNMD, color=sample_label),position = position_jitter(width = 0.02),size=3,alpha = 0.8)+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title = element_text(size=15))+
  ylab("Value")+
  xlab("NNMD")+
  scale_color_brewer(palette = "Spectral")+
  scale_fill_manual(values="slategray1")+
  labs(color='Uveal cell line',fill="Project Score") 

ggsave(filename=paste('nnmd_comparison_plot_gene_means','.png',sep=""),path =nnmd_plot_path, width=5,height=6)

# Time points labelled on the plot cell lines in the key 
ggplot(score_nnmd_results,aes(x= "",y = NNMD,fill="NNMD")) +
  geom_violin(trim=FALSE)+
  #geom_point(gene_mean_lfc_nnmd_results,aes(x="",y=NNMD),position = position_jitter(0.2))+
  geom_jitter(data=gene_mean_lfc_nnmd_results_time_point, aes(x="", y=NNMD, color=cell_line_name),position = position_jitter(width = 0.03),size=3,alpha = 0.8)+
  theme_classic()+
  geom_text_repel(data=gene_mean_lfc_nnmd_results_time_point,aes(label=time_point),size=3,max.overlaps = 30)+
  theme(axis.text=element_text(size=14),axis.title = element_text(size=15),legend.text = element_text(size=11),legend.title = element_text(size=11))+
  ylab("Value")+
  xlab("NNMD")+
  scale_color_brewer(palette = "Spectral")+
  scale_fill_manual(values="slategray1")+
  labs(color='Uveal cell line',fill="Project Score") 

ggsave(filename=paste('nnmd_comparison_plot_gene_means_time_point','.png',sep=""),path =nnmd_plot_path, width=5,height=6)


ggplot(score_ssmd_results,aes(x= "",y = SSMD,fill="SSMD")) +
  geom_violin(trim=FALSE)+
  #geom_point(gene_mean_lfc_nnmd_results,aes(x="",y=NNMD),position = position_jitter(0.2))+
  geom_jitter(data=gene_mean_lfc_ssmd_results, aes(x="", y=SSMD, color=sample_label),position = position_jitter(width = 0.02),size=3,alpha = 0.8)+
  theme_classic()+
  theme(axis.text=element_text(size=14),axis.title = element_text(size=15))+
  ylab("Value")+
  xlab("SSMD")+
  scale_color_brewer(palette = "Spectral")+
  scale_fill_manual(values="slategray1")+
  labs(color='Uveal cell line',fill="Project Score") 

ggsave(filename=paste('ssmd_comparison_plot_gene_means','.png',sep=""),path =nnmd_plot_path, width=5,height=6)

# Time points labelled on the plot cell lines in the key 
ggplot(score_ssmd_results,aes(x= "",y = SSMD,fill="SSMD")) +
  geom_violin(trim=FALSE)+
  #geom_point(gene_mean_lfc_nnmd_results,aes(x="",y=NNMD),position = position_jitter(0.2))+
  geom_jitter(data=gene_mean_lfc_ssmd_results_time_point, aes(x="", y=SSMD, color=cell_line_name),position = position_jitter(width = 0.04,seed=2),size=3,alpha = 0.8)+
  theme_classic()+
  #geom_text(data=gene_mean_lfc_ssmd_results_time_point,position=position_jitter(seed=1,width=0.1),aes(label=time_point),size=3)+
  geom_text_repel(data=gene_mean_lfc_ssmd_results_time_point,aes(label=time_point),size=3,max.overlaps = 30,position=position_jitter(seed=2,width=0.04))+
  theme(axis.text=element_text(size=14),axis.title = element_text(size=15),legend.text = element_text(size=11),legend.title = element_text(size=11))+
  ylab("Value")+
  xlab("SSMD")+
  scale_color_brewer(palette = "Spectral")+
  scale_fill_manual(values="slategray1")+
  labs(color='Uveal cell line',fill="Project Score") 

ggsave(filename=paste('ssmd_comparison_plot_gene_means_time_point','.png',sep=""),path =nnmd_plot_path, width=5,height=6)


