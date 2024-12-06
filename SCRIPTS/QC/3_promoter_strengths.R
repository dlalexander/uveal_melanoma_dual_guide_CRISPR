
suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(ggrepel)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-d", "--dir"), type = "character",
              help = "full path to repository", metavar = "character")
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

#dataset <- "lfcs_unscaled"
dataset <- "lfcs_scaled"

if (dataset == "lfcs_scaled") {
  lfc_matrix <- read.delim(file.path(repo_path,"DATA/preprocessing/lfc_matrix.scaled.all.tsv"), sep = "\t", header = T, check.names = F)
  print(paste("Reading scaled LFC matrix")) 
}
if (dataset == "lfcs_unscaled") {
  lfc_matrix <- read.delim(file.path(repo_path,"DATA/preprocessing/lfc_matrix.unscaled.all.tsv"), sep = "\t", header = T, check.names = F)
  print(paste("Reading unscaled LFC matrix")) 
}

# Pivot longer by sample 
lfc_matrix_longer <- lfc_matrix %>%
  pivot_longer(cols=c(16:72), names_to="sample_label")


############################################################
# Promoter strengths#
############################################################

promoter_plot_path <- file.path(QC_path,'promoter_strengths')

if (!dir.exists(promoter_plot_path)) {
  dir.create(promoter_plot_path)
}

# Classify the guide pairs 
lfc_matrix_promoters <- lfc_matrix_longer %>%
  # Filtering out safe safe combo 
  filter(guide_type %in% c('gene|safe_control', 'safe_control|gene')) %>% 
  # Classify the gene groups for plotting
  mutate(plot_status = case_when(
    sgrna_group == "Essential|Essential" ~ 'Essential',
    sgrna_group == "Non_essential|Non_essential"  ~ 'Nonessential'
  ))%>%
  filter(plot_status%in% c("Nonessential","Essential"))

# Grouping by gene to find the mean guide LFC per gene in a particular orientation
lfc_matrix_promoters_means <- lfc_matrix_promoters %>%
  group_by(singles_target_gene,guide_type,plot_status)%>%
  summarise(avg_gene_orientation = mean(value))
# Pivot wider to have both the orientations by gene
lfc_matrix_promoters_means_short <- lfc_matrix_promoters_means %>%
  pivot_wider(names_from = guide_type,values_from = avg_gene_orientation)%>%
  rename("safe_gene"="safe_control|gene")%>%
  rename("gene_safe"="gene|safe_control")


# Plot separately by Essential and Nonessential 
ggplot(lfc_matrix_promoters_means_short,aes(x=safe_gene,y=gene_safe,col=plot_status))+
  geom_point()+
  theme_classic()+
  geom_hline(yintercept = 0, color = 'gray30', linetype = 'dashed', linewidth = 0.5)+
  geom_vline(xintercept = 0, color = 'gray30', linetype = 'dashed', linewidth = 0.5)+
  geom_abline(intercept = 0, slope = 1)+
  ylim(c(-5.5,1))+
  xlim(c(-5.5,1))
ggsave(filename="promoter_strengths_by_gene.png",path =promoter_plot_path)

# Plot each replicate/cell line separately
lfc_matrix_promoters_per_replicate <- lfc_matrix_promoters %>%
  group_by(singles_target_gene,guide_type,sample_label,plot_status)%>%
  summarise(avg_gene_replicate_orientation = mean(value))

# For plotting 
lfc_matrix_promoters_per_replicate_short <- lfc_matrix_promoters_per_replicate%>%
  pivot_wider(names_from = guide_type,values_from = avg_gene_replicate_orientation)%>%
  rename("safe_gene"="safe_control|gene")%>%
  rename("gene_safe"="gene|safe_control")

#Finding the means to display on the plot
lfc_essentiality_promoter_means <- lfc_matrix_promoters %>%
  group_by(plot_status,guide_type)%>%
  summarise(avg_ess_orientation = mean(value))%>%
  pivot_wider(names_from = guide_type,values_from = avg_ess_orientation)%>%
  rename("safe_gene"="safe_control|gene")%>%
  rename("gene_safe"="gene|safe_control")

ggplot(lfc_matrix_promoters_per_replicate_short)+
  geom_point(aes(x=safe_gene,y=gene_safe,col=plot_status),size=0.25)+
  geom_point(data=lfc_essentiality_promoter_means,aes(x=safe_gene,y=gene_safe,fill=plot_status),size=4,stroke=1.5,shape=22)+
  theme_classic()+
  geom_hline(yintercept = 0, color = 'gray42', linetype = 'dashed')+
  geom_vline(xintercept = 0, color = 'gray42', linetype = 'dashed')+
  geom_abline(intercept = 0, slope = 1,color='grey24',linetype = 'longdash')+
  scale_colour_manual(values=c("lightcoral", "darkcyan"))+
  scale_fill_manual(values=c("lightcoral", "darkcyan"))+
  labs(colour="Gene type",fill="Mean by gene type")+ 
  xlab("LFC for guides in safe_gene orientation")+
  ylab("LFC for guides in gene_safe orientation")+
  theme(axis.title.x = element_text(size=10))
  # Add axis titles
  #ylim(c(-5.5,1))+
  #xlim(c(-5.5,1))
ggsave(filename="promoter_strengths_by_gene_sample.png",path=promoter_plot_path)




