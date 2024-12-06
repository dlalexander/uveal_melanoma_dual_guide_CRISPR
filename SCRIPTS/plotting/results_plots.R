suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(GGally)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(RColorBrewer)))
suppressPackageStartupMessages(suppressWarnings(library(gplots)))
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

setwd(repo_path)

# Check top level directory exists
if (!dir.exists(repo_path)) {
  stop(sprintf("Repository directory not exist: %s", repo_path))
}

# Add helper functions
source(file.path(repo_path, 'SCRIPTS', 'preprocessing', 'subscripts', 'helper.R'))

plot_path <- file.path(repo_path, 'DATA', 'results_plots')

# Creating all output directories
# Create QC output directory
if (!dir.exists(plot_path)) {
  dir.create(plot_path)
}

############################################################
# Sample mapping                                           #
############################################################

check_file_exists(opt$mapping)

print(paste("Reading sample annotations from:", opt$mapping))

# Read in sample mapping
sample_mapping <- read.delim(opt$mapping, header = T, sep = "\t")


############################################################
# Observed and expected LFCs                               #
############################################################

# Plot an example observed vs expected for different cell lines

pred_vs_obvs_with_ngis <- read.delim(file.path(repo_path,"DATA/dual_guide/scaled/pred_vs_obs_with_ngis.tsv"), sep = "\t", header = T)

obs_path <- file.path(plot_path,'obs_vs_exp')

# Creating all output directories
# Create QC output directory
if (!dir.exists(obs_path )) {
  dir.create(obs_path )
}

for(group in unique(pred_vs_obvs_with_ngis$qc_group)){
  obs_vs_exp_plot <- pred_vs_obvs_with_ngis %>% 
    filter(qc_group == group) %>% 
    ggplot(aes(x = pred_y12, y = obs_y12)) +
    geom_point(alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey') +
    geom_hline(yintercept = 0, color = 'darkgrey') +
    geom_vline(xintercept = 0, color = 'darkgrey') +
    theme_bw() +
    theme_classic()+
    facet_wrap(~ sample)
  png(file.path(obs_path,paste('obs_vs_exp_',group,'.png',sep="")),width = 600,height = 300)
  plot(obs_vs_exp_plot)
  dev.off()
}

# Plot observed vs expected gammas

for(group in unique(pred_vs_obvs_with_ngis$qc_group)){
  gamma_gis_plot <- pred_vs_obvs_with_ngis %>% 
    filter(qc_group == group) %>% 
    ggplot(aes(x = predicted_gamma_12, y = gamma_12)) +
    geom_point(alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey') +
    geom_hline(yintercept = 0, color = 'darkgrey') +
    geom_vline(xintercept = 0, color = 'darkgrey') +
    geom_line(aes(x = predicted_gamma_12, y = loesspredicted), color = 'darkred') +
    geom_text(aes)
    theme_bw() +
    facet_wrap(~ sample)
  png(file.path(obs_path,paste('loess_gammas_',group,'.png',sep="")),width = 600,height = 300)
  plot(gamma_gis_plot)
  dev.off()
}


############################################################
# GI scores and t-test results                             #
############################################################

gene_level_results <- read.delim(file.path(repo_path,"DATA/postprocessing/combined_gene_level_results.scaled.tsv"), sep = "\t", header = T)

gene_level_results_binary <- read.delim(file.path(repo_path,"DATA/postprocessing/combined_gene_level_results.binary.scaled.tsv"), sep = "\t", header = T)

gi_fdr_path <- file.path(plot_path,'gi_fdr')

# Creating all output directories
if (!dir.exists(gi_fdr_path)) {
  dir.create(gi_fdr_path)
}

group_hit_counts <- gene_level_results_binary %>%
  select(sorted_gene_pair,bassik__Uveal_melanoma)%>%
  rename(n_group_hits = bassik__Uveal_melanoma)

gene_level_results_counts <- gene_level_results %>%
  left_join(group_hit_counts, by = 'sorted_gene_pair') %>%
  mutate(hit = case_when(
    is_bassik_hit == 1 ~ 'SL_hit',
    is_bassik_hit == 0 ~ 'not_hit'
  ))

# Volcano plot - x axis is mean residual, y axis is neglog fdr, size of point is how many times it is a negative hit
# Filter used was: mean_norm_gi < -0.5 & fdr < 0.01, 1, 0

for(group in unique(gene_level_results_counts$qc_group)){
  volcano_plot <- gene_level_results_counts %>% 
    filter(qc_group == group) %>% 
    ggplot(aes(x = mean_norm_gi, y = -log(fdr, base = 10), size = n_group_hits,label = sorted_gene_pair,color=hit)) +
    geom_point(alpha = 0.3) +
    theme_bw() +
    geom_hline(yintercept = -log(0.01, 10), linetype = 'dashed', col = 'grey') +
    geom_vline(xintercept = -0.5, linetype = 'dashed', col = 'grey')+
    facet_wrap(~qc_group)
  ggsave(filename=paste('volcano_gi_fdr_',group,'.png',sep=""),path =gi_fdr_path,plot=volcano_plot, width=6,height=3)
}

############################################################
# Guide level LFC density across cell lines and replicates #
############################################################

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

lfc_plot_path <- file.path(plot_path,'lfc_densities')

if (!dir.exists(lfc_plot_path)) {
  dir.create(lfc_plot_path)
}

# Pivot longer by sample 
lfc_matrix_longer <- lfc_matrix %>%
  pivot_longer(cols=c(16:72), names_to="sample_label")

# Classify the guide pairs 
lfc_matrix_longer <- lfc_matrix_longer %>%
  filter(guide_type %in% c('gene|safe_control', 'safe_control|gene', 'gene|gene')) %>% 
  # Classify the gene groups for plotting
  mutate(plot_status = case_when(
    sgrna_group == "Essential|Essential" ~ 'Essential control',
    sgrna_group == "Non_essential|Non_essential"  ~ 'Non-essential control',
    guide_type == 'gene|gene'& sgrna_group %in% c("Achilles|Achilles","Jen Uveal|Jen Uveal", "Emanuel|Emanuel","paralog|paralog")~ 'Double knockout',
    ### 
    TRUE ~ 'Single knockout'
))

# LFC distribution 
lfc_distribution <- ggplot(lfc_matrix_longer,aes(x = value, fill = factor(plot_status,levels = c('Essential control','Non-essential control','Single knockout','Double knockout')))) +
    geom_density(alpha = 0.7) +
    theme_bw() +
    theme_classic()+
    xlim(c(-3.1,1.1))+
    xlab("Log2 fold change")+
    ylab("Density")+
    labs(fill="pgRNA type")+
    scale_fill_manual(values=c("#91B7A9","#EBBE62","#D04C65","#4F799D"))+
    theme(legend.position = "inside",
          legend.position.inside = c(0.05,0.95),
          legend.justification = c("left", "top"),
          legend.box.just = "left")
pdf(file.path(lfc_plot_path,paste(dataset,"_density_plot.pdf")), height = 5, width = 5)
print(lfc_distribution)
dev.off()

# Boxplot for all gene groups as defined by the library 
all_group_lfcs_boxplot <- lfc_matrix_longer %>%
  ggplot(aes(x = factor(plot_status,levels = c('Non-essential control','Essential control','Single knockout','Double knockout')), y = value,fill = factor(plot_status,levels = c('Essential control','Non-essential control','Single knockout','Double knockout')))) +
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("#91B7A9","#EBBE62","#D04C65","#4F799D"))+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text=element_text(size=10.5),axis.title=element_text(size=16))+
  ylim(c(-2.8,1.2))+
  xlab("pgRNA type")+
  ylab("Log2 fold change")
#scale_color_manual())
ggsave(filename="all_group_lfcs_library_boxplot.png",path =lfc_plot_path,plot=all_group_lfcs_boxplot, width=6,height=6)

# Only the SKO and DKO plots 
sko_dko_lfcs_boxplot <- lfc_matrix_longer %>% 
  filter(plot_status%in% c("Double knockout","Single knockout"))%>%
  ggplot(aes(x = factor(plot_status,levels = c('Single knockout','Double knockout')), y = value,fill = factor(plot_status,levels = c('Single knockout','Double knockout')))) +
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("#D04C65","#4F799D"))+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text=element_text(size=10.5),axis.title=element_text(size=16))+
  ylim(c(-2.8,1.2))+
  xlab("pgRNA type")+
  ylab("Log2 fold change")
#scale_color_manual())
ggsave(filename="sko_dko_lfcs_library_boxplot.png",path =lfc_plot_path,plot=sko_dko_lfcs_boxplot, width=4,height=6)


lfc_matrix_double_knockout <- lfc_matrix_longer %>%
  filter(plot_status=="Double knockout")
lfc_matrix_single_knockout <- lfc_matrix_longer %>%
  filter(plot_status=="Single knockout")

# Wilcoxon test to compare the groups

# Double knockout vs single knockout 
wilcox.test(lfc_matrix_double_knockout$value, lfc_matrix_single_knockout$value, alternative = 'less', paired = FALSE)

# Essential vs nonessential 
wilcox.test(lfc_matrix_longer_Ess_only$value, lfc_matrix_longer_nonessential_only$value, alternative = 'less', paired = FALSE)

### Plotting to compare BAGEL essentials and non-essentials 
bagel_essentials <- read.delim(file.path(repo_path,"DATA/QC/essentials/CEGv2.txt"), sep = "\t", header = T, check.names = F)
bagel_nonessentials <- read.delim(file.path(repo_path,"DATA/QC/essentials/NEGv1.txt"), sep = "\t", header = T, check.names = F)

lfc_matrix_longer_bagel <- lfc_matrix_longer %>%
  # Filtering out safe safe combo 
  filter(guide_type %in% c('gene|safe_control', 'safe_control|gene', 'gene|gene')) %>% 
  # Classify the gene groups for plotting
  mutate(bagel_status = case_when(
    guide_type %in% c('gene|safe_control', 'safe_control|gene') & singles_target_gene %in% bagel_essentials$GENE~ 'Positive control',
    guide_type %in% c('gene|safe_control', 'safe_control|gene') & singles_target_gene %in% bagel_nonessentials$GENE ~ 'Negative control',
    guide_type == 'gene|gene'& sgrna_group %in% c("Achilles|Achilles","Jen Uveal|Jen Uveal", "Emanuel|Emanuel","paralog|paralog")~ 'Double knockout',
    ### 
    TRUE ~ 'Single knockout'
  ))

lfc_distribution_bagel <- ggplot(lfc_matrix_longer_bagel,aes(x = value, fill = bagel_status)) +
  geom_density(alpha = 0.7) +
  theme_bw() +
  labs(title='LFCs across paired guide RNA type')+
  xlim(c(-3,1))+
pdf(file.path(lfc_plot_path,paste(dataset,"bagel_density_plot.pdf")), height = 5, width = 8)
print(lfc_distribution_bagel)
dev.off()

lfc_matrix_longer_bagel$bagel_status <- factor(lfc_matrix_longer_bagel$bagel_status , levels=c('Negative control','Positive control', 'Single knockout', 'Double knockout'))

all_group_lfcs_bagel <- lfc_matrix_longer_bagel %>% 
  ggplot(aes(x = bagel_status, y = value,fill=bagel_status)) +
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("#f3c16f", "#93bbb0", "#db4642","#6086a4"))+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=17))+
  ylim(c(-2.8,1.2))+
  xlab("pgRNA type")+
  ylab("Log2 fold change")
  #scale_color_manual())
ggsave(filename="all_group_lfcs_bagel_boxplot.png",path =lfc_plot_path,plot=all_group_lfcs_bagel, width=6,height=6)

# Including outliers 
all_group_lfcs_bagel_outliers <- lfc_matrix_longer_bagel %>% 
  ggplot(aes(x = bagel_status, y = value,fill=bagel_status)) +
  geom_boxplot()+
  scale_fill_manual(values=c("#f3c16f", "#93bbb0", "#db4642","#6086a4"))+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=15))+
  xlab("pgRNA type")+
  ylab("Log2 fold change")
ggsave(filename="all_group_lfcs_bagel_boxplot_outliers.png",path =lfc_plot_path,plot=all_group_lfcs_bagel_outliers, width=6,height=6)

# Plot for just the DKO and SKO for comparison 
sko_dko_lfcs_bagel <- lfc_matrix_longer_bagel %>% 
  filter(bagel_status%in% c("Double knockout","Single knockout"))%>%
  ggplot(aes(x = bagel_status, y = value,fill=bagel_status)) +
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("#db4642","#6086a4"))+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=15))+
  ylim(c(-2.8,1.2))+
  xlab("pgRNA type")+
  ylab("Log2 fold change")
#scale_color_manual())
ggsave(filename="sko_dko_lfcs_bagel_boxplot.png",path =lfc_plot_path,plot=sko_dko_lfcs_bagel, width=6,height=6)


############################################################
# Top hits                                                 #
############################################################


top_hits_path <- file.path(plot_path,'top_hits')

if (!dir.exists(top_hits_path)) {
  dir.create(top_hits_path)
}

# Filtering out those that are singly essential in particular cell lines 
# Then re-counting how many cell lines they are hits in
gene_level_results_filtered<- gene_level_results %>% 
  filter(is_bassik_hit == 1 & targetA__is_single_depleted != 1 & targetB__is_single_depleted != 1)

gene_level_results_filtered_counts <- gene_level_results_filtered %>%
  group_by(sorted_gene_pair)%>%
  summarise(n_groups = sum(is_bassik_hit))

gene_hits_and_counts <- gene_level_results_filtered %>%
  left_join(gene_level_results_filtered_counts,by = 'sorted_gene_pair')

counts_to_plot <- gene_level_results_filtered_counts %>%
  filter(n_groups > 4)

counts_plot <- ggplot(filter(counts_to_plot, n_groups>6)) +
  ylab("Gene pair")+
  theme_classic()+
  xlab("Number of cell lines/time points where called as hit")+
  ggtitle("Gene pair hits in 7+ cell lines/time points")+
  geom_col(aes(y = reorder(sorted_gene_pair, n_groups),x=n_groups), width = 0.6)


# Plot those hits present in multiple cell lines 
# Horizontal bar plot to read gene pairs
ggsave(filename='hit_horizontal_bar.png',path =top_hits_path,plot=counts_plot, width=6,height=10)

# Plot GI scores for the hits in the different cell lines and replicates
dot_gi_plot <- ggplot(filter(gene_hits_and_counts,n_groups > 4), aes(x=median_norm_gi,y = reorder(sorted_gene_pair, n_groups))) +
  geom_point(aes(color = qc_group))+
  ylab("Gene pair")
ggsave(filename='top_hit_dot_GIs.png',path =top_hits_path,plot=dot_gi_plot, width=8,height=6)

box_gi_plot <- ggplot(filter(gene_hits_and_counts,n_groups > 7), aes(x=median_norm_gi,y = reorder(sorted_gene_pair, n_groups))) +
  ylab("Gene pair")+
  geom_boxplot(width=0.4,aes(fill = sorted_gene_pair))+
  theme_classic()+
  xlim(-2.8,0)+
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey')+
  ggtitle("GI distribution across hits \nin 8+ cell lines/time points")+
  theme(legend.position = "none",axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=16,face="bold",hjust = 0.5))
ggsave(filename='top_hit_box_GIs.png',path =top_hits_path,plot=box_gi_plot, width=7,height=10)

# Heatmap
# Would be good to include this 
hits_10_groups<-filter(gene_hits_and_counts,n_groups > 10)

top_hits <- gene_level_results %>%
  filter(gene_pair_id %in% hits_10_groups$gene_pair_id)%>%
  select(gene_pair_id,qc_group,median_norm_gi)%>%
  pivot_wider(names_from = qc_group,values_from=median_norm_gi)%>%
  column_to_rownames(var="gene_pair_id")

top_hits <- top_hits %>%
  select(order(colnames(top_hits)))

gene_level_results_filtered_counts <- gene_level_results_filtered_counts %>%
  arrange(desc(n_groups))

jpeg(filename=file.path(top_hits_path,"top_hit_GI_heatmap.jpg"),width = 220, height = 180,units='mm', res = 400)
heatmap.2(as.matrix(top_hits),
          trace = "none",
          Colv=FALSE,
          Rowv=FALSE,
          dendrogram = "none",
          xlab="Cell line and time point",
          ylab="Gene pair",
          col=colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
          key=TRUE,
          keysize=1,
          margins=c(8,10),
          key.title = "GI Score",
          key.xlab = "",
          key.ylab="",
          key.par = list(cex=0.6),
          scale="none")
dev.off()


############################################################
# Key gene pairs                                           #
############################################################

gene_pairs_path <- file.path(plot_path,'gene_pairs')

# Creating all output directories
# Create QC output directory
if (!dir.exists(gene_pairs_path)) {
  dir.create(gene_pairs_path)
}

# Plot CDS1/2 LFCs

pred_vs_obvs_CDS12 <- pred_vs_obvs_with_ngis %>%
  filter(sorted_gene_pair=="CDS1|CDS2")

# Longer matrix 
# Raw LFCs
pred_vs_obvs_CDS12_LFC <- pred_vs_obvs_CDS12 %>%
  rename(CDS1=y1,CDS2=y2,pred_DKO=pred_y12,obs_DKO=obs_y12)%>%
  select(sorted_gene_pair,qc_group,CDS1,CDS2,pred_DKO,obs_DKO)%>%
  pivot_longer(cols=c(CDS1,CDS2,pred_DKO,obs_DKO),names_to = "Category")%>%
  rename(LFC=value)

pred_vs_obvs_CDS12_LFC$Category <- factor(pred_vs_obvs_CDS12_LFC$Category , levels=c("CDS1", "CDS2", "pred_DKO", "obs_DKO"))

# Plotting in individual cell lines/time points 
for(group in unique(pred_vs_obvs_CDS12_LFC$qc_group)){
  CDS12_plot <- pred_vs_obvs_CDS12_LFC %>% 
    filter(qc_group == group) %>% 
    ggplot(aes(x = Category, y = LFC)) +
    geom_point(alpha = 0.5) +
    theme_classic()+
    geom_boxplot(width=0.4, color="grey", alpha=0.2,fill="deepskyblue1")+
    ggtitle(paste("CDS1/CDS2 knockout LFCs in ",group))
  ggsave(filename=paste('CDS12_LFC_',group,'.png',sep=""),path =gene_pairs_path,plot=CDS12_plot, width=6,height=3)
}

CDS12_boxplot_all <- pred_vs_obvs_CDS12_LFC %>% 
    ggplot(aes(x = Category, y = LFC)) +
    geom_point(alpha = 0.5) +
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),axis.text=element_text(size=12))+
    geom_boxplot(width=0.4, color="grey", alpha=0.2,fill="deepskyblue1")+
    ggtitle(paste("CDS1/CDS2 knockout LFCs for individual guides or pairs across \ncell lines, replicates, and time points"))
ggsave(filename=paste('CDS12_LFC_',"all_boxplot",'.png',sep=""),path =gene_pairs_path,plot=CDS12_boxplot_all, width=6,height=6)


CDS12_plot_violin_all <- pred_vs_obvs_CDS12_LFC %>% 
  ggplot(aes(x = Category, y = LFC)) +
  geom_violin(trim = TRUE,fill="lightblue")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),axis.text=element_text(size=12))+
  geom_boxplot(width=0.1,outlier.shape = NA)+
  ylim(-4,2.3)+
  ggtitle(paste("CDS1/CDS2 knockout LFCs for individual guides or pairs across \ncell lines, replicates, and time points"))
ggsave(filename=paste('CDS12_LFC_violin',"all_guides",'.png',sep=""),path =gene_pairs_path,plot=CDS12_plot_violin_all, width=8,height=6)


pred_vs_obvs_CDS12_averages_per_rep <- pred_vs_obvs_CDS12 %>%
  group_by(sample) %>% 
  summarise(CDS1=mean(y1), CDS2=mean(y2), predicted_DKO=mean(pred_y12),observed_DKO=mean(obs_y12))%>%
  pivot_longer(cols=c(CDS1,CDS2,predicted_DKO,observed_DKO),names_to = "Category")%>%
  rename(mean_LFC=value)

pred_vs_obvs_CDS12_averages_per_rep$Category <- factor(pred_vs_obvs_CDS12_averages_per_rep$Category, 
                                                       levels=c("CDS1", "CDS2", "predicted_DKO", "observed_DKO"))

CDS12_violin_plot_means <- pred_vs_obvs_CDS12_averages_per_rep %>% 
  ggplot(aes(x = Category, y = mean_LFC)) +
  geom_violin(trim = TRUE,fill="lightblue")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),axis.text=element_text(size=12))+
  geom_boxplot(width=0.1,outlier.shape = NA)+
  ggtitle(paste("CDS1/CDS2 knockout mean LFCs of guides in different \n cell lines, replicates, and time points"))
ggsave(filename=paste('CDS12_LFC_violin_',"all_means",'.png',sep=""),path =gene_pairs_path,plot=CDS12_violin_plot_means, width=8,height=6)

CDS12_box_plot_means <- pred_vs_obvs_CDS12_averages_per_rep %>% 
  ggplot(aes(x = Category, y = mean_LFC)) +
  geom_boxplot(width=0.4,fill="lightblue")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),axis.text=element_text(size=12))+
  ggtitle(paste("CDS1/CDS2 knockout mean LFCs of guides in different \n cell lines, replicates, and time points"))
ggsave(filename=paste('CDS12_LFC_boxplot_',"all_means",'.png',sep=""),path =gene_pairs_path,plot=CDS12_box_plot_means, width=8,height=6)


CDS12_box_plot_means_no_outliers <- pred_vs_obvs_CDS12_averages_per_rep %>% 
  ggplot(aes(x = Category, y = mean_LFC,fill=Category)) +
  geom_boxplot(width=0.55,outlier.shape = NA)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,size=18),axis.text=element_text(size=14),axis.title = element_text(size=16))+
  scale_fill_manual(values=c("#f3c16f", "#93bbb0", "#db4642","#6086a4"))+
  ylim(-1.6,0)+
  ggtitle(paste("CDS1/CDS2 knockout mean LFCs of guides across \n cell lines, replicates, and time points"))
ggsave(filename=paste('CDS12_LFC_boxplot_',"all_means_no_outliers_updated",'.png',sep=""),path =gene_pairs_path,plot=CDS12_box_plot_means_no_outliers, width=8,height=6)

# Plot ASF1A/ASF1B LFCs

pred_vs_obvs_ASF1AB <- pred_vs_obvs_with_ngis %>%
  filter(sorted_gene_pair=="ASF1A|ASF1B")

pred_vs_obvs_ASF1AB_LFC <- pred_vs_obvs_ASF1AB %>%
  rename(ASF1A=y1,ASF1B=y2,pred_DKO=pred_y12,obs_DKO=obs_y12)%>%
  select(sorted_gene_pair,qc_group,ASF1A,ASF1B,pred_DKO,obs_DKO)%>%
  pivot_longer(cols=c(ASF1A,ASF1B,pred_DKO,obs_DKO),names_to = "Category")%>%
  rename(LFC=value)

pred_vs_obvs_ASF1AB_LFC$Category <- factor(pred_vs_obvs_ASF1AB_LFC$Category , levels=c("ASF1A", "ASF1B", "pred_DKO", "obs_DKO"))

for(group in unique(pred_vs_obvs_ASF1AB_LFC$qc_group)){
  ASF1AB_plot <- pred_vs_obvs_ASF1AB_LFC %>% 
    filter(qc_group == group) %>% 
    ggplot(aes(x = Category, y = LFC)) +
    geom_point(alpha = 0.4) +
    theme_classic()+
    geom_boxplot(width=0.55, color="grey", alpha=0.2,fill="deepskyblue1")+
    ggtitle(paste("ASF1A/ASF1B knockout LFCs in",group))
  ggsave(filename=paste('ASF1AB_LFC_',group,'.png',sep=""),path =gene_pairs_path,plot=ASF1AB_plot, width=6,height=3.5)
}

