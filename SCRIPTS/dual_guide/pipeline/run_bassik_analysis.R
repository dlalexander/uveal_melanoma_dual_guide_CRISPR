suppressPackageStartupMessages(suppressWarnings(library(optparse)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(reshape2)))

############################################################
# OPTIONS                                                  # 
############################################################

option_list = list(
  make_option(c("-f", "--fc"), type = "character",
               help = "fold change matrix file", metavar = "character"),
  make_option(c("-y", "--y12"), type = "character",
               help = "predicted vs observed file", metavar = "character"),
  make_option(c("-m", "--missing"), type = "character",
               help = "missing data file", metavar = "character"),
  make_option(c("-s", "--samples"), type = "character",
              help = "full path to sample mapping", metavar = "character"),
  make_option(c("-c", "--control_guides"), type = "character",
               help = "comma-separated list of control guide identifiers", metavar = "character"),
  make_option(c("-t", "--tsv"), type = "character",
               help = "full path to TSV output directory", metavar = "character"),
  make_option(c("-r", "--rds"), type = "character",
              help = "full path to RDS output directory", metavar = "character")
);

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

############################################################
# GENERAL                                                  #
############################################################

# Check file exists
check_file_exists <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(sprintf("File does not exist: %s", filepath))
  }
}

# Check directory exists
check_dir_exists <- function(dirpath) {
  if (!dir.exists(dirpath)) {
    stop(sprintf("Directory does not exist: %s", dirpath))
  }
}

# Check input files exist
for (opt_file in c('samples', 'fc', 'y12', 'missing')) {
  check_file_exists(opt[[opt_file]])
}

# Check output directories exist
check_dir_exists(opt$tsv)
check_dir_exists(opt$rds)

# Set annotation column names
annotation_colnames <- c(
  'id', 'sgrna_ids', 'sgrna_seqs', 'gene_pair_id',
  'sorted_gene_pair', 'targetA', 'targetB', 
  'sgrna_symbols', 'sgrna_symbol_A', 'sgrna_symbol_B', 
  'sgrna_libraries', 'sgrna_group', 'guide_type', 
  'guide_orientation', 'singles_target_gene')

############################################################
# Prepare control guides                                   #
############################################################

# If no control guides provided, exit
if (is.null(opt$control_guides)) {
  stop("No control guides provided")
}

# Split control guides into vector
control_guides <- str_split(opt$control_guides, ",", simplify = T) %>% as.vector()

############################################################
# Sample mapping                                           #
############################################################

check_file_exists(opt$samples)

print(paste("Reading sample annotations from:", opt$samples))

# Read in sample mapping
sample_mapping <- read.delim(opt$samples, header = T, sep = "\t")

############################################################
# Read in Predicted vs observed                            #
############################################################

print(paste("Reading predicted vs observed:", opt$y12))

# Read in prediced vs observed LFCs per guide pair (all samples)
pred_vs_obs_y12 <- read.delim(opt$y12, sep = "\t", header = T, check.names = F)

# Show first couple of lines
print('Predicted vs observed LFCs:')
print(head(pred_vs_obs_y12, 3))

############################################################
# Read in missing data                                     #
############################################################

print(paste("Reading missing data:", opt$missing))

# Read in missing data - guides removed from analysis, typically because they were filtered during pre-processing (all samples)
missing_data_obs_y12 <- read.delim(opt$missing, sep = "\t", header = T, check.names = F)

# Show first couple of lines
print('Missing data:')
print(head(missing_data_obs_y12, 3))

############################################################
# Read in log fold changes                                 #
############################################################

print(paste("Reading LFCs:", opt$fc))

# Read in log fold changes (all samples)
fc <- read.delim(opt$fc, sep = "\t", header = T, check.names = F)

# Gather fold changes by sample
fc <- fc %>%
  gather(sample, lfc, -all_of(annotation_colnames))

# Show first couple of lines
print('Gathered LFCs:')
print(head(fc, 3))

############################################################
# Check samples between data sets                          #
############################################################

print('Checking samples...')

# Take all samples found in pred_vs_obs_y12
samples <- pred_vs_obs_y12$sample %>% unique()

# Check all samples are in fold change matrix
missing_samples <- setdiff(samples, fc$sample)
if (length(missing_samples) > 0) {
  stop(paste("Samples missing in fold change matrix:", paste(missing_samples, collapse = ', ')))
}

# Filter fold changes for selected samples
if (length(intersect(samples, fc$sample)) != length(samples)) {
  print(paste("Number of samples in fold changes before filtering:", length(unique(fc$sample))))
  fc <- fc %>% filter(sample %in% samples)
  print(paste("Number of samples in fold changes after filtering:", length(unique(fc$sample))))
}

############################################################
# Get samples per cell line                                #
############################################################

print('Identifying samples per cell line...')

# Get selected samples from mapping
sample_mapping.filt <- sample_mapping %>%
  filter(sample_label %in% samples)

# Get samples per cell line
sample_list <- list()
for (i in unique(sample_mapping.filt$qc_group)) {
  sample_list[[i]] <- sample_mapping.filt %>%
    filter(qc_group == i) %>%
    pull(sample_label)
}

############################################################
# Classify fold changes by guide type                      #
############################################################

print("Classifying fold changes by guide type (guide_class)...")

# Add a guide classification based on the guide type
fc.prepared <- fc %>%
  mutate('guide_class' = case_when(
    guide_type == 'safe_control|safe_control' ~ 'control',
    guide_type == 'gene|safe_control' | guide_type == 'safe_control|gene' ~ 'single',
    guide_type == 'gene|gene' ~ 'double',
    TRUE ~ "unknown"))

# Print summary
print(fc.prepared %>% group_by(sample, guide_class) %>% count() %>% spread(guide_class, n) %>% as.data.frame())

############################################################
# Compute category medians by sample                       #
############################################################

print("Calculating median LFC by sample and guide_class...")

# Get median LFC by sample and guide_class
median_fcs_by_sample <- fc.prepared %>%
  group_by(sample, guide_class) %>%
  summarise(median_fc = median(lfc, na.rm = TRUE), .groups = 'keep') %>%
  spread(guide_class, median_fc) %>%
  ungroup() %>%
  rename('median_control_fc' = 'control', 'median_double_fc' = 'double', 'median_single_fc' = 'single')

# Print median LFC by sample and guide class
print(median_fcs_by_sample)

############################################################
# Merge median LFCs with pred_vs_obs_y12                   #
############################################################

print("Merging median LFCs back into pred_vs_obs_y12...")

# Merge fold changes back into pred_vs_obs_y12
pred_vs_obs_y12_with_gene <- pred_vs_obs_y12 %>% left_join(median_fcs_by_sample, by = 'sample')

############################################################
# Calculate gammas                                         #
############################################################

print("Calculating gammas...")

# Calculate gammas by taking the naive fold-change and subtracting off a median fold-change
# Here we use the median of the single guides (gene|safe_control or safe_control|gene)
pred_vs_obs_y12_with_gene <- pred_vs_obs_y12_with_gene %>% 
  rowwise() %>%
  mutate('gamma_12' = obs_y12 - median_single_fc,
         'gamma_1'  = y1 - median_single_fc,
         'gamma_2'  = y2 - median_single_fc,
         'predicted_gamma_12' = gamma_1 + gamma_2)

# Print first few rows
print(head(pred_vs_obs_y12_with_gene, 3))

############################################################
#  Normalised GIs per cell line                            #
############################################################

pred_vs_obs_with_ngis.all <- data.frame()

for (cl in names(sample_list)) {
  
  print(paste('Normalising GIs per cell line for:', paste(sample_list[[cl]], collapse = ", ")))
  
  # Create empty data frame for all samples
  pred_vs_obs_with_ngis <- data.frame()

  # Loop over each sample and calculate normalised GIs
  for (thesample in sample_list[[cl]]) {
    print(paste("Getting normalised GIs for:", thesample))
    
    print("Extracting sample from pred_vs_obs_y12_with_gene (thesubset)...")
    # Note: removed sample == thesample & !(is.na(gamma_12) | is.na(predicted_gamma_12)) as we have no NA values
    thesubset <- pred_vs_obs_y12_with_gene %>% filter(sample == thesample)
    
    print("Adding row numbers to thesubset...")
    thesubset$row_num <- seq(1:nrow(thesubset))
    print(head(thesubset, 3) %>% as.data.frame())
    
    print("Performing LOESS regression...")
    loessMod50 <- loess(gamma_12 ~ predicted_gamma_12, data = thesubset, span = 1.5)
    
    print("Building residuals table...")
    obsvspredresid <- data.frame(loesspredicted = loessMod50$fitted, # loess predicted gamma_12 from observed gamma_12
                                 x = thesubset$predicted_gamma_12,   # Predicted gamma_12 from observed single gamma1 and gamma2
                                 y = thesubset$gamma_12,             # Observed gamma_12
                                 resid = loessMod50$residuals)       # Observed gamma_12 -loess predicted gamma_12
    obsvspredresid$row_num <- seq(1:nrow(obsvspredresid))            # Add row numbers
    print(head(obsvspredresid, 3) %>% as.data.frame())

    print("Add residuals to thesubset...")
    thesubset <- thesubset %>% left_join(obsvspredresid, by="row_num")
    print(head(thesubset, 3) %>% as.data.frame())
    
    print("Sorting thesubset by predicted_gamma_12...")
    thesubset_sorted <- thesubset %>% arrange(predicted_gamma_12)
    print(head(thesubset_sorted, 3) %>% as.data.frame())
    
    print("Calculating number of batches for computing batch variance...")
    bin_size <- 200
    num_batches <- trunc(nrow(thesubset) / bin_size)  + 1
    print(paste('Bin size:', bin_size))
    print(paste('Number of batches:', num_batches))
    
    print("Add batch number to sorted thesubset...")
    thesubset_sorted$resid_variance_batch <- rep(seq(1:num_batches), each = bin_size, length.out = nrow(thesubset_sorted))
    print(head(thesubset_sorted, 3) %>% as.data.frame())
    print(table(thesubset_sorted$resid_variance_batch))
    
    print("Computing batch variance for thesubset...")
    subset_resid_batch_var <- thesubset_sorted %>%
      select(resid_variance_batch, resid) %>%
      group_by(resid_variance_batch) %>%
      summarise(resid_var = var(resid)) %>%
      ungroup()
    print(head(subset_resid_batch_var, 3) %>% as.data.frame())
    
    print("Merging batch variance back into thesubset...")
    thesubset_sorted <- thesubset_sorted %>%
      left_join(subset_resid_batch_var, by = "resid_variance_batch")
    print(head(thesubset_sorted, 3) %>% as.data.frame())
    
    print("Calculating normalised gi...")
    thesubset_sorted <- thesubset_sorted %>%
      rowwise() %>%
      mutate('norm_gi' = resid / sqrt(resid_var))
    print(head(thesubset_sorted, 3) %>% as.data.frame())
    
    print("Expanding annotation of thesubset_sorted...")
    thesubset_sorted <- thesubset_sorted %>%
      left_join(fc %>% select(all_of(annotation_colnames)) %>% unique(), by = c('id', 'gene_pair_id')) %>%
      relocate(all_of(annotation_colnames), .before = 'sample')
    
    print("Adding thesubset_sorted to pred_vs_obs_with_ngis...")
    # accumulate all this together into the return frame
    if (nrow(pred_vs_obs_with_ngis) == 0) {
      pred_vs_obs_with_ngis <- thesubset_sorted
    } else {
      pred_vs_obs_with_ngis <- rbind(pred_vs_obs_with_ngis, thesubset_sorted)
    }
    print(paste("Number of rows in pred_vs_obs_with_ngis:", nrow(pred_vs_obs_with_ngis)))
    print(head(pred_vs_obs_with_ngis, 3) %>% as.data.frame())
  }
  
  # Ungroup pred_vs_obs_with_ngis
  pred_vs_obs_with_ngis <- pred_vs_obs_with_ngis %>% ungroup()
  
  # Add to main data frame
  if (nrow(pred_vs_obs_with_ngis.all) == 0) {
    pred_vs_obs_with_ngis.all <- pred_vs_obs_with_ngis
  } else {
    pred_vs_obs_with_ngis.all <- bind_rows(pred_vs_obs_with_ngis.all, pred_vs_obs_with_ngis)
  }
}

############################################################
# Get t by representation by gene                          #
############################################################

summary_by_gene.all <- data.frame()

for (cl in names(sample_list)) {
  
  print(paste('Getting t score and FDR by gene for:', paste(sample_list[[cl]], collapse = ", ")))

  print('Summarising normalised residuals by gene pair...')
  
  # Set up temporary data frame
  pred_vs_obs_with_ngis.tmp <- data.frame()
  pred_vs_obs_with_ngis.tmp <- pred_vs_obs_with_ngis.all %>% filter(sample %in% sample_list[[cl]])
    
  # Summarise by gene pair
  summary_by_gene <- pred_vs_obs_with_ngis.tmp %>%
    group_by(gene_pair_id) %>%
    summarise('mean_norm_gi' = mean(norm_gi), 
              'median_norm_gi' = median(norm_gi), # Median of normalised GIs of all double-sgRNAs targeting the gene pair (Bassik Uexp)
              'variance_norm_gi' = var(norm_gi), # Variance of normalised GIs of all double-sgRNAs targeting the gene pair (Bassik Vexp)
              'nexp' = n(), .groups = 'keep') # Number of all double-sgRNAs targeting the gene pair (Bassik Nexp)
  print(head(summary_by_gene, 3) %>% as.data.frame())
  
  print('Summarising across all gene pairs...')
  # Summarise corrected residuals across all samples (i.e. for all gene pairs)
  summary_across_all_gene_pairs <- pred_vs_obs_with_ngis.tmp %>%
    ungroup() %>%
    summarise('median_norm_gi_all' = median(norm_gi), # Median of normalised GIs of all double-sgRNAs
              'variance_norm_gi_all' = var(norm_gi)) # Variance of normalised GIs of all double-sgRNAs
  print(summary_across_all_gene_pairs %>% as.data.frame())
  
  print('Calculating number of guide pairs by gene...')
  # Get the number of double-sgRNAs per gene pair
  number_grnas_per_gene <- pred_vs_obs_with_ngis.tmp %>%
    group_by(gene_pair_id) %>%
    summarise(num_grnas_per_gene_pair = n())
  print(head(number_grnas_per_gene, 3) %>% as.data.frame())
  
  print('Calculating mean number of guide pairs across whole dataset...')
  # Get the mean number of double-sgRNAs per gene pair across the whole data set
  mean_grnas_per_gene <- number_grnas_per_gene %>%
    summarise(mean_grna_count_all = mean(num_grnas_per_gene_pair))
  print(mean_grnas_per_gene)

  print('Adding mean guide pairs to data frame...')
  # Add the mean gRNA count onto the data frame
  summary_by_gene <- summary_by_gene %>% cbind(mean_grnas_per_gene)
  print(head(summary_by_gene, 3) %>% as.data.frame())
  
  print('Add broader stats back into gene pair summary...')
  # Add back together
  summary_by_gene <- summary_by_gene %>% cbind(summary_across_all_gene_pairs)
  print(head(summary_by_gene, 3) %>% as.data.frame())
  
  print('Calculate gene level t-statistics...')
  summary_by_gene <- summary_by_gene %>%
    mutate('svar' = variance_norm_gi / (nexp - 1) + variance_norm_gi_all / (mean_grna_count_all - 1)) %>%
    mutate('gi_t_score' =  (median_norm_gi - median_norm_gi_all) / sqrt(svar)) %>%
    mutate('plow' = pnorm(gi_t_score)) %>%
    mutate('phigh' = 1 - pnorm(gi_t_score)) %>%
    mutate('fdr' = p.adjust(plow, method = "BH"))
  print(head(summary_by_gene, 3) %>% as.data.frame())
  
  print('Sort by FDR and median_norm_gi...')
  summary_by_gene <- summary_by_gene %>% arrange(fdr, median_norm_gi)
  print(head(summary_by_gene, 3) %>% as.data.frame())
  
  # Add to main data frame
  if (nrow(summary_by_gene.all) == 0) {
    summary_by_gene.all <- summary_by_gene %>% mutate('qc_group' = cl)
  } else {
    summary_by_gene.all <- bind_rows(summary_by_gene.all, summary_by_gene %>% mutate('qc_group' = cl))
  }
}

############################################################
# Get t by representation by guide                         #
############################################################

summary_by_sgrna.all <- data.frame()

for (cl in names(sample_list)) {
  
  print(paste('Getting t score and FDR by guide for:', paste(sample_list[[cl]], collapse = ", ")))
  
  # Set up temporary data frame
  pred_vs_obs_with_ngis.tmp <- data.frame()
  pred_vs_obs_with_ngis.tmp <- pred_vs_obs_with_ngis.all %>% filter(sample %in% sample_list[[cl]])

  print('Summarising normalised residuals by gene pair...')
  summary_by_sgrna <- pred_vs_obs_with_ngis.tmp %>%
    group_by(id, gene_pair_id) %>%
    summarise('mean_norm_gi' = mean(norm_gi), 
              'median_norm_gi' = median(norm_gi), # Median of normalised GIs of all double-sgRNAs targeting the guide pair (Bassik Uexp)
              'variance_norm_gi' = var(norm_gi), # Variance of normalised GIs of all double-sgRNAs targeting the guide pair (Bassik Vexp)
              'number_of_replicates' = n(), .groups = 'keep') # Number of all double-sgRNAs targeting the guide pair (Bassik Nexp)
  print(head(summary_by_sgrna, 3) %>% as.data.frame())

  print('Summarising across all guide...')
  summary_across_all_guide_pairs <- pred_vs_obs_with_ngis.tmp %>%
    ungroup() %>%
    summarise('median_norm_gi_all' = median(norm_gi), # Median of normalised GIs of all double-sgRNAs
              'variance_norm_gi_all' = var(norm_gi)) # Variance of normalised GIs of all double-sgRNAs
  print(summary_across_all_guide_pairs %>% as.data.frame())
  
  print('Calculating number of replicates per guide pair...')
  number_replicates_per_guide <- pred_vs_obs_with_ngis.tmp %>%
    group_by(id) %>%
    summarise(num_replicates_per_guide_pair = n())
  print(head(number_replicates_per_guide, 3) %>% as.data.frame())
  
  print('Calculating mean number of guide pairs across whole dataset...')
  mean_replicates_per_guide <- number_replicates_per_guide %>%
    summarise(mean_replicate_count_all = mean(num_replicates_per_guide_pair))
  print(mean_replicates_per_guide)

  print('Adding mean replicates to data frame...')
  summary_by_sgrna <- summary_by_sgrna %>% cbind(mean_replicates_per_guide)
  print(head(summary_by_sgrna, 3) %>% as.data.frame())
  
  print('Add broader stats back into guide pair summary...')
  summary_by_sgrna <- summary_by_sgrna %>% cbind(summary_across_all_guide_pairs)
  print(head(summary_by_sgrna, 3) %>% as.data.frame())
  
  print('Calculate gene level t-statistics...')
  summary_by_sgrna <- summary_by_sgrna %>%
    mutate('svar' = variance_norm_gi / (number_of_replicates - 1) + variance_norm_gi_all / (mean_replicate_count_all - 1)) %>%
    mutate('gi_t_score' =  (median_norm_gi - median_norm_gi_all) / sqrt(svar)) %>%
    mutate('plow' = pnorm(gi_t_score)) %>%
    mutate('phigh' = 1 - pnorm(gi_t_score)) %>%
    mutate('fdr' = p.adjust(plow, method = "BH"))
  print(head(summary_by_sgrna, 3) %>% as.data.frame())
  
  print('Sort by FDR and median_norm_gi...')
  summary_by_sgrna <- summary_by_sgrna %>% arrange(fdr, median_norm_gi)
  print(head(summary_by_sgrna, 3) %>% as.data.frame())
  
  # Add to main data frame
  if (nrow(summary_by_sgrna.all) == 0) {
    summary_by_sgrna.all <- summary_by_sgrna %>% mutate('qc_group' = cl)
  } else {
    summary_by_sgrna.all <- bind_rows(summary_by_sgrna.all, summary_by_sgrna %>% mutate('qc_group' = cl))
  }
}

############################################################
# Add number of significant guides to summary by gene      #
############################################################

print('Getting number of guides per gene pair...')
significant_guides_per_gene <- summary_by_sgrna.all %>%
  group_by(gene_pair_id, qc_group) %>%
  summarise('n_guide_pairs' = n(),
            'n_sig_guide_pairs' = sum(fdr < 0.05), .groups = 'keep')

print('Adding guide summaries to gene summary...')
summary_by_gene.all <- summary_by_gene.all %>%
  left_join(significant_guides_per_gene, by = c('gene_pair_id', 'qc_group')) %>%
  relocate(n_guide_pairs, n_sig_guide_pairs, .after = 'nexp') %>%
  mutate('n_replicates' = nexp / n_guide_pairs, .after = 'nexp')

############################################################
# Add sample metadata                                      #
############################################################

print('Adding cancer type and gene pair annotations...')

print('pred_vs_obs_with_ngis...')

pred_vs_obs_with_ngis.all <- pred_vs_obs_with_ngis.all %>%
  left_join(sample_mapping.filt %>% select(sample_label, qc_group, cancer_type) %>% unique(), by = c('sample' = 'sample_label')) %>%
  relocate(sorted_gene_pair, targetA, targetB, .after = 'gene_pair_id') %>%
  relocate(qc_group, cancer_type, .after = 'sample') %>%
  unique()

print('summary_by_gene...')

summary_by_gene.all <- summary_by_gene.all %>%
  left_join(sample_mapping.filt %>% select(qc_group, cancer_type) %>% unique(), by = c('qc_group')) %>%
  left_join(fc %>% select(gene_pair_id, sorted_gene_pair, targetA, targetB) %>% unique(), by = 'gene_pair_id') %>%
  relocate(sorted_gene_pair, targetA, targetB, .after = 'gene_pair_id') %>%
  relocate(qc_group, cancer_type, .after = 'gene_pair_id') %>%
  unique()

print('summary_by_sgrna...')

summary_by_sgrna.all <- summary_by_sgrna.all %>%
  left_join(sample_mapping.filt %>% select(qc_group, cancer_type) %>% unique(), by = c('qc_group')) %>%
  left_join(fc %>% select(all_of(annotation_colnames)) %>% unique(), by = c('id', 'gene_pair_id')) %>%
  relocate(sgrna_ids, sgrna_seqs, sorted_gene_pair, targetA, targetB, 
           sgrna_symbols, sgrna_symbol_A, sgrna_symbol_B, sgrna_libraries, 
           sgrna_group, guide_orientation, .after = 'gene_pair_id') %>%
  relocate(qc_group, cancer_type, .after = 'guide_orientation') %>%
  unique()

############################################################
# TSV OUTPUTS                                              #
############################################################

print('Saving pred_vs_obs_with_ngis (TSV)...')
pred_vs_obs_with_ngis.filepath <- file.path(opt$tsv, 'pred_vs_obs_with_ngis.tsv')
write.table(pred_vs_obs_with_ngis.all, pred_vs_obs_with_ngis.filepath, row.names = F, sep = "\t", quote = F)
print(paste('pred_vs_obs_with_ngis written to:', pred_vs_obs_with_ngis.filepath))

print('Saving summary_by_gene (TSV)...')
summary_by_gene.filepath <- file.path(opt$tsv, 'gene_sample_t_scores.tsv')
write.table(summary_by_gene.all, summary_by_gene.filepath, row.names = F, sep = "\t", quote = F)
print(paste('summary_by_gene written to:', summary_by_gene.filepath))

print('Saving summary_by_sgrna (TSV)...')
summary_by_sgrna.filepath <- file.path(opt$tsv, 'sgrna_sample_t_scores.tsv')
write.table(summary_by_sgrna.all, summary_by_sgrna.filepath, row.names = F, sep = "\t", quote = F)
print(paste('summary_by_sgrna written to:', summary_by_sgrna.filepath))

############################################################
# RDS OUTPUTS                                              #
############################################################

print('Saving pred_vs_obs_with_ngis (RDS)...')
pred_vs_obs_with_ngis.rds.filepath <- file.path(opt$rds, 'pred_vs_obs_with_ngis.rds')
saveRDS(pred_vs_obs_with_ngis.all, pred_vs_obs_with_ngis.rds.filepath)
print(paste('pred_vs_obs_with_ngis RDS written to:', pred_vs_obs_with_ngis.rds.filepath))

print('Saving summary_by_gene (RDS)...')
summary_by_gene.rds.filepath <- file.path(opt$rds, 'summary_by_gene.rds')
saveRDS(summary_by_gene.all, summary_by_gene.rds.filepath)
print(paste('summary_by_gene RDS written to:', summary_by_gene.rds.filepath))

print('Saving summary_by_sgrna (RDS)...')
summary_by_sgrna.rds.filepath <- file.path(opt$rds, 'summary_by_sgrna.rds')
saveRDS(summary_by_sgrna.all, summary_by_sgrna.rds.filepath)
print(paste('summary_by_sgrna RDS written to:', summary_by_sgrna.rds.filepath))

print('Saving workspace...')
workspace.filepath <- file.path(opt$rds, 'bassik_workspace.RData')
save.image(file = workspace.filepath)
print(paste('Workspace written to:', workspace.filepath))