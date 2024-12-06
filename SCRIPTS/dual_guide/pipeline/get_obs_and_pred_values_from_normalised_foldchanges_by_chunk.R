suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))

############################################################
# OPTIONS                                                  #
############################################################

option_list = list(
  make_option(c("-f", "--fc"), type = "character",
              help = "fold change matrix file", metavar = "character"),
  make_option(c("-d", "--doubles_guide_matrix"), type = "character", default = 'norm',
              help = "doubles guide matrix file", metavar = "character"),
  make_option(c("-i", "--chunk_index"), type="integer", default = NULL,
              help="which library chunk to process", metavar="integer"),
  make_option(c("-n", "--chunk_size"), type="integer", default=NULL,
    			    help="library chunk size", metavar="integer"),
  make_option(c("-o", "--out"), type = "character", default = '.',
              help = "output directory [Default: . ]", metavar = "character")
);

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

############################################################
# VALIDATION                                               #
############################################################

# Fold change matrix file
if (is.null(opt$fc)) {
  print_help(opt_parser)
  stop("Please provide a fold change matrix file", call.=FALSE)
}

if (!file.exists(opt$fc)) {
  print_help(opt_parser)
  stop(paste("Fold change matrix file does not exist:", opt$fc), call.=FALSE)
}

# Doubles guide matrix file
if (is.null(opt$doubles_guide_matrix)) {
  print_help(opt_parser)
  stop("Please provide a doubles guide matrix file", call.=FALSE)
}

if (!file.exists(opt$doubles_guide_matrix)) {
  print_help(opt_parser)
  stop(paste("Doubles guide matrix file does not exist:", opt$doubles_guide_matrix), call.=FALSE)
}

############################################################
# FUNCTIONS                                                #
############################################################

get_pred_vs_obs_y12 <- function(sn, gm, fcs) {
  # Create empty data frame for pred vs obs values
  pred_vs_obs_y12 <- data.frame(id = character(), gene_pair_id = character(), sample = character(),
                                obs_y12 = numeric(), pred_y12 = numeric(),
                                y1 = numeric(), y2 = numeric())

  # Create empty data frame for missing guides per sample
  missing_data <- data.frame(id = character(), gene_pair_id = character(), sample = character())

  for (i in 1:nrow(gm)) {
    # Get dual and respective single guide ids
    g12 <- as.vector(gm$id[i])
    g1  <- as.vector(gm$g1[i])
    g2  <- as.vector(gm$g2[i])
    gp  <- as.vector(gm$gene_pair_id[i])

    # Get observed fold changes for left guide (single)
    y1 <- fcs %>%
            filter(id == g1 & gene_pair_id == gp & sample == sn) %>%
            select(fc) %>%
            unlist() %>% as.vector()

    # Get observed fold changes for right guide (single)
    y2 <- fcs %>%
            filter(id == g2 & gene_pair_id == gp & sample == sn) %>%
            select(fc) %>%
            unlist() %>% as.vector()

    # Get observed fold changes for dual guide (dual)
    obs_y12 <- fcs %>%
                filter(id == g12 & gene_pair_id == gp & sample == sn) %>%
                select(fc) %>%
                unlist() %>% as.vector()

    # Check there are fold changes for the dual guide and both related single guides
    if ( length( y1 ) == 1 && length( y2 ) == 1 && length( obs_y12 ) == 1 ) {
      pred_y12 <- y1 + y2
      pred_vs_obs_y12.tmp <- data.frame(  'id' = g12, 'gene_pair_id' = gp, 'sample' = sn,
                                          'obs_y12' = obs_y12, 'pred_y12' = pred_y12,
                                          'y1' = y1, 'y2' = y2 )
      pred_vs_obs_y12 <- rbind(pred_vs_obs_y12, pred_vs_obs_y12.tmp)
    } else {
      # If one of the guides isn't present (filtered), then add to missing data
      missing_data.tmp <- data.frame('id' = g12, 'gene_pair_id' = gp, 'sample' = sn)
      missing_data <- rbind(missing_data, missing_data.tmp)
    }
    if(i %% 200 == 0) { print(i) }
  }
  return(list("pred_vs_obs_y12" = pred_vs_obs_y12, "missing_data" = missing_data))
}

get_pred_vs_obs_y12_for_all_samples <- function(samples, gm, fcs) {
  pred_vs_obs_y12 <- list()
  missing_data    <- list()
  sample.results  <- list()

  for (sn in samples) {
    print(paste("Processing:", sn))
    sample.results <- get_pred_vs_obs_y12(sn, gm, fcs)
    pred_vs_obs_y12[[sn]] <- sample.results[['pred_vs_obs_y12']]
    missing_data[[sn]] <- sample.results[['missing_data']]
  }
  return(list('pred_vs_obs_y12' = pred_vs_obs_y12, 'missing_data' = missing_data))
}

generate_output_filepath <- function( filename, outdir = '.' ) {
  output_filepath <- file.path( outdir, filename )
  return( output_filepath )
}

############################################################
# MAIN SCRIPT                                              #
############################################################

# Read in dual guide matrix (maps dual guide to its single components)
doubles_guide_matrix <- read.delim(file = opt$doubles_guide_matrix , sep = "\t", header = T, check.names = F)
fc <- read.delim(file = opt$fc, sep = "\t", header = T, check.names = F)

# Set annotation column names
annotation_colnames <- c(
  'id', 'sgrna_ids', 'sgrna_seqs', 'gene_pair_id',
  'sorted_gene_pair', 'targetA', 'targetB', 
  'sgrna_symbols', 'sgrna_symbol_A', 'sgrna_symbol_B', 
  'sgrna_libraries', 'sgrna_group', 'guide_type', 
  'guide_orientation', 'singles_target_gene')

# Narrow the FC matrix
fc.narrow <- fc %>% gather(sample, fc, -all_of(annotation_colnames))

# Get sample names
samples <- unique(fc.narrow$sample)

# Get chunks for double guide matrix
chunk_list <- split(c(1:nrow(doubles_guide_matrix)), ceiling(seq_along(c(1:nrow(doubles_guide_matrix))) / opt$chunk_size))
doubles_chunk_indexes <- unlist(chunk_list[ opt$chunk_index ])
chunk_doubles_guide_matrix <- doubles_guide_matrix[doubles_chunk_indexes, ]

# Get observed and predicted fold changes for doubles (per sample lists)
results <- get_pred_vs_obs_y12_for_all_samples(samples, chunk_doubles_guide_matrix, fc.narrow)

# Bring together results into single dataframe
pred_vs_obs_y12 <- data.frame()
missing_data <- data.frame()

for (sn in samples) {
  print(paste("Merging results for:", sn))
  if (nrow(pred_vs_obs_y12) == 0) {
    pred_vs_obs_y12 <- results[['pred_vs_obs_y12']][[sn]]
    missing_data <- results[['missing_data']][[sn]]
  } else {
    pred_vs_obs_y12 <- rbind(pred_vs_obs_y12, as.data.frame(results[['pred_vs_obs_y12']][[sn]]))
    missing_data <- rbind(missing_data, results[['missing_data']][[sn]])
  }
}

pred_vs_obs_y12_filename = generate_output_filepath(paste("pred_vs_obs_y12", opt$chunk_index, "tsv", sep = "."), opt$out)
write.table(pred_vs_obs_y12, file = pred_vs_obs_y12_filename, row.names = F, sep = "\t", quote = F)

missing_data_filename = generate_output_filepath(paste( "missing_data", opt$chunk_index, "tsv", sep = "."), opt$out)
write.table(missing_data, file = missing_data_filename, row.names = F, sep = "\t", quote = F)
