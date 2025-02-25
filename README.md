# Dual gRNA CRISPRko screens in uveal melanoma 

## Overview

CRISPRko screens using a bespoke paired gRNA library were performed in 10 human uveal melanoma cell lines.

This project aims to address the following questions:
1. What gene pairs are essential for cell viability in uveal melanoma and are potential therapeutic targets?
2. Are the gene essentialities observed in uveal melanoma unique or shared with cutaneous melanoma?


## Experimental design 

10 uveal melanoma cell lines were engineered to express Cas9 and had at least 85% Cas9 activity. CRISPRko screens were performed by using a bespoke paired gRNA library. 

Cell lines were transduced in technical triplicate at a MOI of 0.3, puromycin selected for ~30% transduction efficiency and maintained at a representation of 1000x library coverage. Endpoints shown in table below consisting of D28 after transduction +/- timepoint after approximately 20 cell doublings.  Cells were then pelleted, genomic DNA extracted, followed by gRNA PCR amplification and sequencing. Batch B CRISPR screens were performed by the CGaP core facility. 

| Cell line | Batch | First endpoint | Second endpoint |
|-----------|-------|----------------|-----------------|
| Mel202 control | A | D8 | None (non-Cas9 control) |
| Mel202-Cas9 | A | D28 | None (20 cell doublings by D28) | 
| Mel270-Cas9 | A | D28 | D35 |
| Mel285-Cas9 | A | D28 | D49 |
| MP46-Cas9 | A | D28 | D64 | 
| OMM2.5-Cas9 | A | D28 | D42 | 
| OMM2.3 control | B | D8 | None (non-Cas9 control) |
| OMM2.3-Cas9 | B | D28 | D35 | 
| OMM1-Cas9 | B | D28 | D35 | 
| 92.1-Cas9 | B | D28 | D42 | 
| MP38-Cas9 | B | D28 | D65 | 
| MP41-Cas9 | B | D28 | D32 | 


# Analysis steps
## Source

Analysis pipeline is adapted from a pipeline written by [Victoria Offord](https://github.com/vaofford) for the Adams Group at the Sanger Institue [Team 113 Github](https://github.com/team113sanger), available [here](https://github.com/team113sanger/harle_dgCRISPR_paralogs).

## Installation and dependencies

For installation instructions, system and package dependencies see [INSTALL_README.md](https://github.com/dlalexander/uveal_melanoma_dual_guide_CRISPR/blob/master/INSTALL_README.md)

## Repository
For almost all steps, the `REPO_PATH` variable will need to be set. This will depend on where you have cloned the repository.

```
export REPO_PATH='/your/repository/path'
```

## Generating guide abundances (counts)

### Preparing the library

Convert original library into pyCROQUET format:
```
awk -F"\t" 'BEGIN{
        OFS="\t"; 
        print "##library-type: dual"; 
        print "##library-name: dgCRISPR_paralogs"; 
        print "#id\tsgrna_ids\tsgrna_seqs\tgene_pair_id\tsgrna_symbols\tsgrna_libraries\tcustom_annotation";
    } NR > 1 && $6 != "Jen Uveal" && $13 != "Jen Uveal" {
        gsub("_", "|", $15)
        id=$3"|"$10"|"$15
        sgrna_ids=$3"|"$10
        sgrna_seqs=$2"|"$9
        gene_pair_id=$15
        sgrna_symbols=$1"|"$8
        sgrna_libraries=$5"|"$12
        sgrna_sources=$6"|"$13
        print id, sgrna_ids, sgrna_seqs, gene_pair_id, sgrna_symbols, sgrna_libraries, sgrna_sources;
    }' ${REPO_PATH}/METADATA/libraries/paired_lib_original.tsv > ${REPO_PATH}/METADATA/libraries/paired_lib_original.pyCROQUET.tsv
```

### Quantification per CRAM

Run pyCROQUET across all CRAM files in `DATA/CRAM`. Remembering to set the `REPO_PATH` and update the jobscript header to match the number of expected CRAM files and the top level directory for the logs.

```
bsub < ${REPO_PATH}/SCRIPTS/pyCROQUET/jobscript.sh
```

This runs an LSF job array where each CRAM gets quantified independently. Logs for these can be found in `LOGS/pyCROQUET`.

To clean up afterwards:

```
cat ${REPO_PATH}/LOGS/pyCROQUET/pycroquet.[0-9]*.e > ${REPO_PATH}/LOGS/pyCROQUET/pycroquet.e
cat ${REPO_PATH}/LOGS/pyCROQUET/pycroquet.[0-9]*.o > ${REPO_PATH}/LOGS/pyCROQUET/pycroquet.o
find ${REPO_PATH}/LOGS/pyCROQUET -name 'pycroquet.[0-9]*.e' -delete
find ${REPO_PATH}/LOGS/pyCROQUET -name 'pycroquet.[0-9]*.o' -delete
```

## Preprocessing

There are several pre-processing steps which are covered by `SCRIPTS/preprocessing/run_preprocessing_rscripts.sh`. Each step utilises a script found in `SCRIPTS/preprocessing/subscripts` whose outputs and general descriptions are:

| Step | Description | Script | Inputs | Outputs |
| --- | --- | --- | --- | --- |
| 0 | Prepare libraries for downstream analyses | `0_prepare_libraries.R` | Sample metadata: <br/> `METADATA/5429_sample_annotations.tsv` <br/> pyCROQUET library: <br /> `METADATA/libraries/paired_lib_original.pyCROQUET.tsv` | Expanded library annotations: <br/> `METADATA/libraries/expanded_library.tsv` <br/> Dual guide matrix: <br/> `DATA/dual_guide/dual_guide_matrix.tsv'` |
| 1 | Convert per-lane pyCROQUET counts into a sample count matrix with annotations | `1_merge_counts_by_sample.R` | Sample metadata: <br/> `METADATA/5429_sample_annotations.tsv` <br/> Expanded library: <br /> `METADATA/libraries/expanded_library.tsv` | Count matrix: <br/> `DATA/preprocessing/count_matrix.tsv` |
| 2 | Remove user-defined guides by id | `2_remove_user_defined_guides.R` | Count matrix: <br/> `DATA/preprocessing/count_matrix.tsv` <br/> Guides to remove: <br/> `METADATA/guides_ids_to_remove.txt` | Count matrix with user-defined guides removed: <br/> `DATA/preprocessing/count_matrix.guides_removed.tsv` |
| 3 | Normalise counts using BAGEL method | `3_bagel_normalisation.R` | Count matrix with user-defined guides removed: <br/> `DATA/preprocessing/count_matrix.guides_removed.tsv` | Normalised count matrix: <br/> `DATA/preprocessing/count_matrix.norm.tsv` |
| 4 | Filter out guides which have low counts (< 30) in controls  | `4_filter_low_count_guides.R` | Normalised count matrix: <br/> `DATA/preprocessing/count_matrix.norm.tsv` <br/> Sample metadata: <br/> `METADATA/5429_sample_annotations.tsv` | Filtered count matrix: <br/> `DATA/preprocessing/count_matrix.norm.filt.tsv` <br/> Filtered guides: <br/> `DATA/preprocessing/filtered_guides.txt` |
| 5 | Generate LFC matrix from counts  | `5_convert_counts_to_lfcs.R` | Filtered count matrix: <br/> `DATA/preprocessing/count_matrix.norm.filt.tsv` | Unscaled LFC matrix (all samples): <br/> `DATA/preprocessing/lfc_matrix.unscaled.all.tsv` |
| 6 | Remove samples failing QC by sample name | `6_remove_samples.R` | Unscaled LFC matrix (all samples): <br/> `DATA/preprocessing/lfc_matrix.unscaled.all.tsv` | Unscaled LFC matrix (samples removed): <br/> `DATA/preprocessing/lfc_matrix.unscaled.samples_removed.tsv` |
| 7 | Generate scaled LFC matrix | `7_scale_lfcs.R` | Unscaled LFC matrix: <br/> `DATA/preprocessing/lfc_matrix.unscaled.[all\|samples_removed].tsv` <br/> Sample metadata: <br/> `METADATA/5429_sample_annotations.tsv` | Scaled LFC matrix (all samples): <br/> `DATA/preprocessing/lfc_matrix.[scaled\|scaled_pos_zero].[all\|samples_removed].tsv` |
| 8 | Convert scaled LFC matrix into scaled count matrix  | `8_convert_scaled_lfcs_to_scaled_counts.R` | Scaled LFC matrix: <br/> `DATA/preprocessing/lfc_matrix.[scaled\|scaled_pos_zero].[all\|samples_removed].tsv` | Scaled count matrix: <br/> `DATA/preprocessing/count_matrix.[scaled\|scaled_pos_zero].[all\|samples_removed].tsv` |

All steps can be run using `SCRIPTS/preprocessing/run_preprocessing_rscripts.sh` with logs written to `LOGS/preprocessing`, output to `DATA/preprocessing` and `DATA/RDS/preprocessing`.

```
bash ${REPO_PATH}/SCRIPTS/preprocessing/run_preprocessing_rscripts.sh
```

## Single guide analysis

Single guide analyses use [C-SAR](https://github.com/cancerit/C-SAR) version 1.3.6 which converts a count matrix to an LFC matrix, runs MAGeCK and BAGEL2 and summarises the results.

The single guide analysis is run using:

```
SCRIPTS/single_guide/run_single_guide_analyses.sh
```

And when complete, the unrequired intermediate files can be removed with:

```
rm -r DATA/single_guide/*/*/*/work
rm -r DATA/single_guide/*/*/*/.nextflow
```

There are 3 datasets:

* `A` - gene|safe_control and safe_control|safe_control
* `B` - safe_control|gene and safe_control|safe_control
* `combined` - gene|safe_control, safe_control|gene and safe_control|safe_control

Each dataset has been run for all three stages:

* `unscaled`
* `scaled`
* `scaled_pos_zero`

There are several pre-processing steps which are covered by `SCRIPTS/single_guide/run_single_guide_analyses.sh`. Each step utilises a script found in `SCRIPTS/single_guide` whose outputs and general descriptions are:

| Step | Description | Script | Inputs | Outputs |
| --- | --- | --- | --- | --- |
| 0 | Prepare single libraries for singles analyses | `0_prepare_single_guide_libraries.R` | User-defined guides to remove: <br/> `METADATA/guides_ids_to_remove.txt` <br/> Expanded library: <br /> `METADATA/libraries/expanded_library.tsv` | Single library annotations with user-defined guides removed (per stage): <br/> `DATA/single_guide/combined/singles_library.[stage].tsv` |
| 1 | Prepare input files and directories for C-SAR | `1_prepare_files_and_directories_for_analysis.R` | Single guide library with user-defined guides removed: <br /> `DATA/single_guide/[stage]/singles_library.[stage].tsv` <br/> Sample metadata: <br/> `METADATA/5429_sample_annotations.tsv` <br/> Count matrix: <br/> `DATA/preprocessing/count_matrix.[stage].[dataset].tsv` <br/> Log path: <br/> `LOGS/single_guide/C-SAR/[dataset]/[stage]` <br/> Suffix: <br/> `[dataset].[stage]` <br/> Cell lines to exclude: <br/> "A-375 1000x,A-375 500x" <br/> Output directory: <br/> `DATA/single_guide/[dataset]/[stage]` | List of `bsub` commands: <br/> `SCRIPTS/single_guide/bsub_commands.[dataset].[stage].sh` <br/> Filtered singles library (only guides in counts): <br/> `DATA/single_guide/[stage]/singles_library.[stage].filt.tsv` <br/> Sample counts per cell line (per dataset/stage): <br/> `DATA/single_guide/[dataset]/[stage]/[cell line]/sample_counts.tsv` <br/> Sample manifest per cell line (per dataset/stage) <br/>`DATA/single_guide/[dataset]/[stage]/[cell line]/sample_manifest.tsv` <br/> C-SAR results per cell line (per dataset/stage) <br/>`DATA/single_guide/[dataset]/[stage]/[cell line]/results` <br/> Nextflow log per cell line (per dataset/stage) <br/>`DATA/single_guide/[dataset]/[stage]/[cell line]/.nextflow.log` |

The directory structure for C-SAR is as follows:

```
REPO_PATH
    ->DATA
        -> single_guide
            -> [dataset]
                -> [stage]
                    -> [cell_line]
                        -> results
                            -> BAGEL2
                            -> MAGeCK
                            -> ...
```

To clean up all logs and results for a fresh analysis run:

```
rm -r LOGS/single_guide/C-SAR/*/*/*.o
rm -r LOGS/single_guide/C-SAR/*/*/*.e
rm -r DATA/single_guide/*/*/*/results
rm -r DATA/single_guide/*/*/*/work
rm -r DATA/single_guide/*/*/*/.nextflow
rm -r DATA/single_guide/*/*/*/.nextflow.log
```

## Dual guide analyses

Bassik analysis requires an LFC matrix and a dual guide matrix (mapping of dual guide to its corresponding single guides) which are generated during the pre-processing steps.

We generated three LFC matrices:

* Unscaled - raw LFCs
* Scaled - LFCs scaled so that the median of the safe guides are 0 and the median of the essential guides are -1
* Scaled (pos zero) - positive scaled LFCs are set to 0 (i.e. no LFCs will be > 0)

### Generating predicted vs observed LFC matrix

Generation of the observed (dual) and predicted (sum of singles) LFCs is run as an LSF job array due to the time it takes to run as a single job.

To run job arrays to get pred vs obs for each dataset:

```
STAGE=('unscaled' 'scaled' 'scaled_pos_zero')
for stage in "${STAGE[@]}"
do
    bsub < ${REPO_PATH}/SCRIPTS/dual_guide/get_pred_vs_obs.jobarray.${stage}.sh
done    
```

The LSF jobarray generates many results files, one per chunk, which need to be brought back together into a single matrix.

To merge job array results into single logs and matrices per dataset:

```
STAGE=('unscaled' 'scaled' 'scaled_pos_zero')
for stage in "${STAGE[@]}"
do
    bsub < ${REPO_PATH}/SCRIPTS/dual_guide/merge_pred_vs_obs.jobscript.${stage}.sh
done    
```

### Bassik analysis

This script will run the Bassik analysis across all cell lines for all three stages ('unscaled', 'scaled' and 'scaled_pos_zero') by dynamically defining the inputs for `SCRIPTS/dual_guide/pipeline/run_bassik_analysis.R`.

```
bash ${REPO_PATH}/SCRIPTS/dual_guide/run_bassik_analysis.sh
```

Log files can be found in `LOGS/dual_guide/[stage]`. 

TSV files can be found in `DATA/dual_guide/[stage]` and include:

* `DATA/dual_guide/[stage]/pred_vs_obs_with_ngis.tsv` - predicted vs observed values with residuals and normalised GI per guide for all cell lines
* `DATA/dual_guide/[stage]/gene_sample_t_scores.tsv` - gene-level t-statistics for all cell lines
* `DATA/dual_guide/[stage]/sgrna_sample_t_scores.tsv` - guide-level t-statistics for all cell lines

## Combining datasets / results

There are several post-processing steps which are covered by `SCRIPTS/postprocessing/run_postprocessing_rscripts.sh`. Each step utilises a script found in `SCRIPTS/postprocessing/subscripts` whose outputs and general descriptions are:

| Step | Description | Script | Output |
| --- | --- | --- | --- |
| 0 | Combine results of single guide essentiality analyses | `0_combine_single_guide_results.R` | Full results: `DATA/postprocessing/intermediate_tables/combined_singles_results.binary.[stage].tsv`<br/> Binary results: `DATA/postprocessing/intermediate_tables/combined_singles_results.[stage].tsv` |
| 1 | Download and filter DepMap relative copy number | `1_download_depmap_relative_copy_number.R` | `DATA/postprocessing/intermediate_tables/depmap_relative_copy_number.tsv` |
| 2 | Download and filter cell model passport total copy number | `2_download_cmp_total_copy_number.R` | `DATA/postprocessing/intermediate_tables/cmp_total_copy_number.tsv` |
| 3 | Download and filter DepMap expression (log2(TPM + 1)) | `3_download_depmap_expression.R` | `DATA/postprocessing/intermediate_tables/depmap_expression.tsv` |
| 4 | Download and filter DepMap (pre- and post-chronos) and Hart (BAGEL) common essentials and non-essentials | `4_download_common_essentials.R` | `DATA/postprocessing/intermediate_tables/depmap_and_bagel_common_genes.tsv` |
| 5 | Download and filter DepMap depletion probabilities | `5_download_depmap_depletion_probabilities.R` | `DATA/postprocessing/intermediate_tables/depmap_depletion_probability.tsv` | 
| 6 | Download and filter Cell Model Passport cancer driver mutations per cell line | `6_download_ccle_mutations.R` | `DATA/postprocessing/intermediate_tables/cmp_cancer_driver_mutations.tsv` | 
| 7 | Combine intermediate datasets into narrow matrix (1 row = 1 gene pair per cell line) | `7_combine_intermediate_datasets.R` | `DATA/postprocessing/combined_gene_level_results.[stage].tsv` | 
| 8 | Build binary (pass/fail) matrix with data summarised by cancer type (1 row = gene pair) | `8_build_binary_results_matrix.R` | `DATA/postprocessing/combined_gene_level_results.binary.[stage].tsv` | 

Steps 7 and 8 are run for all three stages:

* `unscaled`
* `scaled`
* `scaled_pos_zero`

The combined results matrices (``) are the comprehensive results where each row represents a gene pair in a cell line. 

The binary combined results matrices (``) have the following hard filters applied and are where each row represents a single gene pair and columns are summaries.

Bassik (dual guide):
* `is_bassik_hit` - `mean_norm_gi` < -0.5 and `fdr` < 0.01
* `n_sig_guide_pairs` - number of guide pairs for a gene_pair which have a Bassik guide fdr < 0.05 
* `pct_sig_guide_pairs` - proportion of all guide pairs for a gene pair which are significant (fdr < 0.05) 

Singles analyses (MAGeCK and BAGEL):
* `target[AB]__is_depleted_mageck` - a gene is considered depleted in the singles analysis by MAGeCK when `neg.fdr` < 0.05
* `target[AB]__is_depleted_bagel` - a gene is considered depleted in the singles analysis by BAGEL when BF <= per cell line threshold (i.e. singles analysis binary `is_depleted` = 1)
* `targetA__is_single_depleted` - a gene is considered to be depleted in the singles analyses where it was found depleted in both MAGeCK (`target[AB]__is_depleted_mageck` = 1) and BAGEL (`target[AB]__is_depleted_bagel` = 1)

DepMap expression:
* `target[AB]__is_expressed` - a gene is considered to be expressed when TPM in DepMap log2(TPM + 1) > 1 (i.e. log2(1 + 1) > 1)
* `target[AB]__is_expressed` - a gene is considered to be expressed when TPM in DepMap log2(TPM + 1) > 1 (i.e. log2(1 + 1) > 1)

DepMap dependent:
* `target[AB]__is_depmap_dependent` - a gene is considered DepMap dependent when `depmap_depletion_probability` > 0.5

