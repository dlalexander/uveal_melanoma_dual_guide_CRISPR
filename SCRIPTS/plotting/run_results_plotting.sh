#!/usr/bin/env bash

module load R/4.1.0
module load c-sar/1.2.0

Rscript ${REPO_PATH}/SCRIPTS/plotting/results_plots.R \
	-d ${REPO_PATH} \
	-m ${REPO_PATH}/METADATA/5429_sample_annotations.tsv

