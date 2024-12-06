#!/usr/bin/env bash

module load R/4.1.0
module load c-sar/1.2.0

LOG_PATH="${REPO_PATH}/LOGS/QC"
PLOT_PATH="${REPO_PATH}/DATA/QC"
STAGES=('unscaled' 'scaled' 'scaled_pos_zero')
DATASETS=('all' 'samples_removed')

echo "Repository path: ${REPO_PATH}"
echo "Log path: ${LOG_PATH}"
echo "Plot output directory: ${PLOT_PATH}"

# Check repository path exists or exit
if [[ ! -e ${REPO_PATH} ]]
then
	echo "Repository path does not exist: ${REPO_PATH}"
	exit 1;
fi

# Remove log files
#if [[ -e ${LOG_PATH} ]]
#then
#	echo "Removing logs from ${LOG_PATH}"
#	rm ${LOG_PATH}/*.o
#	rm ${LOG_PATH}/*.e
#else
#	echo "Log path does not exist: ${LOG_PATH}"
#	exit 1;
#fi

# Function for exiting on error
exit_on_error () {
	err_path=$1
	if [[ -e $err_path ]]
	then
		if [[ ! -z $(grep '[^[:space:]]' $err_path) ]]
		then
			echo "Error found in: ${err_path}"
			exit 1;
		fi
	fi
}


# Analysing correlation between controls and between test samples: normalised counts
echo "Analysing correlation between normalised counts of samples..."
Rscript ${REPO_PATH}/SCRIPTS/QC/subscripts/1_correlation.R \
	-d ${REPO_PATH} \
	-m ${REPO_PATH}/METADATA/5429_sample_annotations.tsv \
	-d normalised_counts \
	2> ${LOG_PATH}/1_correlation_normalised_counts.e \
	1> ${LOG_PATH}/1_correlation_normalised_counts.o
exit_on_error ${LOG_PATH}/1_correlation_normalised_counts.e

# Analysing correlation between controls and between test samples: LFCs
echo "Analysing correlation between lfcs of samples..."
Rscript ${REPO_PATH}/SCRIPTS/QC/subscripts/1_correlation.R \
	-d ${REPO_PATH} \
	-m ${REPO_PATH}/METADATA/5429_sample_annotations.tsv \
	-d lfcs \
	2> ${LOG_PATH}/1_correlation_lfcs.e \
	1> ${LOG_PATH}/1_correlation_lfcs.o
exit_on_error ${LOG_PATH}/1_correlation_lfcs.e

echo "Analysing lfc densities for different gene targets ..."
Rscript ${REPO_PATH}/SCRIPTS/QC/subscripts/2_guide_type_lfc_density_singles.R \
	-d ${REPO_PATH} \
	-m ${REPO_PATH}/METADATA/5429_sample_annotations.tsv \
	-d lfcs \
	2> ${LOG_PATH}/2_guide_type_lfc_density_singles.e \
	1> ${LOG_PATH}/2_guide_type_lfc_density_singles.o
exit_on_error ${LOG_PATH}/2_guide_type_lfc_density_singles.e

echo "Analysing promoter strengths..."
Rscript ${REPO_PATH}/SCRIPTS/QC/subscripts/3_promoter_strengths.R \
	-d ${REPO_PATH} \
	2> ${LOG_PATH}/3_promoter_strengths.e \
	1> ${LOG_PATH}/3_promoter_strengths.o
exit_on_error ${LOG_PATH}/3_promoter_strengths.e

echo "Analysing essential and nonessential separation NNMDs and SSMDs ..."
Rscript ${REPO_PATH}/SCRIPTS/QC/subscripts/4_NNMD_SSMD.R \
	-d ${REPO_PATH} \
	-m ${REPO_PATH}/METADATA/5429_sample_annotations.tsv \
	2> ${LOG_PATH}/4_NNMD_SSMD.e \
	1> ${LOG_PATH}/4_NNMD_SSMD.o
exit_on_error ${LOG_PATH}/4_NNMD_SSMD.e
