#!/usr/bin/env bash

module load R/4.1.0
module load C-SAR/1.3.7

# Script to preprocess data prior to single KO and dual KO analysis

LOG_PATH="${REPO_PATH}/LOGS/preprocessing"
TSV_PATH="${REPO_PATH}/DATA/preprocessing"
RDS_PATH="${REPO_PATH}/DATA/RDS"

STAGES=('unscaled' 'scaled' 'scaled_pos_zero')
DATASETS=('all' 'samples_removed')
BATCHES=('batchA' 'batchB')

echo "Repository path: ${REPO_PATH}"
echo "Log path: ${LOG_PATH}"
echo "TSV output directory: ${TSV_PATH}"
echo "RDS output directory: ${RDS_PATH}"

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

# Prepare singles and dual libraries with expanded annotations
echo "Preparing libraries for analysis..."
Rscript ${REPO_PATH}/SCRIPTS/preprocessing/subscripts/0_prepare_libraries.R \
	-d ${REPO_PATH} \
	-m ${REPO_PATH}/METADATA/5429_sample_annotations.tsv \
	-a ${REPO_PATH}/METADATA/libraries/paired_lib_original.pyCROQUET.tsv \
	2> ${LOG_PATH}/0_prepare_libraries.e \
	1> ${LOG_PATH}/0_prepare_libraries.o
exit_on_error ${LOG_PATH}/0_prepare_libraries.e

# Gunzip pyCROQUET count files generated by pyCROQUET 1.5.1 (different from 1.3.0)
gunzip ${REPO_PATH}/DATA/pyCROQUET/COUNTS*counts.tsv.gz

# Merge per-lane pyCROQUET counts into per-sample count matrix with annotations
echo "Converting lane counts into sample matrix..."
Rscript ${REPO_PATH}/SCRIPTS/preprocessing/subscripts/1_merge_counts_by_sample.R \
	-d ${REPO_PATH} \
	-m ${REPO_PATH}/METADATA/5429_sample_annotations.tsv \
	-a ${REPO_PATH}/METADATA/libraries/expanded_library.tsv \
	2>${LOG_PATH}/1_merge_counts_by_sample.e \
	1>${LOG_PATH}/1_merge_counts_by_sample.o
exit_on_error ${LOG_PATH}/1_merge_counts_by_sample.e

# Remove user-defined guides
echo "Removing user-defined guides from count matrix..."
Rscript ${REPO_PATH}/SCRIPTS/preprocessing/subscripts/2_remove_user_defined_guides.R \
	-d ${REPO_PATH} \
	-g ${REPO_PATH}/METADATA/guides_ids_to_remove.txt \
	-c ${REPO_PATH}/DATA/preprocessing/count_matrix.tsv \
	2>${LOG_PATH}/2_remove_user_defined_guides.e \
	1>${LOG_PATH}/2_remove_user_defined_guides.o
exit_on_error ${LOG_PATH}/2_remove_user_defined_guides.e

# Normalise counts using BAGEL method
echo "Normalising counts using BAGEL method..."
Rscript ${REPO_PATH}/SCRIPTS/preprocessing/subscripts/3_bagel_normalisation.R \
	-d ${REPO_PATH} \
	-c ${REPO_PATH}/DATA/preprocessing/count_matrix.guides_removed.tsv \
	2>${LOG_PATH}/3_bagel_normalisation.e \
	1>${LOG_PATH}/3_bagel_normalisation.o
exit_on_error ${LOG_PATH}/3_bagel_normalisation.e

# Split the samples into different batches so can calculate LFCs per batch
echo "Split data into batches..."
Rscript ${REPO_PATH}/SCRIPTS/preprocessing/subscripts/3.5_batch_splitting.R \
	-d ${REPO_PATH} \
	-c ${REPO_PATH}/DATA/preprocessing/count_matrix.norm.tsv \
	-m ${REPO_PATH}/METADATA/5429_sample_annotations.tsv \
	2>${LOG_PATH}/3.5_batch_splitting.e \
	1>${LOG_PATH}/3.5_batch_splitting.o
exit_on_error ${LOG_PATH}/3.5_batch_splitting.e

# Filter out guides with low counts in controls
# Done per batch
echo "Filter guides with low control counts for batch A..."
Rscript ${REPO_PATH}/SCRIPTS/preprocessing/subscripts/4_filter_low_count_guides.R \
	-d ${REPO_PATH} \
	-c ${REPO_PATH}/DATA/preprocessing/batchA/count_matrix.norm.batchA.tsv \
	-m ${REPO_PATH}/METADATA/5429_sample_annotations_batchA.tsv \
	-b batchA \
	2>${LOG_PATH}/4_filter_low_count_guides_batchA.e \
	1>${LOG_PATH}/4_filter_low_count_guides_batchA.o
exit_on_error ${LOG_PATH}/4_filter_low_count_guides_batchA.e

#echo "Filter guides with low control counts for batch B..."
Rscript ${REPO_PATH}/SCRIPTS/preprocessing/subscripts/4_filter_low_count_guides.R \
	-d ${REPO_PATH} \
	-c ${REPO_PATH}/DATA/preprocessing/batchB/count_matrix.norm.batchB.tsv \
	-m ${REPO_PATH}/METADATA/5429_sample_annotations_batchB.tsv \
	-b batchB \
	2>${LOG_PATH}/4_filter_low_count_guides_batchB.e \
	1>${LOG_PATH}/4_filter_low_count_guides_batchB.o
exit_on_error ${LOG_PATH}/4_filter_low_count_guides_batchB.e


# Convert count matrix to LFC matrix for each batch
echo "Convert batch A count matrix to LFC matrix..."
Rscript ${REPO_PATH}/SCRIPTS/preprocessing/subscripts/5_convert_counts_to_lfcs.R \
	-d ${REPO_PATH} \
	-c ${REPO_PATH}/DATA/preprocessing/batchA/count_matrix.norm.filt.tsv \
	-m ${REPO_PATH}/METADATA/5429_sample_annotations_batchA.tsv \
	-b batchA \
	1>${LOG_PATH}/5_convert_counts_to_lfcs_batchA.o
	2>${LOG_PATH}/5_convert_counts_to_lfcs_batchA.e \
exit_on_error ${LOG_PATH}/5_convert_counts_to_lfcs_batchA.e

#echo "Convert batch B count matrix to LFC matrix..."
Rscript ${REPO_PATH}/SCRIPTS/preprocessing/subscripts/5_convert_counts_to_lfcs.R \
	-d ${REPO_PATH} \
	-c ${REPO_PATH}/DATA/preprocessing/batchB/count_matrix.norm.filt.tsv \
	-m ${REPO_PATH}/METADATA/5429_sample_annotations_batchB.tsv \
	-b batchB \
	2>${LOG_PATH}/5_convert_counts_to_lfcs_batchB.e \
	1>${LOG_PATH}/5_convert_counts_to_lfcs_batchB.o
exit_on_error ${LOG_PATH}/5_convert_counts_to_lfcs_batchB.e


# Remove samples failing QC
echo "Remove samples failing QC from batch A..."
samples_to_remove_batchA=""
Rscript ${REPO_PATH}/SCRIPTS/preprocessing/subscripts/6_remove_samples.R \
	-d ${REPO_PATH} \
	-s "${samples_to_remove_batchA}" \
	-f ${REPO_PATH}/DATA/preprocessing/batchA/lfc_matrix.unscaled.all.tsv \
	-b batchA \
	2>${LOG_PATH}/6_remove_samples_batchA.e \
	1>${LOG_PATH}/6_remove_samples_batchA.o
exit_on_error ${LOG_PATH}/6_remove_samples_batchA.e

# Remove samples failing QC
echo "Remove samples failing QC from batch B..."
samples_to_remove_batchB=""
Rscript ${REPO_PATH}/SCRIPTS/preprocessing/subscripts/6_remove_samples.R \
	-d ${REPO_PATH} \
	-s "${samples_to_remove_batchB}" \
	-f ${REPO_PATH}/DATA/preprocessing/batchB/lfc_matrix.unscaled.all.tsv \
	-b batchB \
	2>${LOG_PATH}/6_remove_samples_batchB.e \
	1>${LOG_PATH}/6_remove_samples_batchB.o
exit_on_error ${LOG_PATH}/6_remove_samples_batchB.e


# Scale LFCs
# Done per batch
# Scaling using BAGEL Core essentials list
# Updated from Harle dual guide screen which used library classifications of essential and non-essential
for batch in "${BATCHES[@]}"
do
	for dataset in "${DATASETS[@]}"
	do
		echo "Scale LFC matrix (${dataset} ${batch})..."
		Rscript ${REPO_PATH}/SCRIPTS/preprocessing/subscripts/7_scale_lfcs.R \
			-d ${REPO_PATH} \
			-f ${REPO_PATH}/DATA/preprocessing/${batch}/lfc_matrix.unscaled.${dataset}.tsv \
			-m ${REPO_PATH}/METADATA/5429_sample_annotations_${batch}.tsv \
			-e ${REPO_PATH}/METADATA/essentials_nonessentials/CEGv2.txt \
			-s "${dataset}" \
			-b ${batch} \
			2>${LOG_PATH}/7_scale_lfcs.${dataset}_${batch}.e \
			1>${LOG_PATH}/7_scale_lfcs.${dataset}_${batch}.o
		exit_on_error ${LOG_PATH}/7_scale_lfcs.${dataset}_${batch}.e
	done
done



# Convert scaled LFC matrix to scaled count matrix
# Done per batch
for batch in "${BATCHES[@]}"
do
	for stage in "${STAGES[@]}"
	do
		for dataset in "${DATASETS[@]}"
		do
			echo "Convert scaled LFC matrix to scaled count matrix (${dataset} ${stage} ${batch})..."
			Rscript ${REPO_PATH}/SCRIPTS/preprocessing/subscripts/8_convert_scaled_lfcs_to_scaled_counts.R \
			-d ${REPO_PATH} \
			-f ${REPO_PATH}/DATA/preprocessing/${batch}/lfc_matrix.${stage}.${dataset}.tsv \
			-s "${stage}.${dataset}" \
			-c ${REPO_PATH}/DATA/preprocessing/${batch}/count_matrix.norm.filt.mean_control.tsv \
			-b ${batch} \
			2>${LOG_PATH}/8_convert_scaled_lfcs_to_scaled_counts.${stage}.${dataset}.${batch}.e \
			1>${LOG_PATH}/8_convert_scaled_lfcs_to_scaled_counts.${stage}.${dataset}.${batch}.o
			exit_on_error ${LOG_PATH}/8_convert_scaled_lfcs_to_scaled_counts.${stage}.${dataset}.${batch}.e
		done
	done
done

echo "DONE."
