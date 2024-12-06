#!/usr/bin/env bash

module load R/4.1.0

LOG_PATH="${REPO_PATH}/LOGS/postprocessing"
TSV_PATH="${REPO_PATH}/DATA/postprocessing"
RDS_PATH="${REPO_PATH}/DATA/RDS"

STAGES=('unscaled' 'scaled' 'scaled_pos_zero')

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

echo "Removing intermediate files from single guide analysis..."
rm -R ${REPO_PATH}/DATA/single_guide/*/*/*/*/work

# Combine single guide analysis gene level results from MAGeCK and BAGEL
echo "Combining single guide gene-level results..."
Rscript ${REPO_PATH}/SCRIPTS/postprocessing/subscripts/0_combine_single_guide_results.R \
	-d ${REPO_PATH} \
	-m ${REPO_PATH}/METADATA/5429_sample_annotations.tsv \
	-a ${REPO_PATH}/METADATA/libraries/expanded_library.tsv \
	2> ${LOG_PATH}/0_combine_single_guide_gene_results.e \
	1> ${LOG_PATH}/0_combine_single_guide_gene_results.o
exit_on_error ${LOG_PATH}/0_combine_single_guide_gene_results.e

# Download relative copy number per gene per cell line from DepMap
# Not currently used
#echo "Downloading DepMap relative copy number..."
#Rscript ${REPO_PATH}/SCRIPTS/postprocessing/subscripts/1_download_depmap_relative_copy_number.R \
#	-d ${REPO_PATH} \
#	-m ${REPO_PATH}/METADATA/5429_sample_annotations.tsv \
#	-a ${REPO_PATH}/METADATA/libraries/expanded_library.tsv \
#	2> ${LOG_PATH}/1_download_depmap_relative_copy_number.e \
#	1> ${LOG_PATH}/1_download_depmap_relative_copy_number.o
#exit_on_error ${LOG_PATH}/1_download_depmap_relative_copy_number.e

# Download total copy number per gene per cell line from Cell Model Passport
#echo "Downloading Cell Model Passport total copy number..."
#Rscript ${REPO_PATH}/SCRIPTS/postprocessing/subscripts/2_download_cmp_total_copy_number.R \
#	-d ${REPO_PATH} \
#	-m ${REPO_PATH}/METADATA/5429_sample_annotations.tsv \
#	-a ${REPO_PATH}/METADATA/libraries/expanded_library.tsv \
#	2> ${LOG_PATH}/2_download_cmp_total_copy_number.e \
#	1> ${LOG_PATH}/2_download_cmp_total_copy_number.o
#exit_on_error ${LOG_PATH}/2_download_cmp_total_copy_number.e

# Download DepMap expression per gene per cell line
#echo "Downloading DepMap expression (TPM)..."
#Rscript ${REPO_PATH}/SCRIPTS/postprocessing/subscripts/3_download_depmap_expression.R \
#	-d ${REPO_PATH} \
#	-m ${REPO_PATH}/METADATA/5429_sample_annotations.tsv \
#	-a ${REPO_PATH}/METADATA/libraries/expanded_library.tsv \
#	2> ${LOG_PATH}/3_download_depmap_expression.e \
#	1> ${LOG_PATH}/3_download_depmap_expression.o
#exit_on_error ${LOG_PATH}/3_download_depmap_expression.e

# Download common essentials and non-essentials
#echo "Downloading common essentials and non-essentials..."
#Rscript ${REPO_PATH}/SCRIPTS/postprocessing/subscripts/4_download_common_essentials.R \
#	-d ${REPO_PATH} \
#	-a ${REPO_PATH}/METADATA/libraries/expanded_library.tsv \
#	2> ${LOG_PATH}/4_download_common_essentials.e \
#	1> ${LOG_PATH}/4_download_common_essentials.o
#exit_on_error ${LOG_PATH}/4_download_common_essentials.e

# Download DepMap depletion probabilities per gene per cell line
#echo "Downloading DepMap depletion probabilities..."
#Rscript ${REPO_PATH}/SCRIPTS/postprocessing/subscripts/5_download_depmap_depletion_probabilities.R \
#	-d ${REPO_PATH} \
#	-m ${REPO_PATH}/METADATA/5429_sample_annotations.tsv \
#	-a ${REPO_PATH}/METADATA/libraries/expanded_library.tsv \
#	2> ${LOG_PATH}/5_download_depmap_depletion_probabilities.e \
#	1> ${LOG_PATH}/5_download_depmap_depletion_probabilities.o
#exit_on_error ${LOG_PATH}/5_download_depmap_depletion_probabilities.e

# Download Cell Model Passport cancer driver mutations per cell line
#echo "Downloading Cell Model passport cancer driver mutations..."
#Rscript ${REPO_PATH}/SCRIPTS/postprocessing/subscripts/6_download_ccle_mutations.R \
#	-d ${REPO_PATH} \
#	-m ${REPO_PATH}/METADATA/5429_sample_annotations.tsv \
#	-a ${REPO_PATH}/METADATA/libraries/expanded_library.tsv \
#	2> ${LOG_PATH}/6_download_ccle_mutations.e \
#	1> ${LOG_PATH}/6_download_ccle_mutations.o
#exit_on_error ${LOG_PATH}/6_download_ccle_mutations.e

# Combine intermediate datasets to give gene-level results
for stage in "${STAGES[@]}"
do
	echo "Building gene level results for ${stage}..."
	Rscript ${REPO_PATH}/SCRIPTS/postprocessing/subscripts/7_combine_intermediate_datasets.R \
		-d ${REPO_PATH} \
		-m ${REPO_PATH}/METADATA/5429_sample_annotations.tsv \
		-a ${REPO_PATH}/METADATA/libraries/expanded_library.tsv \
		--dataset "${stage}" \
		2> ${LOG_PATH}/7_combine_intermediate_datasets.${stage}.e \
		1> ${LOG_PATH}/7_combine_intermediate_datasets.${stage}.o
	exit_on_error ${LOG_PATH}/7_combine_intermediate_datasets.e
done

# Build binary gene-level results
for stage in "${STAGES[@]}"
do
	echo "Building binary results matrix for ${stage}..."
	Rscript ${REPO_PATH}/SCRIPTS/postprocessing/subscripts/8_build_binary_results_matrix.R \
		-d ${REPO_PATH} \
		--dataset "${stage}" \
		2> ${LOG_PATH}/8_build_binary_results_matrix.${stage}.e \
		1> ${LOG_PATH}/8_build_binary_results_matrix.${stage}.o
	exit_on_error ${LOG_PATH}/8_build_binary_results_matrix.e
done

echo "DONE."