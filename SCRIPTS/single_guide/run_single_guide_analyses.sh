#!/usr/bin/env bash

module load R/4.1.0
module load C-SAR/1.3.6

# Script for running single guide analysis, with BAGEL core essentials and non-essentials
LOG_PATH="LOGS/single_guide"
DATA_PATH="DATA/single_guide"
DATASETS=('A' 'B' 'combined')
STAGES=('unscaled' 'scaled' 'scaled_pos_zero')
SAMPLES=('all' 'samples_removed')
BATCHES=('batchA' 'batchB')

echo "Repository path: ${REPO_PATH}"
echo "Log path: ${REPO_PATH}/${LOG_PATH}"

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

# Check repository path exists or exit
if [[ ! -e ${REPO_PATH} ]]
then
	echo "Repository path does not exist: ${REPO_PATH}"
	exit 1;
fi

# Create output directories
# Loop over dataset and stages and create directories
for stage in "${STAGES[@]}"
do
	for dataset in "${DATASETS[@]}"
	do	
		for batch in "${BATCHES[@]}"
		do
			if [[ ! -d ${REPO_PATH}/${DATA_PATH}/${dataset}/${stage}/${batch} ]]
			then	
				mkdir -p ${REPO_PATH}/${DATA_PATH}/${dataset}/${stage}/${batch}
			fi
			if [[ ! -d "${REPO_PATH}/${LOG_PATH}/C-SAR/${dataset}/${stage}/${batch}" ]]
			then	
				mkdir -p "${REPO_PATH}/${LOG_PATH}/C-SAR/${dataset}/${stage}/${batch}"
			fi
		done
	done
done

# Prepare singles libraries with user-defined guides removed
echo "Preparing initial singles libraries..."
Rscript ${REPO_PATH}/SCRIPTS/single_guide/0_prepare_single_guide_libraries.R \
	-d ${REPO_PATH} \
	-g ${REPO_PATH}/METADATA/guides_ids_to_remove.txt \
	-a ${REPO_PATH}/METADATA/libraries/expanded_library.tsv \
	2> ${REPO_PATH}/${LOG_PATH}/0_prepare_libraries.e \
	1> ${REPO_PATH}/${LOG_PATH}/0_prepare_libraries.o
exit_on_error ${REPO_PATH}/${LOG_PATH}/0_prepare_libraries.e

# NB: This analysis only uses the 'samples_removed' dataset
# In this case we have not removed any of the samples so this is the same as 'all' but note this for if exclude some samples in future
# Prepare inputs for each stage
echo "Preparing single pipeline inputs for analysis..."

for stage in "${STAGES[@]}"
do
	for dataset in "${DATASETS[@]}"
	do	
		for batch in "${BATCHES[@]}"
		do
			Rscript ${REPO_PATH}/SCRIPTS/single_guide/1_prepare_files_and_directories_for_analysis.R \
				-d ${REPO_PATH} \
				-a ${DATA_PATH}/${dataset}/singles_library.${dataset}.tsv \
				-m ${REPO_PATH}/METADATA/5429_sample_annotations_${batch}.tsv \
				-c ${REPO_PATH}/DATA/preprocessing/${batch}/count_matrix.${stage}.samples_removed.tsv \
				-l ${LOG_PATH}/C-SAR/${dataset}/${stage}/${batch} \
				-b ${batch} \
				-s "${dataset}.${stage}.${batch}" \
				-e "${cell_lines_to_exclude}" \
				-o "${DATA_PATH}/${dataset}/${stage}/${batch}" \
				2> ${REPO_PATH}/${LOG_PATH}/1_prepare_files_and_directories_for_analysis.${stage}.${dataset}.${batch}.e \
				1> ${REPO_PATH}/${LOG_PATH}/1_prepare_files_and_directories_for_analysis.${stage}.${dataset}.${batch}.o
			exit_on_error ${REPO_PATH}/${LOG_PATH}/1_prepare_files_and_directories_for_analysis.${stage}.${dataset}.${batch}.e
		done
	done
done

# Run jobscripts
echo "Running single guide jobscripts..."
for stage in "${STAGES[@]}"
do
	for dataset in "${DATASETS[@]}"
	do	
		for batch in "${BATCHES[@]}"
		do
			bash ${REPO_PATH}/SCRIPTS/single_guide/bsub_commands.${dataset}.${stage}.${batch}.sh
		done
	done
done
cd ${REPO_PATH}

echo "DONE."
