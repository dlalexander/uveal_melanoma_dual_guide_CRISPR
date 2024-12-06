#!/usr/bin/env bash

module load R/4.1.0

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

LOG_PATH="${REPO_PATH}/LOGS/dual_guide"
DATA_PATH="${REPO_PATH}/DATA/dual_guide"
STAGES=('unscaled' 'scaled' 'scaled_pos_zero')
DATASETS=('all' 'samples_removed')

echo "Repository path: ${REPO_PATH}"
echo "Log path: ${REPO_PATH}/${LOG_PATH}"

# Preparing merged LFC matrix
echo "Merging LFC files from batchA and batchB for input to dual pipeline..."
for stage in "${STAGES[@]}"
do
	for dataset in "${DATASETS[@]}"
	do	
		Rscript ${REPO_PATH}/SCRIPTS/dual_guide/pipeline/merge_LFC_batches.R \
			-d ${REPO_PATH} \
			-a ${REPO_PATH}/DATA/preprocessing/batchA/lfc_matrix.${stage}.${dataset}.tsv \
			-b ${REPO_PATH}/DATA/preprocessing/batchB/lfc_matrix.${stage}.${dataset}.tsv \
			-s "${stage}.${dataset}" \
			-o "${REPO_PATH}/DATA/preprocessing" \
			2> ${LOG_PATH}/merge_LFC_batches.${stage}.${dataset}.e \
			1> ${LOG_PATH}/merge_LFC_batches.${stage}.${dataset}.o
		exit_on_error ${LOG_PATH}/merge_LFC_batches.${stage}.${dataset}.e
	done
done

echo "DONE."
