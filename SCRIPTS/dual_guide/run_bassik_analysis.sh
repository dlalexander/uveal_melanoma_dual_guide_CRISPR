#!/usr/bin/env bash

module load R/4.1.0

# Check repository path exists or exit
if [[ ! -e ${REPO_PATH} ]]
then
	echo "Repository path does not exist: ${REPO_PATH}"
	exit 1;
fi

echo "Repository path: ${REPO_PATH}"

STAGES=('unscaled' 'scaled' 'scaled_pos_zero')
TSV_PATH="${REPO_PATH}/DATA/dual_guide"
RDS_PATH="${REPO_PATH}/DATA/RDS/dual_guide"

echo "TSV output directory: ${TSV_PATH}"
echo "RDS output directory: ${RDS_PATH}"

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

# Loop over each stage
for stage in "${STAGES[@]}"
do
	echo "Running Bassik script for ${stage}..."

	LOG_PATH="${REPO_PATH}/LOGS/dual_guide/${stage}"
	echo "Log path: ${LOG_PATH}"

	fc="${REPO_PATH}/DATA/preprocessing/lfc_matrix.${stage}.samples_removed.tsv"
	y12="${REPO_PATH}/DATA/dual_guide/${stage}/pred_vs_obs_y12.tsv"
	missing="${REPO_PATH}/DATA/dual_guide/${stage}/missing_data.tsv"
	samples="${REPO_PATH}/METADATA/5429_sample_annotations.tsv"
	control_guides="F1,F2,F3,F4,F5,F6,F7,F8,F9,F10"
	tsv="${TSV_PATH}/${stage}"
	rds="${RDS_PATH}/${stage}"

	Rscript ${REPO_PATH}/SCRIPTS/dual_guide/pipeline/run_bassik_analysis.R \
		-f ${fc} \
		-y ${y12} \
		-m ${missing} \
		-s ${samples} \
		-c ${control_guides} \
		-t ${tsv} \
		-r ${rds} \
		2>${LOG_PATH}/run_bassik_analysis.e \
		1>${LOG_PATH}/run_bassik_analysis.o
	exit_on_error ${LOG_PATH}/run_bassik_analysis.e
done

echo "DONE."