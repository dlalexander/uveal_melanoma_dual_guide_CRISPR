#!/usr/bin/env bash

module load R/4.1.0

LOG_PATH="${REPO_PATH}/LOGS/preprocessing"
LABELS_PATH="${REPO_PATH}/METADATA"

# Check repository path exists or exit
if [[ ! -e ${REPO_PATH} ]]
then
	echo "Repository path does not exist: ${REPO_PATH}"
	exit 1;
fi

# Remove log files
if [[ -e ${LOG_PATH} ]]
then
	echo "Removing logs from ${LOG_PATH}"
	rm ${LOG_PATH}/*.o
	rm ${LOG_PATH}/*.e
else
	echo "Log path does not exist: ${LOG_PATH}"
	exit 1;
fi

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

# Creating cell line annotations file
echo "Creating cell line annotations file..."
Rscript ${REPO_PATH}/SCRIPTS/preprocessing/subscripts/sample_labelling.R \
	-d ${REPO_PATH} \
	-m ${REPO_PATH}/METADATA/sample_info.csv \
	-s ${REPO_PATH}/METADATA/iRODS/5429_samples.required.tsv \
	2> ${LOG_PATH}/sample_labelling.e \
	1> ${LOG_PATH}/sample_labelling.o
exit_on_error ${LOG_PATH}/sample_labelling.e

echo "DONE."
