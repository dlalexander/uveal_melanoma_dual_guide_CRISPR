#BSUB -q normal
#BSUB -J merge_pred_vs_obs
#BSUB -oo "${REPO_PATH}/LOGS/dual_guide/unscaled/merge_pred_vs_obs.o"
#BSUB -eo "${REPO_PATH}/LOGS/dual_guide/unscaled/merge_pred_vs_obs.e"
#BSUB -R "select[mem>5000] rusage[mem=5000] span[hosts=1]"
#BSUB -M 5000

# Set R library and load R
#export R_LIBS_USER=~vo1/r_libs_dual_crispr
module load R/4.1.0

# Set repository path
repo_path="${REPO_PATH}"

# Set inputs and output paths
stage='unscaled'
output_directory="${repo_path}/DATA/dual_guide/${stage}"
rds_directory="${repo_path}/DATA/RDS/dual_guide/${stage}"
log_directory="${repo_path}/LOGS/dual_guide/${stage}"
script_directory="${repo_path}/SCRIPTS/dual_guide/pipeline"

# Clean out previous results / logs
if [[ "$(ls -A ${log_directory})" ]]
then
	rm ${log_directory}/list_of_pred_vs_obs_y12_files.txt
	rm ${log_directory}/list_of_missing_data_files.txt
	rm ${log_directory}/get_pred_vs_obs.log
	rm ${log_directory}/get_pred_vs_obs.o
	rm ${log_directory}/get_pred_vs_obs.e
fi

if [[ "$(ls -A ${output_directory})" ]]
then
	rm ${output_directory}/pred_vs_obs_y12.tsv
	rm ${output_directory}/missing_data.tsv
fi

# Merge pred_vs_obs_y12 files
find "${output_directory}" -name "pred_vs_obs_y12.[0-9]*" > "${log_directory}/list_of_pred_vs_obs_y12_files.txt"
merge_pred_vs_obs_y12_cmd="Rscript ${script_directory}/merge_chunks.R -f ${log_directory}/list_of_pred_vs_obs_y12_files.txt -p pred_vs_obs_y12 -o ${output_directory} -r ${rds_directory}"
echo "Merge pred_vs_obs_y12 chunks: "$(date +"%T")
echo "${merge_pred_vs_obs_y12_cmd}"
eval "$merge_pred_vs_obs_y12_cmd"

# Merge missing_data files
find "${output_directory}" -name "missing_data.[0-9]*" > "${log_directory}/list_of_missing_data_files.txt"
merge_missing_data_cmd="Rscript ${script_directory}/merge_chunks.R -f ${log_directory}/list_of_missing_data_files.txt -p missing_data -o ${output_directory} -r ${rds_directory}"
echo "Merge missing_data chunks: "$(date +"%T")
echo "${merge_missing_data_cmd}"
eval "$merge_missing_data_cmd"

# Clean up outputs
find ${output_directory} -name 'missing_data.[0-9]*' -delete
find ${output_directory} -name 'pred_vs_obs_y12.[0-9]*' -delete

# Clean up logs
cat ${log_directory}/get_pred_vs_obs_y12.[0-9]*.log > ${log_directory}/get_pred_vs_obs.log
cat ${log_directory}/get_pred_vs_obs.[0-9]*.o > ${log_directory}/get_pred_vs_obs.o
cat ${log_directory}/get_pred_vs_obs.[0-9]*.e > ${log_directory}/get_pred_vs_obs.e
find ${log_directory} -name 'get_pred_vs_obs_y12.[0-9]*.log' -delete
find ${log_directory} -name 'get_pred_vs_obs.[0-9]*.o' -delete
find ${log_directory} -name 'get_pred_vs_obs.[0-9]*.e' -delete