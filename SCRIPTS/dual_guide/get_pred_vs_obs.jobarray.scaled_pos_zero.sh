#BSUB -q normal
#BSUB -J pred_vs_obs[1-1600]
#BSUB -oo "${REPO_PATH}/LOGS/dual_guide/scaled_pos_zero/get_pred_vs_obs.%I.o"
#BSUB -eo "${REPO_PATH}/LOGS/dual_guide/scaled_pos_zero/get_pred_vs_obs.%I.e"
#BSUB -R "select[mem>5000] rusage[mem=5000] span[hosts=1]"
#BSUB -M 5000

# Set R library and load R
#export R_LIBS_USER=~vo1/r_libs_dual_crispr
module load R/4.1.0

# Set repository path
repo_path="${REPO_PATH}"

# Set number of jobs (should mirror array range in LSF header)
njobs=1600

# Set inputs and output paths
stage='scaled_pos_zero'
lfc_matrix="${repo_path}/DATA/preprocessing/lfc_matrix.${stage}.samples_removed.tsv"
output_directory="${repo_path}/DATA/dual_guide/${stage}"
log_directory="${repo_path}/LOGS/dual_guide/${stage}"
script_directory="${repo_path}/SCRIPTS/dual_guide/pipeline"

# Load dual guide matrix (dual guides mapped to their corresponding singles )
doubles_gm_file="${repo_path}/DATA/dual_guide/dual_guide_matrix.tsv"
doubles_gm_rownum=$(wc -l "${doubles_gm_file}" | awk '{print $1-1}')

# Calculate chunk size from length of dual guide matrix and number of jobs
chunk_size=$(Rscript -e "ceiling( ${doubles_gm_rownum} / ${njobs} ) " | awk '{print $2}' ) 

# Prepare command to get predicted and observed LFCs for chunk
obs_pred_cmd="Rscript ${script_directory}/get_obs_and_pred_values_from_normalised_foldchanges_by_chunk.R -f ${lfc_matrix} -d ${doubles_gm_file} -o ${output_directory} -n ${chunk_size} -i ${LSB_JOBINDEX}"

# Clean out previous results / logs
#if [[ "$(ls -A ${log_directory})" ]]
#then
#	rm ${log_directory}/get_pred_vs_obs.[0-9]*.log
#	rm ${log_directory}/get_pred_vs_obs_y12.[0-9]*.o
#	rm ${log_directory}/get_pred_vs_obs_y12.[0-9]*.e
#	rm ${log_directory}/get_pred_vs_obs_y12.log
#	rm ${log_directory}/get_pred_vs_obs_y12.o
#	rm ${log_directory}/get_pred_vs_obs_y12.e
#	rm ${log_directory}/list_of_pred_vs_obs_y12_files.txt
#	rm ${log_directory}/list_of_missing_data_files.txt
#fi

#if [[ "$(ls -A ${output_directory})" ]]
#then
#	rm ${output_directory}/missing_data.[0-9]*.tsv
#	rm ${output_directory}/pred_vs_obs_y12.[0-9]*.tsv
#	rm ${output_directory}/missing_data.tsv
#	rm ${output_directory}/pred_vs_obs_y12.tsv
#fi 

# Run script
echo "Getting observed and predicted values (${LSB_JOBINDEX} of ${njobs}): "$(date +"%T")
echo "${obs_pred_cmd}"
eval "$obs_pred_cmd > ${log_directory}/get_pred_vs_obs_y12.${LSB_JOBINDEX}.log 2>&1 &"