#BSUB -q normal
#BSUB -J pycroquet[1-1008]
#BSUB -oo "${REPO_PATH}/LOGS/pyCROQUET/pycroquet.%I.o"
#BSUB -eo "${REPO_PATH}/LOGS/pyCROQUET/pycroquet.%I.e"
#BSUB -R "select[mem>4000] rusage[mem=4000] span[hosts=1]"
#BSUB -M 4000

#Previously ran with v 1.3.0
#module load pycroquet/1.3.0
module load pycroquet/1.5.1
export SINGULARITYENV_REF_PATH='insert/your/singularity/ref/path'

guides="${REPO_PATH}/METADATA/libraries/paired_lib_original.pyCROQUET.tsv"
cram_directory="${REPO_PATH}/DATA/CRAM"
cram_files=($(ls ${cram_directory}/*.cram))
array_index=$(expr ${LSB_JOBINDEX} - 1)
query_file_name=$(basename ${cram_files[${array_index}]})
query_file_label=${query_file_name::-5}
output_directory="${REPO_PATH}/DATA/pyCROQUET/COUNTS"

pycroquet dual-guide -g ${guides} -q "${cram_directory}/${query_file_name}" -o "${output_directory}/${query_file_label}" --chunks 50000
