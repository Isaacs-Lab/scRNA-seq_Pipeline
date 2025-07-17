##############################################
# MULTIOME ANALYSIS PIPELINE
# Author: Gary Schweickart - Wedemeyer Lab
# Created: 8/27/24      
# Last Updated: 3/21/25
###############################################

setup_to_qc() 
{
    module load R-4.3.3

    local input_dir="${1}" # string: path to input directory
    local output_dir="${2}" # string: path to output directory
    local sample_path="${3}" # string: path to sample.csv file
    #local macs2="${4}" # boolean: run macs2 (in R)
    local scrub="${4}" # boolean: run scrub (in python)
    local project="${5}"
    local gene_dir="${output_dir}"/genes/ # string: path to gene directory


    echo Creating output directory: "$output_dir"
    mkdir "${output_dir}"
    echo Creating scrublet output directory: "${output_dir}"/scrub
    mkdir "${output_dir}"/scrub

    echo Running scrublet...
    source /igm/projects/HD-BAT_Wedemeyer/HD-BAT_Functions/multiome_analysis_pipeline/multiome_analysis_pipeline/code/sh_scripts/scrub.sh; scrub "$input_dir" \
                             "$output_dir"/scrub \
                             "$gene_dir" \
                             "$sample_path"
    echo Scrublet done.
    
    module rm python-2.7.12
    module rm python-3.11.4
    
    echo Running multiome pipeline to QC
    Rscript --vanilla /igm/projects/HD-BAT_Wedemeyer/HD-BAT_Functions/multiome_analysis_pipeline/multiome_analysis_pipeline/code/Rscripts/multome_to_QC.R "${input_dir}" "${output_dir}" "${sample_path}" "${macs2}" "${scrub}" "${project}" 
    echo Done.
}

qc_merge()
{
    module load R-4.3.3

    local input_dir="${1}" # string: path to input directory
    local output_dir="${2}" # string: path to output directory
    local qc_vals="${3}" # string: path to QC_Values.csv file
    local integrate="${4}" # boolean: whether to integrate or not
    
    echo Running multiome pipeline to QC
    Rscript --vanilla /igm/projects/HD-BAT_Wedemeyer/HD-BAT_Functions/multiome_analysis_pipeline/multiome_analysis_pipeline/code/Rscripts/QC_merge.R "${input_dir}" "${output_dir}" "${qc_vals}" "${integrate}"
    echo Done.
}
