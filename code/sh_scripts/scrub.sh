########################
# scrublet function
########################

scrub ()
{

    module rm python-2.7.12
    module load python-3.11.4


    local input_dir="${1}" # string: path to input directory
    local output_dir="${2}" # string: path to output directory
    local gene_dir="${3}" # string: path to gene directory
    local sample_path="${4}" # string: path to sample.csv file


    # make genes folder
    mkdir "$gene_dir"

    # export genes from filtered fb mats
    for d in "$input_dir"/*/; 
    do 
        name=$(basename "$d")
        gunzip -c "$input_dir"/"$name"/outs/filtered_feature_bc_matrix/features.tsv.gz  > "$gene_dir"/"$name"_features.tsv
    done

    # run scrublet
    python3 /igm/projects/HD-BAT_Wedemeyer/HD-BAT_Functions/multiome_analysis_pipeline/scrub.py "$sample_path" \
                    "$input_dir" \
                    "$output_dir"\
                    "$gene_dir"

}