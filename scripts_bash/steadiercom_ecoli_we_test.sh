source ~/opt/anaconda3/etc/profile.d/conda.sh
conda activate steadiercom

run_steadiercom() {
    community_dir="$1"
    samples="$2"
    output_dir="$3"
    growth="$4"
    we="$5"
    
    
    cmd="steadiercom ../models/Guillen2021/top_7_strains/* -c \"$community_dir\" -m M9 --mediadb ../data/Guillen2021/generated/media_db.tsv" 
    
    if [ -n "$growth" ]; then
        cmd+=" --growth $growth"
    fi

    if [ -n "$we" ]; then
        cmd+=" --we $we"
    fi
    
    cmd+=" -o \"$output_dir\" --sample $samples --exchange-detection '^R_EX_(?!nad_c$|thf_c$).*' --unlimited ../data/Guillen2021/generated/unlimited_glc.txt" 
    
    time eval "$cmd"
}

# Run SteadierCom with abundance and growth


input_dir="../data/Guillen2021/generated/communities_smetana.tsv" # No abundance
output_dir1="../output/Guillen2021/results_steadiercom/results_we_test/0_0025_"
output_dir2="../output/Guillen2021/results_steadiercom/results_we_test/0_0035_"

samples=1000
growth=0.221

we1=0.0025
we2=0.0035

run_steadiercom "$input_dir" "$samples" "$output_dir1" "$growth" "$we1"

run_steadiercom "$input_dir" "$samples" "$output_dir2" "$growth" "$we2"




