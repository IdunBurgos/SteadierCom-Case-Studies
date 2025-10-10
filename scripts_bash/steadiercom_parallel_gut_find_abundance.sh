#!/usr/bin/env bash

source ~/opt/anaconda3/etc/profile.d/conda.sh
conda activate steadiercom

run_steadiercom() {
    community_dir="$1"
    samples="$2"
    output_dir="$3"
    growth="$4"
    
    
    community_id=$(echo "$community_dir" | rev | cut -d "/" -f1 | rev)
    output_path="${output_dir}/${community_id%_abundance.tsv}"
    
    echo "${community_id}"
    
    cmd="steadiercom ../models/Venturelli2018/species_used_noO2/* -c \"$community_dir\" -m LB_extend --mediadb ../data/Venturelli2018/generated/media_db.tsv"
    
    if [ -n "$growth" ]; then
        cmd+=" --growth $growth"
    fi
    
    cmd+=" -o \"$output_path\" --sample $samples --exchange-detection 'R_EX_(?!biomass__40__e__41__)' --unlimited ../data/Venturelli2018/generated/unlimited.txt"
    
    time eval "$cmd"
}

samples=2000
growth=0.125

# Run SteadierCom with growth
 
output_dir="../output/Venturelli2018/results_steadiercom/results_growth"
input_dir="../data/Venturelli2018/generated/communities_steadiercom_no_abundance"

for dir in "$input_dir"/*; do
  run_steadiercom "$dir" "$samples" "$output_dir" "$growth"
done



