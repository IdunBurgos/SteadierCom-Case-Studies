source ~/opt/anaconda3/etc/profile.d/conda.sh
conda activate steadiercom

run_steadiercom() {
    community_dir="$1"
    samples="$2"
    output_dir="$3"
    growth="$4"
    
    community_id=$(echo "$community_dir" | rev | cut -d "/" -f1 | rev)
    
    echo "${community_id}"
    
    cmd="steadiercom ../models/Venturelli2018/species_used_noO2/* -c \"$community_dir\" -m LB_extend --mediadb ../data/Venturelli2018/generated/media_db.tsv"
    
    if [ -n "$growth" ]; then
        cmd+=" --growth $growth"
    fi
    
    cmd+=" -o \"$output_dir\" --sample $samples --exchange-detection 'R_EX_(?!biomass__40__e__41__)' --unlimited ../data/Venturelli2018/generated/unlimited.txt"
    
    time eval "$cmd"
}

input_dir1="../data/Venturelli2018/generated/communities_pairwise.tsv"
input_dir2="../data/Venturelli2018/generated/communities_pairwise_no_abundance.tsv" # No abundance
output_dir1="../output/Venturelli2018/results_steadiercom/results_pairwise_growth_abundance/"
output_dir2="../output/Venturelli2018/results_steadiercom/results_pairwise_growth/"

samples1=1000
samples2=2000
growth=0.125

run_steadiercom "$input_dir1" "$samples1" "$output_dir1" "$growth"
run_steadiercom "$input_dir2" "$samples2" "$output_dir2" "$growth"