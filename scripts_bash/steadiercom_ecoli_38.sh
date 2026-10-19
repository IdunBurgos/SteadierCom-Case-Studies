source ~/opt/anaconda3/etc/profile.d/conda.sh
conda activate steadiercom

run_steadiercom() {
    community_dir="$1"
    samples="$2"
    output_dir="$3"
    growth="$4"
    
    
    
    cmd="steadiercom ../models/Guillen2021/strains/* -c \"$community_dir\" -m M9 --mediadb ../data/Guillen2021/generated/media_db.tsv"
    
    if [ -n "$growth" ]; then
        cmd+=" --growth $growth"
    fi
    
    cmd+=" -o \"$output_dir\" --sample $samples --exchange-detection '^R_EX_(?!nad_c$|thf_c$).*' --unlimited ../data/Guillen2021/generated/unlimited.txt"
    
    time eval "$cmd"
}

# Run SteadierCom with abundance and growth

input_dir="../data/Guillen2021/generated/communities_38.tsv"

output_dir="../output/Guillen2021/results_steadiercom/results_38_growth/ecoli_mutants_38_10"


samples=2
growth=0.221

run_steadiercom "$input_dir" "$samples" "$output_dir" "$growth"

