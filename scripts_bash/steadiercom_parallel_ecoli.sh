source ~/opt/anaconda3/etc/profile.d/conda.sh
conda activate steadiercom

run_steadiercom() {
    community_dir="$1"
    samples="$2"
    output_dir="$3"
    growth="$4"
    
    
    
    cmd="steadiercom ../models/Guillen2021/top_7_strains/* -c \"$community_dir\" -m M9 --mediadb ../data/Guillen2021/generated/media_db.tsv"
    
    if [ -n "$growth" ]; then
        cmd+=" --growth $growth"
    fi
    
    cmd+=" -o \"$output_dir\" --sample $samples --exchange-detection '^R_EX_(?!nad_c$|thf_c$).*' --unlimited ../data/Guillen2021/generated/unlimited.txt"
    
    time eval "$cmd"
}

# Run SteadierCom with abundance and growth

input_dir1="../data/Guillen2021/generated/communities_steadiercom.tsv"
input_dir2="../data/Guillen2021/generated/communities_smetana.tsv" # No abundance
output_dir1="../output/Guillen2021/results_steadiercom/results_growth_abundance/ecoli_mutants"
output_dir2="../output/Guillen2021/results_steadiercom/results_abundance/ecoli_mutants"
output_dir3="../output/Guillen2021/results_steadiercom/results_growth/ecoli_mutants"

samples1=2000
samples2=4000

growth=0.221

run_steadiercom "$input_dir1" "$samples1" "$output_dir1" "$growth"&
run_steadiercom "$input_dir1" "$samples1" "$output_dir2"&
wait

run_steadiercom "$input_dir2" "$samples2" "$output_dir3" "$growth" 



