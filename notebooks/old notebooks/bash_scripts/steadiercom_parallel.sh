#!/usr/bin/env bash


run_steadiercom() {
    community_dir="$1"
    samples="$2"
    
    community_id=$(echo "$community_dir" | rev | cut -d "/" -f1 | rev)
    echo "../../output/Venturelli2018/results_steadiercom/${community_id%_abundance.tsv}"
    
    time steadiercom ../../models/species_used/* -c "$community_dir" -m LB_extend --mediadb ../../data/Venturelli2018/generated/media_db_lb_extend.tsv --growth 0.039683 -o "../../output/Venturelli2018/results_steadiercom/${community_id%.tsv}"  --sample $samples --exchange-detection 'R_EX_(?!biomass__40__e__41__)' | tee -a "../../output/Venturelli2018/results_steadiercom/time_log.txt"

}


num_processes=4
input_dir="../../data/Venturelli2018/generated/communities_exp_abundance"
samples=$1
i=0

for dir in "$input_dir"/*; do
    ((i++))
    run_steadiercom "$dir" "$samples" &
    if (( i % num_processes == 0 )); then
        wait  # Wait every 4 processes before continuing
    fi
done

wait



