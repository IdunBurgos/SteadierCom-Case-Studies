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
    
    cmd="steadiercom ../models/Venturelli2018/species_used_test/* -c \"$community_dir\" -m LB_extend --mediadb ../data/Venturelli2018/generated/media_db.tsv"
    
    if [ -n "$growth" ]; then
        cmd+=" --growtimport logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)h $growth"
    fi
    
    cmd+=" -o \"$output_path\" --sample $samples --exchange-detection 'R_EX_(?!biomass__40__e__41__)' --unlimited ../data/Venturelli2018/generated/unlimited.txt"
    
    time eval "$cmd"
}

samples=100
num_processes=2
growth=0.125


# Run SteadierCom with abundance and growth
input_dir="../data/Venturelli2018/generated/communities_steadiercom"
output_dir="../output/Venturelli2018/results_steadiercom/results_growth_abundance_test"
i=0

for dir in "$input_dir"/*; do
    ((i++))
    run_steadiercom "$dir" "$samples" "$output_dir" "$growth"&
    if (( i % num_processes == 0 )); then
        wait  # Wait every 4 processes before continuing
    fi
done

wait




