#!/usr/bin/env bash

set -eux

steadiercom ../models/species_used/* -c ../data/Venturelli2018/generated/communities_12_exp_abundance_part2.tsv -m LB_extend --mediadb ../data/Venturelli2018/generated/media_db_lb_extend.tsv --growth 0.039683 -o ../output/Venturelli2018/results_communities_12_1000s_part2 --sample 1000 --exchange-detection 'R_EX_(?!biomass__40__e__41__)'