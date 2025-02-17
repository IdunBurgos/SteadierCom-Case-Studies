#!/usr/bin/env bash

set -eux

steadiercom ../models/species_used_noO2/* -c ../data/Venturelli2018/generated/communities_NONE_abundance.tsv -m LB_extend --mediadb ../data/Venturelli2018/generated/media_db_lb_extend.tsv --growth 0.039683 -o ../output/Venturelli2018/results_community_NONE_100s/results_communities_12_100s_1 --sample 100 --exchange-detection 'R_EX_(?!biomass__40__e__41__)'


steadiercom ../models/species_used_noO2/* -c ../data/Venturelli2018/generated/communities_NONE_abundance.tsv -m LB_extend --mediadb ../data/Venturelli2018/generated/media_db_lb_extend.tsv --growth 0.039683 -o ../output/Venturelli2018/results_community_NONE_100s/results_communities_12_100s_2 --sample 100 --exchange-detection 'R_EX_(?!biomass__40__e__41__)'

steadiercom ../models/species_used_noO2/* -c ../data/Venturelli2018/generated/communities_NONE_abundance.tsv -m LB_extend --mediadb ../data/Venturelli2018/generated/media_db_lb_extend.tsv --growth 0.039683 -o ../output/Venturelli2018/results_community_NONE_100s/results_communities_12_100s_3 --sample 100 --exchange-detection 'R_EX_(?!biomass__40__e__41__)'


steadiercom ../models/species_used_noO2/* -c ../data/Venturelli2018/generated/communities_NONE_abundance.tsv -m LB_extend --mediadb ../data/Venturelli2018/generated/media_db_lb_extend.tsv --growth 0.039683 -o ../output/Venturelli2018/results_community_NONE_100s/results_communities_12_100s_4 --sample 100 --exchange-detection 'R_EX_(?!biomass__40__e__41__)'


steadiercom ../models/species_used_noO2/* -c ../data/Venturelli2018/generated/communities_NONE_abundance.tsv -m LB_extend --mediadb ../data/Venturelli2018/generated/media_db_lb_extend.tsv --growth 0.039683 -o ../output/Venturelli2018/results_community_NONE_100s/results_communities_12_100s_5 --sample 100 --exchange-detection 'R_EX_(?!biomass__40__e__41__)'