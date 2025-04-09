#!/usr/bin/env bash
# E.coli mutants

echo "Global analysis for E.coli mutants"
smetana "../models/Guillen2021/top_7_strains/*" -c "../data/Guillen2021/generated/communities_smetana.tsv" -m "M9" --mediadb "../data/Guillen2021/generated/media_db.tsv" -g -o "../output/Guillen2021/smetana/" --exchange-detection '^R_EX_(?!nad_c$|thf_c$).*'

echo "Detailed analysis for E.coli mutants"
smetana "../models/Guillen2021/top_7_strains/*" -c "../data/Guillen2021/generated/communities_smetana.tsv" -m "M9" --mediadb "../data/Guillen2021/generated/media_db.tsv" -d -o "../output/Guillen2021/smetana/" --exchange-detection '^R_EX_(?!nad_c$|thf_c$).*'

echo "Global analysis for Gut microbiome"
smetana "../models/Venturelli2018/species_used_noO2/*" -c "../data/Venturelli2018/generated/communities_smetana.tsv" -m "LB_extend" --mediadb "../data/Venturelli2018/generated/media_db.tsv" -g -o "../output/Venturelli2018/smetana/" --exchange-detection '^R_EX_(?!biomass__40__e__41__$).*'

echo "Detailed analysis for Gut microbiome"
smetana "../models/Venturelli2018/species_used_noO2/*" -c "../data/Venturelli2018/generated/communities_smetana.tsv" -m "LB_extend" --mediadb "../data/Venturelli2018/generated/media_db.tsv" -d -o "../output/Venturelli2018/smetana/" --exchange-detection '^R_EX_(?!biomass__40__e__41__$).*' --no-coupling