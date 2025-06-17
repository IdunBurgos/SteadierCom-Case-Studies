sys.path.append("../functions/")
from minimal_medium import find_optimized_medium


## SCRIPT FOR FINDING MINIMAL MEDIA FOR COMMUNITY

# NB: This script needs access to variables from the notebook Step 1.1.

## DEFINE COMPOUNDS WE ARE INTERESTED IN

rxns = ["R_EX_"+met[2:].replace("91","40").replace("93","41") for met in metid2metname.keys() if "__e__" in met]

# Following reactions added based on media composition: http://www.oxoid.com/UK/BLUE/prod_detail/prod_detail.asp?pr=CM0957
rxns_media_composition = [
    "R_EX_malt__40__e__41__",
    "R_EX_cellb__40__e__41__",
    "R_EX_succ__40__e__41__",
    "R_EX_mqn8__40__e__41__",
    "R_EX_arg_L__40__e__41__",
    "R_EX_cys_L__40__e__41__",
    "R_EX_pyr__40__e__41__"
]

# Amino acids in peptone: https://www.usbio.net/media/P3300/peptone
rxns_amino_acids_peptone = [
    "R_EX_asn_L__40__e__41__",
    "R_EX_ala_L__91__e__93__",
    "R_EX_arg_L__91__e__93__",
    "R_EX_asp_L__40__e__41__",
    "R_EX_cys_L__40__e__41__",
    "R_EX_gly__40__e__41__",
    "R_EX_his_L__40__e__41__",
    "R_EX_ile_L__40__e__41__",
    "R_EX_glu_L__40__e__41__",
    "R_EX_leu_L__40__e__41__",
    "R_EX_lys_L__40__e__41__",
    "R_EX_met_L__40__e__41__",
    "R_EX_phe_L__40__e__41__",
    "R_EX_pro_L__40__e__41__",
    "R_EX_ser_L__40__e__41__",
    "R_EX_thr_L__40__e__41__",
    "R_EX_trp_L__40__e__41__",
    "R_EX_tyr_L__40__e__41__",
    "R_EX_val_L__40__e__41__",
    "R_EX_gln_L__40__e__41__"
]

# Nucleotides from yeast extract https://pmc.ncbi.nlm.nih.gov/articles/PMC9998214/
rxns_nucleotides_yeast_extract = [
    "R_EX_ade__40__e__41__",
    "R_EX_csn__40__e__41__",
    "R_EX_gua__40__e__41__",
    "R_EX_thym__40__e__41__",
    "R_EX_ura__40__e__41__"
]

# Other important compounds https://pmc.ncbi.nlm.nih.gov/articles/PMC9998214/
rxns_other_important_compounds = [
    "R_EX_thm__40__e__41__",
    "R_EX_ribflv__40__e__41__",
    "R_EX_nac__40__e__41__",
    "R_EX_ncam__40__e__41__",
    "R_EX_pnto_R__40__e__41__",
    "R_EX_pydxn_R__40__e__41__",
    "R_EX_btn__40__e__41__",
    "R_EX_fol__40__e__41__",
    "R_EX_cbl1__40__e__41__"
]

rxns = rxns+rxns_media_composition+rxns_amino_acids_peptone+rxns_nucleotides_yeast_extract+rxns_other_important_compounds


## FIND MINIMAL MEDIA COMPONENTS FOR EACH MODEL 
# This part of the script iterates through the models and find the minimum number of compounds to add to the medium
# in order for all models to grow. The base of the medium is LB-medium defined in carveme. The compounds are found 
# using a version of minimal_medium from reframed that includes a "preferred_reactions" option.


missing_total = []
sol_count_dfs = {}
rxns_exclude = ["R_EX_h2__40__e__41__","R_EX_h2s__40__e__41__","R_EX_o2__40__e__41__"] # not possible without including h2!!
    

for model_id,model in models.items():
    missing,sol_count_dfs = find_optimized_medium(model,
                                    sol_count_dfs,
                                    initial_env=lb_env,
                                    preferred_reactions=rxns,
                                    excluded_reactions=rxns_exclude,
                                    biomass_reaction='R_EX_biomass__40__e__41__')
    # Add missing reactions
    missing_total.extend(missing)

# This is the new media components
add2lb = list(set(missing_total))