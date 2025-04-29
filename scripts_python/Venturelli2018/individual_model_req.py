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
# in order for all models to grow. The base of the medium is LB-medium defined in carveme The compounds are found 
# using a version of minimal_medium from reframed that includes a "preferred_reactions" option 



missing_total = []
sol_count_dfs = {}
rxns_exclude = ["R_EX_h2__40__e__41__","R_EX_h2s__40__e__41__","R_EX_o2__40__e__41__"] # not possible without including h2!!
    

for model_id,model in models.items():
    print(f"Started finding medium for {model_id}")
    lb_env.apply(model,inplace=True,exclusive=True,warning=False)
    sol = reframed.FBA(model)
    
    if sol.fobj>1e-6:
        print(f"\t No need for extra compounds for {model_id}, growth: {sol.fobj}")
    
        print("\n")
        continue
    else:
        print(f"\t {model_id} cannot grow: {sol.fobj}")
    
    # Define complete medium for model - exclude some reactions
    env_complete = reframed.Environment.complete(model)
    env_complete.apply(model,inplace=True,exclusive=True,warning=False)

    rxns_complete = env_complete.keys()
    rxns_complete = list(set([e for e in rxns_complete if ("R_DM" not in e) and ("_acald_" not in e)])) #("R_DM" not in e) & ("sink" not in e)]
    rxns_complete = [rxn for rxn in rxns_complete if rxn not in rxns_exclude]
    
    # Run minimal medium and gather results
    solutions = minimal_medium(model,exchange_reactions=list(rxns_complete),preferred_reactions=rxns,min_growth=0.1,n_solutions=1000,use_pool=True,min_mass_weight=True)

    for sol in solutions[1]:
        print(f"\t Minimal media fobj: {sol.fobj}")
        break
        
    sol_count = dict()
    for sol in solutions[0]:
        for rxn in sol:
            if rxn in sol_count.keys():
                sol_count[rxn]+=1
            else:
                sol_count[rxn]=1
    
    sol_count_df = pd.Series(sol_count).sort_values(ascending=False)
    sol_count_dfs[model_id]=sol_count_df
    
    # Process results: Use an adjustable cutoff to minimize added compounds + run FBA to confirm growth
    cutoff_values = np.linspace(0,1000,11)[::-1] # going in reverse order from 1000 to 0
    
    for cutoff in cutoff_values:
        missing = list((set(sol_count_df[sol_count_df>=cutoff].index)-set(lb_env.keys())))

        lb_extended = [met2ex[met] for met in compounds if met in met2ex.keys()]
        lb_extended.extend(missing)
        lb_extended.extend(rxns)
        lb_extended = list(set(lb_extended))
        env_lb_extend = reframed.Environment.from_reactions(lb_extended)
        env_lb_extend.apply(model,inplace=True,exclusive=True,warning=False)
        sol = reframed.FBA(model)
        
        if sol.fobj>1e-6:
            print(f"\t {model_id} growth at cutoff {cutoff}, growth: {sol.values['R_EX_biomass__40__e__41__']}\n")
            break
        else:
            print(f"\t No solution for cutoff: {cutoff}")
            if cutoff==0:
                raise Exception("\t There should be a solution when all compounds are included...")
            
    
    # Add missing reactions
    missing_total.extend(missing)

# This is the new media components
add2lb = list(set(missing_total))