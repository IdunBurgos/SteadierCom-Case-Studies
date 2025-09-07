import pandas as pd
import numpy as np
import json

import reframed

import sys
sys.path.append("../functions/")
import general
from minimal_medium import minimal_medium

import logging
logging.basicConfig(
    level=logging.INFO, 
    format="%(levelname)-8s | %(message)s",
)

def compounds_of_interest(ex_rxns_lb):

    ex_rxns_lb = ["R_EX_"+met[2:].replace("91","40").replace("93","41") for met in ex_rxns_lb if "__e__" in met]

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
    rxns = ex_rxns_lb+rxns_media_composition+rxns_amino_acids_peptone+rxns_nucleotides_yeast_extract+rxns_other_important_compounds
    rxns_exclude = ["R_EX_h2__40__e__41__","R_EX_h2s__40__e__41__","R_EX_o2__40__e__41__"] # not possible without including h2!!
    
    return list(set(rxns)),rxns_exclude

def make_bigg2agora_dict(models):
    mets_id_ag = []
    mets_id_bigg = []
    ex_rxns_ag = []
    for model in models.values():
        for rxn in model.get_exchange_reactions():
            if "R_EX_" not in rxn:
                continue
            met_id_ag = model.reactions[rxn].get_substrates()[0]
            mets_id_ag.append(met_id_ag)
            met_id_bigg = model.metabolites[met_id_ag].metadata["bigg.metabolite"][0]
            mets_id_bigg.append(met_id_bigg)
            ex_rxns_ag.append(rxn)

    met2ex = dict(zip(mets_id_bigg,ex_rxns_ag))
    return met2ex

def load_data(models):
    
    logging.info("Loading data and LB medium from carveme media database...")
    
    with open("../data/Venturelli2018/generated/metid2metname.json") as json_file:
        metid2metname = json.load(json_file)
   
    # Load metabolites from carveme dataset
    media_db = pd.read_csv("https://raw.githubusercontent.com/cdanielmachado/carveme/refs/heads/master/carveme/data/input/media_db.tsv",sep="\t")
    media_db_lb = media_db[media_db.medium=="LB[-O2]"]
    mets_lb = media_db_lb.compound.values

    # Make bigg to agora translation dict
    met2ex = make_bigg2agora_dict(models)
    ex_rxns_lb = [met2ex[met] for met in mets_lb if met in met2ex.keys()]

    # Make environment
    lb_env = reframed.Environment.from_reactions(ex_rxns_lb)

    # Load preferred reactions and excluded reactions
    preferred_rxns,rxns_exclude = compounds_of_interest(ex_rxns_lb)
    
    return lb_env,preferred_rxns,rxns_exclude


def find_optimized_media_solutions(model,initial_env,preferred_reactions,excluded_reactions):
    """
    This function finds the minimum number of compounds to add to the medium in order for all models to grow. 
    The base of the medium (initial_env) given in the function and you can define preferred reactions to include and 
    unwanted reactions to exclude.
    """

    initial_env.apply(model,inplace=True,exclusive=True,warning=False)
    sol = reframed.FBA(model)


    log = ""
    if sol.fobj>1e-6:
        log += f"\t No need for extra compounds for {model.id}, growth: {sol.fobj}"

        log += "\n"
        return pd.Series()
    else:
        log += f"\t {model.id} cannot grow: {sol.fobj}"
    
    # Define complete medium for model - exclude some reactions
    env_complete = reframed.Environment.complete(model)
    env_complete.apply(model,inplace=True,exclusive=True,warning=False)

    rxns_complete = env_complete.keys()
    # Some reactions are not exchange reactions even though they are identified as such -> filter these out
    rxns_complete = list(set([e for e in rxns_complete if ("R_DM" not in e) and ("_acald_" not in e)]))
    rxns_complete = [rxn for rxn in rxns_complete if rxn not in excluded_reactions]
    
    # Run minimal medium and gather results
    solutions = minimal_medium(model,
                               exchange_reactions=list(rxns_complete),
                               preferred_reactions=preferred_reactions,
                               min_growth=0.1,
                               n_solutions=1000,
                               min_mass_weight=True,
                               max_uptake=5,
                               milp=True)
    for sol in solutions[1]:
        log += f"\t Minimal media fobj: {sol.fobj}"
        break

    # Count number of solutions where exchange rxn appears        
    sol_count = dict()
    for sol in solutions[0]:
        for rxn in sol:
            if rxn in sol_count.keys():
                sol_count[rxn]+=1
            else:
                sol_count[rxn]=1
    logging.info(log)
    
    sol_count_df = pd.Series(sol_count).sort_values(ascending=False)
    return sol_count_df


def find_optimized_medium(model,sol_count_df,initial_env,min_growth=1e-3):
    """ 
    From the results in sol_count_dfs use adjustable cutoff to minimize added compounds + run FBA to confirm growth.
    """

    if sol_count_df.empty:
        return []

    samples = sol_count_df.max()
    cutoff_frac = np.linspace(0,1,11)[::-1] 
    cutoff_values = cutoff_frac*samples


    log = ""
    for cutoff in cutoff_values:
        missing = list((set(sol_count_df[sol_count_df>=cutoff].index)-set(initial_env.keys())))
        
        initial_env_copy = initial_env.copy()
        env_dict_add = dict(zip(missing,[(-10,np.inf)]*len(missing)))
        initial_env_copy.update(env_dict_add)
        grows, log_growth = general.test_growth_model_venturelli2018(model,
                                                         env=initial_env_copy,
                                                         detailed_log=False,
                                                         min_growth=min_growth) 

        if grows:
            log += f"\t {log_growth} at cutoff {cutoff}\n"
            break
        else:
            log += f"\t No solution for cutoff: {cutoff}\n"
            if cutoff==0:
                raise Exception("\t There should be a solution when all compounds are included...")
    
    logging.info(log)

    return missing

def find_optimized_medium_all(models,lb_env, preferred_rxns,rxns_exclude,min_growth=1e-3):
    """ 
    FIND MINIMAL MEDIA COMPONENTS FOR EACH MODEL 
    This part of the script iterates through the models and find the minimum number of compounds to add to the medium
    in order for all models to grow. The base of the medium is LB-medium defined in carveme. The compounds are found 
    using a version of minimal_medium from reframed that includes a "preferred_reactions" option.
    """
    logging.info("Finding optimized medium for all models...")

    missing_total = []
    sol_count_dfs = {}

    for model_id,model in models.items():
        logging.info(f"Started finding medium for {model.id}")

        sol_count_df = find_optimized_media_solutions(model,
                                                       initial_env=lb_env,
                                                       preferred_reactions=preferred_rxns,
                                                       excluded_reactions=rxns_exclude,
                                                       )
        sol_count_dfs[model.id]=sol_count_df

        missing = find_optimized_medium(model,
                                        sol_count_df,
                                        initial_env=lb_env,
                                        min_growth=min_growth)
       # Add missing reactions
        missing_total.extend(missing)

    # This is the new media components
    add2lb = list(set(missing_total))

    return add2lb, sol_count_dfs


def main():
    
    logging.info("STARTING SCRIPT FOR FINDING MINIMAL MEDIA FOR COMMUNITY")

    models = general.load_models_venturelli2018()

    lb_env,preferred_rxns,rxns_exclude = load_data(models)

    add2lb, sol_count_dfs = find_optimized_medium_all(models, lb_env, preferred_rxns, rxns_exclude)

    rxns_venturelli = list(lb_env.keys())+add2lb


    logging.info("SAVE media_db.tsv in ../data/Venturelli2018/generated/")
    media_db = pd.DataFrame({"medium":["LB_extend"]*len(rxns_venturelli),
                             "compound":[met.replace("_40_","_91_").replace("_41_","_93_") for met in rxns_venturelli]})

    media_db.to_csv("../data/Venturelli2018/generated/media_db.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()