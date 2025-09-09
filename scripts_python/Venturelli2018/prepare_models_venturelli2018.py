import reframed
import pandas as pd

from os import listdir
import copy
import logging

import sys

sys.path.append("../functions/")
from minimal_medium import minimal_medium

from EGC import EGC_identifier
from general import test_growth_models_venturelli2018

import logging
logging.basicConfig(
    level=logging.INFO,  # change to DEBUG for more detail
    format="%(levelname)-8s | %(message)s",
)


def change_models_remove_oxygen(models):
    """
    Models are working in an anaerobe environment so all oxygen-related exchange should be avoided. Relevant reactions and metabolites are therefore removed from the models. 
    """
    
    models_new = {}
    changed_elements = {}

    for name,model in models.items():
        changed_elements[name] = {"added":[], "removed":[]}

        for o2_met in ["M_o2__91__e__93__","M_o2__91__p__93__"]:
            if o2_met in model.metabolites:
                remove_reactions = model.get_metabolite_reactions(o2_met)
                
                
                model.remove_reactions(remove_reactions)
                model.remove_metabolite(o2_met)

                changed_elements[name]["removed"].extend(remove_reactions)
                changed_elements[name]["removed"].append(o2_met)
        
        model.update()
        models_new[name]=model
        
    return models_new,changed_elements
        

def change_models_remove_EGCs(models):
    """
    Standardize formate and fructose transport in models and remove reactions related to bile, antibiotic and other drug transport. Focused on metabolites showing ability to carry flux in empty environemnt. 
    """
    
    R_FORt2r = copy.deepcopy(models["Eggerthella_lenta_DSM_2243"].reactions.R_FORt2r)
    R_FRUt2r = copy.deepcopy(models["Faecalibacterium_prausnitzii_A2_165"].reactions.R_FRUt2r)

    R_FORt2r.set_gpr_association(reframed.GPRAssociation())
    R_FRUt2r.set_gpr_association(reframed.GPRAssociation())

    models_new = {}

    # Mets that are a part of bile, antibiotics or other drugs. 
    mets_cycling = ["M_nzp__91__e__93__","M_anzp__91__e__93__","M_czp__91__e__93__","M_7a_czp__91__e__93__",
                    "M_nchlphncl__91__e__93__","M_cholate__91__e__93__","M_12dhchol__91__e__93__", "M_7dhcdchol__91__e__93__",
                    "M_3dhcdchol__91__e__93__","M_sanilamide__91__e__93__","M_tab__91__e__93__","M_pront__91__e__93__",
                    "M_olsa__91__e__93__","M_neopront__91__e__93__","M_chlphncl__91__e__93__","M_bzd__91__e__93__",
                    "M_5asa__91__e__93__","M_HC02191__91__e__93__","M_alchac__91__e__93__","M_dchac__91__e__93__"]
    
    replace_reactions = {"R_FORt":R_FORt2r,
                         "R_FRUpts":R_FRUt2r,
                         }
    remove_reactions = ["R_FDNADOX_H","R_FDNADOX_Hpp"]
    
    changed_elements = {}

    for member,model in models.items(): 
        changed_elements[member]={"added":[],"removed":[]}
        
        models_new[member]=model.copy()
        
        for rxn_id_rem,rxn_new in replace_reactions.items():
            if rxn_id_rem in models_new[member].reactions:
                try:
                    models_new[member].remove_reaction(rxn_id_rem)
                    models_new[member].add_reaction(copy.deepcopy(rxn_new))
                    changed_elements[member]["removed"].append(rxn_id_rem)
                    changed_elements[member]["added"].append(rxn_new.id)
                except Exception as e:
                    logging.exception("Failed to replace %s with %s for %s: %s", rxn_id_rem, rxn_new.id, member, e)

                
        for rxn_id_rem in remove_reactions:
            if rxn_id_rem in models_new[member].reactions:
                try:
                    models_new[member].remove_reaction(rxn_id_rem)
                    changed_elements[member]["removed"].append(rxn_id_rem)
                except Exception as e:
                    logging.exception("Failed to remove %s for %s: %s", rxn_id_rem, member, e)
 
        for met in mets_cycling:
            met_p = copy.copy(met)
            met_p = met_p.replace("_e_","_p_")

            if met in models_new[member].metabolites:
                rxns_remove = models_new[member].get_metabolite_reactions(met)
                models_new[member].remove_reactions(rxns_remove)
                models_new[member].remove_metabolite(met)
                changed_elements[member]["removed"].extend(rxns_remove)

            if met_p in models_new[member].metabolites:
                rxns_remove = models_new[member].get_metabolite_reactions(met_p)
                models_new[member].remove_reactions(rxns_remove)
                models_new[member].remove_metabolite(met_p)     
                changed_elements[member]["removed"].extend(rxns_remove)

        env_complete = reframed.Environment.complete(models_new[member])
        env_complete.apply(models_new[member],inplace=True)
        reframed.simplify(models_new[member],inplace=True)
        
    
    return models_new, changed_elements

def atp_EGC_community(models,get_EGC_data=False):
    logging.info("Testing for ATP-generating cycles in community model...")
    comm_test = reframed.Community('EGC_test',models.values())
    comm_model = comm_test.merged_model
    env_empty = reframed.Environment.empty(comm_model)
    env_empty.apply(comm_model,inplace=True)


    for member in comm_test.organisms.keys():
        rxn = "R_DM_atp_c__"+member
        comm_model.reactions[rxn].lb=0
        comm_model.reactions[rxn].ub=1000

    rxns_members = {}    
    
    atp_sol = {}
    for member in comm_test.organisms.keys():
        rxn = "R_DM_atp_c__"+member

        objective = {rxn:0 for rxn in comm_model.reactions}
        objective[rxn]=1
        sol = reframed.pFBA(comm_model,objective=objective)
        
        atp_sol[rxn] =sol.values[rxn]
        
        rxns = [key for key,value in sol.values.items() if abs(value)>1e-6]
        rxns_members[member]=rxns
    
    
    atp_string = "".join([f"{rxn} = {value}\n" for rxn,value in atp_sol.items()])
    logging.info("\nATP demand optimum by member:\n%s",atp_string
                 )
    # Safe intersection
    common = set.intersection(*map(set, rxns_members.values())) if rxns_members else set()
    logging.info(
        "Reactions active in ALL ATP-production solutions across species: %s",
        sorted(common)
    )
    if get_EGC_data:
        return rxns_members,common

def load_data_models():
    """
    Specific function
    """
    logging.info("load media data and models from '../models/Venturelli2018/species_used/'...")

    models = {}
    
    for filename in listdir("../models/Venturelli2018/species_used/"):
        if filename.endswith(".xml"):
            models[filename.replace(".xml","")] = reframed.load_cbmodel("../models/Venturelli2018/species_used/"+filename,exchange_detection='^R_EX_(?!biomass__40__e__41__$).*')
            
    return models



def analyze_EGC_models(models,get_EGC_data):
    """
    Find infeasible cycles in single models and ATP-generating cycles in the community.
    Returns (rxns_members, common_reactions).
    """
    rxns_members = {}
    # EGC-single models
    logging.info("Testing infeasible cycles in single models (EGC_identifier)...")
    log_EGC_single = ""
    EGCs_results = {}
    for member, model in models.items():
        try:
            EGCs,log = EGC_identifier(model, agora=True,send_log=True)
            EGCs_results[member] = EGCs
            log_EGC_single = log_EGC_single+f"EGC_identifier on {member}: "+log
        except Exception as e:
            logging.exception("EGC_identifier failed for %s: %s", member, e)
            
    logging.info(f"Results: \n {log_EGC_single}")
    if get_EGC_data:
        return EGCs_results


def analyze_models(models, media_db=False,get_EGC_data=False):
    """
    Analyze the given models using the provided media database or use a complete environment.
    """

    if media_db:
        rxns_or_mets = media_db.compound.values
        is_rxns = rxns_or_mets[0].startswith("R_EX")
        if is_rxns:
            env = reframed.Environment.from_reactions(rxns_or_mets)
        else:
            env = reframed.Environment.from_metabolites(rxns_or_mets)
        complete_env = False
    else:
        complete_env = True
        env = False

    test_growth_models_venturelli2018(models,
                                        env=env,
                                        complete_env=complete_env)

    EGCs_single_results = analyze_EGC_models(models,get_EGC_data=get_EGC_data)
    EGCs_community_results = atp_EGC_community(models,get_EGC_data=get_EGC_data)

    if get_EGC_data:
        return {"community": {"rxns_members": EGCs_community_results[0], 
                              "common": EGCs_community_results[1]},
                "single": EGCs_single_results}


def main():
    logging.info("LOAD DATA AND MODELS")
    models = load_data_models()
        
    logging.info("Make community model")
    comm_test = reframed.Community('EGC_test',models.values())

    logging.info("ANALYZE CURRENT MODELS")
    
    analyze_models(models,get_EGC_data=False)
    
    logging.info("CHANGE MODELS: Remove reactions+metabolites in EGCs")
    models_new,changed_elements = change_models_remove_EGCs(models)
    
    logging.info("CHANGE MODELS: Remove reactions+metabolites in oxygen exchange")
    models_new2,changed_elements2 = change_models_remove_oxygen(models_new)
    
    
    logging.info("ANALYZE NEW MODELS")
    analyze_models(models_new2,get_EGC_data=False)

    logging.info("SAVE MODELS AND DATA")

    for name,model in models_new2.items():

        reframed.save_cbmodel(model,"../models/Venturelli2018/species_used_noO2/"+name+".xml")


if __name__ == "__main__":
    
    main()