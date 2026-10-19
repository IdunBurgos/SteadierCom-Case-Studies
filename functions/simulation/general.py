import pandas as pd
import reframed
from os import listdir
import logging
import warnings

logging.basicConfig(level=logging.INFO, 
    format="%(levelname)-8s | %(message)s",
)

def load_models_venturelli2018():
    """
    Load the models from the Venturelli et al (2018) study.
    """
    logging.info("Loading models for Venturelli et al (2018)...")
    
    models = {}
    for filename in listdir("../models/Venturelli2018/species_used_noO2/"):
        if filename.endswith(".xml"):
            name = filename[:-4]  # Remove .xml extension
            model = reframed.load_cbmodel(f"../models/Venturelli2018/species_used_noO2/{filename}", 
                                          exchange_detection="R_EX_(?!biomass__40__e__41__)")
            models[name] = model
    return models


def load_models_guillen2021(all_38_strains=True):
    """
    Load models from the Guillen et al (2021) study.
    """
    logging.info("Loading models for Guillen et al (2021)...")

    if all_38_strains:
        sub_dir="strains"
    else:
        sub_dir="top_7_strains"

    models = {}
    for filename in listdir(f"../models/Guillen2021/{sub_dir}/"):
        if filename.endswith(".xml"):
            name = filename[:-4]  # Remove .xml extension
            model = reframed.load_cbmodel(f"../models/Guillen2021/{sub_dir}/{filename}",
                                          exchange_detection="^R_EX_(?!nad_c$|thf_c$).*")
            models[name] = model
    return models

def load_ecoli_model():
    """
    Load E.coli model saved in ../models/Guillen2021/iML1515.xml
    """

    model_ecoli = reframed.load_cbmodel("../models/Guillen2021/iML1515.xml")
    return model_ecoli


def test_growth_model(model, biomass_reaction,env=False,only_return_log=False,min_growth=1e-3,detailed_log=True):
    if env is False:
        env = reframed.Environment.complete(model)
        if detailed_log:
            logging.info("No environment provided, using complete environment for the model.")

    log_growth = ""
    model_copy = model.copy()
    assert isinstance(model_copy, reframed.CBModel), f"Expected model to be a CBModel, but got {type(model_copy)}"
    assert biomass_reaction in model_copy.reactions, f"Expected biomass reaction '{biomass_reaction}' to be in model reactions"

    objective={rxn:0 for rxn in model_copy.reactions.keys()}
    objective[biomass_reaction] = 1

    try: 
        env.apply(model_copy,inplace=True,exclusive=True,warning=False)
    except Exception as e:
        logging.exception("Setting environment failed for %s: %s", model.id, e)
    try:
        sol = reframed.FBA(model_copy,objective=objective)
        log_growth += f"{model.id}: {biomass_reaction} = {sol.fobj} \n"   
    except Exception as e:
        logging.exception("FBA failed for %s: %s", model.id, e)
    
    if only_return_log:
        return log_growth

    return sol.fobj>min_growth,log_growth

def test_growth_model_venturelli2018(model,env=False,complete_env=False,only_return_log=False,min_growth=1e-3,detailed_log=True):
    if detailed_log:
        logging.info(f"Testing growth for model {model.id}...")
    
    if (not env) and (not complete_env):
        # Use env for the venturelli models
        media_db = pd.read_csv("../data/Venturelli2018/generated/media_db.tsv",sep="\t")
        env_rxns = media_db.compound.map(lambda x:x.replace("91","40").replace("93","41")).values
        env = reframed.Environment.from_reactions(env_rxns)
        if detailed_log:
            logging.info("No environment provided, using LB_extend media for case study.")

    if complete_env:
        if env:
            warnings.warn("Using complete_env and overwriting provided env.")
        env = False


    return test_growth_model(model, 
                             biomass_reaction='R_EX_biomass__40__e__41__', 
                             env=env, 
                             only_return_log=only_return_log,
                             min_growth=min_growth,
                             detailed_log=detailed_log)

def test_growth_models_venturelli2018(models,env=False,complete_env=False,min_growth=1e-3):
    logging.info("Check growth of single models")
    if (not env) and (not complete_env):
        logging.info("No environment provided, using LB_extend media for case study.")
    log_growth = ""
    for model in models.values():

        log_growth += test_growth_model_venturelli2018(model,
                      env=env,
                      complete_env=complete_env,
                      only_return_log=True,
                      min_growth=min_growth,
                      detailed_log=False)
    logging.info(f"Results: \n {log_growth}")


