import pandas as pd
import sys
sys.path.append("../functions/")
from simulation.general import load_models_venturelli2018

import logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

def load_chebi_data():

    """
    Loads the ChEBI class data for metabolites and returns dataframe with BIGG id as index

    Returns
    -------
    pandas.DataFrame
        DataFrame containing ChEBI class information for metabolites.
    """
    met_chebi_class= pd.read_csv("../data/met_chebi_class.tsv",sep="\t",index_col=0)
    extract = lambda s: s[s.index("_")+1:s.rindex("_")] #extract bigg id from the format M_bigg_id_e

    met_chebi_class.index = met_chebi_class.index.map(extract)
    met_chebi_class.rename_axis("bigg_id",axis=0,inplace=True)

    return met_chebi_class


def find_mets_bigg_id(models):
    """
    Finds the bigg IDs for metabolites in the metadata of the given models.

    Returns
    -------
    met_bigg_dict: A dictionary mapping metabolite IDs to their BIGG IDs
    no_bigg_annotation: A list of metabolite IDs without BIGG annotations
    """

    # Create a mapping of metabolite IDs to their BIGG IDs (met_bigg_dict)
    met_bigg_dict = {} # dict for storing mapping
    no_bigg_annotation = [] # list for storing metabolites without BIGG annotation

    for model in models.values(): #for each model
        for met in model.get_external_metabolites(): #for each metabolite

            if met in met_bigg_dict:
                continue 
            
            if "bigg.metabolite" in model.metabolites[met].metadata:

                if model.metabolites[met].metadata["bigg.metabolite"][0]!="NaN":
                    met_bigg_dict[met] = model.metabolites[met].metadata["bigg.metabolite"][0]
                else:
                    no_bigg_annotation.append(met)
            else:
                no_bigg_annotation.append(met)

    return met_bigg_dict, list(set(no_bigg_annotation))


def make_met_chebi_class_agora_df(met_bigg_dict, met_chebi_class):
    """
    Map AGORA metabolite IDs to BiGG IDs and join ChEBI class information.

    Parameters
    ----------
    met_bigg_dict : dict
        Mapping from AGORA metabolite ID -> BiGG metabolite ID.
    met_chebi_class : pandas.DataFrame
        DataFrame indexed by BiGG ID containing ChEBI class information. 
        
    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by AGORA met IDs. Contains mapping for all external metabolites with bigg identifier mapping in the AGORA models.

        Columns: 'bigg id', 'chebi id', 'chebi class', 'self defined super class'

    """
    met_bigg_df = (
        pd.Series(met_bigg_dict, name="bigg_id")
        .rename_axis("agora_id")
        .reset_index()
        .set_index("bigg_id")
    )

    chebi_agora_df = met_bigg_df.join(met_chebi_class, how="left")
    chebi_agora_df = chebi_agora_df.reset_index().set_index("agora_id")

    return chebi_agora_df

def find_annotation_level(chebi_agora_df, no_bigg_annotation):
    """
    Find the annotation level for each metabolite based on its BIGG and ChEBI annotations.

    Annotation level is of three levels:
        1. No annotation
        2. BiGG ID only
        3. BiGG ID + ChEBI class

    Parameters
    ----------
    chebi_agora_df : pandas.DataFrame
        DataFrame containing ChEBI and BiGG metabolite information.
    no_bigg_annotation : list
        List of metabolite IDs without BIGG annotations.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by AGORA ID. Contains mapping for all external metabolites in the AGORA models.
        Columns: 'annotation_level','bigg_id','chebi_id','chebi class','self defined super class'
    """

    ## Define annotation level
    # Define sets of mets based on annotation levels
    mets_no_bigg = set(no_bigg_annotation) # mets without bigg id
    mets_w_chebi = set(chebi_agora_df[chebi_agora_df["chebi id"].notna()].index) # mets with chebi class map
    mets_w_bigg = set(chebi_agora_df.index)-mets_w_chebi # mets with only bigg id but no chebi

    # Collect results in dict
    mets_annot_level_dict = {met:"0. no annotation" for met in mets_no_bigg}
    mets_annot_level_dict.update({met:"1. bigg id only" for met in mets_w_bigg})
    mets_annot_level_dict.update({met:"2. bigg id + chebi class" for met in mets_w_chebi})

    mets_annot_level_ser = pd.Series(mets_annot_level_dict)
    mets_annot_level_df = pd.DataFrame(mets_annot_level_ser.rename("annotation_level").rename_axis("agora_id"))

    chebi_agora_annot_df = pd.merge(mets_annot_level_df,chebi_agora_df, left_index=True, right_index=True, how="left")

    chebi_agora_annot_df["self defined super class"] = chebi_agora_annot_df["self defined super class"].fillna("other")
    return chebi_agora_annot_df


def main():

    # Load data
    logging.info("Loading ChEBI class data and AGORA models...")
    chebi_df = load_chebi_data()
    models = load_models_venturelli2018()


    logging.info("Create mapping for metabolites in AGORA models to ChEBI classes...")
    # Find BIGG IDs for metabolites in each model
    met_bigg_dict, no_bigg_annotation = find_mets_bigg_id(models)

    # Map ChEBI classes to AGORA metabolite IDs through BIGG ids
    chebi_agora_df = make_met_chebi_class_agora_df(met_bigg_dict, chebi_df)

    # Add column to chebi_agora_df with annotation levels for metabolites in the AGORA models
    chebi_agora_annot_df = find_annotation_level(chebi_agora_df, no_bigg_annotation)

    assert chebi_agora_annot_df.shape[0] == len(no_bigg_annotation) + len(met_bigg_dict), f"sizes: {len(chebi_agora_annot_df)}, {len(no_bigg_annotation)}, {len(met_bigg_dict)}"

    logging.info("Saving annotated DataFrame to ../data/Venturelli2018/generated/met_chebi_class_agora.tsv")
    chebi_agora_annot_df.to_csv("../data/Venturelli2018/generated/met_chebi_class_agora.tsv",sep="\t")

    return chebi_agora_annot_df
if __name__ == "__main__":
    main()

