from itertools import product
import pandas as pd
import json
import numpy as np

def merge_smetana_steadiercom(smetana_results,steadiercom_results):
    """
    Function for creating donor-compound-receiver triplets for SMETANA and SteadierCom results. 
    
    Also assigns super_class

    Parameters
    -----------
    smetana_results (pd.DataFrame): Dataframe with SMETANA results
    steadiercom_results (pd.DataFrame): Dataframe with SteadierCom results

    Returns
    ----------
    merged_sim (pd.DataFrame): Dataframe with merged data from simulations.
    """

    # Remove environment
    steadiercom_subset = steadiercom_results.query("receiver != 'environment' & donor != 'environment'")[["receiver","compound","donor","community","frequency","mass_rate*frequency"]]

    merged_sim = pd.merge(left=smetana_results, right=steadiercom_subset, on=["receiver", "compound","donor","community"], how="outer")
    merged_sim.loc[:,["donor","compound","receiver","community","frequency","mass_rate*frequency","smetana"]]
    return merged_sim


def _make_metabolomics_triplets(metabolomics, exp_thresh=1):

    """
    Function for creating donor-compound-receiver triplets from metabolomics data.

    NB! This function shows all inferred interactions (even when the compound name cannot be translated to the AGORA database).
    """

    # Make df with all possible combination
    receiver_list = metabolomics.index.get_level_values("member").unique()
    donor_list = metabolomics.index.get_level_values("member").unique()
    compound_list = metabolomics.index.get_level_values("compound").unique()
    key_cols = ["receiver", "donor", "compound_name"]

    universe = pd.DataFrame(
        product(receiver_list, donor_list, compound_list),
        columns=key_cols
        )
    universe = universe.query("receiver != donor") # Remove self interactions

    ## Create triplets for metabolomics data
    # Filter by log change for finding consumers and producers of compounds 
    consumer = (
        metabolomics[metabolomics < -exp_thresh]
                .rename("consumer")
                .reset_index()[["compound","member"]]
                )
    producer = (
        metabolomics[metabolomics > exp_thresh]
                .rename("producer")
                .reset_index()[["compound","member"]]
                )


    # Find triplets from merging compound producer and consumer...
    exp_inter_df = (
        pd.merge(left=producer,right=consumer,on="compound",how="inner")
        .rename(columns={"member_x":"donor","member_y":"receiver"})
        )
    
    cols = ["compound","receiver","donor"]
    exp_inter_df = exp_inter_df[cols].rename(columns={"compound":"compound_name"})

    # combine with universal data
    exp_class_df = exp_inter_df.rename(columns={"compound":"compound_name"}).merge(universe, on=key_cols, how="outer", indicator=True)

    exp_class_df["exp_interact"] = exp_class_df["_merge"].map({
        "both": True,
        "right_only": False
    })
    exp_class_df = exp_class_df.drop(columns="_merge")

    return exp_class_df

def add_metabolomics(merged_sim,metabolomics,metid2metname,sim_data_columns,exp_thresh=1):
    """
    Function for labeling interactions predicted simulated data (steadiercom and smetana) based on inferred interactions from metabolomics data. 

    NB! This function keeps all simulated interactions (even when the compound has not been meassured in the metabolomics). 

    Parameters
    ----------
    merged_sim (pd.DataFrame): Dataframe with merged data from simulations.
    metabolomics (pd.DataFrame): Dataframe with metabolomics data
    metids_mapped (list): met IDs that could be mapped to metabolomics dataset.
    sim_data_columns (list): List of columns with simulated data to be filled in with 0 when no simulated data is available.

    Returns
    ----------
    merged_sim_exp (pd.DataFrame): Dataframe with merged data from simulations and metabolomics.
    """

    exp_class_df = _make_metabolomics_triplets(metabolomics, exp_thresh=exp_thresh)

    ### NB! filter data###
    # Experimental data: only include compounds that could be mapped to simulated data. 
    exp_class_df = exp_class_df[exp_class_df.compound_name.isin(metid2metname.values())].copy() # Only keep compounds that could be mapped to AGORA database.

    # Simulated data: Compound ID is mapped to compound name
    merged_sim = merged_sim.assign(compound_name=lambda x: x["compound"].map(metid2metname))

    # Merge exp triplets with sim triplets   
    merged_sim_exp = merged_sim.merge(exp_class_df, on=["compound_name","receiver","donor"], how="outer", indicator=True)

    ### Clean up columns
    sim_missing = merged_sim_exp["_merge"].eq("right_only") # Identify missing simulations

    merged_sim_exp["compound"] = ( # fill in missing compound ids
        merged_sim_exp["compound"]
        .fillna(merged_sim_exp["compound_name"])
    )

    merged_sim_exp.loc[sim_missing, sim_data_columns] = 0 # fill in missing simulated data

    merged_sim_exp["experimental_support"] = np.select( # create experimental support column
        [
            merged_sim_exp["_merge"].eq("left_only"),
            merged_sim_exp["exp_interact"].eq(True),
            merged_sim_exp["exp_interact"].eq(False),
        ],
        [
            "met not measured",
            "supported",
            "not supported",
        ],
        default=pd.NA,
    )
    
    return merged_sim_exp.drop(["compound_name","exp_interact","_merge"], axis=1) # drop unnecessary columns


def add_quantiles(merged_sim):

    # Prepare for calculating statistics: fill in missing values in simulated data (ensure all communities are represented with the same interactions)
    sim_wide = merged_sim.pivot_table(
        index=["donor", "compound", "receiver"],
        columns="community",
        values=["frequency", "smetana"],
        fill_value=0,
    )

    # Calculate statistics around median value
    df_stats = (
        pd.DataFrame({
            "steadiercom_median": sim_wide["frequency"].median(axis=1),
            "steadiercom_q25":    sim_wide["frequency"].quantile(0.25, axis=1),
            "steadiercom_q75":    sim_wide["frequency"].quantile(0.75, axis=1),

            "smetana_median":      sim_wide["smetana"].median(axis=1),
            "smetana_q25":         sim_wide["smetana"].quantile(0.25, axis=1),
            "smetana_q75":         sim_wide["smetana"].quantile(0.75, axis=1),
        })
        .reset_index()
    )

    return df_stats