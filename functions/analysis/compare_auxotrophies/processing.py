import pandas as pd

def merge_smetana_steadiercom_triplets(smetana_results,steadiercom_results,super_class_map):
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
    steadiercom_subset = steadiercom_results.query("receiver != 'environment' & donor != 'environment'")[["receiver","compound","donor","frequency"]]

    merged_sim = pd.merge(left=smetana_results, right=steadiercom_subset, on=["receiver", "compound","donor"], how="outer")
    
    merged_sim["super_class"] = merged_sim.compound.map(super_class_map)
    merged_sim.super_class = merged_sim.super_class.fillna("other")
    return merged_sim

def add_auxotrophies_to_df(merged_sim, expected_aux=[("mutant_panZ","M_pnto__R_e"),
 ("mutant_panD","M_pnto__R_e"),
("mutant_bioC" ,"M_btn_e"),
("mutant_bioF" ,"M_btn_e"),
("mutant_bioH" ,"M_btn_e"),
("mutant_cysW","M_cys__L_e"),
("mutant_icd","M_glu__L_e"),
]):
    
    mask = pd.Series(zip(merged_sim.receiver, merged_sim.compound))
    merged_sim.loc[:,"expected_auxotrophy"] = mask.isin(expected_aux)
    return merged_sim