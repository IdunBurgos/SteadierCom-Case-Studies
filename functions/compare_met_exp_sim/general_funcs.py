import pandas as pd
import operator

import logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)



def merge_sim_exp(simulated_results,signal_exp,metric="frequency",receiver_or_donor="receiver",agg_method="sum",met_column = "compound",verbose=True): #combine_exp_sim_met
    """
    Function to combine simulation results (steadycom,steadiercom) with experimental result.
    
    Input:
    simulated_results (pd.DataFrame): The raw data from steadiercom, steadycom, or smetana simulation.
    signal_exp (pd.Series): The experimental data, where with two indeces ("member","group"). This also has to be the names of the indeces.
    metric (str): Determines what metric in the simulation result to study (ex. frequency, rate, mass_rate)
    receiver_or_donor (str): 'receiver' or 'donor
    agg_method (str): how to aggregate signals in simulation results (ex. sum, mean, median)
    met_column (str): The column containing mets being compared.
    Output:
    exp_sim_combined (pd.DataFrame): Dataframe with experimental and simulated results combined.
    
    """
    
    if verbose:
        logging.info(f"merge_sim_exp removes nan values in 'experiment' column of the return exp_sim_combined object")
        
    
    valid_methods = ['sum', 'mean', 'median']
    if agg_method not in valid_methods:
        raise Exception(f"agg_method must be one of {valid_methods}")
        
        
    if receiver_or_donor!="receiver" and receiver_or_donor!="donor":
        raise Exception("receiver_or_donor must be either 'receiver' or 'donor'")
    
    exp_idx = signal_exp.index.names
    if len(exp_idx)!=2:
        raise Exception("signal_exp must have a double index (member,'met_column')")
    
    if (exp_idx[1]!=met_column):
        raise Exception(f"Concatenation will fail because 'met_column': {met_column}, does not match the index in 'signal_exp': {exp_idx[1]}. Make sure you have given the correct name for index in simulated results (default met_column='compound').")
        
    simulated_res_reduced = simulated_results.copy()
    # get mean receiving/donating quantity (or frequency) - combine from different donors/communities

    simulated_res_grouped = simulated_res_reduced[[receiver_or_donor,met_column,metric]]\
        .groupby([receiver_or_donor,met_column])\
        .agg(agg_method)\
        .reset_index()
    
    # Process results
    simulated_res_grouped = simulated_res_grouped.set_index([receiver_or_donor,met_column])[metric]
    simulated_res_grouped.index.names = ["member",met_column]
    simulated_res_grouped.name = "simulated"

    # Combine results with experiments
    exp_sim_combined = pd.concat([simulated_res_grouped,signal_exp],axis=1)
    
    # NB - removes results for compounds that are not in the dataset
    exp_sim_combined = exp_sim_combined[~exp_sim_combined.experimental.isna()].copy()
    return exp_sim_combined


