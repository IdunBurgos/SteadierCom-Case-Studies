import pandas as pd
import numpy as np
import json
import csv
import operator

from sklearn.metrics import classification_report
import matplotlib.pyplot as plt

def prepare_exp_met_venturelli(dataset4,metname_groups,metname2group):
    """
    Function for processing metabolomics data which gives output 'signal_exp'
    """
    
    dataset4.index = dataset4.index.map(lambda met: met.replace("'", "").replace(" ", "").lower())

    
    ## Assign metabolites to groups and join signal
    
    dataset4_stack = dataset4.stack().reset_index()
    dataset4_stack.columns = ["metname","member","change"]

    # Add group name to df
    dataset4_stack = dataset4_stack[dataset4_stack.metname.isin(metname2group.keys())].copy()
    dataset4_stack["met_group"] = dataset4_stack.metname.map(lambda x: metname2group[x])

    # Add met weight to df
    dataset4_stack["weight"]=dataset4_stack.apply(lambda x: metname_groups[x["met_group"]][x["metname"]],axis=1)
    dataset4_stack["change*weight"]=dataset4_stack.change*dataset4_stack.weight

    # Find the sum of the change and weights within a group
    combined_signal = dataset4_stack[["member","met_group","change*weight"]].groupby(["member","met_group"]).sum()["change*weight"]
    combined_weights = dataset4_stack[["member","met_group","weight"]].groupby(["member","met_group"]).sum()["weight"]

    # Find the geometric average of the signal
    signal_exp = (combined_signal/combined_weights).unstack()

    # Adjust name of member
    signal_exp.index = signal_exp.index.map(lambda x: x.replace("'","").replace(" ",""))

    # Prepare data for joining with steadiercom data
    signal_exp = signal_exp.stack()
    signal_exp.name = "experimental"
    return signal_exp


def prepare_sim_met_venturelli(simulated_results,metid2group,receiver_or_donor="receiver"):
    
    """
    Function to create new column in simulation results with metabolite group "met_group".
    """
    # Map met id to met group
    if receiver_or_donor=="receiver":
        simulated_results = simulated_results[simulated_results.receiver!="environment"].copy()
        simulated_results = simulated_results[~simulated_results.receiver.isna()].copy()
        simulated_results.receiver = simulated_results.receiver.map(lambda x: x.split("_")[0][0].upper()+x.split("_")[1][0].upper())
        
    elif receiver_or_donor=="donor":
        
        simulated_results = simulated_results[simulated_results.donor!="environment"].copy()
        simulated_results = simulated_results[~simulated_results.donor.isna()].copy()
        simulated_results.donor = simulated_results.donor.map(lambda x: x.split("_")[0][0].upper()+x.split("_")[1][0].upper())
    else:
        raise Exception("receiver_or_donor must be either 'receiver' or 'donor'")

    
    simulated_res_reduced = simulated_results[simulated_results.compound.isin(metid2group.keys())].copy()
    simulated_res_reduced["met_group"] = simulated_res_reduced.compound.map(lambda x: metid2group[x])
    
    return simulated_res_reduced


def filter_exp_by(exp_sim_combined, condition_str):

    op_str, thresh_str = condition_str.strip().split()
    thresh = float(thresh_str)

    # Map string to operator function
    ops = {
        '<': operator.lt,
        '>': operator.gt,
        '<=': operator.le,
        '>=': operator.ge,
        '==': operator.eq,
        '!=': operator.ne
    }

    compare_func = ops[op_str]
    exp_sim_combined.experimental = compare_func(exp_sim_combined.experimental, thresh)

    return exp_sim_combined


def combine_exp_sim_met(simulated_results,signal_exp,metric="frequency",receiver_or_donor="receiver",agg_method="sum",met_column = "compound"):
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



def scoring_met(simulation_res_dict,signal_exp,exp_condition_str="< -1",sim_thresh=0.5,sim_metric="frequency",receiver_or_donor="receiver",agg_method="sum",met_column="compound"):
    
    """
    Input:
    simulation_res_dict (dict(str,pd.DataFrame)): A dict with raw data from steadiercom, steadycom, or smetana simulation.
    signal_exp (pd.Series): The experimental data, where with two indeces ("member",met_column). This also has to be the names of the indeces.
    exp_thresh (float): The threshold for experimental results. Anything above/below will be True/False.
    sim_thresh (float): The threshold for simulation results. Anything above will be True.
    sim_metric (str): Determines what metric in the simulation result to study (ex. frequency, rate, mass_rate)
    
    receiver_or_donor (str): 'receiver' or 'donor
    agg_method (str): how to aggregate signals in simulation results (ex. sum, mean, median)
    
    Output:
    scores_df_all (pd.DataFrame): Dataframe with results for all elements in simulation_res_dict. Columns correspond to keys in simulation_res_dict.
    """
    
    scores_list = []
    
    # For each result
    for name,simulation_res in simulation_res_dict.items():
        
        # Find average receiving/donating quantity
        exp_sim_combined = combine_exp_sim_met(simulation_res,
                                           signal_exp,
                                           metric=sim_metric,
                                           receiver_or_donor=receiver_or_donor,
                                           agg_method=agg_method,
                                           met_column=met_column)
        # Make boolean output
        ## Determine positive (True) or negative (False) results
        exp_sim_combined = filter_exp_by(exp_sim_combined, exp_condition_str)
        exp_sim_combined.simulated = exp_sim_combined.simulated>sim_thresh
            
        # Score experiment vs simulation overlap
        scores = classification_report(exp_sim_combined.experimental,exp_sim_combined.simulated,output_dict=True)
        scores_df = pd.DataFrame(scores)
        
        # Extract True class (positive class)
        scores_df_copy1 = scores_df["True"].copy()
        scores_df_copy1 = scores_df_copy1.drop("support")
        scores_df_copy1.index = ["precision","sensitivity","f1-score"]
        
        # Extract False class (negative class)
        scores_df_copy2 = scores_df["False"].copy()
        scores_df_copy2 = scores_df_copy2.drop(["support","precision","f1-score"])
        scores_df_copy2.index = ["specificity"]
        
        # Combine results
        scores_df_copy = pd.concat([scores_df_copy1,scores_df_copy2])
        scores_df_copy["accuracy"] = scores["accuracy"]
        scores_df_copy.name = name

        scores_list.append(scores_df_copy)

    scores_df_all = pd.concat(scores_list,axis=1)
    scores_df_all.sort_index(inplace=True)
    return scores_df_all



def plot_radar_chart(scores_df,ax=False):
    """
    This functions plots the results from scores_df as a radar chart. 
    
    """

    fig, ax = plt.subplots(figsize=(5, 5), subplot_kw={'projection': 'polar'})

    # Set up radar chart
    metrics = scores_df.index.values
    N = len(metrics)
    theta = np.linspace(0, 2 * np.pi, N, endpoint=False)
    theta = np.concatenate([theta, [theta[0]]])
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rlabel_position(90)
    ax.spines['polar'].set_zorder(1)
    ax.spines['polar'].set_color('lightgrey')
    ax.set_ylim(0, 0.5)


    for idx, (i, row) in enumerate(scores_df.T.iterrows()):
        values = row[metrics].values.flatten().tolist()
        values = values + [values[0]] # circle back
        ax.plot(theta, values, linewidth=1.75, linestyle='solid', label=i, marker='o', markersize=10)


    # Use only the non-duplicated angles for labeling
    ax.set_xticks(theta[:-1])
    ax.set_xticklabels(metrics, color='black', size=12)

    # Optional: customize radial ticks
    ax.set_yticks([0.0, 0.25, 0.5,1.0])
    plt.legend(loc='upper right', bbox_to_anchor=(1.6, 1.1))
    
    return fig,ax