import pandas as pd
import numpy as np
import json
import csv
import operator
import seaborn as sns

import scipy as sp

from sklearn.metrics import confusion_matrix,recall_score,precision_score,classification_report,average_precision_score, roc_auc_score,precision_recall_curve,roc_curve

import matplotlib.pyplot as plt

from collections import defaultdict

import logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

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


def prepare_sim_met_venturelli(simulated_results,metid2group,receiver_or_donor="receiver",verbose=True):
    
    """
    Function to create new column in simulation results with metabolite group "met_group" and to filter out results we are not interested in.
    """
    
    if verbose:
        logging.info(f"prepare_sim_met_venturelli removes environment from {receiver_or_donor}s and all compounds that are not found in metabolomics and map to AGORA identifiers")
    
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



def prepare_all_simulations_venturelli(sim_dicts, metid2group,verbose=True):
    """
    Function to prepare all simulation results for analysis with prepare_sim_met_venturelli.
    """
    
    simulation_dict = {"receiver": {}, "donor": {}}
    for role in simulation_dict:
        for method, sim_data in sim_dicts.items():
            simulation_dict[role][method] = {
                name: prepare_sim_met_venturelli(
                    simulated_results=res,
                    metid2group=metid2group,
                    receiver_or_donor=role,
                    verbose=verbose
                )
                for name, res in sim_data.items()
            }
    return simulation_dict




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


def binarize_sim_exp(exp_sim_combined,sim_thresh=0.5,exp_condition_str="< -1"):
    
    exp_sim_combined = filter_exp_by(exp_sim_combined, exp_condition_str)
    exp_sim_combined.simulated = exp_sim_combined.simulated>sim_thresh
    
    return exp_sim_combined




def eval_sim_dict(simulation_res_dict,signal_exp,exp_condition_str="< -1",sim_thresh=0.5,sim_metric="frequency",receiver_or_donor="receiver",agg_method="sum",met_column="compound",verbose=True):
    
    exp_sim_combined_bool_dict = {}
    
    # For each result
    
    for name,simulation_res in simulation_res_dict.items():

        # Find average receiving/donating quantity
        exp_sim_combined = merge_sim_exp(simulation_res,
                                           signal_exp,
                                           metric=sim_metric,
                                           receiver_or_donor=receiver_or_donor,
                                           agg_method=agg_method,
                                           met_column=met_column,
                                           verbose=verbose)
        # Make boolean output
        ## Determine positive (True) or negative (False) results

        exp_sim_combined_bool_dict[name]=binarize_sim_exp(exp_sim_combined,sim_thresh=sim_thresh,exp_condition_str=exp_condition_str)

    return exp_sim_combined_bool_dict


def run_all_evaluations(simulation_dict, signal_exp, role_config,verbose=True):
    exp_sim_bool_dict = {}

    for role, cfg in role_config.items():
        exp_sim_bool_dict[role] = {}

        eval = lambda sim_dict, metric, thresh: eval_sim_dict(
            simulation_res_dict=sim_dict,
            signal_exp=signal_exp,
            sim_metric=metric,
            exp_condition_str=cfg["exp_cond"],
            sim_thresh=thresh,
            agg_method="sum",
            receiver_or_donor=role,
            met_column="met_group",
            verbose=verbose
        )

        exp_sim_bool_dict[role]["steadiercom_frequency"] = eval(simulation_dict[role]["steadiercom"], "frequency", cfg["freq_thresh"])
        exp_sim_bool_dict[role]["steadiercom_rate"]      = eval(simulation_dict[role]["steadiercom"], "rate",      cfg["rate_thresh"])
        exp_sim_bool_dict[role]["steadycom_rate"]        = eval(simulation_dict[role]["steadycom"],   "rate?",      cfg["rate_thresh"])
        exp_sim_bool_dict[role]["smetana"]               = eval(simulation_dict[role]["smetana"],     "mus",       cfg["freq_thresh"])

    return exp_sim_bool_dict

def merge_eval_results(exp_sim_bool_dict, role, keys):
    return {
        key: val
        for metric in keys
        for key, val in exp_sim_bool_dict[role][metric].items()
    }


def score_match(exp_sim_combined_bool,name):
    
    """
    Input:
    exp_sim_combined_bool (pd.DataFrame): One column for experimental results, and one for simulation.
    
    """
    # Score experiment vs simulation overlap
    scores = classification_report(exp_sim_combined_bool.experimental,exp_sim_combined_bool.simulated,output_dict=True)
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
    
    return scores_df_copy

def score_dict_matches(exp_sim_combined_bool_dict):
    
    """
    Input:
    exp_sim_combined_bool_dict (dict): keys show the  One column for experimental results, and one for simulation.
    
    Output:
    scores_df_all (pd.DataFrame): Dataframe with results for all elements in simulation_res_dict. Columns correspond to keys in simulation_res_dict.
    """
    scores_list = []
    
    # For each result
    for name,exp_sim_combined_bool in exp_sim_combined_bool_dict.items():
            
        # Score experiment vs simulation overlap
        scores = score_match(exp_sim_combined_bool,name)

        scores_list.append(scores)

    scores_df_all = pd.concat(scores_list,axis=1)
    scores_df_all.sort_index(inplace=True)
    return scores_df_all

def plot_scatter_plots(simulation_dict,signal_exp,simulation_method="steadiercom",receiver_or_donor="receiver",agg_method="sum",metric="frequency",verbose=True):
    
    if verbose:
        logging.info(f"plot_scatter_plots uses fillna(0) to for missing entries in simulations that exist in experimental results.")
        
    # Process simulated result
    collected_result = defaultdict(dict)
    for receiver_or_donor in ["donor","receiver"]:
        for name, simulated_res in simulation_dict[receiver_or_donor][simulation_method].items():

            exp_sim_combined = merge_sim_exp(
                simulated_res,
                signal_exp,
                metric=metric,
                agg_method=agg_method,
                receiver_or_donor=receiver_or_donor,
                met_column="met_group",
                verbose=verbose
            ).fillna(0)

            collected_result[name][receiver_or_donor] = exp_sim_combined


    # Plot results    
    fig, axs = plt.subplots(2,3,figsize=(15,10))

    i = 0
    for name, subdict in collected_result.items():
        j = 0
        
        for receiver_or_donor,color in {"receiver":"tab:blue","donor":"tab:orange"}.items():
            ax = axs[j][i]

            merged = subdict[receiver_or_donor]
            sns.scatterplot(ax=ax,x="experimental", y="simulated", data=merged,alpha=0.5,color=color)
            
            if receiver_or_donor=="receiver":
                ax.set_title(f"{name}")

            ax.set_ylabel(f"Simulated - {agg_method} of {metric}")
            ax.set_xlabel(f"Experimental - log2 conc. change")
            if metric!="frequency":
                ax.set_yscale('log')

            ax.axvspan(-1, 1, alpha=0.2, color='grey')

            r, p = sp.stats.spearmanr(subdict[receiver_or_donor])
            plt.text(.1, .75, "Spearman's \n r={:.2f}, p-value={:.4f}".format(r,p), transform=ax.transAxes)
            j+=1

        i +=1
        
    fig.text(0.01, 0.70, "Receivers", va='center', ha='left', fontsize=12)
    
    fig.text(0.01, 0.30, "Donors", va='center', ha='left', fontsize=12)
    return fig,axs

def plot_diagnostic_curve(simulation_dict,signal_exp,PR_or_ROC="PR", receiver_or_donor="receiver",simulation_method="steadiercom",metric="frequency",met_column="met_group",agg_method="sum",verbose=True):
    fig,axs =  plt.subplots(1,3,figsize=(20, 6))
    
    if verbose:
        logging.info(f"plot_diagnostic_curve uses fillna(0) to for missing entries in simulations that exist in experimental results.")
        

    i=0
    for name,sim_df in simulation_dict[receiver_or_donor][simulation_method].items():


        exp_sim_combined = merge_sim_exp(sim_df,
                                         signal_exp=signal_exp,metric=metric,
                                         receiver_or_donor=receiver_or_donor,
                                         agg_method=agg_method,
                                         met_column = met_column,
                                         verbose=verbose).fillna(0)
        if receiver_or_donor=="receiver":
            exp_sim_combined["experimental"]=(exp_sim_combined["experimental"]<-1).map(int)
        elif receiver_or_donor=="donor":
            exp_sim_combined["experimental"]=(exp_sim_combined["experimental"]>1).map(int)
            
        y_true=exp_sim_combined["experimental"].values
        y_scores=exp_sim_combined["simulated"]

        
        
        if PR_or_ROC=="PR":
            baseline = sum(y_true)/len(y_true)
            auc = average_precision_score(y_true, y_scores)
            precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
            # Plot the curve
            ax=axs[i]
            ax.plot(recall, precision, marker='.', label='PR curve')
            ax.set_xlabel('Recall')
            ax.set_ylabel('Precision')
            plt.text(.1, .75, f"Area Under Precision-Recall curve (PR-AUC): {auc:.3f}",transform=ax.transAxes)
            ax.axhline(y=baseline, color='black', linestyle='--',label=f"baseline={baseline:.3f}")
        elif PR_or_ROC=="ROC":
            auc = roc_auc_score(y_true, y_scores)
            fpr, tpr, thresholds = roc_curve(y_true, y_scores)
            # Plot the curve
            ax=axs[i]
            ax.plot(fpr, tpr, marker='.', label='ROC curve')
            ax.set_xlabel('FPR')
            ax.set_ylabel('TPR')

            plt.text(.1, .75, f"Area Under the Receiver Operating Characteristic curve \n (ROC-AUC): {auc:.3f}",transform=ax.transAxes)
        
        ax.set_title(f"{name}")
        ax.legend()

        i+=1
    fig.suptitle(f"{receiver_or_donor.upper()}")
    return fig,axs


def plot_confusion_matrix_venturelli(exp_sim_combined_bool_dict,suptitle="No title"):

    i = 0
    fig,axs=plt.subplots(2,3,figsize=(15,10))
    axs_flat = axs.flatten()

    for name,exp_sim_combined_bool in exp_sim_combined_bool_dict.items():
            
        cf_matrix = confusion_matrix(exp_sim_combined_bool.experimental, exp_sim_combined_bool.simulated)

        ax = sns.heatmap(cf_matrix, 
                         annot = cf_matrix,
                         fmt="", cmap='winter_r',vmin=0,vmax=200, ax=axs_flat[i])
        ax.set_ylabel("Experimental")
        ax.set_xlabel("Simulated")
        ax.set_yticklabels(["False","True"])
        ax.set_xticklabels(["False","True"])
        ax.set_title(name)
        i +=1
    fig.suptitle(suptitle)
    return fig,axs



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




