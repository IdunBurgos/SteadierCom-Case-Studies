import pandas as pd
import numpy as np
from os import listdir

import operator
import seaborn as sns

import scipy as sp

from sklearn.metrics import confusion_matrix,recall_score,precision_score,classification_report,average_precision_score, roc_auc_score,precision_recall_curve,roc_curve,PrecisionRecallDisplay,RocCurveDisplay

import matplotlib.pyplot as plt

from collections import defaultdict

import logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


from compare_met_exp_sim.general_funcs import merge_sim_exp

def _eval_exp(exp_ser,receiver_or_donor="receiver",exp_thresh=1):
    """
    Function to classify experimental interaction as True or False.
    """
    
    if receiver_or_donor!="receiver" and receiver_or_donor!="donor":
        raise Exception("receiver_or_donor must be either 'receiver' or 'donor'")
    
    exp_ser = exp_ser.copy()
    
    if receiver_or_donor=="receiver":
        exp_ser = (exp_ser<-exp_thresh).map(int)
    elif receiver_or_donor=="donor":
        exp_ser = (exp_ser>exp_thresh).map(int)
        
    return exp_ser

def identify_support(signal_exp, exp_thresh_list=[0.25,0.5,1]):
    """
    Function to identify number of interactions that are supported by experimental data at different thresholds.
    """
    
    support_dict = {}

    for receiver_or_donor in ["receiver", "donor"]:
        for exp_thresh in exp_thresh_list:

            signal_exp_bool = signal_exp.copy()
            if receiver_or_donor=="receiver":
                signal_exp_bool = (signal_exp_bool<-exp_thresh).map(int)
            elif receiver_or_donor=="donor":
                signal_exp_bool = (signal_exp_bool>exp_thresh).map(int)
            

            positive_support = signal_exp_bool.sum()
            negative_support = (signal_exp_bool==0).sum()
            total = len(signal_exp_bool)

            support_dict[(receiver_or_donor, exp_thresh)] = {"positive_support":positive_support,
                                                            "negative_support":negative_support,
                                                            "total":total,
                                                            "no-skill line": positive_support/total if total > 0 else 0,
                                                            }

    return pd.DataFrame.from_dict(support_dict, orient="index")



def _eval_sim(sim_ser,sim_thresh=0.5):
    """
    Function to classify simulated interaction as True or False.
    """
    
    sim_ser = sim_ser.copy()
    sim_ser = (sim_ser>sim_thresh).map(int)
    
    return sim_ser

def _eval_sim_exp(exp_sim_combined,sim_thresh=0.5,exp_thresh=1,receiver_or_donor="receiver"):

    """
    Function to classify both experimental and simulated interactions as True or False.
    """

    exp_sim_combined.experimental= _eval_exp(exp_sim_combined=exp_sim_combined.experimental,receiver_or_donor=receiver_or_donor,exp_thresh=exp_thresh)
    exp_sim_combined.simulated= _eval_sim(exp_sim_combined=exp_sim_combined.simulated,sim_thresh=sim_thresh)
    
    return exp_sim_combined

def run_all_evaluations(simulation_dict, signal_exp, role_config=None,verbose=True):
    """
    Function to merge simulated and experimental results for all simulation results in simulation_dict_method, 
    and classify both experimental and simulated interactions as True or False. simulation_dict_method is a dictionary 
    with keys as the mode for the simulation (with fixed abundance, fixed growth...), and values as pd.DataFrame 
    with simulation results.
    """

    def _eval_sim_exp_dict(simulation_dict_method,signal_exp,exp_thresh=1,sim_thresh=0.5,sim_metric="frequency",receiver_or_donor="receiver",agg_method="sum",met_column="compound",verbose=True):
    
        """
        Function to merge simulated and experimental results for all simulation results in simulation_dict_method, 
        and classify both experimental and simulated interactions as True or False. simulation_dict_method is a dictionary 
        with keys as the mode for the simulation (with fixed abundance, fixed growth...), and values as pd.DataFrame 
        with simulation results.
        """

        exp_sim_combined_method_bool_dict = {}
        
        # For each result

        i=0
        for name,simulation_res in simulation_dict_method.items():
            # Find average receiving/donating quantity
            exp_sim_combined = merge_sim_exp(simulation_res,
                                            signal_exp,
                                            metric=sim_metric,
                                            receiver_or_donor=receiver_or_donor,
                                            agg_method=agg_method,
                                            met_column=met_column,
                                            verbose=(i==0 and verbose),
                                            )
            # Make boolean output
            ## Determine positive (True) or negative (False) results

            exp_sim_combined_method_bool_dict[name]=_eval_sim_exp(exp_sim_combined,
                                                            sim_thresh=sim_thresh,
                                                            exp_thresh=exp_thresh,
                                                            receiver_or_donor=receiver_or_donor)
            i+=1

        return exp_sim_combined_method_bool_dict
    


    if not role_config:
        role_config = {
            "receiver": {"exp_thresh":1,"freq_thresh":0.5,"rate_thresh":0.1},
            "donor":    {"exp_thresh":1,"freq_thresh":0.5,"rate_thresh":0.1}
        }

    exp_sim_bool_dict = {}
    
    i=0
    for role, cfg in role_config.items():
        exp_sim_bool_dict[role] = {}

        eval = lambda sim_dict_method, metric, thresh: _eval_sim_exp_dict(
            simulation_dict_method=sim_dict_method,
            signal_exp=signal_exp,
            sim_metric=metric,
            exp_thresh=cfg["exp_thresh"],
            sim_thresh=thresh,
            agg_method="sum",
            receiver_or_donor=role,
            met_column="met_group",
            verbose=(i==0 and verbose)
        )

        exp_sim_bool_dict[role]["steadiercom_frequency"] = eval(simulation_dict[role]["steadiercom"], "frequency", cfg["freq_thresh"])
        exp_sim_bool_dict[role]["steadiercom_rate"]      = eval(simulation_dict[role]["steadiercom"], "rate",      cfg["rate_thresh"])
        exp_sim_bool_dict[role]["steadycom_rate"]        = eval(simulation_dict[role]["steadycom"],   "rate?",      cfg["rate_thresh"])
        exp_sim_bool_dict[role]["smetana"]               = eval(simulation_dict[role]["smetana"],     "mus",       cfg["freq_thresh"])
        i+=1

    return exp_sim_bool_dict


def _create_classification_report(exp_sim_combined_bool,name):
    
    """
    Input:
    exp_sim_combined_bool (pd.DataFrame): One column for experimental results, and one for simulation.
    
    """
    # Score experiment vs simulation overlap
    report = classification_report(exp_sim_combined_bool.experimental,exp_sim_combined_bool.simulated,output_dict=True)
    report_df = pd.DataFrame(report)

    # Extract True class (positive class)
    report_df_copy1 = report_df['1'].copy()
#    report_df_copy1 = report_df_copy1.drop("support")
    report_df_copy1.index = ["precision","sensitivity","f1-score","support positive"]

    # Extract False class (negative class)
    report_df_copy2 = report_df['0'].copy()
    report_df_copy2 = report_df_copy2.drop(["precision","f1-score"])
    report_df_copy2.index = ["specificity","support negative"]

    # Combine results
    report_df_copy = pd.concat([report_df_copy1,report_df_copy2])
    report_df_copy["accuracy"] = report["accuracy"]
    report_df_copy.name = name

    return report_df_copy.sort_index()

def create_classification_reports(exp_sim_combined_bool_dict):
    """
    Input:
    exp_sim_combined_bool_dict (dict): Dictionary with keys as names, and values as pd.DataFrame with one column for experimental results, and one for simulation.
    
    Output:
    scores_df (pd.DataFrame): Dataframe with classification report scores for each name in exp_sim_combined_bool_dict.
    
    """
    
    reports_dict = {}
    
    for receiver_or_donor, exp_sim_combined_bool_dict_role in exp_sim_combined_bool_dict.items():
            for method, exp_sim_combined_bool_dict_role_method in exp_sim_combined_bool_dict_role.items():
                for name,exp_sim_combined_bool in exp_sim_combined_bool_dict_role_method.items():

                    report_df_copy = _create_classification_report(exp_sim_combined_bool,
                                                                name,
                                                                )
                    reports_dict[(receiver_or_donor, method, name)] = report_df_copy

    reports_df = pd.concat(reports_dict)

    return reports_df


def _permutation_test(y_true,y_scores,N):

    observed_auc = average_precision_score(y_true,y_scores)

    auc_values = []

    for _ in range(N):
        y_true_i = np.random.permutation(y_true)
        auc_i = average_precision_score(y_true_i, y_scores)
        auc_values.append(auc_i)

    p_value = np.mean(auc_values >= observed_auc)

    return p_value

def pr_auc_pvalue(simulation_dict,signal_exp,exp_thresh=1,receiver_or_donor="receiver",name="SteadierCom - fixed abundance",simulation_method="steadiercom",metric="frequency",met_column="met_group",agg_method="sum",N=1000,verbose=True):


    sim_df = simulation_dict[receiver_or_donor][simulation_method][name]
    exp_sim_combined = merge_sim_exp(sim_df,
                                    signal_exp=signal_exp,
                                    metric=metric,
                                    receiver_or_donor=receiver_or_donor,
                                    agg_method=agg_method,
                                    met_column = met_column,
                                    verbose=verbose).fillna(0)
    
    exp_sim_combined.experimental=_eval_exp(exp_sim_combined.experimental,
                                            receiver_or_donor=receiver_or_donor,
                                            exp_thresh=exp_thresh)

    y_true=exp_sim_combined["experimental"].values
    y_scores=exp_sim_combined["simulated"].values

    p_value = _permutation_test(y_true,y_scores,N=N)

    return {"p_value":p_value}


def pr_auc_pvalue_simulation_method(simulation_dict,signal_exp,exp_thresh=1,simulation_method="steadiercom",metric="frequency",met_column="met_group",agg_method="sum",N=1000,verbose=True):
    
    # Process simulated result
    pvalues_dict = defaultdict(dict)
    
    for receiver_or_donor in ["receiver", "donor"]:

        for i,name in enumerate(simulation_dict[receiver_or_donor][simulation_method].keys()):
            print(f"{simulation_method}, {name}, {receiver_or_donor}, {metric}")
            verbose_spec = (i==0 and verbose)
            pvalue = pr_auc_pvalue(simulation_dict,
                                          signal_exp,
                                          exp_thresh=exp_thresh,
                                          name=name,
                                          simulation_method=simulation_method,
                                          receiver_or_donor=receiver_or_donor,
                                          metric=metric,
                                          agg_method=agg_method,
                                          met_column = met_column,
                                          N=N,
                                          verbose=verbose_spec)

            pvalues_dict[(receiver_or_donor,simulation_method,name,agg_method,metric,exp_thresh)] = pvalue

    pvalues_df = pd.DataFrame(pvalues_dict).T 
    pvalues_df.index.names = ["role","simulation_method","condition","agg_method","metric","exp_thresh"]
    pvalues_df = pvalues_df.sort_index()
    return pvalues_df


def plot_diagnostic_curve(simulation_dict,signal_exp,PR_or_ROC="PR", exp_thresh=1,receiver_or_donor="receiver",name="SteadierCom - fixed abundance",simulation_method="steadiercom",metric="frequency",met_column="met_group",agg_method="sum",ax=None,verbose=True):
    if verbose:
            logging.info(f"plot_diagnostic_curve uses fillna(0) to for missing entries in simulations that exist in experimental results.")
    
    if not ax:
        fig, ax = plt.subplots(figsize=(6, 6))
    else:
        fig = ax.get_figure()

    sim_df = simulation_dict[receiver_or_donor][simulation_method][name]

    exp_sim_combined = merge_sim_exp(sim_df,
                                    signal_exp=signal_exp,metric=metric,
                                    receiver_or_donor=receiver_or_donor,
                                    agg_method=agg_method,
                                    met_column = met_column,
                                    verbose=verbose).fillna(0)
    
    exp_sim_combined.experimental=_eval_exp(exp_sim_combined.experimental,receiver_or_donor=receiver_or_donor,exp_thresh=exp_thresh)

    y_true=exp_sim_combined["experimental"].values
    y_scores=exp_sim_combined["simulated"]

    if PR_or_ROC=="PR":
        baseline = sum(y_true)/len(y_true)
        auc = average_precision_score(y_true, y_scores)
        precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
        # Plot the curve
        
        PrecisionRecallDisplay(precision,recall).plot(ax, marker='.', label='PR curve')
        ax.set_xlabel('Recall TP/(FN+TP)')
        ax.set_ylabel('Precision TP/(FP+TP)')
        plt.text(.1, .75, f"Area Under Precision-Recall curve (PR-AUC): {auc:.3f}",transform=ax.transAxes)
        ax.axhline(y=baseline, color='black', linestyle='--',label=f"baseline={baseline:.3f}")
    
    elif PR_or_ROC=="ROC":
        auc = roc_auc_score(y_true, y_scores)
        fpr, tpr, thresholds = roc_curve(y_true, y_scores)
        # Plot the curve
        
        RocCurveDisplay(fpr=fpr, tpr=tpr).plot(ax=ax, marker='.', label='ROC curve')
        ax.set_xlabel('FPR (FP/(FP+TN))')
        ax.set_ylabel('TPR (TP/(FN+TP))')

        plt.text(.1, .75, f"Area Under the Receiver Operating Characteristic curve \n (ROC-AUC): {auc:.3f}",transform=ax.transAxes)
    
    else:
        raise Exception("PR_or_ROC must be either 'PR' or 'ROC'")
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.set_title(f"{name}")
    ax.legend()

    return fig,ax


def plot_diagnostic_curves(simulation_dict,signal_exp,PR_or_ROC="PR",exp_thresh=1,simulation_method="steadiercom",metric="frequency",met_column="met_group",agg_method="sum",verbose=True):
    fig,axs =  plt.subplots(2, 3, figsize=(18, 10))
    
    if verbose:
        logging.info(f"plot_diagnostic_curve uses fillna(0) to for missing entries in simulations that exist in experimental results.")
        

    i=0

    for receiver_or_donor in ["receiver","donor"]:
        j=0
        for name in simulation_dict[receiver_or_donor][simulation_method].keys():

            _ , _ = plot_diagnostic_curve(simulation_dict,
                                signal_exp,
                                exp_thresh=exp_thresh,
                                PR_or_ROC=PR_or_ROC,
                                receiver_or_donor=receiver_or_donor,
                                name=name,
                                simulation_method=simulation_method,
                                metric=metric,
                                met_column=met_column,
                                agg_method=agg_method,
                                ax=axs[i][j],
                                verbose=(i==0 and j==0 and verbose))

            j+=1
    
        i+=1
        

    fig.text(0.05, 0.75, "Receiver",fontsize=14)
    fig.text(0.05, 0.25, "Donor",fontsize=14)
    fig.tight_layout()

    #fig.suptitle(f"{receiver_or_donor.upper()}")
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




