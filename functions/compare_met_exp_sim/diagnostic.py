import pandas as pd
import numpy as np

from sklearn.metrics import confusion_matrix,recall_score,precision_score,classification_report,average_precision_score, roc_auc_score,precision_recall_curve,roc_curve,PrecisionRecallDisplay,RocCurveDisplay

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

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



def _pr_permutation_test(simulation_dict,signal_exp,exp_thresh=1,receiver_or_donor="receiver",name="SteadierCom - fixed abundance",simulation_method="steadiercom",metric="frequency",met_column="met_group",agg_method="sum",N=1000,verbose=True):

    def _permutation_test(y_true,y_scores,N):

        observed_auc = average_precision_score(y_true,y_scores)
        auc_values = []

        for _ in range(N):
            y_true_i = np.random.permutation(y_true)
            auc_i = average_precision_score(y_true_i, y_scores)
            auc_values.append(auc_i)

        p_value = np.mean(auc_values >= observed_auc)

        return p_value

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


def _pr_auc_pvalue_simulation_method(simulation_dict,signal_exp,exp_thresh=1,simulation_method="steadiercom",metric="frequency",met_column="met_group",agg_method="sum",N=1000,verbose=True):
    
    # Process simulated result
    pvalues_dict = defaultdict(dict)
    
    for receiver_or_donor in ["receiver", "donor"]:

        for i,name in enumerate(simulation_dict[receiver_or_donor][simulation_method].keys()):
            verbose_spec = (i==0 and verbose)
            pvalue = _pr_permutation_test(simulation_dict,
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


def find_pvalues_of_PR(simulation_dict,signal_exp):

    N=10000
    results_pr_auc = []
    sim_method_metric = [("steadycom","rate?"),("steadiercom","frequency*mass_rate"),("steadiercom","frequency"),("smetana","mus")]

    verbose=True
    for exp_thresh in [0.5,1]:

        agg_method ="sum"

        for simulation_method,metric in sim_method_metric:
            results_pr_auc.append(_pr_auc_pvalue_simulation_method(simulation_dict,
                                                                            signal_exp,
                                                                            exp_thresh=exp_thresh,
                                                                            simulation_method=simulation_method,
                                                                            metric=metric,
                                                                            agg_method=agg_method,
                                                                            N=N,
                                                                            verbose=verbose)
            )
            verbose=False

    pr_auc_pvalue_df = pd.concat(results_pr_auc)

    return pr_auc_pvalue_df


def plot_diagnostic_curve(simulated_res,signal_exp,PR_or_ROC="PR", exp_thresh=1,receiver_or_donor="receiver",metric="frequency",met_column="met_group",agg_method="sum",ax=None,title="",plot_cbar=True, verbose=True):
    if verbose:
            logging.info(f"plot_diagnostic_curve uses fillna(0) to for missing entries in simulations that exist in experimental results.")
    
    if not ax:
        fig, ax = plt.subplots(figsize=(6, 6))
    else:
        fig = ax.get_figure()

    exp_sim_combined = merge_sim_exp(simulated_res,
                                    signal_exp=signal_exp,metric=metric,
                                    receiver_or_donor=receiver_or_donor,
                                    agg_method=agg_method,
                                    met_column=met_column,
                                    verbose=verbose).fillna(0)
    
    exp_sim_combined.experimental=_eval_exp(exp_sim_combined.experimental,receiver_or_donor=receiver_or_donor,exp_thresh=exp_thresh)

    y_true=exp_sim_combined["experimental"].values
    y_scores=exp_sim_combined["simulated"]

    if PR_or_ROC=="PR":
        baseline = sum(y_true)/len(y_true)
        auc = average_precision_score(y_true, y_scores)
        precision, recall, thresholds = precision_recall_curve(y_true, y_scores)

        # Plot the curve
        prec_pts   = precision[:-1]
        recall_pts = recall[:-1]
        thr        = np.ravel(thresholds)

        # Avoid log(0) issues
        thr = np.clip(thr, a_min=1e-8, a_max=None)

        # Plot base PR curve
        disp = PrecisionRecallDisplay(precision=precision, recall=recall)
        disp.plot(ax=ax, marker='', label='PR curve')

        sc = ax.scatter(
            recall_pts, prec_pts,
            c=thr,
            cmap='viridis',
            norm=LogNorm(vmin=1e-6, vmax=1),
            s=30,
            alpha=0.9,
            edgecolor='none'
        )
        if plot_cbar:
            cbar = plt.colorbar(sc, ax=ax,shrink=0.7)
            cbar.set_label('Threshold (log scale)')

        ax.set_xlabel('Recall TP/(FN+TP)')
        ax.set_ylabel('Precision TP/(FP+TP)')
        
        # No skill line and annotation
        ax.axhline(y=baseline, color='black', linestyle='--',label=f"no-skill={baseline:.3f}")
        plt.text(.1, .75, f"Area Under PR-curve \n (PR-AUC): {auc:.3f}",transform=ax.transAxes)
    
    elif PR_or_ROC=="ROC":
        auc = roc_auc_score(y_true, y_scores)
        fpr, tpr, thresholds = roc_curve(y_true, y_scores)

        # Drop the first (infinite) threshold to avoid inf and align arrays
        fpr_pts = fpr[1:]
        tpr_pts = tpr[1:]
        thr = thresholds[1:]

        # Define fixed log normalization range (since thresholds > 0)
        vmin, vmax = 1e-6, 1
        norm = LogNorm(vmin=vmin, vmax=vmax)

        # Plot ROC curve
        RocCurveDisplay(fpr=fpr, tpr=tpr).plot(ax=ax, marker='', label='ROC curve')
        ax.set_xlabel('FPR (FP/(FP+TN))')
        ax.set_ylabel('TPR (TP/(FN+TP))')

        # Threshold-colored scatter
        sc = ax.scatter(
            fpr_pts, tpr_pts,
            c=thr,
            cmap='viridis',
            norm=norm,
            s=40,
            edgecolors='none',
            alpha=1.0
        )
        if plot_cbar:
            cbar = plt.colorbar(sc, ax=ax,shrink=0.7)
            cbar.set_label('Threshold (log scale)')

        # Diagonal no-skill line and annotation
        ax.plot([0, 1], [0, 1], '--', color='black', label='no-skill')
        ax.text(.1, .75, f"Area Under the ROC-curve\n(ROC-AUC): {auc:.3f}", transform=ax.transAxes)

    else:
        raise Exception("PR_or_ROC must be either 'PR' or 'ROC'")
    
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.set_title(title)
    ax.legend()

    return fig,ax


def plot_diagnostic_curves_method(simulation_dict,signal_exp,PR_or_ROC="PR",exp_thresh=1,simulation_method="steadiercom",metric="frequency",met_column="met_group",agg_method="sum",verbose=True):
    
    names = simulation_dict["receiver"][simulation_method].keys()
    fig,axs =  plt.subplots(2, len(names), figsize=(7*len(names), 12))

    if verbose:
        logging.info(f"plot_diagnostic_curve uses fillna(0) to for missing entries in simulations that exist in experimental results.")
        
    i=0

    for receiver_or_donor in ["receiver","donor"]:
        j=0

        for name in names:
            simulated_res = simulation_dict[receiver_or_donor][simulation_method][name]
            
            _ , _ = plot_diagnostic_curve(simulated_res,
                                signal_exp,
                                exp_thresh=exp_thresh,
                                PR_or_ROC=PR_or_ROC,
                                receiver_or_donor=receiver_or_donor,
                                metric=metric,
                                met_column=met_column,
                                agg_method=agg_method,
                                ax=axs[i][j],
                                verbose=(i==0 and j==0 and verbose),
                                plot_cbar=False)
            j+=1
        i+=1

    col_titles = [f"{name} \n {agg_method} of {metric}, exp. thresh={exp_thresh} " for name in names]
    for col_idx, col_title in enumerate(col_titles):
        pos = axs[0][col_idx].get_position()  # top row’s subplot in this column
        x_center = pos.x0 + pos.width / 2
        y_top = pos.y1 + 0.02  
        fig.text(x_center, y_top, col_title, ha='center', va='bottom', fontsize=16)

    fig.text(0.05, 0.75, "Receiver",fontsize=14)
    fig.text(0.05, 0.25, "Donor",fontsize=14)

    cmap = plt.cm.viridis
    vmin = 1e-6
    vmax = 1
    norm = LogNorm(vmin=vmin, vmax=vmax)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])  # required by mpl
    cbar = fig.colorbar(sm, ax=axs, location='right', fraction=0.03, pad=0.02,shrink=0.7)
    cbar.set_label('Threshold (log scale)')

    return fig,axs


def plot_diagnostic_curves_pvalues(pr_auc_pvalue_df,simulation_dict,signal_exp,receiver_or_donor="receiver"):

    elements = pr_auc_pvalue_df.index

    fig, axs = plt.subplots(len(elements),2,figsize=(10,4*len(elements)))

    i = 0

    verbose=True

    row_titles = []
    for  index in elements:
        dict_ = dict(zip(elements.names,index))

        agg_method = dict_["agg_method"]
        metric = dict_["metric"]
        simulation_method = dict_["simulation_method"]
        name = dict_["condition"]
        exp_thresh = dict_["exp_thresh"]
        if len(elements)>1:
            ax1 = axs[i,0]
            ax2 = axs[i,1]
        else:
            ax1 = axs[0]
            ax2 = axs[1]



        simulated_res = simulation_dict[receiver_or_donor][simulation_method][name]
        _,_ = plot_diagnostic_curve(simulated_res,
                                                        signal_exp,
                                                        PR_or_ROC="PR", 
                                                        exp_thresh=exp_thresh,
                                                        receiver_or_donor=receiver_or_donor,
                                                        metric=metric,
                                                        met_column="met_group",
                                                        agg_method=agg_method,
                                                        ax=ax1,
                                                        verbose=verbose,
                                                        plot_cbar=False,
                                                        )

        _,_ = plot_diagnostic_curve(simulated_res,
                                                        signal_exp,
                                                        PR_or_ROC="ROC", 
                                                        exp_thresh=exp_thresh,
                                                        receiver_or_donor=receiver_or_donor,
                                                        metric=metric,
                                                        met_column="met_group",
                                                        agg_method=agg_method,
                                                        ax=ax2,
                                                        verbose=verbose,
                                                        plot_cbar=False,
                                                        )
        verbose=False

        row_title = f"{name}, \n {agg_method} of {metric}, exp. thresh={exp_thresh}"
        row_titles.append(row_title)
        
        if i<len(elements)-1:
            ax1.set_xlabel("")
            ax2.set_xlabel("")
            ax1.set_xticklabels([])
            ax2.set_xticklabels([])
        ax2.set_yticklabels([])

        i+=1


    for row_axes, row_title in zip(axs, row_titles):
        # Find the vertical center of this row by averaging y-lims of axes positions
        pos = row_axes[0].get_position()  
        y_center = pos.y0 + pos.height / 2  
        fig.text(0.1, y_center, row_title, va='center', ha='right', fontsize=14)

    fig.text(0.3, 0.9, "PR",fontsize=16)
    fig.text(0.7, 0.9, "ROC",fontsize=16)

    cmap = plt.cm.viridis
    vmin = 1e-6
    vmax = 1
    norm = LogNorm(vmin=vmin, vmax=vmax)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])  # required by mpl
    cbar = fig.colorbar(sm, ax=axs, location='right', fraction=0.03, pad=0.02,shrink=0.5)
    cbar.set_label('Threshold (log scale)')

    return fig,axs
