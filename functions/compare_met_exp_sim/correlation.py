import pandas as pd
import numpy as np
import seaborn as sns

import scipy as sp
import matplotlib.pyplot as plt
import compare_met_exp_sim.general_funcs as general_funcs

from collections import defaultdict

import logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def find_correlation(simulated_res,signal_exp,simulation_method="steadiercom",receiver_or_donor="receiver",metric="frequency",agg_method="sum",met_column="met_group",verbose=True):
    if verbose:
        logging.info(f"find_correlation uses fillna(0) for missing entries in simulations that exist in experimental results.")
       

    exp_sim_combined = general_funcs.merge_sim_exp(
                simulated_res,
                signal_exp,
                simulation_method=simulation_method,
                metric=metric,
                agg_method=agg_method,
                receiver_or_donor=receiver_or_donor,
                met_column=met_column,
                verbose=verbose
            ).fillna(0)


    r, p = sp.stats.kendalltau(x=exp_sim_combined.experimental,y=exp_sim_combined.simulated)


    correlation = {"r":r,"p-value":p,"datapoints":exp_sim_combined.shape[0]}

    return correlation


def find_correlation_simulation_method(simulation_dict,signal_exp,simulation_method="steadiercom",metric="frequency",agg_method="sum",verbose=True):
    """
    Function to calculate correlation between experimental and simulated results.
    """
 
    # Process simulated result
    correlations_dict = defaultdict(dict)

    
    for receiver_or_donor in ["receiver", "donor"]:

        for i,name in enumerate(simulation_dict[receiver_or_donor][simulation_method].keys()):

            verbose_spec = (i==0 and verbose)
            simulated_res = simulation_dict[receiver_or_donor][simulation_method][name]

            correlation = find_correlation(simulated_res,
                                          signal_exp,
                                          simulation_method=simulation_method,
                                          receiver_or_donor=receiver_or_donor,
                                          metric=metric,
                                          agg_method=agg_method,
                                          verbose=verbose_spec)

            correlations_dict[(receiver_or_donor,simulation_method,name,agg_method,metric)] = correlation


    correlations_df = pd.DataFrame(correlations_dict).T
    correlations_df.index.names = ["role","simulation_method","condition","agg_method","metric"]
    correlations_df = correlations_df.sort_index()
    return correlations_df

def find_correlation_simulation_method_detailed(simulation_dict,signal_exp,member_or_met="member",simulation_method="steadiercom",metric="frequency",agg_method="sum",met_column="met_group",verbose=True):
    """
    Function to calculate correlation between experimental and simulated results.
    """
 
    # Process simulated result
    correlations_dict = defaultdict(dict)

    if member_or_met=="member":
        classes = signal_exp.index.get_level_values(0).unique()
    elif member_or_met==met_column:
        classes = signal_exp.index.get_level_values(1).unique()

    for receiver_or_donor in ["receiver", "donor"]:

        for i,name in enumerate(simulation_dict[receiver_or_donor][simulation_method].keys()):

            verbose_spec = (i==0 and verbose)
            
            for class_ in classes:
                simulated_res = simulation_dict[receiver_or_donor][simulation_method][name]
                
                if member_or_met=="member":
                    simulated_res_copy = simulated_res[simulated_res[receiver_or_donor]==class_].copy()
                    signal_exp_copy = signal_exp[signal_exp.index.get_level_values(0)==class_].copy()
                elif member_or_met==met_column:
                    simulated_res_copy = simulated_res[simulated_res[met_column]==class_].copy()
                    signal_exp_copy = signal_exp[signal_exp.index.get_level_values(1)==class_].copy()

                correlation = find_correlation(simulated_res_copy,
                                                signal_exp_copy,
                                                simulation_method=simulation_method,
                                                receiver_or_donor=receiver_or_donor,
                                                metric=metric,
                                                agg_method=agg_method,
                                                verbose=verbose_spec)

                correlations_dict[(receiver_or_donor,simulation_method,name,agg_method,metric,class_)] = correlation

    correlations_df = pd.DataFrame(correlations_dict).T

    correlations_df.index.names = ["role","simulation_method","condition","agg_method","metric",member_or_met]
    correlations_df = correlations_df.sort_index()

    return correlations_df


def make_scatterplot(simulated_res, signal_exp,simulation_method="steadiercom",receiver_or_donor="receiver",metric="frequency",agg_method="sum", met_column="met_group",ax=None,xlabel="",ylabel="",title="",legend=True,verbose=True):
    color = {"receiver":"black","donor":"black"}  
    color_palettes = {"receiver":"RdGy_r","donor":"RdGy_r"}
    
    if verbose:
            logging.info(f"make_scatterplot removes nan values in 'simulated' column for scatterplot and uses all values for density plot.")

    if not ax:
        fig, ax = plt.subplots(figsize=(5,5))
    else:
        fig = ax.get_figure()


    # Combine experimental and simulated results
    exp_sim_combined = general_funcs.merge_sim_exp(
        simulated_res,
        signal_exp,
        metric=metric,
        agg_method=agg_method,
        receiver_or_donor=receiver_or_donor,
        met_column=met_column,
        verbose=verbose
    )

    # Find entries with and without simulated values
    exp_sim_combined_fillna = exp_sim_combined.fillna(0) # For scatterplot (to represent non-simulated values in the plot)
    exp_sim_combined_is_nan = exp_sim_combined.copy() # For density plot and to quantify number of missing values - to show how many metabolites were not simulated and potential bias
    exp_sim_combined_is_nan.simulated = ~exp_sim_combined.simulated.isna()

    # Create twin axis for density plot
    ax_twin = ax.twinx()
    ax_twin.set_ylim([0,1])
    ax.set_zorder(2)
    ax_twin.set_zorder(1)
    ax.patch.set_visible(False)  # hide the 'canvas'

    # Plot density of simulated and not simulated
    sns.kdeplot(data=exp_sim_combined_is_nan,
                x="experimental",
                hue="simulated",
                hue_order=[True,False],
                cut=0,
                fill=True,
                alpha=0.6,
                linewidth=0.5,
                palette=color_palettes[receiver_or_donor],
                common_norm=False,
                ax=ax_twin,
                legend=True
                )
    # Plot scatter of experimental vs simulated

    sns.scatterplot(data=exp_sim_combined_fillna,
                    x="experimental",
                    y="simulated",
                    alpha=0.3,
                    color=color[receiver_or_donor],
                    ax=ax,
                    label="Experimental vs Simulated",
                    legend=legend
                    )

    ax.set_yscale('log')

    ax.set_ylabel(f"{ylabel}")
    ax.set_xlabel(f"{xlabel}")
    ax.set_title(f"{title}\n #simulated={(~exp_sim_combined_is_nan.simulated).sum()}, #not simulated={exp_sim_combined_is_nan.simulated.sum()}")

    correlation = find_correlation(simulated_res,
                                   signal_exp,
                                   simulation_method=simulation_method,
                                   receiver_or_donor=receiver_or_donor,
                                   metric=metric,
                                   agg_method=agg_method,
                                   met_column=met_column,
                                   verbose=verbose)
    plt.text(.1, .65, "Correlation (spearman) \n rho={:.2f}, p-value={:.3f}".format(correlation["r"],correlation["p-value"]), transform=ax.transAxes)

    return fig, ax, ax_twin

def make_scatterplots_simulation_method(simulation_dict, signal_exp, simulation_method="steadiercom", metric="frequency", agg_method="sum",  met_column="met_group",xlabel="Experimental - log2 conc. change",verbose=True):

    
    names = simulation_dict["receiver"][simulation_method].keys()
    fig, axs = plt.subplots(2, len(names), figsize=(8*len(names), 12))

    twin_axes = []
    for i, receiver_or_donor in enumerate(["receiver", "donor"]):
        
        for j, name in enumerate(names):
            
            verbose_spec = (j==0 and verbose)
            legend = (i==0 and j==len(names)-1)

            ylabel = f"Simulated - {agg_method} of {metric}"
            
            simulated_res = simulation_dict[receiver_or_donor][simulation_method][name]
            _,_,ax_twin = make_scatterplot(simulated_res, 
                             signal_exp, 
                             simulation_method=simulation_method,
                             receiver_or_donor=receiver_or_donor, 
                             metric=metric, 
                             agg_method=agg_method, 
                             met_column=met_column, 
                             ax=axs[i,j], 
                             verbose=verbose_spec,
                             legend=legend,
                             xlabel=xlabel,
                             ylabel=ylabel,
                             )

            if legend:
                sns.move_legend(ax_twin,loc="upper left",bbox_to_anchor=(1.1, 1))
                axs[i,j].legend(loc="upper left",bbox_to_anchor=(1.1, 0.7))
            else:
                ax_twin.get_legend().remove()

            twin_axes.append(ax_twin)

    col_titles = [f"{name}" for name in names]

    for col_idx, col_title in enumerate(col_titles):
        pos = axs[0][col_idx].get_position()  # top row’s subplot in this column
        x_center = pos.x0 + pos.width / 2
        y_top = pos.y1 + 0.02  
        fig.text(x_center, y_top, col_title, ha='center', va='bottom', fontsize=16)


    fig.text(0.05, 0.75, "Receiver",fontsize=14)
    fig.text(0.05, 0.25, "Donor",fontsize=14)

    return fig, axs, twin_axes


def plot_scatterplots_from_corr(correlations_df, simulation_dict,signal_exp,member_or_met="member",met_column="met_group",xlabel="Experimental - log2 conc. change"):

    parameters = correlations_df.index
    parameters_list = [dict(zip(parameters.names,param)) for param in parameters]

    rows = int(np.ceil(len(parameters_list)/3))
    cols = min(len(parameters_list),3)

    fig,axs = plt.subplots(rows,cols,figsize=(5*cols,5*rows))

    axs_twin = []
    for i,param_dict in enumerate(parameters_list):

        if rows==1:
            ax = axs[i]
        else:
            row = int(np.floor(i/3))
            col = i%3
            
            ax = axs[row,col]

        role = param_dict["role"]
        simulation_method = param_dict["simulation_method"]
        condition = param_dict["condition"]
        agg_method = "sum"
        metric = param_dict["metric"]

        simulated_res = simulation_dict[role][simulation_method][condition]
        

        if member_or_met=="member":
            member = param_dict["member"]
            signal_exp_member_or_met = signal_exp[signal_exp.index.get_level_values(0)==member].copy()
            simulated_res_member_or_met = simulated_res[simulated_res[role]==member].copy()
            title = f"{condition} \n -  member={member} - {role}"
        elif member_or_met=="met":
            met = param_dict[met_column]
            signal_exp_member_or_met = signal_exp[signal_exp.index.get_level_values(1)==met].copy()
            simulated_res_member_or_met = simulated_res[simulated_res[met_column]==met].copy()
            title = f"{condition} \n -  met={met} - {role}"

        elif not member_or_met:
            signal_exp_member_or_met = signal_exp
            simulated_res_member_or_met = simulated_res
            title = f"{condition}"
        else:
            raise Exception("member_or_met must be either 'met', 'member' or None")

        _, _, ax_twin = make_scatterplot(simulated_res_member_or_met,
                                    signal_exp_member_or_met,
                                    simulation_method=simulation_method,
                                    receiver_or_donor=role,
                                    metric=metric,
                                    agg_method=agg_method,
                                    xlabel=xlabel,
                                    ylabel=f"Simulated - {agg_method} of {metric}",
                                    title=title,
                                    ax=ax,
                                    legend=True,
                                    verbose=False)
        if rows==1 and i!=(cols-1):
            ax_twin.set_ylabel("")
            ax_twin.set_yticklabels([])


        elif rows>1 and col!=(cols-1):
            ax_twin.set_ylabel("")
            ax_twin.set_yticklabels([])



        axs_twin.append(ax_twin)





    return fig, axs, axs_twin
        



