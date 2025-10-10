import pandas as pd
import seaborn as sns

import scipy as sp
import matplotlib.pyplot as plt
import compare_met_exp_sim.general_funcs as general_funcs

from collections import defaultdict

import logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def find_correlation(simulation_dict,signal_exp,name="SteadierCom - fixed abundance",simulation_method="steadiercom",receiver_or_donor="receiver",metric="frequency",agg_method="sum",verbose=True):
    if verbose:
        logging.info(f"find_correlation uses dropna() to for missing entries in simulations that exist in experimental results.")
       
    simulated_res = simulation_dict[receiver_or_donor][simulation_method][name]

    exp_sim_combined = general_funcs.merge_sim_exp(
                simulated_res,
                signal_exp,
                metric=metric,
                agg_method=agg_method,
                receiver_or_donor=receiver_or_donor,
                met_column="met_group",
                verbose=verbose
            ).dropna()

    
    r, p = sp.stats.kendalltau(x=exp_sim_combined.experimental,y=exp_sim_combined.simulated)

    correlation = {"r":r,"p-value":p}

    return correlation


def find_correlation_simulation_method(simulation_dict,signal_exp,simulation_method="steadiercom",metric="frequency",agg_method="sum",verbose=True):
    """
    Function to calculate spearman correlation between experimental and simulated results.
    """
 
    # Process simulated result
    correlations_dict = defaultdict(dict)

    
    for receiver_or_donor in ["receiver", "donor"]:

        for i,name in enumerate(simulation_dict[receiver_or_donor][simulation_method].keys()):

            verbose_spec = (i==0 and verbose)

            correlation = find_correlation(simulation_dict,
                                          signal_exp,
                                          name=name,
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


def make_scatterplot(simulation_dict, signal_exp,name="SteadierCom - fixed abundance",simulation_method="steadiercom",receiver_or_donor="receiver",metric="frequency",agg_method="sum", met_column="met_group",ax=None,verbose=True,legend=True):
    color = {"receiver":"black","donor":"black"}  
    color_palettes = {"receiver":"RdGy_r","donor":"RdGy_r"}
    
    if verbose:
            logging.info(f"make_scatterplot removes nan values in 'simulated' column for scatterplot and uses all values for density plot.")

    if not ax:
        fig, ax = plt.subplots(figsize=(5,5))
    else:
        fig = ax.get_figure()


    # Combine experimental and simulated results
    simulated_res = simulation_dict[receiver_or_donor][simulation_method][name]
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
    exp_sim_combined_no_nan = exp_sim_combined.dropna()
    exp_sim_combined_is_nan = exp_sim_combined.copy()
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
    sns.scatterplot(data=exp_sim_combined_no_nan,
                    x="experimental",
                    y="simulated",
                    alpha=0.3,
                    color=color[receiver_or_donor],
                    ax=ax,
                    label="Experimental vs Simulated",
                    legend=legend
                    )



    ax.set_yscale('log')

    ax.set_title(f"{name}, {agg_method} of {metric} \n #simulated={exp_sim_combined_no_nan.shape[0]}, #not simulated={exp_sim_combined.simulated.isna().sum()}")

    ax.set_ylabel(f"Simulated - {agg_method} of {metric}")
    ax.set_xlabel(f"Experimental - log2 conc. change")

    #ax.axvspan(-1, 1, alpha=0.1, color='grey')
    r, p = sp.stats.kendalltau(x=exp_sim_combined_no_nan.experimental,y=exp_sim_combined_no_nan.simulated)
    plt.text(.1, .65, "Correlation (Kendall's Tau) \n r={:.2f}, p-value={:.3f}".format(r,p), transform=ax.transAxes)

    return fig, ax, ax_twin

def make_scatterplots_simulation_method(simulation_dict, signal_exp, simulation_method="steadiercom", metric="frequency", agg_method="sum",  met_column="met_group",verbose=True):

    
    names = simulation_dict["receiver"][simulation_method].keys()
    fig, axs = plt.subplots(2, max(len(names),2), figsize=(6*max(len(names),2), 10))

    twin_axes = []
    for i, receiver_or_donor in enumerate(["receiver", "donor"]):
        
        for j, name in enumerate(names):
            
            verbose_spec = (j==0 and verbose)

            legend = (i==0 and j==len(names)-1)

            _,_,ax_twin = make_scatterplot(simulation_dict, 
                             signal_exp, 
                             name=name,
                             simulation_method=simulation_method, 
                             receiver_or_donor=receiver_or_donor, 
                             metric=metric, 
                             agg_method=agg_method, 
                             met_column=met_column, 
                             ax=axs[i,j], 
                             verbose=verbose_spec,
                             legend=legend)

            
            if legend:
                sns.move_legend(ax_twin,loc="upper left",bbox_to_anchor=(1.1, 1))
                axs[i,j].legend(loc="upper left",bbox_to_anchor=(1.1, 0.7))
            else:
                ax_twin.get_legend().remove()

            twin_axes.append(ax_twin)


    fig.text(-0.06, 0.75, "Receiver",fontsize=14)
    fig.text(-0.06, 0.25, "Donor",fontsize=14)
    fig.tight_layout()

    return fig, axs, twin_axes
