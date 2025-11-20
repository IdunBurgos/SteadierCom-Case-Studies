import pandas as pd
import numpy as np
import scipy as sp
import seaborn as sns

import matplotlib.pyplot as plt
import warnings

import logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

from .general_funcs import merge_exp_sim


def _find_correlation(exp_sim_combined_isna,verbose=True):

    exp_sim_combined = exp_sim_combined_isna[exp_sim_combined_isna["Was simulated"]]
    not_simulated_count = (~exp_sim_combined_isna["Was simulated"]).sum()
    simulated_count = exp_sim_combined_isna["Was simulated"].sum()


    r, p = sp.stats.kendalltau(x=exp_sim_combined.experimental,y=exp_sim_combined.simulated)

    correlation = {"r":r,"p-value":p}
    correlation["#simulated"]=simulated_count
    correlation["#not simulated"]=not_simulated_count
    return correlation

def find_correlation_simulation_method(simulation_dict,signal_exp,simulation_method="steadiercom",metric="frequency"):
    
    correlations_dict = {}
    for coop_or_comp in ["coop","comp","coop-comp"]:

        for name, simulation_res in simulation_dict[simulation_method].items():

            assert isinstance(simulation_res,pd.DataFrame), "simulation_res is not a pd.DataFrame"
            assert not simulation_res.empty, "simulation_res is empty"

            exp_sim_combined_isna = merge_exp_sim(simulation_res,
                                        signal_exp,
                                        coop_or_comp=coop_or_comp,
                                        metric=metric)
            assert signal_exp.shape[0]==exp_sim_combined_isna.shape[0], f"The size of signal_exp (experimental signal) and merged simulation and experimental results do not match {signal_exp.shape[0]}!={exp_sim_combined_isna.shape[0]}."

            correlation = _find_correlation(exp_sim_combined_isna,
                                          signal_exp,
                                          )

            correlations_dict[(coop_or_comp, simulation_method,name, metric)]=correlation

    correlations_df = pd.DataFrame(correlations_dict).T
    correlations_df.index.names = ["coop_or_comp","simulation_method","condition","metric"]
    correlations_df = correlations_df.sort_index()

    return correlations_df
        
def find_correlation_simulation_method_detailed(simulation_dict,signal_exp,simulation_method="steadiercom",metric="frequency"):
    
    correlations_dict = {}

    main_members = signal_exp.index.get_level_values(1).unique()
    for coop_or_comp in ["coop","comp","coop-comp"]:

        for name, simulation_res in simulation_dict[simulation_method].items():
            for main_member in main_members:

                assert isinstance(simulation_res,pd.DataFrame), "simulation_res is not a pd.DataFrame"
                assert not simulation_res.empty, "simulation_res is empty"

                
                exp_sim_combined_isna = merge_exp_sim(simulation_res,
                                            signal_exp,
                                            coop_or_comp=coop_or_comp,
                                            metric=metric)
                
                # Only look at main member (the one being competed or cooperated with)
                exp_sim_combined_isna_member = exp_sim_combined_isna[exp_sim_combined_isna.index.get_level_values(1)==main_member]
                signal_exp_member = signal_exp.xs(main_member,level="member_main",drop_level=False)
                assert signal_exp_member.shape[0]==exp_sim_combined_isna_member.shape[0], f"The size of signal_exp (experimental signal) and merged simulation and experimental results do not match {signal_exp_member.shape[0]}!={exp_sim_combined_isna_member.shape[0]}."

                correlation = _find_correlation(exp_sim_combined_isna_member,verbose=True)
                correlations_dict[(coop_or_comp, simulation_method,name, metric, main_member)]=correlation

    correlations_df = pd.DataFrame(correlations_dict).T
    correlations_df.index.names = ["coop_or_comp","simulation_method","condition","metric","member"]
    correlations_df = correlations_df.sort_index()

    return correlations_df


def make_scatterplot(simulated_res,signal_exp,coop_or_comp="coop",metric="frequency",ax=None,xlabel="",ylabel="",title="",include_legend=True):

    if not ax:
        fig,ax = plt.subplots(1,1)

    else:
        fig = ax.get_figure()
        
    exp_sim_combined_isna = merge_exp_sim(simulated_res,
                                                     signal_exp,
                                                     coop_or_comp=coop_or_comp,
                                                     metric=metric)
    
    assert "Was simulated" in exp_sim_combined_isna.columns, "exp_sim_combined_isna does not have the column 'Was simulated'"
    exp_sim_combined = exp_sim_combined_isna[exp_sim_combined_isna["Was simulated"]]

    ax_twin = ax.twinx()
    ax.set_zorder(2)
    ax.patch.set_alpha(0)
    ax_twin.set_zorder(1)

    sns.scatterplot(data=exp_sim_combined,x="experimental", y=f"simulated",alpha=0.5,ax=ax,color="black")
    
    # Plot density of simulated and not simulated
    g = sns.kdeplot(data=exp_sim_combined_isna,
                x="experimental",
                hue="Was simulated",
                hue_order=[True,False],
                cut=0,
                fill=True,
                alpha=0.6,
                linewidth=0.5,
                palette="RdGy_r",
                common_norm=False,
                ax=ax_twin,
                legend=include_legend
                )
    
    if include_legend:
        sns.move_legend(g, "upper right", bbox_to_anchor=(1.4, 1.0), title='Was simulated')
    
    not_simulated_count = (~exp_sim_combined_isna["Was simulated"]).sum()
    simulated_count = exp_sim_combined_isna["Was simulated"].sum()
    assert not_simulated_count+simulated_count==exp_sim_combined_isna.shape[0], "Number of simulated+non-simulated items do not match the shape of exp_sim_combined_isna"

    ax.set_title(f"{title} \n  #simulated={simulated_count}, #not simulated={not_simulated_count}")
    ax.set_ylabel(f"{ylabel}")
    ax.set_xlabel(f"{xlabel}")

    correlation = _find_correlation(exp_sim_combined)
    plt.text(.1, .85, "Kendall's Tau r={:.2f}, p-value={:.4f}".format(correlation["r"],correlation["p-value"]), transform=ax.transAxes)
    not_simulated = exp_sim_combined_isna[~exp_sim_combined_isna["Was simulated"]].index.values
    return fig,ax,ax_twin,not_simulated


def make_scatterplots_simulation_method(simulation_dict,signal_exp,simulation_method="steadiercom",metric="frequency"):

    names = simulation_dict[simulation_method].keys()

    fig, axs = plt.subplots(3,len(names),figsize=(len(names)*6,15))

    axs_twin = []
    not_simulated_list = []
    for i,coop_or_comp in enumerate(["coop","comp","coop-comp"]):

        j=0
        
        for name in names:
            
            if len(names)==1:
                ax = axs[i]
            else:
                ax = axs[i,j]
            
            include_legend= (i==0) & (j==len(names)-1)
            simulated_res = simulation_dict[simulation_method][name]
            _,_,ax_twin,not_simulated = make_scatterplot(simulated_res,
                                    signal_exp,
                                    coop_or_comp=coop_or_comp,
                                    metric=metric,
                                    ax=ax,
                                    xlabel="gLV",
                                    ylabel=f"scaled {metric})",
                                    include_legend=include_legend,
                                    )
            
            axs_twin.append(ax_twin)
            not_simulated_list.extend(not_simulated)
            if j>0:
                ax.set_ylabel("")
            if j<len(names)-1:
                ax_twin.set_ylabel("")

            j+=1

    for col_idx, col_title in enumerate(names):
        pos = axs[0][col_idx].get_position() 
        x_center = pos.x0 + pos.width / 2
        y_top = pos.y1 + 0.02  
        fig.text(x_center, y_top, col_title, ha='center', va='bottom', fontsize=16)


    fig.text(0.0, 0.75, "Cooperation",fontsize=14)
    fig.text(0.0, 0.45, "Competition",fontsize=14)
    fig.text(0.0, 0.2, "Cooperation-\nCompetition",fontsize=14)
    
    return fig,axs,axs_twin,not_simulated_list
