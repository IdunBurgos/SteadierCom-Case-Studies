import pandas as pd
import numpy as np
import scipy as sp
import seaborn as sns

import matplotlib.pyplot as plt

import logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

def _cooperation_agg(simulation_res, metric="rate"):
    """

    Finds fraction of flux going to member_main coming from member_second. 

    Details: Finds the norm of flux for each member and calculate the sum of normalized flux coming to member_main
    from member_second. 
    """

    df = simulation_res.copy()

    # total inflow per receiver (guard zeros)
    total = df.groupby("receiver")[metric].transform("sum")
    total_safe = total.replace(0, np.nan)

    # normalized contribution to A
    df["norm_flux"] = df[metric] / total_safe
    df["norm_flux"] = df["norm_flux"].fillna(0.0)

    # drop environment as a 'cooperating member'
    df = df[df["donor"] != "environment"]

    # sums all fluxes in donor-receiver relationship
    coop = df.groupby(["donor", "receiver"], as_index=True)["norm_flux"].sum()
    coop.index.names = ["member_second", "member_main"]
    coop.name = "simulated"
    return coop


def _competition_agg(simulation_res, metric="rate", met_column="compound"):
    
    """
    Finds the fraction of flux to member_main that member_second is competing for. The magnitude
    of flux going to member_second is not 
    """
    df = simulation_res.copy()

    # normalize inflow per receiver (guard zeros)
    total = df.groupby("receiver")[metric].transform("sum")
    total_safe = total.replace(0, np.nan)
    df["norm_flux"] = df[metric] / total_safe
    df["norm_flux"] = df["norm_flux"].fillna(0.0)

    # receiver x compound matrix of normalized inflow
    M = (df.groupby(["receiver", met_column], as_index=True)["norm_flux"]
            .sum()
            .unstack(fill_value=0.0))

    U = (M > 1e-6).astype(int)

    C = M.values @ U.values.T
    np.fill_diagonal(C, np.nan)

    comp = (pd.DataFrame(C, index=M.index, columns=M.index)
              .stack()
              .rename("simulated"))

    # index levels: (member_second=B competitor, member_main=A affected)
    comp.index.names = ["member_main", "member_second"]

    comp = comp.swaplevel(0, 1).sort_index()
    return comp


def _merge_exp_sim(simulation_res,exp_ser, coop_or_comp="coop",metric="frequency"):

    # Get score for cooperation by sim and exp
    if coop_or_comp=="coop":
        simulated_score = _cooperation_agg(simulation_res,
                                                metric=metric,
                                                )
    elif coop_or_comp=="comp":
        simulated_score = _competition_agg(simulation_res,
                                    metric=metric,
                                    )
    else:
        raise Exception(f"coop_or_comp={coop_or_comp}, but must be either 'coop' or 'comp'")

    # Combine simulated score with experimental values
    exp_sim_combined = pd.concat([exp_ser,simulated_score],axis=1)   

    return exp_sim_combined


def _find_correlation(exp_sim_combined,verbose=True):

    r, p = sp.stats.kendalltau(x=exp_sim_combined.experimental,y=exp_sim_combined.simulated)

    correlation = {"r":r,"p-value":p,"datapoints":exp_sim_combined.shape[0]}
    return correlation


def make_scatterplot(simulation_dict,signal_exp,coop_or_comp="coop",simulation_method="steadiercom",name="SteadierCom - Fixed growth",metric="frequency",ax=None,xlabel="",ylabel="",title=""):

    if not ax:
        fig,ax = plt.subplots(1,1)

    else:
        fig = ax.get_figure()
    
    if coop_or_comp=="coop-comp":
        simulated_res = simulation_dict[simulation_method][name]
        exp_sim_combined_coop =  _merge_exp_sim(simulated_res,
                                      signal_exp,
                                      coop_or_comp="coop",
                                      metric=metric,
                                      ).fillna(0)
        exp_sim_combined_comp =  _merge_exp_sim(simulated_res,
                                      signal_exp,
                                      coop_or_comp="comp",
                                      metric=metric,
                                      ).fillna(0)

        
        exp_sim_combined = exp_sim_combined_coop.copy()
        exp_sim_combined.simulated= exp_sim_combined_coop.simulated-exp_sim_combined_comp.simulated

        # Remove simulated data that was 0 - means essentially now prediction of neitehr competition nor cooperation (fillna(0))
        isna = (exp_sim_combined_coop.simulated<1e-6) & (exp_sim_combined_comp.simulated<1e-6)
        

    
    else:
        simulated_res = simulation_dict[simulation_method][name]
        exp_sim_combined = _merge_exp_sim(simulated_res,
                                        signal_exp,
                                        coop_or_comp=coop_or_comp,
                                        metric=metric,
                                        )

        # Remove entries that do not have a match...
        isna = exp_sim_combined.simulated.isna()

    exp_sim_combined_isna = exp_sim_combined.copy()

    exp_sim_combined_isna["Was simulated"] = ~isna
    exp_sim_combined = exp_sim_combined[~isna]

    ax_twin = ax.twinx()
    ax.set_zorder(2)
    ax.patch.set_alpha(0)
    ax_twin.set_zorder(1)

    sns.scatterplot(data=exp_sim_combined,x="experimental", y=f"simulated",alpha=0.5,ax=ax,color="black")
    # Plot density of simulated and not simulated
    
    sns.kdeplot(data=exp_sim_combined_isna,
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
                legend=True
                )

    ax.set_title(f"{title} \n  #simulated={(~isna).sum()}, #not simulated={isna.sum()}")
    ax.set_ylabel(f"{ylabel}")
    ax.set_xlabel(f"{xlabel}")

    correlation = _find_correlation(exp_sim_combined)
    plt.text(.1, .85, "Kendall's Tau r={:.2f}, p-value={:.4f}".format(correlation["r"],correlation["p-value"]), transform=ax.transAxes)
    not_simulated = isna[isna].index.values
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
            
            _,_,ax_twin,not_simulated = make_scatterplot(simulation_dict,
                                    signal_exp,
                                    coop_or_comp=coop_or_comp,
                                    simulation_method=simulation_method,
                                    name=name,
                                    metric=metric,
                                    ax=ax,
                                    xlabel="gLV",
                                    ylabel=f"scaled {metric})"
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
