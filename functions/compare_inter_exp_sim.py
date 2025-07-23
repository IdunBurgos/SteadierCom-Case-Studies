import pandas as pd
import numpy as np
import scipy as sp
import seaborn as sns

import matplotlib.pyplot as plt

import logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def cooperation_data(simulation_res,metric="frequency",count_or_sum = "count",sim_thresh =0.1):
    """
    Function to score level of pairwise cooperation between members based on simulation results.
    
    There are two main options: 
        - count cooperative interactions above a threshold >sim_thresh (donor-receiver)
        - sum up cooperative interactions (donor-receiver)
    """
    
    if count_or_sum=="sum" and sim_thresh:
        logging.warning(f"sim_thresh>{sim_thresh} is not active when using 'sum' alternative for count_or_sum")
        
    if count_or_sum!="sum" and count_or_sum!="count":
        raise Exception(f"count_or_sum={count_or_sum}, but must be either 'sum' or 'count'")
        
    # Process simulated result
    simulated = simulation_res.copy()
    
    # NB: Here we look at donor-receiver-pairs
    if count_or_sum=="count":
        # If 'count' use threshold and count instances above threshold for each receiver-donor-combination
        simulated_red = simulated[simulated[metric]>sim_thresh]
        simulated_grouped = simulated_red[["receiver","donor"]].value_counts() 
        simulated_grouped.name=f"Cooperation met {count_or_sum} - {metric}"
        
    elif count_or_sum=="sum":
        # If 'sum' sum up interaction for each receiver-donor-combination
        simulated_grouped = simulated.groupby(["receiver","donor"]).sum()[metric]
        simulated_grouped.name=f"Cooperation met {count_or_sum} - {metric}"
    return simulated_grouped

def competition_data(simulation_res,metric="frequency",count_or_sum = "count",sim_thresh =0.1):
    """
    Function to score level of pairwise competition based on simulation results.
    
    There are two main options: 
        - count number of metabolite both are receiveing above a threshold >sim_thresh (receiver-metabolite)
        - sum up rate/frequency of compounds consumed by both members (receiver-metabolite)
    
    NB: both multiply the receiver-compound matrix by each other and sum the results.
    
    """
    if count_or_sum=="sum" and sim_thresh:
        logging.warning(f"sim_thresh>{sim_thresh} is not active when using 'sum' alternative for count_or_sum")
        
    if count_or_sum!="sum" and count_or_sum!="count":
        raise Exception(f"count_or_sum={count_or_sum}, but must be either 'sum' or 'count'")
    
    if "environment" not in simulation_res.donor.values:
        logging.warning(f"'environment' is not among donors. You are therefore currently only studying competition for compounds donated by other members.")
        
    ## Process simulated result
    simulated = simulation_res.copy()

    # Filter out the less frequently shared compounds
    simulated_red = simulated[simulated[metric]>sim_thresh]
    if count_or_sum=="count":
        counts = (simulated_red[["receiver","compound"]].value_counts().unstack().fillna(0)>0).map(int)
        
    elif count_or_sum=="sum":
        counts = simulated_red[["receiver","compound",metric]].groupby(["receiver","compound"]).sum()[metric].unstack().fillna(0)
    
    # Calculate overlapping compound usage
    simulated_grouped = counts.dot(counts.T).stack()
    simulated_grouped.name=f"Competition met {count_or_sum} - {metric}"
    
    simulated_grouped.index.names = ["receiver","donor"]

    return simulated_grouped


def coop_or_comp_exp_communities(simulated_dict,exp_ser,coop_or_comp="coop",metric="frequency",count_or_sum = "count",sim_thresh =0.1):
    """
    Function to get score for competition/cooperation for all communities that also joins simulation data for 
    """
    
    if count_or_sum=="sum" and sim_thresh:
        logging.warning(f"sim_thresh>{sim_thresh} is not active when using 'sum' alternative for count_or_sum")
        
    if count_or_sum!="sum" and count_or_sum!="count":
        raise Exception(f"count_or_sum={count_or_sum}, but must be either 'sum' or 'count'")
    
    if coop_or_comp!="coop" and coop_or_comp!="comp":
        raise Exception(f"coop_or_comp={coop_or_comp}, but must be either 'coop' or 'comp'")
    
    i = 0
    
    data_dict = {}
    
    for name, simulation_res in simulated_dict.items():
        
        # Get score for cooperation by sim and exp
        if coop_or_comp=="coop":
            kind="Cooperation"
            simulated_score = cooperation_data(simulation_res,
                                                 metric=metric,
                                                 count_or_sum=count_or_sum,
                                                 sim_thresh=sim_thresh)
        elif coop_or_comp=="comp":
            kind="Competition"
            simulated_score = competition_data(simulation_res,
                                     metric=metric,
                                     count_or_sum=count_or_sum,
                                     sim_thresh=sim_thresh)

        # Combine simulated score with experimental values
        exp_sim_combined = pd.concat([exp_ser,simulated_score],axis=1)       

        data_dict[name]=exp_sim_combined

        i+=1
    
    exp_name = exp_sim_combined.columns.values[0]
    
    label_dict = {"coop_or_comp":coop_or_comp,"metric":metric,"count_or_sum":count_or_sum,"sim_thresh":sim_thresh,"exp_name":exp_name}

    return data_dict,label_dict

def plot_coop_or_comp(exp_sim_combined,name,ax=False,label_dict = {"coop_or_comp":np.nan,"metric":np.nan,"count_or_sum": np.nan,"sim_thresh":np.nan},verbose=True):
    """
    Function for making scatterplots comparing 'experimental'-interaction score with simulated-interaction score.
    """
    if verbose:
        logging.info("plot_coop_comp_venturelli filters exp_sim_combined (receiver!=donor + exp_sim_combined.dropna())")
    
    if ax:
        has_ax = True
    else:
        has_ax = False
        fig,ax = plt.subplots(1,1)
    
    coop_or_comp = label_dict["coop_or_comp"]
    if coop_or_comp=="coop":
        kind="Cooperation"
    elif coop_or_comp=="comp":
        kind="Competition"
    
    metric = label_dict["metric"]
    count_or_sum = label_dict["count_or_sum"]
    sim_thresh = label_dict["sim_thresh"]
    exp_name = label_dict["exp_name"]

    # Remove entries that do not have a match...
    exp_sim_combined = exp_sim_combined.dropna()
    # Remove entries where receiver==donor
    exp_sim_combined = exp_sim_combined.loc[exp_sim_combined.index.get_level_values(0)!=exp_sim_combined.index.get_level_values(1)]


    sns.scatterplot(ax=ax, x=exp_name, y=f"{kind} met {count_or_sum} - {metric}", data=exp_sim_combined,alpha=0.5)

    ax.set_title(name)
    r, p = sp.stats.spearmanr(exp_sim_combined)
    plt.text(.1, .85, "Spearman's r={:.2f}, p-value={:.4f}".format(r,p), transform=ax.transAxes)
        
    if has_ax:
        return ax
    else:
        return fig,ax
    
    

    
def plot_coop_comp_combined(simulated_dict_rec,simulated_dict_inter,exp_ser,metric="frequency",count_or_sum = "count",sim_thresh =0.1,verbose=True):
    
    """
    Function for plotting scatterplots of the output of simulated (coop-comp) interactions vs. the 'experimental' interaction score. 
    """
    
    if verbose:
        logging.info("'plot_coop_comp_combined' filters data: receiver!=donor. Before combining coop and comp, we use fillna(0)")

    data_comp,label_dict = coop_or_comp_exp_communities(simulated_dict_rec,
                                        exp_ser,
                                        coop_or_comp="comp",
                                        metric=metric,
                                        count_or_sum=count_or_sum,
                                        sim_thresh=sim_thresh)

    data_coop,label_dict = coop_or_comp_exp_communities(simulated_dict_inter,
                                        exp_ser,
                                        coop_or_comp="coop",
                                        metric=metric,
                                        count_or_sum=count_or_sum,
                                        sim_thresh=sim_thresh)
    
    comp_df = pd.concat(data_comp)
    coop_df = pd.concat(data_coop)
    comp_df.index.names = [None, "receiver","donor"]
    diff_df = pd.concat([coop_df.iloc[:,-1],comp_df.iloc[:,-1]],axis=1).fillna(0)
    diff = diff_df.iloc[:,0]-diff_df.iloc[:,1]
    
    fig, axs = plt.subplots(1,3,figsize=(18,5))

    diff_new = diff.unstack(level=0)
    
    exp_name = label_dict["exp_name"]
    
    i = 0
    for name,series in diff_new.T.iterrows():

        # Remove entries that do not have a match...
        combined = pd.concat([series,exp_ser],axis=1)
        combined = combined[combined.index.get_level_values(0)!=combined.index.get_level_values(1)]
        combined.plot.scatter(y=name,x=exp_name,ax=axs[i],alpha=0.5)
        axs[i].set_ylabel("Count simulated (coop-comp)")
        axs[i].set_title(name)

        r, p = sp.stats.spearmanr(combined)
        plt.text(.1, .85, "Spearman's r={:.2f}, p-value={:.4f}".format(r,p), transform=axs[i].transAxes)
        i+=1
        
    plt.tight_layout()
    return fig,axs