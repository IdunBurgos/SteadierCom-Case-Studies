import pandas as pd
import numpy as np
import seaborn as sns

from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import matplotlib as mpl
from collections import OrderedDict

import logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

def make_co_growth_count_matrix(abundance_filtered,community,verbose=True):
    """
    This function creates a matrix that shows the number of samples where member A and member B are growing together. It sets diagonal entries (memberA==memberB) to 0.
    """
    if verbose:
        logger.info("Make_co_growth_count_matrix sets diagonal entries (memberA==memberB) to 0.")
    
    # Filter out data not related to community

    abundance_filtered_subdf = abundance_filtered.query(f"community=='{community}'")
        
    # Create data series with 'sample','member' as index and 'relative_abundance' as values
    member_abund = abundance_filtered_subdf.set_index(["sample","member"])["relative_abundance"].unstack().fillna(0)

    co_growth_count_df = (member_abund>0).map(int)
    
    # Find cases of co-growth
    co_growth_count_ser = co_growth_count_df.T.dot(co_growth_count_df).stack()
    
    # Annotate data series and revert to dataframe
    co_growth_count_ser.index.names = ["member1","member2"]
    co_growth_count_ser.name="co-growth"
    co_growth_count = co_growth_count_ser.unstack().fillna(0)
    
    # Set diagonal entries (memberA==memberB) to 0
    np.fill_diagonal(co_growth_count.values, 0)
    
    return co_growth_count,member_abund

def order_co_growth_count(co_growth_count):
    """
    This function orders members based on how similar they are in showing up in co-growth with other members.
    """
    condensed = squareform(co_growth_count.values)
    linkage_matrix = linkage(condensed, method="average")
    order = leaves_list(linkage_matrix)
    order_of_members = co_growth_count.index[order]
    
    return list(order_of_members)


def plot_co_growth_count(co_growth_count,title,ax=False,order_of_members=False,**kwargs):
    """
    This function plots the co-growth matrix as a heatmap for a specific condition. It's optional to include the order of members to order the matrix to be able to study the differences between conditions. 
    """
    
    if not ax:
        fig,ax = plt.subplots(1,1)
    
    if not order_of_members:
        order_of_members = order_co_growth_count(co_growth_count)
    
    co_growth_count =  co_growth_count.reindex(index=order_of_members, columns=order_of_members, fill_value=0) 
    
    co_growth_count_ordered = co_growth_count.loc[order_of_members, order_of_members]
    
    g = sns.heatmap(co_growth_count_ordered,
                    cmap="binary",
                    ax=ax,
                    square=True,
                    vmin=0, 
                    **kwargs)
    g.set_title(f"{title}", y=1.02)
    
    if ax:
        return ax
    else:
        return fig,ax
    

def plot_co_growth_count_venturelli(abundance_filtered,vmax=40):
    """
    This function plots all the co-growth heatmaps for the different conditions in the Venturelli case study.
    """
    # For NONE community: Find co-growth matrix and order
    co_growth_count_none,member_abund_none = make_co_growth_count_matrix(abundance_filtered,"NONE",verbose=True)
    order_of_members_none = order_co_growth_count(co_growth_count_none)
    
    # START PLOTTING
    fig,axs = plt.subplots(5,3, figsize=(12, 23))
    axs_flat = [ax for axs_sub in axs for ax in axs_sub]
    
    ax = plot_co_growth_count(co_growth_count_none,
                              title="Full",
                              ax=axs_flat[0],
                              order_of_members=order_of_members_none,
                              vmax=vmax,
                              cbar=False)

    # Prepare for loop
    communities = list(abundance_filtered.community.sort_values().unique())
    communities.remove("NONE")

    member_abund_dict = {"Full":member_abund_none}

    for i,community in enumerate(communities):
        # For each community: Find co-growth matrix and order
        co_growth_count,member_abund=make_co_growth_count_matrix(abundance_filtered,community)
        
        ax = plot_co_growth_count(co_growth_count,
                              title=community,
                              ax=axs_flat[i+1],
                              order_of_members=order_of_members_none,
                              vmax=vmax,
                              cbar=False)
        
        member_abund_dict[community]=member_abund
        
    # Add common 
    norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
    sm = mpl.cm.ScalarMappable(cmap="binary", norm=norm)
    sm.set_array([])  
    cbar_ax = fig.add_axes([0.95, 0.6, 0.015, 0.2])  # [left, bottom, width, height]
    fig.colorbar(sm, cax=cbar_ax, label="# co-growth samples")
    
    axs_flat[-1].set_axis_off()
    axs_flat[-2].set_axis_off()
    
    
def plot_co_growth_count_guillen(abundance_filtered):
    co_growth_count,member_abund_df = make_co_growth_count_matrix(abundance_filtered,"ecoli_auxotrophs")
    order_of_members = order_co_growth_count(co_growth_count)
    
    # START PLOTTING
    fig,ax = plt.subplots(1,1, figsize=(10, 10))

    ax = plot_co_growth_count(co_growth_count,
                          title="E.coli auxotrophs",
                          ax=ax,
                          order_of_members=order_of_members,
                          )
    
def plot_hist_helper(abundance_df,community,members,co_or_ind_str,palette,ax,title=False):
    """
    Helper functiom for plotting histogram of predicted community sizes in plot_co_ind_growers functions.
    
    """
    abundance_co_growers = abundance_df[abundance_df.member.isin(members)]
    counts  = abundance_df.value_counts("sample").reset_index()
    abundance_counts = pd.merge(abundance_co_growers,counts,on="sample")
    abundance_counts = abundance_counts.sort_values("member")

    g = sns.histplot(abundance_counts,
                     x="count",
                     hue="member",
                     ax=ax,
                     palette=palette,
                     element="bars",
                     multiple='dodge',
                     shrink=.7,
                     discrete=True)
    if title:
        ax.set_title(f"{title}: {co_or_ind_str}")
    else:
        ax.set_title(f"{community}: {co_or_ind_str}")
    ax.set_ylabel("# samples")
    ax.set_xlabel("# members in same sample")
    ax.set_xticks([1,2,3,4,5,6,7])
    return ax
    
def plot_density_helper(abundance_df,community,members,co_or_ind_str,palette,ax,title=False):
    """
    Helper functiom for plotting relative abundances in plot_co_ind_growers functions.
    """
    
    abundance_df = abundance_df.sort_values("member")

    g = sns.kdeplot(abundance_df[abundance_df.member.isin(members)],
                hue="member",
                x="relative_abundance",
                cut=0,
                ax=ax,
                common_norm=False,
                palette=palette)
    ax.set_ylim([0,5])
    
    if title:
        ax.set_title(f"{title}: {co_or_ind_str}")
    else:
        ax.set_title(f"{community}: {co_or_ind_str}")
    ax.set_xlabel("Relative abundance")
    return ax
    

def plot_co_ind_growers_venturelli(abundances_df,co_growers):
    """
    Function for plotting histograms and density plots for number of members in each sample and the distribution of the relative abundance across the samples.  
    """
    all_members = set(abundances_df.member.unique())
    non_co_growers = all_members-set(co_growers)
    
    communities = list(abundances_df.community.sort_values().unique())
    communities.remove("NONE")

    colors_members =["#2e3ef3", "#8f6039", "#86a95d", "#df1645", "#24a5f7", "#f76015", "#7c28a4", "#f250a3", "#096013", "#cb66ea", "#e6d0c2", "#ffd65b"]
    palette = OrderedDict(zip(all_members,colors_members))
    
    fig,axs = plt.subplots(13,4,figsize=(15,40))
    
    for community,(ax1,ax2,ax3,ax4) in zip(communities,axs[1:]):
    
        abundance_comm = abundances_df[abundances_df.community==community].copy()

        plot_hist_helper(abundance_comm,community,co_growers,"Co-growers",palette,ax1)
        plot_hist_helper(abundance_comm,community,non_co_growers,"Independent growers",palette,ax3)

        plot_density_helper(abundance_comm,community,co_growers,"Co-growers",palette,ax2)
        plot_density_helper(abundance_comm,community,non_co_growers,"Independent growers",palette,ax4)
        
    ax1,ax2,ax3,ax4 = axs[0]
    plot_hist_helper(abundance_comm,community,co_growers,"Co-growers",palette,ax1,title="Full")
    plot_hist_helper(abundance_comm,community,non_co_growers,"Independent growers",palette,ax3,title="Full")

    plot_density_helper(abundance_comm,community,co_growers,"Co-growers",palette,ax2,title="Full")
    plot_density_helper(abundance_comm,community,non_co_growers,"Independent growers",palette,ax4,title="Full")

    plt.tight_layout()
    
    return fig,axs


def plot_co_ind_growers_guillen(abundances_df):
    
    all_members = set(abundances_df.member.unique())
    colors_members =["#2e3ef3", "#8f6039", "#86a95d", "#df1645", "#24a5f7", "#f76015", "#7c28a4", "#f250a3", "#096013", "#cb66ea", "#e6d0c2", "#ffd65b"]
    palette = OrderedDict(zip(all_members,colors_members))
    
    
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
    
    abundance_comm = abundances_df.copy()
    
    plot_hist_helper(abundance_comm,"E.coli auxotrophs",all_members,"All members",palette,ax1)

    plot_density_helper(abundance_comm,"E.coli auxotrophs",all_members,"All members",palette,ax2)
    plt.tight_layout()
  







