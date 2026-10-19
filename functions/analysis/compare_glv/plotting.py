import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib as mpl
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec
from matplotlib.collections import LineCollection

def plot_interactions(df_merged,signal_exp,title=""):

    # Normalize colors so that 0 is the center
    norm = TwoSlopeNorm(
        vmin=-3,
        vcenter=0,
        vmax=3
    )
    palette=mpl.colormaps["coolwarm"]

    fig = plt.figure(figsize=(10,8))

    gs = GridSpec(2, 2, width_ratios=[1, 2], height_ratios=[2, 1],)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])
    ax4 = fig.add_subplot(gs[3])

    sns.scatterplot(ax=ax2,
                    data=df_merged,
                    x="simulated_coop",
                    y="simulated_comp",
                    hue="experimental",
                    palette=palette,
                    hue_norm=norm,
                    s=30
    )

    sns.move_legend(ax2, "upper left", bbox_to_anchor=(1.0, 1.0),title="gLV coefficient")
    sns.violinplot(df_merged,x="exp_categories",y="simulated_comp",hue="exp_categories",ax=ax1,palette="coolwarm",cut=0,alpha=0.8,linecolor="gray",inner_kws=dict(box_width=8, whis_width=1),zorder=0)
    sns.stripplot(df_merged,x="exp_categories",y="simulated_comp",hue="exp_categories",palette="coolwarm",ax=ax1,alpha=0.5,linewidth=0.5,edgecolor="gray",size=4,jitter=0.2,zorder=1)

    sns.violinplot(df_merged,y="exp_categories",x="simulated_coop",hue="exp_categories",ax=ax4,palette="coolwarm",cut=0,alpha=0.8,linecolor="gray",inner_kws=dict(box_width=8, whis_width=1),legend=True,zorder=0)
    sns.stripplot(df_merged,y="exp_categories",x="simulated_coop",hue="exp_categories",palette="coolwarm",ax=ax4,alpha=0.5,linewidth=0.5,edgecolor="gray",size=4,jitter=0.2,zorder=1)
    sns.move_legend(ax4, "upper left", bbox_to_anchor=(1.0, 1.0),title="gLV coefficient (categorical)")


    sns.kdeplot(signal_exp,
                ax=ax3
    )
    x, y = ax3.get_lines()[0].get_data()
    segments = np.array([x[:-1], y[:-1], x[1:], y[1:]]).T.reshape(-1, 2, 2)
    lc = LineCollection(segments, cmap='coolwarm', norm=norm)
    lc.set_array(x[:-1])
    lc.set_linewidth(2)
    ax3.get_lines()[0].remove()
    line = ax3.add_collection(lc)
    ax3.set_xlabel("gLV coefficient")

    
    ax2.set_xlabel("")
    ax2.set_ylabel("")

    ax1.set_xlabel("")
    ax4.set_ylabel("")

    fig.suptitle(title, fontsize=14)
    
    plt.show()

def plot_interactions_median(df_merged,title=""):

    # Column used for the color gradient
    values = df_merged["experimental"]

    # Normalize colors so that 0 is the center
    norm = TwoSlopeNorm(
        vmin=-3,
        vcenter=0,
        vmax=3
    )
    palette=mpl.colormaps["coolwarm"]

    fig = plt.figure(figsize=(11,8))

    gs = GridSpec(2, 2, width_ratios=[1, 2], height_ratios=[2, 1],)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])
    ax4 = fig.add_subplot(gs[3])

    sns.scatterplot(data=df_merged, 
                    x="coop_median", 
                    y="comp_median", 
                    hue="experimental", 
                    palette=palette,
                    hue_norm=norm,
                    s=30,
                    alpha=1,
                    ax=ax2)
    y_errors = [df_merged["comp_median"] - df_merged["comp_q25"], df_merged["comp_q75"] - df_merged["comp_median"]]
    x_errors = [df_merged["coop_median"] - df_merged["coop_q25"], df_merged["coop_q75"] - df_merged["coop_median"]]
    ax2.errorbar(
        x=df_merged["coop_median"], 
        y=df_merged["comp_median"], 
        yerr=y_errors, 
        xerr=x_errors,
        fmt="none",
        color="gray",
        alpha=0.3,
    )

    sns.move_legend(ax2, "upper left", bbox_to_anchor=(1.0, 1.0),title="gLV coefficient")

    for ax,comp_coop in zip([ax1, ax4], ["comp_median", "coop_median"]):

        if comp_coop=="comp_median":
            x,y ="exp_categories",comp_coop
        else:
            x,y =comp_coop,"exp_categories"

        sns.violinplot(df_merged,x=x,y=y,hue="exp_categories",ax=ax,palette="coolwarm",cut=0,alpha=0.5,linecolor="gray",inner_kws=dict(box_width=8, whis_width=1,alpha=0.8),zorder=0,legend=(comp_coop=="coop_median"))
        sns.stripplot(df_merged,x=x,y=y,hue="exp_categories",palette="coolwarm",ax=ax,alpha=0.5,linewidth=0.5,edgecolor="gray",size=4,jitter=0.2,zorder=1)

    sns.move_legend(ax4, "upper left", bbox_to_anchor=(1.0, 1.0),title="gLV coefficient (categorical)")


    sns.kdeplot(df_merged,
                x="experimental",
                ax=ax3
    )
    x, y = ax3.get_lines()[0].get_data()
    segments = np.array([x[:-1], y[:-1], x[1:], y[1:]]).T.reshape(-1, 2, 2)
    lc = LineCollection(segments, cmap='coolwarm', norm=norm)
    lc.set_array(x[:-1])
    lc.set_linewidth(2)
    ax3.get_lines()[0].remove()
    line = ax3.add_collection(lc)
    ax3.set_xlabel("gLV coefficient")
    
    ax2.set_xlabel("")
    ax2.set_ylabel("")

    ax1.set_xlabel("")
    ax4.set_ylabel("")

    fig.suptitle(title, fontsize=14)

    plt.show()