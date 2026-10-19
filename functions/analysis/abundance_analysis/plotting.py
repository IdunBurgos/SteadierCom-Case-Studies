import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kendalltau

def plot_scatter_steadiercom(joined_spec,ax,title,member_colors,legend=None):
    joined_mean = joined_spec.groupby(["community","member"]).agg({'relative_abundance_sim':'mean','relative_abundance_exp':'mean'})\
                            .rename(columns={'relative_abundance_sim':'SteadierCom relative abundance','relative_abundance_exp':'Experimental relative abundance'})

    stat = kendalltau(joined_mean["Experimental relative abundance"],joined_mean["SteadierCom relative abundance"])
    g = sns.lineplot(joined_spec,
                x="relative_abundance_exp",
                y="relative_abundance_sim",
                hue="member",
                marker="o",
                markersize=10,
                linestyle="",
                err_style="bars",
                errorbar="ci",
                palette=member_colors,
                legend=legend,
                ax=ax)
    
    plt.text(0.1, 0.8, f"Kendalls tau={stat.correlation:.2f} \n p-value={stat.pvalue:.2e}", transform=ax.transAxes)
    ax.set_title(title)
    return ax

def plot_scatter_steadycom(joined_spec,ax,title,member_colors,legend=None):
    
    stat = kendalltau(joined_spec["relative_abundance_exp"],joined_spec["relative_abundance_sim"])
    g = sns.scatterplot(joined_spec,
                x="relative_abundance_exp", 
                y="relative_abundance_sim",
                hue="member",
                marker="o",
                s=80,
                linestyle="",
                palette=member_colors,
                legend=legend,
                ax=ax)
    
    plt.text(0.1, 0.8, f"Kendalls tau={stat.correlation:.2f} \n p-value={stat.pvalue:.2e}", transform=ax.transAxes)
    ax.set_title(title)
    return ax
