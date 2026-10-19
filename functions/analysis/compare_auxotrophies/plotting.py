import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D


def plot_triplet_violin(merged,metclass2color_dict,threshold=1e-6):

    threshold_main_classes = 0.15 # Min value of highest reported interaction score for compound

    main_metclasses_steadier = (
        merged.groupby("super_class")["frequency"]
        .max()
        .loc[lambda s: s > threshold_main_classes]
        .index
        .to_list()
    )

   # Find relevant met classes for smetana
    main_metclasse_smetana = (
        merged.groupby("super_class")["smetana"]
        .max()
        .loc[lambda s: s > threshold_main_classes]
        .index
        .to_list()
    )
   # Combine sets
    super_classes=list(set(main_metclasses_steadier+main_metclasse_smetana))

    # Get results for relevant classes only
    merged_copy = (
        merged[merged.super_class.isin(super_classes) ]
        ).copy()

    # Filter steadiercom and smetana results based on threshold
    merged_steadier = merged_copy.query("frequency > @threshold")
    merged_smetana = (
        merged_copy.query("smetana > @threshold")
        )

    # For each dataset -> count number of  in each categroy
    class_count = pd.DataFrame(
        {"freq_count": (
            merged_steadier.groupby("super_class")["frequency"]
            .count()
            ),
        "mus_count":(
            merged_smetana.groupby("super_class")["smetana"]
            .count()
            ) 
            }
    ).fillna(0).sort_values("freq_count",ascending=False)

    class_unique_count = pd.DataFrame(
        {"freq_count": (
            merged_steadier.groupby("super_class")["compound"]
            .nunique()
            ),
        "mus_count":(
            merged_smetana.groupby("super_class")["compound"]
            .nunique()
            ) 
            }
    ).fillna(0).sort_values("freq_count",ascending=False)


    hue_order = class_count.index
    
    fig = plt.figure(figsize =([16, 8]))
    gs = gridspec.GridSpec(2, 6)
    ax1a = plt.subplot(gs[:1, :2])
    ax1b = plt.subplot(gs[1, :2])
    ax2 = plt.subplot(gs[0, 2:])
    ax3 = plt.subplot(gs[1, 2:])
    alpha=0.6

    ax1a = class_count.T.plot.bar(color=metclass2color_dict,
                ax=ax1a,
                legend=None,
                alpha=alpha
                )
    ax1a.set_ylabel("Counts")
    ax1a.set_xticklabels(["frequency","SMETANA"])


    ax1b = class_unique_count.T.plot.bar(color=metclass2color_dict,
                ax=ax1b,
                legend=None,
                alpha=alpha
                )
    ax1b.set_ylabel("Unique Counts")
    ax1b.set_xticklabels(["frequency","SMETANA"])

    axs = [ax2,ax3]
    merged_list = [merged_steadier,merged_smetana]
    metric_list = ["frequency","smetana"]
    legend_list = [True, None]

    for ax,merged_spec,metric,legend in zip(axs,merged_list,metric_list,legend_list):
        # stripplot without edge color
        ax = sns.stripplot(merged_spec[~merged_spec.expected_auxotrophy],
                    y=metric,
                    ax=ax,
                    x="super_class",
                    hue="super_class",
                    order=hue_order,
                    linewidth=0.,
                    legend=legend,
                    s=6,
                    alpha=0.8,
                    zorder=2,
                    palette=metclass2color_dict,
                    )
        
        handles, labels = ax.get_legend_handles_labels()
        
        # stripplot with edge color
        ax = sns.stripplot(
                    data=merged_spec[merged_spec.expected_auxotrophy],
                    x="super_class",
                    y=metric,
                    hue="super_class",
                    palette=metclass2color_dict,
                    linewidth=0.8,
                    edgecolor="red",
                    zorder=2,
                    ax=ax,
                    legend=False,
                    )
    
        ax = sns.boxplot(merged_spec,
                    y=metric,
                    ax=ax,
                    x="super_class",
                    hue="super_class",
                    order=hue_order,
                    linewidth=0.5,
                    legend=False,
                    zorder=1,
                    palette='dark:black',
                    fill=False,
                    fliersize=0,
                    medianprops={"color": "r", "linewidth": 0.5},
                    width=0.7,
                    )
        ax.set_ylim([0,1])


        if legend:

            color_legend = ax.legend(
                handles=handles,
                labels=labels,
                title="Super class",
                loc="upper left",
                frameon=True,
                bbox_to_anchor=(1.02, 0.9),
            )

            edge_handle = Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                markerfacecolor="white",
                markeredgecolor="red",
                markeredgewidth=0.9,
                markersize=7,
                label="Matches auxotrophy = True"
            )

            edge_legend = ax.legend(
                handles=[edge_handle],
                title="Edge color",
                loc="upper left",
                frameon=True,
                bbox_to_anchor=(1.02, -0.5),
            )

            ax.add_artist(color_legend)
            ax.add_artist(edge_legend)
        
    ax2.set_title("SteadierCom - frequency")
    ax2.set_xlabel("")
    ax2.set_xticklabels([])
    
    ax3.set_title("SMETANA")
    ax3.set_xlabel("Metabolite class")
    ax3.set_ylabel("SMETANA")
    ax3.tick_params(axis='x', rotation=90)
    ax3.set_xticklabels([])


    fig.tight_layout(rect=[0, 0, 0.60, 1])
    return fig,(ax1a,ax1b,ax2,ax3)