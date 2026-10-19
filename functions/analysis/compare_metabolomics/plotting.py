from itertools import product
import pandas as pd
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns


cmap_dict_overlap = {"met not measured":"#ffe599","not supported":"#b05ab4","supported":"#35af00"}


def plot_triplets_scatter(df_merge):

    # prepare figure with 3 subplots
    fig = plt.figure(constrained_layout=True,figsize=(8,6))
    spec = gridspec.GridSpec(ncols=4, nrows=4, figure=fig)

    ax1 = fig.add_subplot(spec[1:, 3])
    ax2 = fig.add_subplot(spec[0, :3])
    ax3 = fig.add_subplot(spec[1:, :3])
    ax_legend = fig.add_subplot(spec[0,3])

    # gather data for plotting
    plot_df = df_merge.sort_values("experimental_support")
    hue_order = plot_df.experimental_support.sort_values(ascending=False).unique()


    # Main figure: scatter plot with error bars
    ax3 = sns.scatterplot(
        plot_df,
        x="smetana_median",#"smetana",
        y="steadiercom_median",#"frequency",
        hue="experimental_support", #"overlap",
        s=50,
        palette=cmap_dict_overlap,
        hue_order=hue_order,
        alpha=0.7,
        legend=True,
        ax=ax3,
        zorder=1,
    )

    plot_df_reduced = plot_df.query("experimental_support != 'met not measured'").copy()
    y_errors = [plot_df_reduced["steadiercom_median"] - plot_df_reduced["steadiercom_q25"], plot_df_reduced["steadiercom_q75"] - plot_df_reduced["steadiercom_median"]]
    x_errors = [plot_df_reduced["smetana_median"] - plot_df_reduced["smetana_q25"], plot_df_reduced["smetana_q75"] - plot_df_reduced["smetana_median"]]

    ax3.errorbar(
        x=plot_df_reduced["smetana_median"], 
        y=plot_df_reduced["steadiercom_median"], 
        yerr=y_errors, 
        xerr=x_errors,
        fmt="none",
        color="gray",
        alpha=0.3,
        capsize=2,
        zorder=2,
    )

    ax3.set_ylabel("steadiercom: frequency")
    ax3.set_xlabel("smetana")


    # Histograms for steadiercom and smetana
    ax1 = sns.histplot(plot_df, 
                    y="steadiercom_median",
                        hue="experimental_support",
                        hue_order=hue_order,
                        palette=cmap_dict_overlap,
                        ax=ax1, 
                        bins=25, 
                        alpha=0.4,
                        legend=False)
    ax1.set_xscale("log")
    ax1.set_ylabel("")

    ax2 = sns.histplot(plot_df,
                    x="smetana_median",
                    hue="experimental_support",
                    hue_order=hue_order,
                    palette=cmap_dict_overlap,
                    ax=ax2,
                    bins=25,
                    alpha=0.4,
                    legend=False,
                    )
    ax2.set_yscale("log")
    ax2.set_xlabel("")

    handles, labels = ax3.get_legend_handles_labels()

    ax_legend.legend(handles=handles, loc='center')
    ax_legend.set_axis_off()
    if ax3.get_legend() is not None:
        ax3.get_legend().remove()

    return fig,(ax1,ax2,ax3,ax_legend)


def plot_triplets_violin(merged,metclass2color_dict,threshold=1e-6):

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
    
    assert "experimental_support" in merged_steadier.columns, f"experimental_support not in merged_steadier: {merged_steadier.columns}"
    assert "experimental_support" in merged_smetana.columns, f"experimental_support not in merged_smetana: {merged_smetana.columns}"


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
    
    fig = plt.figure(figsize =([12, 8]))
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
        merged_spec = merged_spec.sort_values("experimental_support")
        ax = sns.stripplot(merged_spec,
                    y=metric,
                    ax=ax,
                    x="super_class",
                    hue="experimental_support",
                    order=hue_order,
                    linewidth=0.,
                    legend=legend,
                    s=6,
                    alpha=0.5,
                    zorder=2,
                    palette=cmap_dict_overlap,
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
    ax2.set_title("SteadierCom - frequency")
    ax2.set_xlabel("")
    ax2.set_xticklabels([])
    ax2.legend(loc='upper left', bbox_to_anchor=(1.0, 1),title="Metabolite class")  
    
    ax3.set_title("SMETANA")
    ax3.set_xlabel("Metabolite class")
    ax3.set_ylabel("SMETANA")
    ax3.tick_params(axis='x', rotation=90)
    plt.tight_layout()

    return fig,(ax1a,ax1b,ax2,ax3)

