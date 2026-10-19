
import numpy as np
import pandas as pd

from matplotlib.patches import Rectangle
from matplotlib.patches import Patch

import matplotlib.pyplot as plt
from pycirclize import Circos

from .steadiercom_samples_processing import data_ReceiverOrDonor,data_for_links 

import logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


metclass2color_dict = {"carbohydrate derivatives":"#DEFD3C",  # carbohydrates and derivatives   e5fa79
            "amino acids and derivatives":"#ddaa7b",  # amino acids
            "nucleotides and derivatives":"#fbdc5d",  # nucleosides
            "oligosaccharides":"#8BBDAA",  # oligosaccharides #"#a7dfb6"
            "oligopeptides":"#FFD19B",  # oligopeptides 
            "simple sugars":"#9BF276",  # simple sugars   bfffa4
            "urea and urea derivatives":"#FFFFBA",  # urea 
            "gases":"#b1fbf6",  # gasses
            "B-vitamins":"#6F7FB0",  # B-vitamins #adb0e8
            "fatty acids":"#ffd1fd",   # fatty acids
            "carboxylic acids and anions":"#EB857E",  # carboxylic acids FFB3BA
            "alcohols and aldehydes":"#C3BBF3",  # alcohols and aldehydes 
            "lipids":"#8FCFFF", # phospholipids    
            "cofactors":"#FF85A1",
            "other":"#eaeaea"
            }


def validate_simulated_data(steadiercom_samples_preprocessed,metric):
    # Check that the columns that are referenced in steadiercom_samples_preprocessed exist
    
    required_columns = ["donor", "receiver", "super_class","frequency", metric]
    for col in required_columns:
        if col not in steadiercom_samples_preprocessed.columns:
            raise ValueError(f"Column {col} not found in steadiercom_samples_preprocessed columns.")


    super_classes = steadiercom_samples_preprocessed.super_class.unique()
    for class_ in super_classes:
        if class_ not in metclass2color_dict.keys():
            raise ValueError(f"Super class {class_} not found in metabolite class color dictionary.")


def validate_relative_abundance_data(experimental_abundance):
    # Check that the columns that are referenced in experimental_abundance exist
    required_columns_abundance = ["community", "member", "relative_abundance"]
    for col in required_columns_abundance:
        if col not in experimental_abundance.columns:
            raise ValueError(f"Column {col} not found in experimental_abundance columns.")
        
def validate_data(steadiercom_samples_preprocessed, experimental_abundance, metric):

    validate_simulated_data(steadiercom_samples_preprocessed, metric)
    validate_relative_abundance_data(experimental_abundance)

    # Check that the members in experimental_abundance and steadiercom_samples_preprocessed match
    exp_members = set(experimental_abundance.member.unique())
    sim_members = set(list(steadiercom_samples_preprocessed.donor.unique())+list(steadiercom_samples_preprocessed.receiver.unique())) - {"environment"}

    if not sim_members.issubset(exp_members):
        raise ValueError(f"Members in steadiercom_samples_preprocessed is not a subset of members in experimental_abundance. Members in experimental_abundance: {exp_members}. Members in steadiercom_samples_preprocessed: {sim_members}.")
    
    if metric!="mass_rate*frequency":
        logging.warning(f"{metric} is not yet supported. Links in the circos plot are scaled assuming that the metric is based on mass rate.")
            

def circos_plot_process_data(steadiercom_samples_preprocessed,experimental_abundance,metric="mass_rate*frequency"):
    """
    This function processes the data for the circos plot and creates a list of links.

    Input:
    steadiercom_samples_preprocessed: Simulation results DataFrame (SteadierCom assumed) with additional columns for analysis: super_class (for compound grouping), and a mass rate metric
    experimental_abundance: Experimental abundance DataFrame with the following columns: (community, member, relative_abundance). 
    metric: Metric to use for analysis which should be based on mass rate for links to be scaled correctly. (default: "mass_rate*frequency")

    Output:
    links: List of links for the circos plot
    members: List of members (sorted according to relative abundance)

    """
    
    # sort the members according to relative abundance
    members = experimental_abundance.sort_values("relative_abundance",ascending=False).member.values

    # remove all data related to environment
    steadiercom_crossfeeding = steadiercom_samples_preprocessed.query("donor != 'environment' & receiver != 'environment'")
   
    # Aggregate results by donor, receiver, and super_class
    steadiercom_crossfeeding = steadiercom_crossfeeding.groupby(["donor","receiver","super_class"]).sum()[metric].reset_index()
    steadiercom_crossfeeding.donor = steadiercom_crossfeeding.donor.astype("category")
    steadiercom_crossfeeding.donor = steadiercom_crossfeeding.donor.cat.set_categories(members)
    steadiercom_crossfeeding.receiver = steadiercom_crossfeeding.receiver.astype("category")
    steadiercom_crossfeeding.receiver = steadiercom_crossfeeding.receiver.cat.set_categories(members)
    steadiercom_crossfeeding = steadiercom_crossfeeding.sort_values(["donor","receiver"],ascending=[False,True]).copy()

    links = data_for_links(steadiercom_crossfeeding,members,metric=metric)
    return links,members


def data_uptake_prod(steadiercom_samples_preprocessed,experimental_abundance,community_id=False,only_media=False,metric="mass_rate*frequency"):
    """
    Prepare data for uptake/production plotting.

    Parameters
    ----------
    steadiercom_samples_preprocessed: pd.DataFrame
        Preprocessed SteadierCom samples data.
    experimental_abundance: pd.DataFrame
        Experimental abundance data.
    community_id: bool
        Whether to filter by community ID.
    only_media: bool
        Whether to include only media data.
    metric: str
        The metric to use for analysis.

    Returns
    ---------
    data_receiver_df: pd.DataFrame
        Data for the receivers.
    data_donor_df: pd.DataFrame
        Data for the donors.
    community_abundance: pd.DataFrame
        Abundance data for the specified community.
    """

    if community_id:
        steadiercom_samples_preprocessed_copy = steadiercom_samples_preprocessed.query("community == @community_id")
        community_abundance = experimental_abundance.query("community == @community_id").copy()
    else:
        steadiercom_samples_preprocessed_copy = steadiercom_samples_preprocessed.copy()
        community_abundance = experimental_abundance.copy()

    # DATA CONSUMERS
    data_receiver_df = data_ReceiverOrDonor(steadiercom_samples_preprocessed_copy,"receiver",only_media=only_media,metric=metric)

    # DATA PRODUCERS
    data_donor_df = data_ReceiverOrDonor(steadiercom_samples_preprocessed_copy,"donor",only_media=only_media,metric=metric)

    # DATA FOR ABUNDANCE "PLOT
    community_abundance = community_abundance.sort_values(by="relative_abundance",ascending=False).set_index("member")

    ## Sort dataframes
    data_receiver_df = data_receiver_df.reindex(columns=community_abundance.index, fill_value=0).copy()
    data_receiver_df = data_receiver_df.sort_index().copy()

    data_donor_df = data_donor_df.reindex(columns=community_abundance.index, fill_value=0).copy()
    data_donor_df = data_donor_df.sort_index().copy()
    
    return data_receiver_df,data_donor_df,community_abundance


def plot_circos_plot(links,members,member2color_dict,fontsize=15,move_legend=0):
    """
    Function for plotting circos plots.
    
    Parameters
    ----------
    links: List of links between members as given from running 'circos_plot_process_data'
    members: List of members in the community
    member2color_dict: Dictionary mapping members to colors
    fontsize: Font size for the plot
    move_legend: Amount to move the legend

    Returns
    -------
    fig: The figure object for the plot
    line_legend: The legend for the metabolite classes
    member_legend: The legend for the members
    """
    
    ### Create Circos Plot ###

    sectors = {member:1 for member in members}
    circos = Circos(sectors,space=5)

    i = 0
    members_specific = []
    
    for sector in circos.sectors:
        
        # Add color for group
        track = sector.add_track((95, 100))
        
        track.axis(fc=member2color_dict[sector.name])
        members_specific.append(sector.name)

        i+=1
    classes = []

    # Create links in the circos plot from data
    for link in links:
        source = link[0]
        target = link[1]
        class_ = link[2] # Class of compound

        if class_ not in metclass2color_dict.keys():
            raise ValueError(f"Class {class_} not found in metabolite class color dictionary.")
        
        color = metclass2color_dict[class_]

        circos.link(source,target,direction=1,color=color,ec=color,arrow_length_ratio=0.15,allow_twist=False)
        classes.append(class_)

    fig = circos.plotfig()

    ### Create Metabolite Legend ###
    classes = sorted(list(set(classes)))
    
    line_handles = [Patch(color=metclass2color_dict[class_], label=class_,alpha=0.8) for class_ in classes]
    line_legend = circos.ax.legend(
        handles=line_handles,
        bbox_to_anchor=(1.0, 1.2+move_legend),  # Positioning to the right of the plot
        loc='upper left', 
        fontsize=fontsize,
        title="Compound",
        handlelength=2,
    )
    
    line_legend.get_title().set_fontsize('15')
    line_legend._legend_box.align = "left"
    
    ### Create Member Legend ###
    member_handles = [Patch(color=color, label=member) for member, color in member2color_dict.items() if member in members_specific]
    member_legend = circos.ax.legend(
        handles=member_handles,
        bbox_to_anchor=(1,0.5+move_legend),  # Positioning to the lower right
        loc='upper left',
        fontsize=fontsize,
        title="Member",
        handlelength=2,
    )
    
    member_legend.get_title().set_fontsize('15')
    member_legend._legend_box.align = "left"
    

    ### Add Legends and Plot ###
    circos.ax.add_artist(line_legend)
    circos.ax.add_artist(member_legend)

    #plt.tight_layout(rect=[0, 0, 0.8, 0.8])
    plt.close()

    return fig,line_legend,member_legend


def plot_uptake_prod(data_receiver_df,data_donor_df,community_abundance,member2color_dict,title=None,legend=True):
    
    fig, axs = plt.subplots(4,gridspec_kw={'height_ratios': [1,3,0.5, 3]},figsize=(20,15),sharex=True)

    # Relative abundance
    community_abundance = community_abundance.reset_index()

    community_abundance.plot(ax=axs[0],kind="line",x="member",y="relative_abundance",ylabel="",xlabel="",color="red",style='o-',markersize=15,lw=4)
    axs[0].set_xticks([])
    axs[0].tick_params(labelsize=22)
    axs[0].legend('')
    axs[0].set_title("Relative \n abundance",x=-.1,y=0.0,fontsize=28,rotation=90)


    # Producing bars
    data_donor_df = data_donor_df.loc[data_donor_df.sum(axis=1)>1e-20,:].copy()
    data_donor_df.transpose().plot(kind="bar",ax=axs[1],color = metclass2color_dict,legend=None,stacked=True,alpha=0.8)
    axs[1].set_xticks([])
    axs[1].tick_params(labelsize=22)
    axs[1].set_xlabel('')
    axs[1].set_title("Production \n [mg/(h*g$_{bio}$)]",x=-.1,y=0.2,fontsize=28,rotation=90)

    # Member colors
    xticks_labels = data_receiver_df.columns
    xticks = np.arange(0,len(xticks_labels))
    
    axs[2].set_yticks([])
    axs[2].set_xticks(xticks)

    members_specific = []

    for i,x in enumerate(xticks_labels):
        patch=Rectangle((i-0.5, 0),1, 1, color=member2color_dict[x])
        members_specific.append(x)
        axs[2].add_patch(patch)
    axs[2].set_xticks([])

    plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.1, 
                        hspace=0.05)
                        
    axs[2].set_xticks([])

    # Consuming bars
    data_receiver_df = data_receiver_df.loc[data_receiver_df.sum(axis=1)>1e-20,:].copy() # TO-DO: Figure out why this is necessary
    data_receiver_df.transpose().apply(lambda x: (-1)*x).plot(kind="bar",ax=axs[3],color=metclass2color_dict,legend=None,stacked=True,alpha=0.8)
    axs[3].set_xticks([])
    axs[3].tick_params(labelsize=22)
    axs[3].set_xlabel('')
    axs[3].set_title("Consumption \n [mg/(h*g$_{bio}$)]",x=-.1,y=0.2,fontsize=28,rotation=90)

    if legend:    
        # Custom legend for compounds
        custom_lines = [Patch(color=metclass2color_dict[super_class], label=super_class,alpha=0.8) for super_class in  pd.concat([data_donor_df,data_receiver_df]).index.unique()]
        fig.legend(custom_lines,pd.concat([data_donor_df,data_receiver_df]).index.unique(), fontsize=15, loc='upper left', bbox_to_anchor=(0.9, 0.8), title="COMPOUNDS")

        members_labels = [x for x in xticks_labels if x in member2color_dict.keys()]
        members_handles = [Patch(color=member2color_dict[member], label=member) for member in members_labels]
        members_legend = fig.legend(members_handles, members_labels, fontsize=15, loc='upper left', bbox_to_anchor=(0.9, 0.5), title="MEMBERS")

        fig.suptitle(title,fontsize=30)
        plt.close()
    else:
        
        extra_compounds = list(pd.concat([data_donor_df,data_receiver_df]).index.unique())
        plt.close()
        return fig,extra_compounds
    return fig



def combined_figure(steadier_sample,experimental_abundance,member2color_dict,community_id,path = "Figures/circos_plots_2_1/",move_legend=0):

    data_receiver_df,data_donor_df,community_abundance = data_uptake_prod(steadier_sample,
                                                                          experimental_abundance,
                                                                          community_id=community_id,
                                                                          only_media=False)
    
    fig,extra_compounds = plot_uptake_prod(data_receiver_df,data_donor_df,community_abundance,member2color_dict,legend=False)
    fig.savefig(path+"bar_"+community_id+".png",bbox_inches='tight',dpi=300)    
    
    links,members = circos_plot_process_data(steadier_sample,experimental_abundance)
    fig,line_legend,genus_legend = plot_circos_plot(links,members,member2color_dict,move_legend=move_legend)
    fig.savefig(path+"circos_"+community_id+".png",bbox_extra_artists=(line_legend,genus_legend),bbox_inches='tight',dpi=300)
    
    # create figure 
    fig = plt.figure(figsize=(20, 10)) 

    # setting values to rows and column variables 
    rows = 1
    columns = 2

    # reading images 
    Image1 = plt.imread(path+"bar_"+community_id+".png") 


    Image2 = plt.imread(path+"circos_"+community_id+".png")

    # Adds a subplot at the 1st position 
    fig.add_subplot(rows, columns, 1) 

    # showing image 
    plt.imshow(Image1) 
    plt.axis('off') 
    plt.title("A",loc="left",fontdict={"fontsize": 23, "weight": "bold"}) 

    # Adds a subplot at the 2nd position 
    fig.add_subplot(rows, columns, 2) 

    # showing image 
    plt.imshow(Image2) 
    plt.axis('off') 
    plt.title("B",loc="left",fontdict={"fontsize": 23, "weight": "bold"})

    plt.tight_layout()
    plt.savefig(path+"combined_"+community_id+".png",dpi=300)
    plt.show()

