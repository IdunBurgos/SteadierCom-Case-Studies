import pandas as pd
from os import listdir
import re

def _change_name(x):
    names =x.split("_")
    name = names[0][0].upper()+names[1][0].upper()
    return name


def load_experimental_results():

    exp_abundance = []

    for file in listdir("../data/Venturelli2018/generated/communities_steadiercom/"):
        if file.endswith(".tsv"):
            
            exp_abund = pd.read_csv("../data/Venturelli2018/generated/communities_steadiercom/"+file,sep="\t",header=None)
            name = file.split("_")[0]
            exp_abund.name=f"{name}"
            exp_abundance.append(exp_abund)

    exp_abundance_df = pd.concat(exp_abundance)
    exp_abundance_df.columns = ["community","member","abundance"]
    exp_abundance_df = exp_abundance_df.set_index(["community","member"])["abundance"].unstack().fillna(0)
    exp_abundance_df.columns = exp_abundance_df.columns.map(_change_name)

    exp_abundance_df.sort_index(axis=1,inplace=True)
    return exp_abundance_df

def load_simulation_results():

    abundances = {}
    for file in listdir("../output/Venturelli2018/results_steadiercom/results_growth/"):
        if file.endswith("_abundance.tsv"):

            abundances[file.replace("_abundance.tsv","")] = pd.read_csv(f"../output/Venturelli2018/results_steadiercom/results_growth/{file}",sep="\t")

    abundances_df = pd.concat(abundances)
    abundances_df.reset_index(drop=True)
    abundances_df.member = abundances_df.member.map(_change_name)

    samples = abundances_df["sample"].max()
    steadiercom_abundance_df = abundances_df.groupby(['community', 'member'], as_index=False).agg({"relative_abundance": (lambda x: sum(x)/samples)})
    steadiercom_abundance_df = steadiercom_abundance_df.pivot_table(index='member', columns='community', values='relative_abundance', fill_value=0)
    steadiercom_abundance_df = steadiercom_abundance_df.T.sort_index(axis=1)

    # Filter out communities where there is only one present member
    count_members_sample = abundances_df[["sample","community"]].value_counts()
    count_above_1 = count_members_sample[count_members_sample>1].index
    abundances_df_filtered = abundances_df[abundances_df[["sample", "community"]].apply(tuple, axis=1).isin(count_above_1)].copy()
    # Count number of samples per community
    samples_in_community = (abundances_df_filtered[["community", "sample"]]
                            .drop_duplicates()
                            .groupby("community")
                            .size())
    # Calculate average relative abundances using precomputed sample counts per community
    steadiercom_abundance_filtered_df = (abundances_df_filtered
                                         .groupby(['community', 'member'], as_index=False)['relative_abundance']
                                         .sum())
    
    # Normalize by number of samples per community
    steadiercom_abundance_filtered_df["relative_abundance"] = (
        steadiercom_abundance_filtered_df["relative_abundance"] /
        steadiercom_abundance_filtered_df["community"].map(samples_in_community)
    )
    steadiercom_abundance_filtered_df = steadiercom_abundance_filtered_df.pivot_table(index='member', columns='community', values='relative_abundance', fill_value=0)
    steadiercom_abundance_filtered_df = steadiercom_abundance_filtered_df.T.sort_index(axis=1)

    steadycom_abundance_df = pd.read_csv("../output/Venturelli2018/steadycom/no_abundance/relative_abundance.tsv",sep="\t",index_col=0)
    steadycom_abundance_df.index = steadycom_abundance_df.index.map(lambda x:x.split("_")[1][0].upper()+x.split("_")[2][0].upper())
    steadycom_abundance_df = steadycom_abundance_df.T
    steadycom_abundance_df.sort_index(inplace=True)
    steadycom_abundance_df.sort_index(inplace=True,axis=1)
    return steadiercom_abundance_df, steadiercom_abundance_filtered_df, steadycom_abundance_df


def load_experimental_results_pairwise():

    exp_abundance_df = pd.read_csv("../data/Venturelli2018/generated/communities_pairwise.tsv",sep="\t",header=None)

    exp_abundance_df.columns = ["community","member","abundance"]
    exp_abundance_df = exp_abundance_df.set_index(["community","member"])["abundance"].unstack()
    exp_abundance_df.columns = exp_abundance_df.columns.map(_change_name)
    exp_abundance_df

    exp_abundance_df.sort_index(axis=1,inplace=True)
    exp_abundance_ser = exp_abundance_df.stack()
    exp_abundance_ser.name="exp"
    return exp_abundance_ser


def load_simulation_results_pairwise():
    ## STEADIERCOM - without filtering
    abundances_df = pd.read_csv("../output/Venturelli2018/results_steadiercom/results_pairwise_growth/_abundance.tsv",sep="\t")
    abundances_df.member = abundances_df.member.map(_change_name)
    steadiercom_abundance = abundances_df.groupby(['community', 'member']).sum()["relative_abundance"]/(abundances_df["sample"].max()+1)
    
    ## STEADIERCOM - with filtering
    # abundances_df = pd.read_csv("../output/Venturelli2018/results_steadiercom/results_pairwise_growth/_abundance.tsv",sep="\t")
    # abundances_df.member = abundances_df.member.map(_change_name)

    # # Filter out communities where there is only one present member
    # count_members_sample = abundances_df[["sample","community"]].value_counts()
    # count_above_1 = count_members_sample[count_members_sample>1].index
    # abundance_filtered = abundances_df[abundances_df[["sample", "community"]].apply(tuple, axis=1).isin(count_above_1)].copy()
    # # Count number of samples per community
    # sample_counts = abundance_filtered["community"].value_counts()/2 # Divide by 2 because we have two entries per member
    # steadiercom_abundance = abundance_filtered.groupby(['community', 'member']).sum()["relative_abundance"]/sample_counts
    # steadiercom_abundance = steadiercom_abundance.sort_index()
    # steadiercom_abundance.name = "sim"

    ## STEADYCOM
    steadycom_abundance = pd.read_csv("../output/Venturelli2018/steadycom/pairwise_no_abundance/relative_abundance.tsv",sep="\t",index_col=0)

    steadycom_abundance.index = steadycom_abundance.index.map(lambda x: re.sub(r"^M_", "", x))
    steadycom_abundance.index = steadycom_abundance.index.map(_change_name)
    steadycom_abundance = (steadycom_abundance.T
                           .stack()
                           .sort_index())
    steadycom_abundance.index.names = ["community","member"]
    steadycom_abundance.name = "sim"

    return steadiercom_abundance,steadycom_abundance











