import pandas as pd

def _change_name(x):
    names =x.split("_")
    name = names[0][0].upper()+names[1][0].upper()
    return name


def load_experimental_results():

    exp_abundance_df  = pd.read_csv("../data/Guillen2021/generated/communities_steadiercom.tsv",sep="\t",header=None)
    exp_abundance_df.columns = ["community","member","abundance"]
    exp_abundance_df = exp_abundance_df.set_index(["community","member"])["abundance"].unstack().fillna(0)

    exp_abundance_df.sort_index(axis=1,inplace=True)
    return exp_abundance_df

def load_simulation_results():

    abundances_df = pd.read_csv(f"../output/Guillen2021/results_steadiercom/results_growth/ecoli_mutants_abundance.tsv",sep="\t")
    abundances_df.reset_index(drop=True)
    #abundances_df.member = abundances_df.member.map(_change_name)

    samples = abundances_df["sample"].max()
    steadiercom_abundance_df = abundances_df.groupby(['community', 'member'], as_index=False).agg({"relative_abundance": (lambda x: sum(x)/samples)})
    steadiercom_abundance_df = steadiercom_abundance_df.pivot_table(index='member', columns='community', values='relative_abundance', fill_value=0)
    steadiercom_abundance_df = steadiercom_abundance_df.T.sort_index(axis=1)

    steadycom_abundance_df = pd.read_csv("../output/Guillen2021/steadycom/no_abundance/relative_abundance.tsv",sep="\t",index_col=0)
    steadycom_abundance_df = steadycom_abundance_df.T
    steadycom_abundance_df.sort_index(inplace=True)
    steadycom_abundance_df.sort_index(inplace=True,axis=1)
    return steadiercom_abundance_df, steadycom_abundance_df

















