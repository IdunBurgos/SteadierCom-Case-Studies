import reframed
from reframed import SteadyCom
import pandas as pd
import numpy as np
from os import listdir
import warnings
import seaborn as sns
from collections import OrderedDict
import json
import sys
sys.path.append("../functions/")

import general

## Load data and models


def load_data_models_guillen2021():
    # NB! yhhK in data == panZ in paper
    relative_abundance_ecoli = pd.read_csv("../data/Guillen2021/relative_abundance.csv",index_col=0)
    # NB! yhhK in data == panZ in paper
    relative_abundance_ecoli.index = relative_abundance_ecoli.index.map(lambda x:x if x!="yhhK" else "panZ")
    relative_abundance_ecoli = relative_abundance_ecoli.T.mean().sort_values(ascending=False).head(7)
    
    models_ecoli = general.load_models_guillen2021(all_38_strains=False)

    return models_ecoli, relative_abundance_ecoli

def load_data_models_venturelli2018():
    # Load media
    media_db_gut = pd.read_csv("../data/Venturelli2018/generated/media_db.tsv",sep="\t")
    reactions = media_db_gut[media_db_gut.medium=="LB_extend"].compound.values
    env_gut = reframed.Environment.from_reactions(reactions,max_uptake=10)

    datasetev3 = pd.read_excel("../data/Venturelli2018/msb178157-sup-0004-datasetev3.xlsx",sheet_name="Dataset EV3")
    
    models_gut = general.load_models_venturelli2018()

    return models_gut, env_gut, datasetev3

def make_community_smetana_guillen2021(models_ecoli):
    model_ids = [model.id for model in models_ecoli.values()]
    community_ids = ["ecoli_auxotrophs" for model in models_ecoli.values()]

    return pd.DataFrame([community_ids,model_ids])

def make_community_steadiercom_guillen2021(relative_abundance_ecoli):

    community_steadiercom_ecoli = relative_abundance_ecoli.reset_index()
    community_steadiercom_ecoli["community"] = community_steadiercom_ecoli["index"].map(lambda x:"Ecoli_mutants")
    community_steadiercom_ecoli["index"] = community_steadiercom_ecoli["index"].map(lambda x:"mutant_"+x)
    community_steadiercom_ecoli.columns=["member","abundance","community"]
    community_steadiercom_ecoli = community_steadiercom_ecoli[["community","member","abundance"]]

    return community_steadiercom_ecoli

def make_community_smetana_venturelli2018(models_gut):
    df_comms = []

    # Full community
    model_ids = [model_id for model_id in models_gut.keys()]
    community_ids = ["Full" for model_id in model_ids]
    df_comm = pd.DataFrame([community_ids,model_ids])
    df_comms.append(df_comm)

    # Reduced communities
    for model_id in models_gut.keys():
        model_ids = list(models_gut.keys())
        model_ids.remove(model_id)
        community_id = model_id.split("_")[0][0].upper() + model_id.split("_")[1][0].upper()+"-"
        community_ids = [community_id for id_ in model_ids]
        
        df_comm = pd.DataFrame([community_ids,model_ids])
        df_comms.append(df_comm)
    communities_gut = pd.concat(df_comms,axis=1).T

    return communities_gut

def make_community_steadiercom_venturelli2018(datasetev3, models_gut):
    member_id2member_name = {x.split("_")[0][0]+x.split("_")[1][0].upper():x for x in models_gut.keys()}
    abundance_dict = {}
    last_item = "'ER-'"
    abundance_dict[last_item]=[]

    for i,row in datasetev3.iterrows():
        list_ = row.iloc[0].split("    ")
        
        if len(list_)==1:
            last_item =list_[0]
            abundance_dict[last_item]=[]
        else:
            abundance_dict[last_item].append(list_[1:])

    abundance_dataframes = {community:pd.DataFrame(abundances).applymap(lambda x: float(x)) for community,abundances in abundance_dict.items()}

    for community,dataframe in abundance_dataframes.items(): 
        dataframe.columns=["BH","CA","BU","PC","BO","BV","BT","EL","FP","CH","DP","ER"] 
    for community,dataframe in abundance_dataframes.items(): 
        dataframe.columns = dataframe.columns.map(lambda x:member_id2member_name[x])

    abundance_series = {community_id.replace("'",""):dataframe.iloc[-1,:] for community_id,dataframe in abundance_dataframes.items()}
    for community,series in abundance_series.items(): 
        if community!="NONE": 
            abundance_series[community]=series.drop(member_id2member_name[community[:-1]])

    communities_steadiercom_gut = pd.concat(abundance_series).reset_index()
    communities_steadiercom_gut.columns = ["community","member","abundance"]

    return communities_steadiercom_gut

def make_media_smetana_venturelli2018(env_gut):

    mets_gut = []
    for rxn in env_gut.keys():
        mets_gut.append(rxn.replace("R_EX_","").replace("__91__e__93__",""))

    media_name = ["LB_extend" for met in mets_gut]
    media_db_gut_smetana = pd.DataFrame({"medium":media_name,"compound":mets_gut})
    
    return media_db_gut_smetana

def main():

    models_ecoli,relative_abundance_ecoli = load_data_models_guillen2021()
    models_gut, env_gut, datasetev3 = load_data_models_venturelli2018()

    community_smetana_ecoli = make_community_smetana_guillen2021(models_ecoli)
    communities_smetana_gut = make_community_smetana_venturelli2018(models_gut)

    media_db_gut_smetana = make_media_smetana_venturelli2018(env_gut)
    communities_steadiercom_gut = make_community_steadiercom_venturelli2018(datasetev3, models_gut)

    
    ## SAVE DATA

    media_db_gut_smetana.to_csv("../data/Venturelli2018/generated/media_db_smetana.tsv",
                                index=False,
                                sep="\t")

    communities_smetana_gut.to_csv("../data/Venturelli2018/generated/communities_smetana.tsv",
                                            index=None,
                                            header=None,
                                            sep="\t")
    
    community_smetana_ecoli.to_csv("../data/Guillen2021/generated/communities_smetana.tsv",
                                            index=None,
                                            header=None,
                                            sep="\t")
    


    community_steadiercom_ecoli = make_community_steadiercom_guillen2021(relative_abundance_ecoli)
    community_steadiercom_ecoli\
        .set_index("community")\
        .reset_index()\
        .to_csv("../data/Guillen2021/generated/communities_steadiercom.tsv",
                                                    sep="\t",
                                                    header=None,
                                                    index=None)
    
    for community_id in communities_steadiercom_gut.community.unique():
        
        new_id = community_id.replace("-","_minus")
        communities_steadiercom_gut.query(f"community=='{community_id}'").to_csv(f"../data/Venturelli2018/generated/communities_steadiercom/{new_id}_abundance.tsv",
                                                        sep="\t",
                                                        header=None)
        
        communities_steadiercom_gut.query(f"community=='{community_id}'").to_csv(f"../data/Venturelli2018/generated/communities_steadiercom_no_abundance/{new_id}_abundance.tsv",
                                                        sep="\t",
                                                        header=None,
                                                        index=None)
        


        
if __name__ == "__main__":
    main()