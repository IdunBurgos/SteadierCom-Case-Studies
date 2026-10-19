import reframed
from reframed import SteadyCom
import pandas as pd
import numpy as np
import sys
sys.path.append("../functions/")

import simulation.general as general

import logging
logging.basicConfig(
    level=logging.INFO, 
    format="%(levelname)-8s | %(message)s",
)

## Load data and models


def load_data_guillen2021():
    logging.info("Loading relative abundance data from Guillen et al (2021)...")
    relative_abundance_ecoli = pd.read_csv("../data/Guillen2021/relative_abundance.csv",index_col=0)
    # NB! yhhK in data == panZ in paper
    relative_abundance_ecoli.index = relative_abundance_ecoli.index.map(lambda x:x if x!="yhhK" else "panZ")
    relative_abundance_ecoli = relative_abundance_ecoli.T.mean().sort_values(ascending=False).head(7)
    
    return relative_abundance_ecoli

def load_data_models_venturelli2018():

    logging.info("Loading models and data from Venturelli et al (2018)...")
    # Load media
    media_db_gut = pd.read_csv("../data/Venturelli2018/generated/media_db.tsv",sep="\t")
    reactions = media_db_gut[media_db_gut.medium=="LB_extend"].compound.values
    env_gut = reframed.Environment.from_reactions(reactions,max_uptake=10)
    
    datasetev1 = pd.read_excel("../data/Venturelli2018/msb178157-sup-0002-datasetev1.xlsx",sheet_name="PW1",header=None,index_col=0)
    datasetev3 = pd.read_excel("../data/Venturelli2018/msb178157-sup-0004-datasetev3.xlsx",sheet_name="Dataset EV3")
    
    models_gut = general.load_models_venturelli2018()

    return models_gut, env_gut, datasetev1, datasetev3

def make_media_smetana_venturelli2018(env_gut):
    """
    SMETANA requires a difference prefix for the media_db compared to SteadierCom. 
    """
    logging.info("Creating media db for SMETANA from environment used in Venturelli et al (2018)...")
    mets_gut = []
    for rxn in env_gut.keys():
        mets_gut.append(rxn.replace("R_EX_","").replace("__91__e__93__",""))

    media_name = ["LB_extend" for met in mets_gut]
    media_db_gut_smetana = pd.DataFrame({"medium":media_name,"compound":mets_gut})
    
    return media_db_gut_smetana



def make_community_df_guillen2021(relative_abundance_ecoli):

    logging.info("Creating community df for Guillen et al (2021) for SteadierCom, SteadyCom, and SMETANA ...")
    community_df_ecoli = relative_abundance_ecoli.reset_index()
    community_df_ecoli["community"] = community_df_ecoli["index"].map(lambda x:"Ecoli_mutants")
    community_df_ecoli["index"] = community_df_ecoli["index"].map(lambda x:"mutant_"+x)
    community_df_ecoli.columns=["member","abundance","community"]
    community_df_ecoli = community_df_ecoli[["community","member","abundance"]]

    return community_df_ecoli


def make_community_df_venturelli2018(datasetev3, models_gut):

    logging.info("Creating community df for Venturelli et al (2018) for SteadierCom, SteadyCom, and SMETANA ...")
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

    communities_df_gut = pd.concat(abundance_series).reset_index()
    communities_df_gut.columns = ["community","member","abundance"]

    return communities_df_gut


def make_pairwise_df_venturelli2018(datasetev1,models_gut):
    """
    Make df from pairwise data from Venturelli 2018 study for SteadierCom, SteadyCom and SMETANA.
    """
    logging.info("Creating pairwise df for Venturelli et al (2018) for SteadierCom, SteadyCom, and SMETANA ...")    

    # Clean input data
    datasetev1.index = datasetev1.index.map(lambda x: x.replace("'","").replace(" ",""))
    new_index = [(string[0:2],string[2:]) for string in datasetev1.index.values]
    datasetev1.index = pd.MultiIndex.from_tuples(new_index, names=["member1","member2"])
    relative_abundance_pairwise = datasetev1.iloc[:,-1].fillna(0) # Get only final timepoint

    # Make community pairwise gut
    id2name_member = {name.split("_")[0][0].upper()+name.split("_")[1][0].upper():name for name in models_gut.keys()}
    community_pairwise_gut=[]

    for index,value in relative_abundance_pairwise.items():
        member1=index[0]
        member2=index[1]
        community_pairwise_gut.append({"community":"-".join(index),"member":member1,"abundance":value})
        community_pairwise_gut.append({"community":"-".join(index),"member":member2,"abundance":1-value})

    community_df_pairwise_gut = pd.DataFrame(community_pairwise_gut)
    community_df_pairwise_gut["member"] = community_df_pairwise_gut["member"].map(id2name_member)
    return community_df_pairwise_gut

def make_pairwise_df_guillen2021(dataset_pairwise):
    """
    Make df from pairwise data from Guillen 2021 study for SteadierCom, SteadyCom and SMETANA.
    """
    #return pairwise_df


def main():

    logging.info("STARTING SCRIPT FOR PREPARING FILES FOR SIMULATION")
    logging.info("LOADING DATA...")
    relative_abundance_ecoli = load_data_guillen2021()
    models_gut, env_gut, datasetev1,datasetev3 = load_data_models_venturelli2018()

    logging.info("CREATING DATA FILES FOR SIMULATION...")
    media_db_gut_smetana = make_media_smetana_venturelli2018(env_gut)
    communities_df_gut = make_community_df_venturelli2018(datasetev3, models_gut)
    community_df_ecoli = make_community_df_guillen2021(relative_abundance_ecoli)
    community_df_pairwise_gut = make_pairwise_df_venturelli2018(datasetev1,models_gut)    
    # pairwise_ecoli_df = make_pairwise_df_guillen2021(dataset_pairwise)

    ## SAVE DATA
    logging.info("SAVING DATA...")

    logging.info("Saving media database for SMETANA in ../data/Venturelli2018/generated/media_db_smetana.tsv \n" \
        " Saving SMETANA community data for Venturelli et al (2018) in ../data/Venturelli2018/generated/communities_smetana.tsv \n" \
        " Saving SMETANA community data for Guillen et al (2021) in ../data/Guillen2021/generated/communities_smetana.tsv \n" \
        " Saving STEADIERCOM community data for Venturelli et al (2018) in ../data/Venturelli2018/generated/communities_steadiercom/ and communities_steadiercom_no_abundance/ \n" \
        " Saving STEADIERCOM community data for Guillen et al (2021) in ../data/Guillen2021/generated/communities_steadiercom.tsv \n" \
        " Saving STEADIERCOM community data for PAIRWISE COMMUNITIES Venturelli et al (2018) in ../data/Venturelli2018/generated/communities_pairwise.tsv and communities_pairwise_no_abundance.tsv \n"
    )
    media_db_gut_smetana.to_csv("../data/Venturelli2018/generated/media_db_smetana.tsv",
                                index=False,
                                sep="\t")

    communities_df_gut.iloc[:,:-1].to_csv("../data/Venturelli2018/generated/communities_smetana.tsv",
                                            index=None,
                                            header=None,
                                            sep="\t")
    
    community_df_ecoli.iloc[:,:-1].to_csv("../data/Guillen2021/generated/communities_smetana.tsv",
                                            index=None,
                                            header=None,
                                            sep="\t")
    community_df_ecoli\
        .set_index("community")\
        .reset_index()\
        .to_csv("../data/Guillen2021/generated/communities_steadiercom.tsv",
                                                    sep="\t",
                                                    header=None,
                                                    index=None)
    
    for community_id in communities_df_gut.community.unique():
        
        new_id = community_id.replace("-","_minus")
        communities_df_gut.query(f"community=='{community_id}'").to_csv(f"../data/Venturelli2018/generated/communities_steadiercom/{new_id}_abundance.tsv",
                                                        sep="\t",
                                                        header=None,
                                                        index=None)

        communities_df_gut.query(f"community=='{community_id}'").iloc[:,:-1].to_csv(f"../data/Venturelli2018/generated/communities_steadiercom_no_abundance/{new_id}_abundance.tsv",
                                                        sep="\t",
                                                        header=None,
                                                        index=None)
    community_df_pairwise_gut.to_csv("../data/Venturelli2018/generated/communities_pairwise.tsv",
                                     sep="\t",
                                     index=None,
                                     header=None)

    community_df_pairwise_gut.iloc[:,:-1].to_csv("../data/Venturelli2018/generated/communities_pairwise_no_abundance.tsv",
                                    sep="\t",
                                    index=None,
                                    header=None)
    # pairwise_ecoli_df.to_csv("../data/Guillen2021/generated/pairwise_df.tsv",sep="\t",index=False)

if __name__ == "__main__":
    main()