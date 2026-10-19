
import pandas as pd
import libsbml
import logging
from os import listdir
import re


logging.basicConfig(level=logging.INFO, 
    format="%(levelname)-8s | %(message)s",
)


def _abbreviate_member(name):
    if name == "environment":
        return name

    names =name.split("_")
    new_name = names[0][0].upper()+names[1][0].upper()
    return new_name


def abbreviate_member_columns(df):
    df = df.copy()

    for col in ["donor", "receiver"]:
        df[col] = df[col].map(_abbreviate_member)

    return df

def load_glv_data():
    """
    Returns a pandas Series with the gLV coefficients from the Venturelli2018 model. 
        The index is a MultiIndex with the first level being the second member and the second level being the main member. 
        The values are the gLV coefficients.

    Loading data from gLV-model 
    """

    logging.info("Loading gLV coefficients from ../data/Venturelli2018/msb178157-sup-0006-codeev1/twelveSpeciesT4.xml")
    
    glv_model_file = libsbml.readSBMLFromFile("../data/Venturelli2018/msb178157-sup-0006-codeev1/twelveSpeciesT4.xml")
    glv_model = glv_model_file.getModel()
    glv_model.getNumSpecies()
    glv_parameters = {}

    for param in glv_model.parameters:
        if param.id.startswith("a"):
            key = param.id
            glv_parameters[(key[1:3],key[3:])]=param.value

    signal_exp = pd.Series(glv_parameters)
    signal_exp.index.names = ["member_second","member_main"]
    signal_exp.name="experimental"

    return signal_exp

def load_metabolomics_data():

    """
    Returns a pandas Series with the log-change in metabolite concentration from the Venturelli2018 study. 
    """
    logging.info("Loading metabolomics data from ../data/Venturelli2018/msb178157-sup-0005-datasetev4.xlsx"
                 )
    # Load dataset for log-change in concentration
    metabolomics = pd.read_excel("../data/Venturelli2018/msb178157-sup-0005-datasetev4.xlsx",sheet_name="Dataset EV4",index_col=0)

    # Adjust name of metabolites
    metabolomics.index = metabolomics.index.map(lambda met: met.replace("'", "").replace(" ", "").lower())

    # Adjust name of member
    metabolomics.columns = metabolomics.columns.map(lambda x: x.replace("'","").replace(" ",""))
    metabolomics.columns.name="member"
    metabolomics.index.name="compound"
    metabolomics.name="log-change"
    metabolomics = metabolomics.stack()
    return metabolomics

def load_experimental_abundance(all_comms=True):

    """
    Returns a pandas DataFrame with the experimental abundances from the Venturelli2018 study.
        The index is the community, the columns are the members and the values are the abundances.
    """
    if all_comms:
        logging.info("Loading experimental abundances for all communities from ../data/Venturelli2018/generated/communities_steadiercom/")

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
        exp_abundance_df.columns = exp_abundance_df.columns.map(_abbreviate_member)

        exp_abundance_df.sort_index(axis=1,inplace=True)

    else:
        logging.info("Loading experimental abundances for FULL community from ../data/Venturelli2018/generated/communities_steadiercom/NONE_abundance.tsv")
        exp_abundance_df = pd.read_csv("../data/Venturelli2018/generated/communities_steadiercom/NONE_abundance.tsv",sep="\t",header=None)
        exp_abundance_df.columns = ["community", "member", "relative_abundance"]
        exp_abundance_df.member = exp_abundance_df.member.map(_abbreviate_member)

    return exp_abundance_df


def load_experimental_abundance_pairwise():

    """
    Returns a pandas Series with the experimental abundances from the Venturelli2018 study for pairwise communities.
        The index is a MultiIndex with the first level being the community and the second level being the member.
    """

    exp_abundance_df = pd.read_csv("../data/Venturelli2018/generated/communities_pairwise.tsv",sep="\t",header=None)

    exp_abundance_df.columns = ["community","member","abundance"]
    exp_abundance_df = exp_abundance_df.set_index(["community","member"])["abundance"].unstack()
    exp_abundance_df.columns = exp_abundance_df.columns.map(_abbreviate_member)
    exp_abundance_df

    exp_abundance_df.sort_index(axis=1,inplace=True)
    exp_abundance_ser = exp_abundance_df.stack()
    exp_abundance_ser.name="exp"
    return exp_abundance_ser



def _process_steadiercom_results(steadiercom_df):

    assert ("rate" in steadiercom_df.columns) and ("mass_rate" in steadiercom_df.columns) and ("frequency" in steadiercom_df.columns), f"This is not standard SteadierCom results. Columns: {steadiercom_df.columns}"
                    
    steadiercom_df["mass_rate*frequency"] = steadiercom_df.frequency*steadiercom_df.mass_rate 
    steadiercom_df["rate*frequency"] = steadiercom_df.frequency*steadiercom_df.rate

    steadiercom_df = abbreviate_member_columns(steadiercom_df)

    return steadiercom_df

def _process_steadycom_results(simulated):

    # Fill missing values for absent members with "environment" (SteadyCom results have NaN for absent members)
    simulated = simulated.fillna("environment")

    # SteadyCom results have "M_" prefix for member names, remove it
    simulated.donor = simulated.donor.map(lambda x: x.replace("M_","") if isinstance(x,str) else x)
    simulated.receiver = simulated.receiver.map(lambda x: x.replace("M_","") if isinstance(x,str) else x)
    
    simulated = abbreviate_member_columns(simulated)

    return simulated


def _process_steadystate_result_dict(steadiercom_dict,steadycom_dict,smetana_dict, verbose=True):

    # SAME AS FOR GUILLEN - CONSIDER MOVING TO SHARED SPACE
    if verbose:
        logging.info(f"_process_steadystate_result_dict cleans up namespace for the simulation results.")
    simulation_dict = {}

    for simulation_method, sim_dict in [("steadiercom",steadiercom_dict),("steadycom",steadycom_dict),("smetana",smetana_dict)]:
        
        simulation_dict[simulation_method]={}

        for name, simulated in sim_dict.items():
            if simulation_method=="steadiercom" or simulation_method=="steadycom":

                # steadycom-specific processing
                if simulation_method=="steadycom":
                    simulated = _process_steadycom_results(simulated)
                # steadiercom-specific processing
                if simulation_method=="steadiercom":
                    assert ("rate" in simulated.columns) and ("mass_rate" in simulated.columns) and ("frequency" in simulated.columns), f"This is not standard SteadierCom results. Columns: {simulated.columns}"
                    simulated = _process_steadiercom_results(simulated)

        # CONTINUED
            simulation_dict[simulation_method][name]=simulated

    return simulation_dict

def load_simulated_steadystate_results(verbose=True):
    """
    Loads experimental data from all simulation methods (But only for the full community 'NONE')
    """

    if verbose:
        logging.info(f"load_simulation_results loads simulation results from \n\t../output/Venturelli2018/results_steadiercom, \n\t../output/Venturelli2018/steadycom and \n\t../output/Venturelli2018/smetana"\
                     "\n Only loads results for community 'NONE' (full community)")

    steadiercom_abundance= pd.read_csv("../output/Venturelli2018/results_steadiercom/results_abundance/NONE.tsv",sep="\t")
    steadiercom_growth = pd.read_csv("../output/Venturelli2018/results_steadiercom/results_growth/NONE.tsv",sep="\t")
    steadiercom_growth_abundance = pd.read_csv("../output/Venturelli2018/results_steadiercom/results_growth_abundance/NONE.tsv",sep="\t")

    smetana = pd.read_csv("../output/Venturelli2018/smetana/no_coupling_detailed.tsv",sep="\t")
    smetana = smetana[smetana.community=="NONE"]


    steadycom_fixed_abundance = pd.read_csv("../output/Venturelli2018/steadycom/fixed_abundance/cross_feeding.tsv",sep="\t")
    steadycom = pd.read_csv("../output/Venturelli2018/steadycom/no_abundance/cross_feeding.tsv",sep="\t")

    steadycom = steadycom[steadycom.community=="NONE"]
    steadycom_fixed_abundance = steadycom_fixed_abundance[steadycom_fixed_abundance.community=="NONE"]


    steadiercom_dict = {
        "SteadierCom - Fixed growth":steadiercom_growth,
        "SteadierCom - Fixed growth & abundance":steadiercom_growth_abundance,
        "SteadierCom - Fixed abundance":steadiercom_abundance
    }

    steadycom_dict = {
        "SteadyCom - Nothing fixed":steadycom,
        "SteadyCom - Fixed abundance":steadycom_fixed_abundance
    }

    smetana_dict = {
        "SMETANA": smetana
    }

    simulation_dict = _process_steadystate_result_dict(steadiercom_dict,steadycom_dict,smetana_dict)

    return simulation_dict


def load_all_steadiercom_steadystate_results(verbose=True):
    """
    Loads experimental data from all simulation methods (But only for the full community 'NONE')
    """
    if verbose:
        logging.info(f"load_simulation_results loads simulation results from \n\t../output/Venturelli2018/results_steadiercom,"\
                     "\n Loads for ALL communities")

    df_res = []
    for file in listdir("../output/Venturelli2018/results_steadiercom/results_growth_abundance/"):
        if file.endswith(".tsv"):
            df = pd.read_csv(f"../output/Venturelli2018/results_steadiercom/results_growth_abundance/{file}", sep="\t")
            df_res.append(df)

    steadiercom_results_all = pd.concat(df_res, axis=0)

    steadiercom_results_all = _process_steadiercom_results(steadiercom_results_all)

    return steadiercom_results_all



def load_simulated_abundance():

    # Load SteadierCom results for all communities
    abundance_list = []
    for file in listdir("../output/Venturelli2018/results_steadiercom/results_growth/"):
        if file.endswith("_abundance.tsv"):
            abundance_list.append(pd.read_csv(f"../output/Venturelli2018/results_steadiercom/results_growth/{file}",sep="\t"))
    steadiercom_abundance_df = pd.concat(abundance_list) # Join all dataframes
    # Map member names to their initials

    steadiercom_abundance_df.member = steadiercom_abundance_df.member.map(_abbreviate_member)
    # Fill missing values for absent members with zero abundance
    steadiercom_abundance_df = (
        steadiercom_abundance_df.pivot_table(index=["community","sample"], columns='member', values='relative_abundance', fill_value=0)
        .stack()
        .reset_index()
    )
    steadiercom_abundance_df.rename(columns={0:"relative_abundance"}, inplace=True)

    # Load SteadyCom abundance
    steadycom_abundance_df = pd.read_csv("../output/Venturelli2018/steadycom/no_abundance/relative_abundance.tsv",sep="\t",index_col=0)

    steadycom_abundance_df.index = steadycom_abundance_df.index.map(lambda x: x.replace("M_","") if isinstance(x,str) else x)
    steadycom_abundance_df.index = steadycom_abundance_df.index.map(_abbreviate_member)
    steadycom_abundance_df = steadycom_abundance_df.T

    steadycom_abundance_df.sort_index(inplace=True)
    steadycom_abundance_df.sort_index(inplace=True,axis=1)

    return steadiercom_abundance_df, steadycom_abundance_df



def load_simulated_abundance_pairwise():
    ## STEADIERCOM - without filtering
    abundances_df = pd.read_csv("../output/Venturelli2018/results_steadiercom/results_pairwise_growth/_abundance.tsv",sep="\t")
    abundances_df.member = abundances_df.member.map(_change_name)
    steadiercom_abundance = abundances_df.groupby(['community', 'member']).sum()["relative_abundance"]/(abundances_df["sample"].max()+1)

    ## STEADYCOM
    steadycom_abundance = pd.read_csv("../output/Venturelli2018/steadycom/pairwise_no_abundance/relative_abundance.tsv",sep="\t",index_col=0)

    steadycom_abundance.index = steadycom_abundance.index.map(lambda x: re.sub(r"^M_", "", x))
    steadycom_abundance.index = steadycom_abundance.index.map(_abbreviate_member)
    steadycom_abundance = (steadycom_abundance.T
                           .stack()
                           .sort_index())
    steadycom_abundance.index.names = ["community","member"]
    steadycom_abundance.name = "sim"

    return steadiercom_abundance,steadycom_abundance











