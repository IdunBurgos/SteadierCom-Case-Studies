import pandas as pd
import numpy as np
import json
import os

import seaborn as sns
import matplotlib.pyplot as plt


from sklearn.metrics import confusion_matrix

import sys
sys.path.append("../functions/")
import general
import compare_met_exp_sim.correlation as correlation
import logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def _assign_groups_and_filter_by_models(dataset4_stack,metname_groups,verbose=True):
    """
    This is a helper function for load_experimental_results_venturelli. It assigns metabolite groups to the metabolites 
    and filters out metabolites that do not exist in any model. It returns the filtered dataframe,
    the set of metabolite groups found in the models and the list of metabolite names not included in any model.

    Inputs:
    ------------
    dataset4_stack: pd.DataFrame
        Dataframe with columns: metname, member, change.
    metname_groups: dict
        Dictionary with metabolite groups as keys and dictionaries as values. The inner dictionaries have metabolite names as keys and weights as values.
    verbose: bool
        Whether to print information during the function execution. Default is True.  
    Outputs:
    ------------
    dataset4_stack_filtered: pd.DataFrame
        Dataframe with metabolite groups assigned and filtered by model metabolites.
    met_groups_model: set of tuples
        Set of tuples with (member_id, met_group) found in the models.
    metnames_not_included: list
        List of metabolite names not included in any model.
    ------------
    """
    
    # VALIDATE INPUT
    assert isinstance(metname_groups,dict), f"metname_groups is not a dict: {type(metname_groups)}"
    assert len(metname_groups)!=0, f"metname_groups is empty: {type(metname_groups)}"
    assert isinstance(dataset4_stack,pd.DataFrame), f"dataset4_stack is not a pd.DataFrame: {type(dataset4_stack)}"
    assert not dataset4_stack.empty, f"dataset4_stack is empty"

    if verbose:
        logging.info(f"_assign_groups_and_filter_by_models assigns metabolite groups to the metabolites \n"\
                     " and filters out metabolites that do not exist in any model.")      

    dataset4_stack = dataset4_stack.copy()
    ## LOAD MODELS
    models = general.load_models_venturelli2018()
    
    ## LOAD DATA
    with open("../data/Venturelli2018/generated/metid2metname.json") as file:
        metid2metname = json.load(file)

    ## PREPARE ADDITIONAL DATA
    metname2group = {name:[] for names in metname_groups.values() for name in names}
    for group,group_dict in metname_groups.items():
        for name in group_dict.keys():
            metname2group[name].append(group)

    metid2group = {met_id:metname2group[met_name] for met_id,met_name in metid2metname.items()}

    mapping_df = (
        pd.Series(metname2group, name='met_group')
        .explode()
        .rename_axis('metname')
        .reset_index()
    )

    assert len(metid2metname)!=0, "metid2metname is empty"
    assert len(metid2group)!=0, "metid2group is empty"
    assert not mapping_df.empty, "mapping_df is empty"

    # Find mets and met_groups matching with model
    met_ids_keys = metid2metname.keys()
    met_groups_model  = []
    groups_included = []
    mets_included = []
    for member_name,model in models.items():
        # For each model check all metabolites 
        for met_id in model.metabolites:
            if ("__e__" in met_id) and (met_id in met_ids_keys):
                # If metabilite is external and in curated dict based on experimental data.
                for group in metid2group[met_id]:
                    member_id = member_name.split("_")[0][0].upper()+member_name.split("_")[1][0].upper()
                    # Add groups that the met is assigned to list
                    met_groups_model.append((member_id,group))
                    groups_included.append(group)
                mets_included.append(met_id)
                    
    met_groups_model = set(met_groups_model)

    groups_included = set(groups_included)
    mets_not_included = set(met_ids_keys) - set(mets_included)
    metnames_not_included = [metid2metname[metid] for metid in mets_not_included]

    assert len(met_groups_model)!=0, "met_groups_model is empty meaning that no metabolites from experimental data are found in any model"

    ## ASSIGN GROUPS
    # Add group name to df
    dataset4_stack = dataset4_stack.merge(mapping_df,on="metname",how="left",)
    assert not dataset4_stack.empty, "dataset4_stack is empty after merging with mapping_df"

    ## FILTER METABOLITES BASED ON MODEL COMPONENTS
    #mask = dataset4_stack[["member","met_group"]].apply(tuple,axis=1).isin(met_groups_model)
    mask = dataset4_stack["met_group"].isin(groups_included)
    dataset4_stack_filtered = dataset4_stack[mask].copy() 
    assert not dataset4_stack_filtered.empty, "dataset4_stack is empty after filtering based on model metabolites"
    
    return dataset4_stack_filtered,met_groups_model,metnames_not_included


def load_experimental_results(verbose=True):
    """
    Function for processing metabolomics data which gives output 'signal_exp' to be compared with simulation results.
    """

    ## LOAD DATA
    if verbose:
        logging.info(f"load_experimental_results_venturelli loads and processes metabolomics data from Venturelli et al. 2018 (Dataset EV4) in"\
                     "\n\t ../data/Venturelli2018/msb178157-sup-0005-datasetev4.xlsx"\
                     "\n\n\t Results are filtered by metabolites present in models.")

    with open("../data/Venturelli2018/generated/metname_groups.json") as file:
        metname_groups = json.load(file)
    
    dataset4 = pd.read_excel("../data/Venturelli2018/msb178157-sup-0005-datasetev4.xlsx",sheet_name="Dataset EV4",index_col=0)
    
    assert len(metname_groups)!=0, "metname_groups is emppty"
    assert not dataset4.empty, "dataset4 is empty"
    assert len(dataset4.index)==97, f"The number of rows (metabolites) in dataset4 should be 97 but is: {len(dataset4.index)}"
    assert len(dataset4.columns)==12, f"The number of columns (members) in dataset4 should be 12 but is: {len(dataset4.columns)}"

    # Adjust name of metabolites
    dataset4.index = dataset4.index.map(lambda met: met.replace("'", "").replace(" ", "").lower())
    
    # Adjust name of member
    dataset4.columns = dataset4.columns.map(lambda x: x.replace("'","").replace(" ",""))
    
    ## RESHAPE DATA AND ASSIGN GROUP
    dataset4_stack = dataset4.stack().reset_index()
    assert len(dataset4_stack.columns)==3, "number of columns in dataset4_stack should be 3"
    dataset4_stack.columns = ["metname","member","change"]

    # Filter out metabolites not in any model and assign groups
    dataset4_stack_filtered,met_groups_model,metnames_not_included = _assign_groups_and_filter_by_models(dataset4_stack,metname_groups)

    ## ASSIGN WEIGHTS AND AGGREGATE
    # Add met weight to df
    dataset4_stack_filtered["weight"]=dataset4_stack_filtered.apply(lambda x: metname_groups[x["met_group"]][x["metname"]],axis=1)
    dataset4_stack_filtered["change*weight"]=dataset4_stack_filtered.change*dataset4_stack_filtered.weight

    # Find the sum of the change and weights within a group

    combined_signal = dataset4_stack_filtered[["member","met_group","change*weight"]].groupby(["member","met_group"]).sum()["change*weight"]
    combined_weights = dataset4_stack_filtered[["member","met_group","weight"]].groupby(["member","met_group"]).sum()["weight"]

    # Find the geometric average of the signal
    signal_exp = (combined_signal/combined_weights).unstack()

    ## RESHAPE RESULTS
    # Prepare data for joining with steadiercom data
    signal_exp = signal_exp.stack()
    signal_exp.name = "experimental"

    raw_data_exp = dataset4.copy()
 
    return raw_data_exp,signal_exp,met_groups_model,metnames_not_included


def _prepare_sim_met(simulated_results,metid2group,receiver_or_donor="receiver",verbose=True):
    
    """
    Function to create new column in simulation results with metabolite group "met_group" and to filter out results we are not interested in.
    """
    
    if verbose:
        logging.info(f"_prepare_sim_met_venturelli removes environment from {receiver_or_donor}s and all compounds that are not found in metabolomics and map to AGORA identifiers")
    
    # Map met id to met group
    if receiver_or_donor=="receiver":
        simulated_results = simulated_results[simulated_results.receiver!="environment"].copy()
        simulated_results = simulated_results[~simulated_results.receiver.isna()].copy()
        simulated_results.receiver = simulated_results.receiver.map(lambda x: x.split("_")[0][0].upper()+x.split("_")[1][0].upper())

    elif receiver_or_donor=="donor":
        
        simulated_results = simulated_results[simulated_results.donor!="environment"].copy()
        simulated_results = simulated_results[~simulated_results.donor.isna()].copy()
        simulated_results.donor = simulated_results.donor.map(lambda x: x.split("_")[0][0].upper()+x.split("_")[1][0].upper())
    else:
        raise Exception("receiver_or_donor must be either 'receiver' or 'donor'")

    # Filter metabolites based on those present in experimental data
    simulated_res_reduced = simulated_results[simulated_results.compound.isin(metid2group.keys())].copy()
    simulated_res_reduced["met_group"] = simulated_res_reduced.compound.map(lambda x: metid2group[x])
    
    return simulated_res_reduced

def _prepare_all_simulations(sim_dicts, metid2group,verbose=True):
    """
    Function to prepare all simulation results for analysis with prepare_sim_met_venturelli.
    """
    
    simulation_dict = {"receiver": {}, "donor": {}}
    for role in simulation_dict:
        i=0
        for method, sim_data in sim_dicts.items():
            verbose_spec= (verbose and i==0)

            simulation_dict[role][method] = {
                name: _prepare_sim_met(
                    simulated_results=res,
                    metid2group=metid2group,
                    receiver_or_donor=role,
                    verbose=verbose_spec
                )
                for name, res in sim_data.items()
            }

            i+=1
    return simulation_dict


def load_simulation_results(metid2group,all_communities=False,verbose=True):

    if verbose:
        logging.info(f"load_simulation_results loads simulation results from \n\t../output/Venturelli2018/results_steadiercom, \n\t../output/Venturelli2018/steadycom and \n\t../output/Venturelli2018/smetana")

    # SteadierCom results

    if not all_communities:
        if verbose:
            logging.info(f"Only loads SteadierCom results for full community ('NONE')")
        steadiercom_abundance = pd.read_csv("../output/Venturelli2018/results_steadiercom/results_abundance/NONE.tsv",sep="\t")
        steadiercom_growth = pd.read_csv("../output/Venturelli2018/results_steadiercom/results_growth/NONE.tsv",sep="\t")
        steadiercom_growth_abundance = pd.read_csv("../output/Venturelli2018/results_steadiercom/results_growth_abundance/NONE.tsv",sep="\t")
    else:
        steadiercom_abundance_list = []
        for file in os.listdir("../output/Venturelli2018/results_steadiercom/results_abundance/"):
            if file.endswith(".tsv") and "growth" not in file:
                df = pd.read_csv(f"../output/Venturelli2018/results_steadiercom/results_abundance/{file}",sep="\t")
                steadiercom_abundance_list.append(df)
        steadiercom_abundance = pd.concat(steadiercom_abundance_list)

        steadiercom_growth_list = []
        for file in os.listdir("../output/Venturelli2018/results_steadiercom/results_growth/"):
            if file.endswith(".tsv") and "abundance" not in file:
                df = pd.read_csv(f"../output/Venturelli2018/results_steadiercom/results_growth/{file}",sep="\t")
                steadiercom_growth_list.append(df)
        steadiercom_growth = pd.concat(steadiercom_growth_list)

        steadiercom_growth_abundance_list = []
        for file in os.listdir("../output/Venturelli2018/results_steadiercom/results_growth_abundance/"):
            if file.endswith(".tsv"):
                df = pd.read_csv(f"../output/Venturelli2018/results_steadiercom/results_growth_abundance/{file}",sep="\t")
                steadiercom_growth_abundance_list.append(df)
        steadiercom_growth_abundance = pd.concat(steadiercom_growth_abundance_list) 


    steadiercom_abundance["frequency*mass_rate"] = steadiercom_abundance.frequency * steadiercom_abundance.mass_rate
    steadiercom_growth["frequency*mass_rate"] = steadiercom_growth.frequency * steadiercom_growth.mass_rate
    steadiercom_growth_abundance["frequency*mass_rate"] = steadiercom_growth_abundance.frequency * steadiercom_growth_abundance.mass_rate

    # SMETANA results
    smetana = pd.read_csv("../output/Venturelli2018/smetana/no_coupling_detailed.tsv",sep="\t")

    if not all_communities:
        smetana = smetana[smetana.community=="NONE"]

    # SteadyCom results
    steadycom_fixed_abundance = pd.read_csv("../output/Venturelli2018/steadycom/fixed_abundance/cross_feeding.tsv",sep="\t")
    steadycom_fixed_abundance = steadycom_fixed_abundance[steadycom_fixed_abundance.community=="NONE"]
    steadycom = pd.read_csv("../output/Venturelli2018/steadycom/no_abundance/cross_feeding.tsv",sep="\t")
    if not all_communities:
        steadycom = steadycom[steadycom.community=="NONE"]

    # Adjust SteadyCom results for easier merging
    steadycom.donor = steadycom.donor.map(lambda x: x.replace("M_","") if isinstance(x,str) else x)
    steadycom.receiver = steadycom.receiver.map(lambda x: x.replace("M_","") if isinstance(x,str) else x)

    steadycom_fixed_abundance.donor = steadycom_fixed_abundance.donor.map(lambda x: x.replace("M_","") if isinstance(x,str) else x)
    steadycom_fixed_abundance.receiver = steadycom_fixed_abundance.receiver.map(lambda x: x.replace("M_","") if isinstance(x,str) else x)

    # Combine results
    steadiercom_dict = {
        "SteadierCom - fixed growth":steadiercom_growth,
        "SteadierCom - fixed growth & abundance":steadiercom_growth_abundance,
        "SteadierCom - fixed abundance":steadiercom_abundance,
    }

    steadycom_dict = {
        "SteadyCom - Nothing fixed":steadycom,
        "SteadyCom - Fixed  abundance":steadycom_fixed_abundance
    }

    smetana_dict = {"Smetana":smetana}

    simulation_sources = {
        "steadiercom": steadiercom_dict,
        "steadycom": steadycom_dict,
        "smetana": smetana_dict
    }

    # Prepare data for analysis
    simulation_dict = _prepare_all_simulations(simulation_sources, metid2group,verbose=verbose)

    return simulation_dict


def load_simulation_results_pairwise(metid2group,verbose=True):

    """
    Function to load pairwise simulation results for Venturelli2018 case study.
    """

    if verbose:
        logging.info(f"load_simulation_results loads simulation results from \n\t../output/Venturelli2018/results_steadiercom, \n\t../output/Venturelli2018/steadycom and \n\t../output/Venturelli2018/smetana")

    # SteadierCom results
    steadiercom_growth = pd.read_csv("../output/Venturelli2018/results_steadiercom/results_pairwise_growth/.tsv",sep="\t")
    steadiercom_growth_abundance = pd.read_csv("../output/Venturelli2018/results_steadiercom/results_pairwise_growth_abundance/.tsv",sep="\t")
    
    steadiercom_growth["frequency*mass_rate"] = steadiercom_growth.frequency * steadiercom_growth.mass_rate
    steadiercom_growth_abundance["frequency*mass_rate"] = steadiercom_growth_abundance.frequency * steadiercom_growth_abundance.mass_rate
    steadiercom_growth["frequency*rate"] = steadiercom_growth.frequency * steadiercom_growth.rate
    steadiercom_growth_abundance["frequency*rate"] = steadiercom_growth_abundance.frequency * steadiercom_growth_abundance.rate

    # SMETANA results
    smetana = pd.read_csv("../output/Venturelli2018/smetana/pairwise_detailed.tsv",sep="\t")

    # SteadyCom results
    steadycom_fixed_abundance = pd.read_csv("../output/Venturelli2018/steadycom/pairwise_no_abundance/cross_feeding.tsv",sep="\t")
    steadycom = pd.read_csv("../output/Venturelli2018/steadycom/pairwise_no_abundance/cross_feeding.tsv",sep="\t")


    # Adjust SteadyCom results for easier merging
    steadycom.donor = steadycom.donor.map(lambda x: x.replace("M_","") if isinstance(x,str) else x)
    steadycom.receiver = steadycom.receiver.map(lambda x: x.replace("M_","") if isinstance(x,str) else x)

    steadycom_fixed_abundance.donor = steadycom_fixed_abundance.donor.map(lambda x: x.replace("M_","") if isinstance(x,str) else x)
    steadycom_fixed_abundance.receiver = steadycom_fixed_abundance.receiver.map(lambda x: x.replace("M_","") if isinstance(x,str) else x)

    # Combine results
    steadiercom_dict = {
        "SteadierCom - fixed growth":steadiercom_growth,
        "SteadierCom - fixed growth & abundance":steadiercom_growth_abundance,
    }

    steadycom_dict = {
        "SteadyCom - Nothing fixed":steadycom,
        "SteadyCom - Fixed  abundance":steadycom_fixed_abundance
    }

    smetana_dict = {"Smetana":smetana}

    simulation_sources = {
        "steadiercom": steadiercom_dict,
        "steadycom": steadycom_dict,
        "smetana": smetana_dict
    }

    # Prepare data for analysis
    simulation_dict = _prepare_all_simulations(simulation_sources, metid2group,verbose=verbose)

    return simulation_dict


def find_correlations(simulation_dict,signal_exp):

    agg_method ="sum"
    simulation_method="steadycom"
    metric = "rate?"

    results_corr = []
    results_corr_member = []
    results_corr_mets = []

    sim_method_metric = [("steadycom","rate?"),("steadiercom","frequency*mass_rate"),("steadiercom","frequency"),("smetana","mus")]

    for simulation_method,metric in sim_method_metric:
        results_corr.append(correlation.find_correlation_simulation_method(simulation_dict,
                                                    signal_exp,
                                                    simulation_method=simulation_method,
                                                    agg_method=agg_method,
                                                    metric=metric,
                                                    verbose=False)
        )

        results_corr_member.append(correlation.find_correlation_simulation_method_detailed(simulation_dict,
                                                    signal_exp,
                                                    member_or_met="member",
                                                    simulation_method=simulation_method,
                                                    agg_method=agg_method,
                                                    metric=metric,
                                                    verbose=False)
        )

        results_corr_mets.append(correlation.find_correlation_simulation_method_detailed(simulation_dict,
                                                    signal_exp,
                                                    member_or_met="met_group",
                                                    simulation_method=simulation_method,
                                                    agg_method=agg_method,
                                                    metric=metric,
                                                    verbose=False)
        )

    results_corr_df = pd.concat(results_corr)
    results_corr_member_df = pd.concat(results_corr_member)
    results_corr_met_df = pd.concat(results_corr_mets)

    return results_corr_df,results_corr_member_df,results_corr_met_df


def plot_confusion_matrix_venturelli(exp_sim_combined_bool_dict,suptitle="No title"):

    i = 0
    fig,axs=plt.subplots(2,3,figsize=(15,10))
    axs_flat = axs.flatten()

    for name,exp_sim_combined_bool in exp_sim_combined_bool_dict.items():
            
        cf_matrix = confusion_matrix(exp_sim_combined_bool.experimental, exp_sim_combined_bool.simulated)

        ax = sns.heatmap(cf_matrix, 
                         annot = cf_matrix,
                         fmt="", cmap='winter_r',vmin=0,vmax=200, ax=axs_flat[i])
        ax.set_ylabel("Experimental")
        ax.set_xlabel("Simulated")
        ax.set_yticklabels(["False","True"])
        ax.set_xticklabels(["False","True"])
        ax.set_title(name)
        i +=1
    fig.suptitle(suptitle)
    return fig,axs


