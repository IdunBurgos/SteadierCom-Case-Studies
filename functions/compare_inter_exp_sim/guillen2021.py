
import pandas as pd
import libsbml

import sys
sys.path.append("../functions/")
import compare_inter_exp_sim.correlation as correlation


def load_experimental_results():
    """
    Loading data from experimental results from Guillen et al. 2021 
    """


    exp_interaction = pd.read_csv("../data/Guillen2021/Guillen_interaction.csv",index_col=0)

    mutants = ["icd","cysW","panZ","panD","bioH","bioF","bioC"]

    exp_interaction = exp_interaction.loc[exp_interaction.index.isin(mutants),exp_interaction.columns.isin(mutants)]
    exp_interaction.index = exp_interaction.index.map(lambda x:"mutant_"+x)
    exp_interaction.columns = exp_interaction.columns.map(lambda x:"mutant_"+x)
    signal_exp = exp_interaction.stack()
    signal_exp.name = "experimental" 
    signal_exp.index.names = ['member_second', 'member_main']

    return signal_exp



def _process_simulation_results(steadiercom_dict,steadycom_dict,smetana_dict):

    simulation_dict = {}

    for simulation_method, sim_dict in [("steadiercom",steadiercom_dict),("steadycom",steadycom_dict),("smetana",smetana_dict)]:
        
        simulation_dict[simulation_method]={}

        for name, simulated in sim_dict.items():
            if simulation_method=="steadiercom" or simulation_method=="steadycom":

                # steadycom-specific processing
                if simulation_method=="steadycom":
                    simulated = simulated.fillna("environment")   
                    #simulated.donor = simulated.donor.map(lambda x: x.replace("M_","") if isinstance(x,str) else x)
                    #simulated.receiver = simulated.receiver.map(lambda x: x.replace("M_","") if isinstance(x,str) else x)

                # steadiercom-specific processing
                if simulation_method=="steadiercom":
                    assert ("rate" in simulated.columns) and ("mass_rate" in simulated.columns) and ("frequency" in simulated.columns), f"This is not standard SteadierCom results. Columns: {simulated.columns}"
                    simulated["frequency*mass_rate"] = simulated.frequency*simulated.mass_rate 
                    simulated["frequency*rate"] = simulated.frequency*simulated.rate

                simulated = simulated[ (simulated.receiver!="environment")].copy()

            #simulated.donor = simulated.donor.map(lambda x: x.split("_")[0][0].upper()+x.split("_")[1][0].upper() if x!="environment" else x)
            #simulated.receiver = simulated.receiver.map(lambda x: x.split("_")[0][0].upper()+x.split("_")[1][0].upper())
        
            simulation_dict[simulation_method][name]=simulated

    return simulation_dict

    
def load_simulation_results():
    """
    Loads experimental data from all simulation methods and 
    """

    steadiercom_abundance = pd.read_csv("../output/Guillen2021/results_steadiercom/results_abundance/ecoli_mutants.tsv",sep="\t")

    steadiercom_growth = pd.read_csv("../output/Guillen2021/results_steadiercom/results_growth/ecoli_mutants.tsv",sep="\t")

    steadiercom_growth_abundance = pd.read_csv("../output/Guillen2021/results_steadiercom/results_growth_abundance/ecoli_mutants.tsv",sep="\t")

    smetana = pd.read_csv("../output/Guillen2021/smetana/_detailed.tsv",sep="\t")

    steadycom_fixed_abundance = pd.read_csv("../output/Guillen2021/steadycom/fixed_abundance/cross_feeding.tsv",sep="\t")
    steadycom = pd.read_csv("../output/Guillen2021/steadycom/no_abundance/cross_feeding.tsv",sep="\t")
    
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

    simulation_dict = _process_simulation_results(steadiercom_dict,steadycom_dict,smetana_dict)

    return simulation_dict




def find_correlations(simulation_dict,signal_exp):

    results_corr = []
    results_corr_member = []

    sim_method_metric = [("steadycom","rate?"),("steadiercom","frequency*mass_rate"),("steadiercom","frequency"),("smetana","scs")]

    for simulation_method,metric in sim_method_metric:
        print(simulation_method)
        results_corr.append(correlation.find_correlation_simulation_method(simulation_dict,
                                                            signal_exp,
                                                            simulation_method=simulation_method,
                                                            metric=metric,
                                                            coop_or_comp_list=["coop","coop-comp"]
                                                            )) 

        results_corr_member.append(correlation.find_correlation_simulation_method_detailed(simulation_dict,
                                                            signal_exp,
                                                            simulation_method=simulation_method,
                                                            metric=metric,
                                                            coop_or_comp_list=["coop","coop-comp"]
                                                            ))

    results_corr_df = pd.concat(results_corr)
    results_corr_member_df = pd.concat(results_corr_member)


    return results_corr_df,results_corr_member_df

