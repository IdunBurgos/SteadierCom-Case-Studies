import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, 
    format="%(levelname)-8s | %(message)s",
)

def load_experimental_abundance():
    """
    Returns a pandas DataFrame with the experimental abundances from the Guillen2021 study
        The columns are: community, member, relative_abundance
    """
    experimental_abundance = pd.read_csv("../data/Guillen2021/generated/communities_steadiercom.tsv",sep="\t",header=None)
    experimental_abundance.columns = ["community", "member", "relative_abundance"]
    return experimental_abundance

def load_experimental_abundance_replicas(experimental_abundance):
    exp_abundance_replicas = pd.read_csv("../data/Guillen2021/relative_abundance.csv",index_col=0)
    # NB! yhhK in data == panZ in paper
    exp_abundance_replicas.index =exp_abundance_replicas.index.map(lambda x:"mutant_"+x if x!="yhhK" else "mutant_panZ")
    # Reduce to only top 7 members
    exp_abundance_replicas = exp_abundance_replicas[exp_abundance_replicas.index.isin(experimental_abundance.member)]
    exp_abundance_replicas = (exp_abundance_replicas
                            .reset_index(names=["member"])
                            .melt(id_vars=["member"],
                                var_name="community",
                                    value_name="relative_abundance",
                                    ))
    
    return exp_abundance_replicas

def load_simulated_abundance():
    # Load SteadierCom abundances
    steadiercom_abundance_df = pd.read_csv(f"../output/Guillen2021/results_steadiercom/results_growth/ecoli_mutants_abundance.tsv",sep="\t")
    steadiercom_abundance_df.reset_index(drop=True)

    steadiercom_abundance_df = (steadiercom_abundance_df.pivot_table(index=["community","sample"], columns='member', values='relative_abundance', fill_value=0)
        .stack()
        .reset_index())
    steadiercom_abundance_df.rename(columns={0:"relative_abundance"}, inplace=True)

    # Load SteadyCom abundances
    steadycom_abundance_df = pd.read_csv("../output/Guillen2021/steadycom/no_abundance/relative_abundance.tsv",sep="\t",index_col=0)
    steadycom_abundance_df = steadycom_abundance_df.T
    steadycom_abundance_df.sort_index(inplace=True)
    steadycom_abundance_df.sort_index(inplace=True,axis=1)
    steadycom_abundance_reshaped = (
        pd.melt((steadycom_abundance_df.fillna(0).rename_axis("community").reset_index()),
            id_vars="community",
            var_name="member",
            value_name="relative_abundance")
            )

    return steadiercom_abundance_df, steadycom_abundance_reshaped

def _process_simulation_results(steadiercom_dict,steadycom_dict,smetana_dict, verbose=True):
    # almost SAME AS FOR VENTURELLI - CONSIDER MOVING TO SHARED SPACE
    if verbose:
        logging.info(f"_process_simulation_results cleans up namespace for the simulation results.")
    simulation_dict = {}

    for simulation_method, sim_dict in [("steadiercom",steadiercom_dict),("steadycom",steadycom_dict),("smetana",smetana_dict)]:
        
        simulation_dict[simulation_method]={}

        for name, simulated in sim_dict.items():
            if simulation_method=="steadiercom" or simulation_method=="steadycom":

                # steadycom-specific processing
                if simulation_method=="steadycom":
                    simulated = simulated.fillna("environment")   
                    simulated.donor = simulated.donor.map(lambda x: x.replace("M_","") if isinstance(x,str) else x)
                    simulated.receiver = simulated.receiver.map(lambda x: x.replace("M_","") if isinstance(x,str) else x)

                # steadiercom-specific processing
                if simulation_method=="steadiercom":
                    assert ("rate" in simulated.columns) and ("mass_rate" in simulated.columns) and ("frequency" in simulated.columns), f"This is not standard SteadierCom results. Columns: {simulated.columns}"
                    simulated["mass_rate*frequency"] = simulated.frequency*simulated.mass_rate 
                    simulated["mass_rate*frequency"] = simulated.frequency*simulated.rate
            
            simulation_dict[simulation_method][name]=simulated

    return simulation_dict

def load_steadystate_simulation_results():
    logging.info(f"Loading simulation results for Guillen2021 study for smetana, steadycom, and steadiercom.")
    # SteadierCom
    steadiercom_growth = pd.read_csv("../output/Guillen2021/results_steadiercom/results_growth/ecoli_mutants.tsv",sep="\t")
    steadiercom_abundance = pd.read_csv("../output/Guillen2021/results_steadiercom/results_abundance/ecoli_mutants.tsv",sep="\t")
    steadiercom_growth_abundance = pd.read_csv("../output/Guillen2021/results_steadiercom/results_growth_abundance/ecoli_mutants.tsv",sep="\t")

    # SMETANA
    smetana = pd.read_csv("../output/Guillen2021/smetana/_detailed.tsv",sep="\t")

    # SteadyCom
    steadycom_fixed_abundance = pd.read_csv("../output/Guillen2021/steadycom/fixed_abundance/cross_feeding.tsv",sep="\t")
    steadycom = pd.read_csv("../output/Guillen2021/steadycom/no_abundance/cross_feeding.tsv",sep="\t")

    # Organize data
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