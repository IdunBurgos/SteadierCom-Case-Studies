import pandas as pd


def load_experimental_data():
    """Load experimental data from auxostrophies analysis. The models have been tailored to match the data in Guillen et al. 2021 based on 
    the results from the auxotrophies and minimal medium analysis.
    """
    auxotrophies_ecoli_df = pd.read_csv("../output/Guillen2021/auxotrophies.tsv",sep="\t")

    auxotrophies_ecoli_df.compound = auxotrophies_ecoli_df.compound.map(lambda x: x.replace("R_EX_","M_"))

    mutants = ['mutant_icd', 'mutant_panZ', 'mutant_bioF',
            'mutant_bioH', 'mutant_cysW', 'mutant_bioC', 'mutant_panD']
    
    auxotrophies_ecoli_df = auxotrophies_ecoli_df.loc[auxotrophies_ecoli_df.member.isin(mutants)]
    auxotrophies_ecoli_df = auxotrophies_ecoli_df[auxotrophies_ecoli_df.compound!="M_o2_e"]

    signal_exp = auxotrophies_ecoli_df.set_index(["member","compound"])['0']
    signal_exp  = signal_exp.unstack().fillna(0).stack()
    signal_exp.name="experimental"
    signal_exp = signal_exp*-1 #invert signal to match venturelli output convention (negative uptake in experiments)
    return signal_exp

def load_simulated_data():

    # SteadierCom
    steadiercom_growth = pd.read_csv("../output/Guillen2021/results_steadiercom/results_growth/ecoli_mutants.tsv",sep="\t")
    steadiercom_abundance = pd.read_csv("../output/Guillen2021/results_steadiercom/results_abundance/ecoli_mutants.tsv",sep="\t")
    steadiercom_growth_abundance = pd.read_csv("../output/Guillen2021/results_steadiercom/results_growth_abundance/ecoli_mutants.tsv",sep="\t")

    steadiercom_growth["frequency*mass_rate"] = steadiercom_growth.frequency * steadiercom_growth.mass_rate
    steadiercom_growth["frequency*rate"] = steadiercom_growth.frequency * steadiercom_growth.rate
    steadiercom_growth_abundance["frequency*mass_rate"] = steadiercom_growth_abundance.frequency * steadiercom_growth_abundance.mass_rate
    steadiercom_growth_abundance["frequency*rate"] = steadiercom_growth_abundance.frequency * steadiercom_growth_abundance.rate
    steadiercom_abundance["frequency*rate"] = steadiercom_abundance.frequency * steadiercom_abundance.rate
    steadiercom_abundance["frequency*mass_rate"] = steadiercom_abundance.frequency * steadiercom_abundance.mass_rate


    # SMETANA
    smetana = pd.read_csv("../output/Guillen2021/smetana/_detailed.tsv",sep="\t")

    # SteadyCom
    steadycom_fixed_abundance = pd.read_csv("../output/Guillen2021/steadycom/fixed_abundance/cross_feeding.tsv",sep="\t")
    steadycom = pd.read_csv("../output/Guillen2021/steadycom/no_abundance/cross_feeding.tsv",sep="\t")


    # Organize data
    steadiercom_dict = {
        "SteadierCom - fixed growth":steadiercom_growth,
        "SteadierCom - fixed growth & abundance":steadiercom_growth_abundance,
        "SteadierCom - fixed abundance":steadiercom_abundance
    }

    steadycom_dict = {
        "SteadyCom - Nothing fixed":steadycom,
        "SteadyCom - Fixed  abundance":steadycom_fixed_abundance
    }

    smetana_dict = {"Smetana":smetana}


    # Create simulation dict
    simulation_dict = {"receiver" : {"steadiercom": {key: value[value["receiver"]!="environment"].copy()
                                            for key, value in steadiercom_dict.items()},
                                "steadycom":{key: value[~value["receiver"].isna()].copy()
                                            for key, value in steadycom_dict.items()},
                                    "smetana":{key: value
                                            for key, value in smetana_dict.items()}
                                    },
                        "donor": {"steadiercom": {key: value[value["donor"]!="environment"].copy()
                                            for key, value in steadiercom_dict.items()},
                                "steadycom":{key: value[~value["donor"].isna()].copy()
                                            for key, value in steadycom_dict.items()},
                                    "smetana":{key: value
                                            for key, value in smetana_dict.items()}
                                    },
                        }
    return simulation_dict