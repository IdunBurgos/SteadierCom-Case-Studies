
import pandas as pd
import libsbml


def load_experimental_results():
    """
    Loading data from gLV-model 
    """
    
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



def _process_simulation_results(steadiercom_dict,steadycom_dict,smetana_dict):

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
                    simulated["frequency*mass_rate"] = simulated.frequency*simulated.mass_rate 
                    simulated["frequency*rate"] = simulated.frequency*simulated.rate

                simulated = simulated[ (simulated.receiver!="environment")].copy()

                simulated.donor = simulated.donor.map(lambda x: x.split("_")[0][0].upper()+x.split("_")[1][0].upper() if x!="environment" else x)
                simulated.receiver = simulated.receiver.map(lambda x: x.split("_")[0][0].upper()+x.split("_")[1][0].upper())
            
            simulation_dict[simulation_method][name]=simulated

    return simulation_dict

    
def load_simulation_results():
    """
    Loads experimental data from all simulation methods and 
    """

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

    simulation_dict = _process_simulation_results(steadiercom_dict,steadycom_dict,smetana_dict)

    return simulation_dict

