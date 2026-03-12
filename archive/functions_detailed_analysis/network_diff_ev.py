import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import functions.old_functions.networks as networks

def process_results(pop,energies):
    energies = np.array(energies)

    pop_up = [pop[0]]
    energies_up = [energies[0]]
    
    for i in range(1, len(pop)):
        new_energies = np.copy(energies_up[-1])
        new_pop = np.copy(pop_up[-1])
        pos = energies[i] < new_energies
        new_energies[pos] = energies[i][pos]
        new_pop[pos] = pop[i][pos]
        pop_up.append(new_pop)
        energies_up.append(new_energies)

    best_evo = pd.DataFrame(pop_up[-1])

    best_evo["penalty"]=energies_up[-1]

    best_evo_sorted = best_evo.sort_values("penalty")

    return best_evo_sorted


def competition_matthew_data(results_env,dataset4,freq_cutoff=0,mass_cutoff=0,cutoff_overlap=0,exp_cutoff1=-1,exp_cutoff2=0):
    
    # Find consumed metabolites - experiments
    consumption = dataset4<=exp_cutoff1
    consumption = consumption.map(int)
    network_exp = consumption.T.dot(consumption)
    network_exp = network_exp.mask(np.tril(np.ones(network_exp.shape)).astype(np.bool)).stack() 
    network_exp = network_exp[network_exp>exp_cutoff2]

    # Cutoff mass_rate*frequency
    results_env = results_env[results_env["mass_rate*frequency"]>mass_cutoff].copy()
    
    # Cutoff frequency and generate network
    network_data = networks.competition_network_data(results_env,freq_cutoff=freq_cutoff)
    
    if (not isinstance(network_data,pd.Series)):
        return np.nan,np.nan,np.nan,np.nan
        
    network_data = network_data[network_data>cutoff_overlap]
    
    # Edges_steadiercom
    edges_steadiercom = network_data.size
    
    # Edges experiment
    edges_exp = network_exp.size
    # Edges matching
    
    network_data.index.names=["r1","r2"]
    match_p1 = set(network_exp.index).intersection(network_data.index)
    network_data2 =network_data.reorder_levels(["r2","r1"])
    
    match_p2 = set(network_exp.index).intersection(network_data2.index)
    
    # True positive = the ones that match (in whichever orde)
    TP = len(set(list(match_p1)+list(match_p2)))
    
    # False negative = The ones that were negative in steadiercom but shouldn't have been
    FN = edges_exp-TP
    
    # False positive = 
    FP = edges_steadiercom-TP
    
    # True negative
    members = dataset4.shape[1]
    TN = (members*(members-1))/2 - TP - FN - FP
    
    assert (TN>=0) and (TP>=0) and (FN>=0) and (FP>=0)
    
    return TP,TN,FP,FN
    

def inner_problem(combination,results_env,dataset4):
    freq_cutoff,mass_cutoff,cutoff_overlap,exp_cutoff1,exp_cutoff2 = combination
    freq_cutoff = 10**freq_cutoff
    mass_cutoff = 10**mass_cutoff

    TP,TN,FP,FN = competition_matthew_data(results_env,dataset4,
                             freq_cutoff=freq_cutoff,
                             mass_cutoff=mass_cutoff,
                             cutoff_overlap=cutoff_overlap,
                             exp_cutoff1=exp_cutoff1,
                             exp_cutoff2=exp_cutoff2)


    denominator = (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)

    if (denominator<=0) or (np.isnan(denominator)):
        return np.inf
    
   
    MCC = ((TP*TN)-(FP*FN))/np.sqrt(denominator)
    if (MCC>1) or (MCC<-1):
        raise Exception(f"Matthew correlation coefficient is out of bounds {MCC}")

    else:

        return 1-MCC