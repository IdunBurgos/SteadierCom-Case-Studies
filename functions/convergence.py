import cobra
import pandas as pd
import numpy as np
import math
from numpy.random import lognormal
import argparse
import textwrap

# Reframed
from reframed.community.solution import CommunitySolution
from reframed.community.model import Community
from reframed.solvers.solution import Status

# Steadiercom
import steadiercom.steadiercom as steadiercom
from steadiercom.steadiercom import allocation_constraints
from steadiercom.cli import load_communities,load_media_db,precompute_exchange_map
from steadiercom.steadiercom import build_problem 
import time

def SteadierSampleConvergence(community,output_dir, objective=None,n=100, growth=0.1, abundance=None, allocation=False, constraints=None,
             w_e=0.002, w_r=0.2, solver=None):

    if abundance:
        norm = sum(abundance.values())
        abundance = {org_id: abundance.get(org_id, 0) / norm for org_id in community.organisms}

    if not solver:
        solver = build_problem(community, growth=growth, abundance=abundance)

    if not allocation:
        w_e, w_r = 0, 0

    random_vars = allocation_constraints(community, solver, w_e=w_e, w_r=w_r, abundance=abundance)

    sols = []
    if abundance:
        w1 = {org_id: 1 / abundance[org_id] if abundance[org_id] > 0 else 0 for org_id in community.organisms}


    #if growth is None and objective is None:
    objective = community.merged_model.biomass_reaction

    external_metabolites, exchange_map = precompute_exchange_map(community)
    
    dicts={}
    j=0
    
    time1 = time.time()
    for i in range(n):
        
        if not abundance:
            w1 = {org_id: lognormal(0, 1) for org_id in community.organisms}

        objective = {vi: w1[org_id] * lognormal(0, 1) for org_id, v_org in random_vars.items() for vi in v_org}

        sol = solver.solve(objective, minimize=True, constraints=constraints)

        if sol.status == Status.OPTIMAL:
            if sol.values[community.merged_model.biomass_reaction] > 0:
                sol = CommunitySolution(community, sol.values)
                dicts[i]=sol.compute_exchanges(external_metabolites, exchange_map)
                
        if np.mod(i,200)==0 and i!=0:
            sol_df=pd.DataFrame.from_dict(dicts)
            sol_df.to_csv(output_dir+"sol_df"+str(j)+".tsv",sep="\t")
            j+=1
            del dicts
            dicts = {} 
            
            time2 = time.time()
            print(f"Round {i}. Time passed: {(time2-time1)/60:.2f} minutes")
            time1 = time2
            

            
            

def std_sample_size(dfs_tot_abs,sample_sizes,cutoff=1e-6 ):
    std = lambda x:np.std(x,ddof=1)
    cv = lambda x: std(x)/ abs(np.mean(x)) * 100 

    datas_std = {}
    datas_cv = {}
    for n in sample_sizes:
        sample_size = int(n)

        number = math.floor(dfs_tot_abs.shape[1]/sample_size)
        data = pd.DataFrame.from_dict({f"set{i}":dfs_tot_abs[np.arange(0+i*sample_size,(i+1)*sample_size)].T.mean() for i in np.arange(0,number)})
    
        # Cutoff data that is actually zero
        equal_zero= data.sum(axis=1).apply(abs).sort_values() <cutoff 
        data = data.loc[~equal_zero]

        data_std = data.T.apply(std)
        datas_std[sample_size]=data_std 
        data_cv = data.T.apply(cv)
        datas_cv[sample_size]=data_cv 
    return datas_std,datas_cv

def std_sampling_grouped(dfs_tot_abs,sample_sizes):
    std = lambda x:np.std(x,ddof=1)
    cv = lambda x: std(x)/ abs(np.mean(x)) * 100 

    statistics = {}
    statistics_cv = {}
    for n in sample_sizes:
        b = -6

        statistics[n]={}
        statistics_cv[n]={}
        for a in np.linspace(-5.5,1,14):

            sample_size = int(n)
            
            # Number of sets it's possible to make with same number of samples
            number = math.floor(dfs_tot_abs.shape[1]/sample_size)

            # For sample sizes of x calculate the mean value in each group of samples
            data = pd.DataFrame.from_dict({f"set{i}":dfs_tot_abs[np.arange(0+i*sample_size,(i+1)*sample_size)].T.mean() for i in np.arange(0,number)})
            
            # Remove variables with average frequency below or above thresholds
            data_set_lower= data.mean(axis=1).apply(abs).sort_values() <10**(a)
            data_set_upper = data.mean(axis=1).apply(abs).sort_values() >10**(b)
            data_bound = data_set_lower * data_set_upper
            data = data.loc[data_bound]
            
            data_std = data.T.apply(std)
            mean_std = data_std.mean()

            data_cv = data.T.apply(cv)
            mean_cv = data_cv.mean()

            statistics[n][(round(10**b,6),round(10**a,6))] =mean_std
            statistics_cv[n][(round(10**b,6),round(10**a,6))] =mean_cv
            b=a
    return statistics,statistics_cv


def std_sampling_grouped_freq(dfs_tot_abs,sample_sizes):
    std = lambda x:np.std(x,ddof=1)
    cv = lambda x: std(x)/ abs(np.mean(x)) * 100 
    
    statistics = {}
    statistics_cv = {}
    for n in sample_sizes:
        b = 0

        statistics[n]={}
        statistics_cv[n]={}
        for a in np.linspace(0.1,1,10):

            sample_size = int(n)
            number = math.floor(dfs_tot_abs.shape[1]/sample_size)

            # For sample sizes of x calculate the mean value in each group of samples
            data = pd.DataFrame.from_dict({f"set{i}":dfs_tot_abs[np.arange(0+i*sample_size,(i+1)*sample_size)].T.mean() for i in np.arange(0,number)})
            
            # For sample sizes of x calculate the frequency of the observation of values above 1e-6
            data_freq = pd.DataFrame.from_dict({f"set{i}":(dfs_tot_abs[np.arange(0+i*sample_size,(i+1)*sample_size)].T>1e-6).sum()/sample_size for i in np.arange(0,number)})
            
            # Remove variables with average frequency below or above thresholds
            data_set_lower= data_freq.mean(axis=1).sort_values() <a
            data_set_upper = data_freq.mean(axis=1).sort_values() >b
            data_bound = data_set_lower * data_set_upper
            data = data.loc[data_bound]

            # Calculate std and cv
            data_std = data.T.apply(std)
            mean_std = data_std.mean()

            data_cv = data.T.apply(cv)
            mean_cv = data_cv.mean()
            
            statistics[n][(round(b,2),round(a,2))] =mean_std
            statistics_cv[n][(round(b,2),round(a,2))] =mean_cv
            b=a
    return statistics,statistics_cv
                          
                          