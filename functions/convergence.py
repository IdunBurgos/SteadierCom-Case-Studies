import reframed
import cobra
import pandas as pd
import numpy as np
from os import listdir
import warnings
import seaborn as sns
from collections import OrderedDict
import json
import os
import matplotlib.pyplot as plt

import math

from reframed.core.model import ReactionType
from reframed.solvers import solver_instance
from math import inf, isinf
from reframed.community.solution import CommunitySolution
from reframed.solvers.solution import Status
from numpy.random import lognormal
from reframed.core.elements import molecular_weight

import argparse
import textwrap
from reframed.community.model import Community
from reframed.solvers.solution import Status

import glob
import pandas as pd
from reframed.io.cache import ModelCache
from reframed import Environment, set_default_solver
from math import inf

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
                          
                          