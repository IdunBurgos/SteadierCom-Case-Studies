import json
import reframed
import pandas as pd
import json

import logging
logging.basicConfig(
    level=logging.INFO,  # change to DEBUG for more detail
    format="%(levelname)-8s | %(message)s",
)

# Import reaction dataframe
EGC_rxns_bigg = pd.read_excel("../data/EGC.xlsx",
                         sheet_name="Sheet2",
                         usecols="C:E")

EGC_rxns_agora = pd.read_excel("../data/EGC AGORA.xlsx",
                         sheet_name="Sheet2",
                         usecols="C:E")

# Make reframed reactions

        
        

def EGC_identifier(model,print_results=False,agora=False,send_log=False):
    if agora:
        EGC_rxns = EGC_rxns_agora.copy()
    else:
        EGC_rxns = EGC_rxns_bigg.copy()
        
        
    # Make reframed reactions
    rxns={}
    for index, row in EGC_rxns.iterrows():
        if not row["Exists"]:
            reaction_id = row['rxn_ID ']
            reversible = True
            stoichiometry = json.loads(row['stoichiometry'])
            rxn = reframed.CBReaction(reaction_id=reaction_id, reversible=reversible, stoichiometry=stoichiometry)
            rxns[row['rxn_ID ']]=rxn
    
    EGCs = {}
    model_copy= model.copy()
    
    # Set lower bound to 0 to avoid infeasible solution because of NGAM
    if agora:
        model_copy.reactions["R_DM_atp_c_"].lb=0
    else:
        model_copy.reactions["R_ATPM"].lb=0
        
    log = ""
    
    for index, row in EGC_rxns.iterrows():

        if not row["Exists"]:
            
            reaction_compounds = rxns[row["rxn_ID "]].get_substrates()+rxns[row["rxn_ID "]].get_products()
            has_compound = []

            for met in reaction_compounds:
                has_compound.append(met in model_copy.metabolites)
            
            if len(reaction_compounds)==sum(has_compound):
                model_copy.add_reaction(rxns[row["rxn_ID "]])
            else:
                continue
                
        for rxn in model_copy.reactions:
            if model_copy.reactions[rxn].lb<0:
                model_copy.reactions[rxn].lb=-1000
            if model_copy.reactions[rxn].ub>0:
                model_copy.reactions[rxn].ub=1000
                
        # Set boundary to avoid infeasible problem
        model_copy.reactions[row["rxn_ID "]].lb=0
        model_copy.reactions[row["rxn_ID "]].ub=1000
        

        # Set environment to empty
        env_empty = reframed.Environment.empty(model_copy)

        # Set objective as the reaction in question
        objective= {rxn:0 for rxn in model_copy.reactions}
        objective[row["rxn_ID "]]=1

        # Solve pFBA (not FBA to avoid finding energy-balanced cycles)
        sol_pfba = reframed.pFBA(model_copy,objective=objective,constraints=env_empty)

        # If flux through objective function is 0 -> no infeasible cycle for this compound
        # -> continue to next cycle
        
        if abs(sol_pfba.fobj)<1e-7 :
            if print_results:
                log = log+'There are NO energy producing cycles in the model for '+str(row["rxn_ID "])+"\n"
            continue
        
        if print_results: 
            log = log+'There are energy producing cycles in the model for '+str(row["rxn_ID "]) + "\n"
        
        EGCs[row["rxn_ID "]]=[]

        # For each reaction in the model I'm finding the associated rxn with flux.
        for rxn,value in sol_pfba.values.items():
            if abs(value)<1e-7:
                continue
            if print_results: print("\t" + str(rxn)+": " + str(value))
            
            EGCs[row["rxn_ID "]].append((rxn,value))
        
        # Set upper bound to 0 to avoid flux in the nex cycle. 
        model_copy.reactions[row["rxn_ID "]].ub=0
    
    message = "\033[92mThere are NO energy producing cycles in the model\033[0m  \n" if len(EGCs) == 0 else "\033[91mThere ARE energy producing cycles in the model\033[0m \n"
    
    log = log+ message
    
    if send_log:
        return EGCs,log
    else:
        return EGCs