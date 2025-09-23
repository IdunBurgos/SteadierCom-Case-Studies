import pandas as pd
import time
from reframed import SteadyCom,Community

def build_constraints_from_abundance(abundance_df,prefix="x_M_"):
    
    constraints_gut = {
    cid: dict(zip((prefix + g["member"].astype(str)).values, g["abundance"].values))
    for cid, g in abundance_df.groupby("community_id", sort=False)
    }
    return constraints_gut


def SteadyComCommunities(models,constraints,env,fixed_relative_abundance=True):    
    cross_feeding_comms = {}
    relative_abundance_comms = {}
    growth_comms = {}

    for community_id,relative_abundance_dict in constraints.items():
        print(f"Started SteadyCom analysis for {community_id}")
        # Prepare for SteadyCom
        model_ids = [id_.replace("x_M_","").replace("x_","") for id_ in relative_abundance_dict.keys()]
        models_specific = [model for id_, model in models.items() if id_ in model_ids]

        # Make community and apply environment
        comm = Community(community_id=community_id,models=models_specific)
        env.apply(comm.merged_model,inplace=True,exclusive=True,warning=False)
        
        # Run SteadyCom
        time_1 = time.time()
        if fixed_relative_abundance:

            # Normalize relative abundance
            s = pd.Series(relative_abundance_dict)
            relative_abundance_dict = (s / s.sum()).to_dict()

            # Solve SteadyCom problem
            sol_steadycom = SteadyCom(comm,constraints=relative_abundance_dict)
        else:
            sol_steadycom = SteadyCom(comm)
        time_2 = time.time()
        print(f"\tRunning SteadyCom took: {(time_2 - time_1)/60:.2f} minutes")
        
        # Gather data
        cross_feeding = pd.DataFrame(sol_steadycom.cross_feeding())
        cross_feeding.columns = ["donor","receiver","compound","rate?","mass_rate?"]
        relative_abundance = sol_steadycom.abundance
        growth_gut = sol_steadycom.growth
        
        # Add data to dicts
        cross_feeding_comms[community_id] = cross_feeding
        relative_abundance_comms[community_id] = relative_abundance
        growth_comms[community_id] = growth_gut
        
        time_3 = time.time()
        print(f"\tCollecting results took: {(time_3 - time_2)/60:.2f} minutes")

    # Turn into DataFrames
    cross_feeding_comms_df = pd.concat(cross_feeding_comms).reset_index(level=0,names="community")
    relative_abundance_comms_df = pd.DataFrame.from_dict(relative_abundance_comms)
    growth_comms_df = pd.DataFrame.from_dict(data=growth_comms,orient='index',columns=["growth"])

    return cross_feeding_comms_df,relative_abundance_comms_df,growth_comms_df