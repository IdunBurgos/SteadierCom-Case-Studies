def combine_exp_sim_venturelli_met(simulated_results,signal_exp,metid2group,metric="frequency",receiver_or_donor="receiver",agg_method="sum"):
    """
    Input:
    simulated_results (pd.DataFrame): The raw data from steadiercom, steadycom, or smetana simulation.
    signal_exp (pd.Series): The experimental data, where with two indeces ("member","group"). This also has to be the names of the indeces.
    metid2group (dict): Dict. for mapping between met. id and met. group
    metric (str): Determines what metric in the simulation result to study (ex. frequency, rate, mass_rate)
    receiver_or_donor (str): 'receiver' or 'donor
    agg_method (str): how to aggregate signals in simulation results (ex. sum, mean, median)
    
    Output:
    exp_sim_combined (pd.DataFrame): Dataframe with experimental and simulated results combined.
    
    """
    
    if len(signal_exp.index.names)!=2:
        raise Exception("signal_exp must have a double index (member,group)")
        
    if signal_exp.index.names[0]!="member" or signal_exp.index.names[1]!="group":
        raise Exception("The names of the double index in signal_exp must be (member,group)")
    
    if receiver_or_donor!="receiver" and receiver_or_donor!="donor":
        raise Exception("receiver_or_donor must be either 'receiver' or 'donor'")
    
    valid_methods = ['sum', 'mean', 'median']
    if agg_method not in valid_methods:
        raise Exception(f"agg_method must be one of {valid_methods}")
    
    # Map met id to met group
    simulated_res_reduced = simulated_results[simulated_results.compound.isin(metid2group.keys())].copy()
    simulated_res_reduced["metname_group"] = simulated_res_reduced.compound.map(lambda x: metid2group[x])
    
    # Group by receiver/donor and metname_group -> sum
    
    simulated_res_grouped = simulated_res_reduced[[receiver_or_donor,"metname_group",metric]]\
                                                        .groupby([receiver_or_donor,"metname_group"])\
                                                        .agg(agg_method)\
                                                        .reset_index()
    simulated_res_grouped = simulated_res_grouped[simulated_res_grouped[receiver_or_donor]!="environment"]
    simulated_res_grouped[receiver_or_donor] = simulated_res_grouped[receiver_or_donor]\
                                                        .map(lambda x: x.split("_")[0][0].upper()+x.split("_")[1][0].upper())
    simulated_res_grouped = simulated_res_grouped.set_index([receiver_or_donor,"metname_group"])[metric]
    simulated_res_grouped.index.names = ["member","group"]
    simulated_res_grouped.name = "simulated"

    # Combine results    
    exp_sim_combined = pd.concat([simulated_res_grouped,signal_exp],axis=1)
    
    return exp_sim_combined


def scoring_met(simulation_res_dict,signal_exp,exp_thresh=0.5,sim_thresh=0.5,sim_metric="frequency",receiver_or_donor="frequency"):
    
    """
    Input:
    simulation_res_dict (dict(str,pd.DataFrame)): A dict with raw data from steadiercom, steadycom, or smetana simulation.
    signal_exp (pd.Series): The experimental data, where with two indeces ("member","group"). This also has to be the names of the indeces.
    exp_thresh (float): The threshold for experimental results. Anything above/below will be True/False.
    sim_thresh (float): The threshold for simulation results. Anything above will be True.
    sim_metric (str): Determines what metric in the simulation result to study (ex. frequency, rate, mass_rate)
    
    receiver_or_donor (str): 'receiver' or 'donor
    agg_method (str): how to aggregate signals in simulation results (ex. sum, mean, median)
    
    Output:
    scores_df_all (pd.DataFrame): Dataframe with results for all elements in simulation_res_dict. Columns correspond to keys in simulation_res_dict.
    """
    
    
    scores_list = []
    
    # For each result
    for key,simulation_res in simulation_res_dict.items():
        
        # Find average receiving/donating quantity
        exp_sim_combined = combine_exp_sim_venturelli_met(simulation_res,
                                           signal_exp,
                                           metric=sim_metric,
                                           receiver_or_donor=receiver_or_donor,
                                           agg_method=agg_method)
        # Make boolean output
        ## Determine positive (True) or negative (False) results
        
        if receiver_or_donor=="receiver":
            exp_sim_combined.experimental = exp_sim_combined.experimental<exp_thresh # Consumption is negative
        elif receiver_or_donor=="donor":
            exp_sim_combined.experimental = exp_sim_combined.experimental>exp_thresh # Production is positive
        else:
            raise Exception("receiver_or_donor must be either 'receiver' or 'donor'")
            
            
        # Score experiment vs simulation overlap
        scores = classification_report(exp_sim_combined.experimental,exp_sim_combined.simulated,output_dict=True)
        scores_df = pd.DataFrame(scores)
        
        # Extract True class (positive class)
        scores_df_copy1 = scores_df["True"].copy()
        scores_df_copy1 = scores_df_copy1.drop("support")
        scores_df_copy1.index = ["precision","sensitivity","f1-score"]
        
        # Extract False class (negative class)
        scores_df_copy2 = scores_df["False"].copy()
        scores_df_copy2 = scores_df_copy2.drop(["support","precision","f1-score"])
        scores_df_copy2.index = ["specificity"]
        
        # Combine results
        scores_df_copy = pd.concat([scores_df_copy1,scores_df_copy2])
        scores_df_copy["accuracy"] = scores["accuracy"]
        scores_df_copy.name = key

        scores_list.append(scores_df_copy)

    scores_df_all = pd.concat(scores_list,axis=1)
    scores_df_all.sort_index(inplace=True)
    return scores_df_all