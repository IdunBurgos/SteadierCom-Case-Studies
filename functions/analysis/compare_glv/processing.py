import pandas as pd
import numpy as np
import scipy as sp
import seaborn as sns

import warnings

import logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

def _cooperation_agg(simulation_res, metric="rate"):
    """

    Finds fraction of flux going to member_main coming from member_second. 

    Details: Finds the norm of flux for each member and calculate the sum of normalized flux coming to member_main
    from member_second. 

    Returns a Series with the cooperation scores.
    """

    df = simulation_res.copy()

    # total inflow per receiver (guard zeros)
    total = df.groupby("receiver")[metric].transform("sum")
    total_safe = total.replace(0, np.nan)

    # normalized contribution to A
    df["norm_flux"] = df[metric] / total_safe
    df["norm_flux"] = df["norm_flux"].fillna(0.0)

    # drop environment as a 'cooperating member'
    df = df[df["donor"] != "environment"]

    # sums all fluxes in donor-receiver relationship
    # TO-DO: consider if aggregating at this point is necessary (I guess it is because it's showing the sum and not the average.)
    coop = df.groupby(["donor", "receiver"], as_index=True)["norm_flux"].sum()
    coop.index.names = ["member_second", "member_main"]
    coop.name = "simulated"
    return coop


def _competition_agg(simulation_res, metric="rate", met_column="compound", method=1):

    """
    Finds the competition scores between members.

    Returns a Series with the competition scores.
    """
    df = simulation_res.copy()

    # absolute consumption: member x compound
    R = (
        df.groupby(["receiver", met_column])[metric]
          .sum()
          .unstack(fill_value=0.0)
    )

    # normalized consumption within each member
    M = R.div(R.sum(axis=1), axis=0)

    X = M.values

    if method == 1:
        # weighted binary overlap
        U = (M > 1e-6).astype(int)
        C = X @ U.values.T

    elif method == 2:
        # cosine similarity
        dot = X @ X.T
        norm = np.linalg.norm(X, axis=1)
        C = dot / np.outer(norm, norm)

    elif method == 3:
        # community-turnover-weighted overlap
        turnover = R.sum(axis=0)

        C = np.zeros((len(R), len(R)))

        for i in range(len(R)):       # A
            for j in range(len(R)):   # B

                A = R.iloc[i]
                B = R.iloc[j]

                shared = (A > 1e-6) & (B > 1e-6)

                compound_scores = (A + B) / turnover
                compound_scores[~shared] = 0

                n_B = (B > 1e-6).sum()

                C[i, j] = compound_scores.sum() / n_B

    np.fill_diagonal(C, np.nan)

    comp = (
        pd.DataFrame(C, index=M.index, columns=M.index)
          .stack()
          .rename("simulated")
    )

    comp.index.names = ["member_main", "member_second"]
    comp = comp.swaplevel(0, 1).sort_index()

    return comp


def merge_exp_sim(simulated_res,signal_exp,coop_or_comp="coop",metric="frequency",verbose=True,method_comp=1):

    """
   coop_or_comp is either 'coop', 'comp' or 'coop-comp'

   Returns a DataFrame with the merged experimental and simulated results.
    """
    def _merge_exp_sim_helper(simulation_res,signal_exp, coop_or_comp="coop",metric="frequency"):

        assert metric in simulation_res.columns, f"{metric} is not in the columns of simulation_res"
        
        # Get score for cooperation by sim and exp
        if coop_or_comp=="coop":
            simulated_score = _cooperation_agg(simulation_res,
                                                    metric=metric,
                                                    )
        elif coop_or_comp=="comp":
            simulated_score = _competition_agg(simulation_res,
                                        metric=metric,
                                        method=method_comp
                                        )
        else:
            raise Exception(f"coop_or_comp={coop_or_comp}, but must be either 'coop' or 'comp'")

        # Combine simulated score with experimental values
        exp_sim_combined = pd.concat([signal_exp,simulated_score],axis=1)   

        exp_sim_combined = exp_sim_combined[~exp_sim_combined.experimental.isna()].copy()

        return exp_sim_combined

    if coop_or_comp=="coop-comp":
        exp_sim_combined_coop =  _merge_exp_sim_helper(simulated_res,
                                    signal_exp,
                                    coop_or_comp="coop",
                                    metric=metric,
                                    ).fillna(0)
        exp_sim_combined_comp =  _merge_exp_sim_helper(simulated_res,
                                    signal_exp,
                                    coop_or_comp="comp",
                                    metric=metric,
                                    method=method_comp
                                    ).fillna(0)

        
        exp_sim_combined = exp_sim_combined_coop.copy()
        exp_sim_combined.simulated= exp_sim_combined_coop.simulated-exp_sim_combined_comp.simulated

        # Remove simulated data that was 0 - means essentially now prediction of neither competition nor cooperation (fillna(0))
        isna = (exp_sim_combined_coop.simulated<1e-6) & (exp_sim_combined_comp.simulated<1e-6)
        
    else:
        exp_sim_combined = _merge_exp_sim_helper(simulated_res,
                                        signal_exp,
                                        coop_or_comp=coop_or_comp,
                                        metric=metric,
                                        )

        # Remove entries that do not have a match...
        isna = exp_sim_combined.simulated.isna()

    exp_sim_combined_isna = exp_sim_combined.copy()
    exp_sim_combined_isna["Was simulated"] = ~isna

    return exp_sim_combined_isna



def process_simulation_results(simulation_res_all, glv_data,method_comp=1):
    """
    Combine experimental and simulated results across all communities. 
    """

    df_merged_list = []

    # Find the competition and cooperation results for each community
    for community, simulated_res in simulation_res_all.groupby("community"):

        # Find the cooperation results
        df_coop = merge_exp_sim(simulated_res,
                                                glv_data,
                                                coop_or_comp="coop",
                                                metric="mass_rate*frequency",
                                                verbose=True
                                                )
        # Find the competition results
        df_comp = merge_exp_sim(simulated_res,
                                                glv_data,
                                                coop_or_comp="comp",
                                                metric="mass_rate*frequency",
                                                verbose=True,
                                                method_comp=method_comp)
        df_merged = pd.merge(df_coop, df_comp.drop("experimental",axis=1),left_index=True,right_index=True, suffixes=("_coop", "_comp"))
        df_merged.loc[:,"community"] = community
        df_merged = df_merged.set_index("community",append=True)

        df_merged_list.append(df_merged)

    df_merged_unfiltered = pd.concat(df_merged_list,axis=0)
    return df_merged_unfiltered

def make_glv_categories(df_merged_unfiltered, glv_cutoff=1):
    df_merged = df_merged_unfiltered.copy().dropna() # NaN values is only for the A-A combinations and interactions with PC (not growing!)
    df_merged["exp_categories"] = pd.cut(df_merged["experimental"], bins=[-np.inf, -glv_cutoff, glv_cutoff, np.inf], labels=["negative", "neutral", "positive"])
    return df_merged
