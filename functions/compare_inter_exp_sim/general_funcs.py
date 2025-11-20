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
    coop = df.groupby(["donor", "receiver"], as_index=True)["norm_flux"].sum()
    coop.index.names = ["member_second", "member_main"]
    coop.name = "simulated"
    return coop


def _competition_agg(simulation_res, metric="rate", met_column="compound"):
    
    """
    Finds the fraction of flux to member_main that member_second is competing for. The magnitude
    of flux going to member_second is not 
    """
    df = simulation_res.copy()

    # normalize inflow per receiver (guard zeros)
    total = df.groupby("receiver")[metric].transform("sum")
    total_safe = total.replace(0, np.nan)
    df["norm_flux"] = df[metric] / total_safe
    df["norm_flux"] = df["norm_flux"].fillna(0.0)

    # receiver x compound matrix of normalized inflow
    M = (df.groupby(["receiver", met_column], as_index=True)["norm_flux"]
            .sum()
            .unstack(fill_value=0.0))

    U = (M > 1e-6).astype(int)

    C = M.values @ U.values.T
    np.fill_diagonal(C, np.nan)

    comp = (pd.DataFrame(C, index=M.index, columns=M.index)
              .stack()
              .rename("simulated"))

    # index levels: (member_second=B competitor, member_main=A affected)
    comp.index.names = ["member_main", "member_second"]

    comp = comp.swaplevel(0, 1).sort_index()
    return comp


def merge_exp_sim(simulated_res,signal_exp,coop_or_comp="coop",metric="frequency",verbose=True):

    """
   coop_or_comp is either 'coop', 'comp' or 'coop-comp'
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
                                    ).fillna(0)

        
        exp_sim_combined = exp_sim_combined_coop.copy()
        exp_sim_combined.simulated= exp_sim_combined_coop.simulated-exp_sim_combined_comp.simulated

        # Remove simulated data that was 0 - means essentially now prediction of neitehr competition nor cooperation (fillna(0))
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
