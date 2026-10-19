import pandas as pd
import numpy as np
import sys
sys.path.append('../functions')

#import translation_dicts
#import colors_MAGs as color_func
import logging

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)



def data_ReceiverOrDonor(steadiercom_samples_preprocessed,ReceiverOrDonor,metric="mass_rate*frequency",only_media=False):
    """
    Aggregate production/consumption data for each receiver/donor based on super_class for bar plots.

    Parameters
    ----------
    steadiercom_samples_preprocessed : pd.DataFrame
        Preprocessed SteadierCom samples DataFrame.
    ReceiverOrDonor : str
        Either 'receiver' or 'donor'.
    metric : str, optional
        The metric to use for the analysis (default is "mass_rate*frequency").
    only_media : bool, optional
        Whether to consider only media data (default is False).
    """
    if metric not in steadiercom_samples_preprocessed.columns:
        raise ValueError(f"Metric {metric} not found in steadiercom_samples_preprocessed columns: {steadiercom_samples_preprocessed.columns}")

    if ReceiverOrDonor not in {"receiver", "donor"}:
        raise ValueError(
            f"ReceiverOrDonor must be 'receiver' or 'donor', got {ReceiverOrDonor}"
        )

    # find members
    steadiercom_only_producers = steadiercom_samples_preprocessed[(steadiercom_samples_preprocessed.donor!="environment")]
    steadiercom_only_consumers = steadiercom_samples_preprocessed[(steadiercom_samples_preprocessed.receiver!="environment")]
    members =  set(list(steadiercom_only_producers.donor.unique())+list(steadiercom_only_consumers.receiver.unique()))
    
    if only_media:
        if ReceiverOrDonor=="receiver":
            steadiercom_crossfeeding_from_media = steadiercom_samples_preprocessed[steadiercom_samples_preprocessed.donor=="environment"].copy()
        elif ReceiverOrDonor=="donor":
            steadiercom_crossfeeding_from_media = steadiercom_samples_preprocessed[steadiercom_samples_preprocessed.receiver=="environment"].copy()
    else:
        steadiercom_crossfeeding_from_media = steadiercom_samples_preprocessed.copy()
    
    # Aggregate and pivot
    aggregated = (
        steadiercom_crossfeeding_from_media.groupby(["super_class", ReceiverOrDonor])[metric]
        .sum()
        .unstack(fill_value=0)
    )

    # Remove environment column if present
    aggregated = aggregated.drop(columns="environment")#, errors="ignore")

    # Ensure all members are present as columns (this is only relevant for filtering)
    aggregated = aggregated.reindex(columns=members, fill_value=0)

    return aggregated

        
def data_for_links(steadiercom_crossfeeding_only,members, metric="mass_rate*frequency"):
    """
    Prepare data for links in the circos plot.

    Parameters
    ----------
    steadiercom_crossfeeding_only : pd.DataFrame
        DataFrame containing the cross-feeding data.
    members : set
        Set of member IDs.


    metric : str, optional
        The metric to use for the links (default is "mass_rate*frequency").
        The minimum flux value to consider for the links (default is 0).

    Returns
    -------
    list
        A list of links for the circos plot.
    """

    # Check that the columns that are referenced in steadiercom_samples_preprocessed exist
    required_columns = ["donor", "receiver", "super_class", metric]
    for col in required_columns:
        if col not in steadiercom_crossfeeding_only.columns:
            raise ValueError(f"Column {col} not found in steadiercom_crossfeeding_only. Preparation of links for circos plots will fail")

    if metric!="mass_rate*frequency":
        logging.warning(f"{metric} is not yet supported. Links in the circos plot are scaled assuming that the metric is based on mass rate.")

    links = []
    member_range = {member:{"upper_donor":[],"upper_receiver":[]} for member in members}

    # Calculates the sum of flux for normalisation of the link thickness. 
    sum_fluxes = steadiercom_crossfeeding_only[metric].sum()

    for index,row in steadiercom_crossfeeding_only.iterrows():
        link = []

        # Donor
        donor =row.donor
        # Receiver
        receiver =row.receiver
        flux_g = row[metric]


        lower_receiver_i = 0 if len(member_range[receiver]["upper_receiver"]) == 0 else member_range[receiver]["upper_receiver"][-1]
        upper_receiver_i = lower_receiver_i + flux_g/sum_fluxes
        
        lower_donor_i = 0 if len(member_range[donor]["upper_donor"]) == 0 else member_range[donor]["upper_donor"][-1]
        upper_donor_i = lower_donor_i + flux_g/sum_fluxes
        

        member_range[receiver]["upper_receiver"].append(upper_receiver_i)

        member_range[donor]["upper_donor"].append(upper_donor_i)
        
        link.append((donor,1-lower_donor_i,1-upper_donor_i))   
        link.append((receiver,lower_receiver_i,upper_receiver_i))
     
        link.append(row.super_class)
        links.append(link)
        
    return links