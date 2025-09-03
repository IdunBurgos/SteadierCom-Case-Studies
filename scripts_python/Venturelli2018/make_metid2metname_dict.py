import numpy as np
import json
import requests
import csv
import warnings
from os import listdir
import reframed

import logging
logging.basicConfig(
    level=logging.INFO, 
    format="%(levelname)-8s | %(message)s",
)


def load_data_models():
    logging.info("Loading data and models for mapping process...")

    with open('../data/Venturelli2018/generated/metname2metanetx.json') as json_file:
        metname2metanetx = json.load(json_file)

    models = {}
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        for path in listdir("../models/Venturelli2018/species_used/"):
            models[path[:-4]] = reframed.load_cbmodel("../models/Venturelli2018/species_used/"+path,exchange_detection="R_EX_(?!biomass__40__e__41__)")

    return models, metname2metanetx

def create_metanetx2other(metname2metanetx):
    logging.info("Finding matches of metanetx identifiers in other databases...")
    databases = ['bigg.metabolite','biocyc','chebi','hmdb','kegg.compound','metanetx.chemical','pubchem.compound', 'seed.compound']

    url = "https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv"

    # Stream the file using requests
    response = requests.get(url, stream=True)
    response.raise_for_status()

    # Use csv.reader to process the file and skip lines
    reader = csv.reader((line.decode("utf-8") for line in response.iter_lines()), delimiter="\t")


    # Skip the first 350 lines
    for _ in range(351):
        next(reader, None)  # Move the iterator forward
        
    metanetx2other = {metanetx_id:{} for metanetx_id in metname2metanetx.values()}
    # Process the remaining lines
    i = 0

    
    for row in reader:
        if np.mod(i,50000)==0:
            logging.info(f"Processed {i} rows")

        for met_name, metanetx_id in metname2metanetx.items():
            if metanetx_id==row[1]:
                database_entry = row[0].split(":")
                if len(database_entry)<2:
                    continue
                metanetx2other[metanetx_id][database_entry[0]]=database_entry[1]
                
                
        i+=1    
    return metanetx2other

def map_metanetx_to_model(metanetx2other,metname2metanetx, models):
    logging.info("Searching for identifiers in model metabolites...")
    final_mapping = {}
    for metanetx_id,database_ids in metanetx2other.items():
        for database_id,id_ in database_ids.items():
            
            database_id = database_id.lower()
            
            if database_id=="chebi":
                id_="CHEBI:"+id_
                
            if database_id=="metacyc.compound":
                database_id = "biocyc"

                
    # 
            for model in models.values():
                for met_id, met_obj in model.metabolites.items():
                    # Skip if the target db_name isn't in the metadata
                    if database_id not in met_obj.metadata:
                        continue

                    # Compare metadata value(s) to id_
                    ids = met_obj.metadata[database_id]
                    if isinstance(ids, list):
                        if id_ in ids:
                            final_mapping[met_id] = metanetx_id
                    elif id_ == ids:
                        final_mapping[met_id] = metanetx_id

    metanetx2metname = {metanetx_id:name for name,metanetx_id in metname2metanetx.items()}

    metid2metname = {}
    for met_id,metanetx_id in final_mapping.items():
        metid2metname[met_id]=metanetx2metname[metanetx_id]

    return metid2metname


def main():
    logging.info("STARTING SCRIPT FOR METID2METNAME")

    models, metname2metanetx = load_data_models()
    metanetx2other = create_metanetx2other(metname2metanetx)
    metid2metname = map_metanetx_to_model(metanetx2other,metname2metanetx, models)

    logging.info("Save data ../data/Venturelli2018/generated/metid2metname.json")
    with open("../data/Venturelli2018/generated/metid2metname.json","w") as outfile:
        json.dump(metid2metname,outfile)

if __name__ == "__main__":
    main()
