import pandas as pd
import numpy as np
import json
import requests
import csv
import warnings
from os import listdir
import reframed


print("Loading data and models for mapping process...")

with open('../../data/generated/metname2metanetx.json') as json_file:
    metname2metanetx = json.load(json_file)

models = {}
with warnings.catch_warnings():
    warnings.filterwarnings('ignore')
    for path in listdir("../../models/species_used/"):
        models[path[:-4]] = reframed.load_cbmodel("../../models/species_used/"+path,exchange_detection="R_EX_(?!biomass__40__e__41__)")
    

    
print("Finding matches of metanetx identifiers in other databases...")
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
        print(i)
        
    for met_name, metanetx_id in metname2metanetx.items():
        if metanetx_id==row[1]:
            database_entry = row[0].split(":")
            if len(database_entry)<2:
                continue
            metanetx2other[metanetx_id][database_entry[0]]=database_entry[1]
            
            
    i+=1    
    
    

print("Searching for identifiers in model metabolites...")
final_mapping = {}
for metanetx_id,database_ids in metanetx2other.items():
    for database_id,id_ in database_ids.items():
        
        database_id = database_id.lower()
        
        if database_id=="chebi":
            id_="CHEBI:"+id_
            
        if database_id=="metacyc.compound":
            database_id = "biocyc"
            
        for model in models.values():
            for met in model.metabolites:
                metadata = model.metabolites[met].metadata

                if database_id not in metadata.keys():
                    continue
                else:

                    if metadata[database_id][0]==id_:
                        final_mapping[met] = metanetx_id


metanetx2metname = {metanetx_id:name for name,metanetx_id in metname2metanetx.items()}

metid2metname = {}

for met_id,metanetx_id in final_mapping.items():
    metid2metname[met_id]=metanetx2metname[metanetx_id]

print("Save data ../../data/generated/metid2metname.json")
with open("../../data/generated/metid2metname.json","w") as outfile:
    json.dump(metid2metname,outfile)