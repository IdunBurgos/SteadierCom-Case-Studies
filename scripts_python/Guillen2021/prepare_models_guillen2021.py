
import requests
import pandas as pd
import reframed
import json

import logging

logging.basicConfig(
    level=logging.INFO, 
    format="%(levelname)-8s | %(message)s",
    )

def load_iM1515_model():
    """
    The model from Bernstein et al (2023) represent the latest changes to the E.coli model.
    The study focused on improving genotype-phenotype relationships and incorporating new experimental data.
    """
    
   # Load curated model from Bernstein et al (2023)
    logging.info("Downloading and loading iML1515 model from Bernstein et al (2023) github repository...")
    url = "https://github.com/dbernste/E_coli_GEM_validation/raw/refs/heads/main/Analysis/iML1515_model_adjusted_all_corrections.xml"
    response = requests.get(url)

    if response.status_code == 200:
        with open("../models/Guillen2021/iML1515.xml", "wb") as f:
            f.write(response.content)
    else:
        logging.error("Failed to download the file. Check the URL.")

    model_ecoli = reframed.load_cbmodel("../models/Guillen2021/iML1515.xml")

    return model_ecoli

def adjust_btn_biomass(model_ecoli):
    """
    Changing stoichiometry of biotin in biomass (Because value was below solver tolerance)
    - Biotin in biomass originally: 2e-06 -> 2e-06(mmol)*244,31(g/mol) = 4.8862e-7 (g)
    - Increase to: 2e-04 -> 2e-04(mmol)*244,31(g/mol) = 4.8862e-3 (g)
    """
    
    but_prev= model_ecoli.reactions.R_BIOMASS_Ec_iML1515_core_75p37M.stoichiometry["M_btn_c"]
    but_new = -2e-04
    
    logging.info(f"Adjusting biomass stoichiometry for butanoate from {but_prev} to {but_new}")
    model_ecoli.reactions.R_BIOMASS_Ec_iML1515_core_75p37M.stoichiometry["M_btn_c"] = but_new

    return model_ecoli

def add_transport_reactions(model_ecoli):
    """
    Adding necessary reactions to the model.
    """
    logging.info("Adding reactions to the model...")
    model_ecoli.add_reaction_from_str("R_PNTOtex2: M_pnto__R_c --> M_pnto__R_e")
    model_ecoli.add_reaction_from_str("R_BTNtex2: M_btn_c --> M_btn_e")

    return model_ecoli

def remove_molybdate_genes(model_ecoli):
    """
    Removing molybdate genes from the model.
    """
    logging.info("Removing molybdate genes from the model...")
    model_ecoli.reactions.R_SULabcpp.gpr.proteins.pop(1) # Remove molybdate genes

    return model_ecoli

def map_geneid2name(model_ecoli):
    logging.info("Mapping gene IDs to gene names...")
    id2name = {}
    for gene in model_ecoli.genes:
        id2name[gene] = model_ecoli.genes[gene].name

    id2name_ser = pd.Series(id2name)
    id2name_ser.name="name"
    id2name_ser.index.name="id"
    return id2name_ser

def load_relative_abundance_data():
    """
    Loading relative abundance data for E.coli mutants.

    NB! Removing the gene for panD instead because panZ is a regulatory protein. 
    """

    logging.info("Loading relative abundance data...")
    relative_abundance = pd.read_csv("../data/Guillen2021/relative_abundance.csv",index_col=0)
    # NB! yhhK in data == panZ in paper
    relative_abundance.index =relative_abundance.index.map(lambda x:x if x!="yhhK" else "panZ")
    return relative_abundance


def create_mutant_models(model_ecoli, id2name_ser,relative_abundance):
    """
    Creating models for each E.coli mutant in the Guillen 2021 study.
    NB: remove cysZ for cysW mutant to get desired behavior, and panD instead of panZ (panZ is a regulatory protein).
    """

    logging.info("Creating single knockout models for each gene in relative abundance data...") 

    name2id_38 = id2name_ser[id2name_ser.isin(relative_abundance.index)]
    name2id_38 = dict(zip(name2id_38.values, name2id_38.index))

    # additional details for special cases
    name2id_38["panZ"]=name2id_38["panD"] # panZ is a regulatory protein for panD

    assert len(name2id_38)==38, f"Expected 38 genes in name2id mapping, but got {len(name2id_38)}. Check if all genes in relative abundance data are present in the model."

    # make models
    ecoli_models = {}

    for name,id_ in name2id_38.items():
        if name=="cysW": # Special case: cysW - also remove cysZ
            ecoli_models[name] = reframed.hard_knockout(model_ecoli,[id_,"G_b2413"],
                                                        rename="mutant_"+name) #remove cysZ gene
        else: 
            ecoli_models[name] = reframed.hard_knockout(model_ecoli,[id_],rename="mutant_"+name)
    return ecoli_models,name2id_38

def main():
    logging.info("STARTING SCRIPT FOR MODEL PREPARATION")
    model_ecoli = load_iM1515_model()
    model_ecoli = adjust_btn_biomass(model_ecoli)
    model_ecoli = add_transport_reactions(model_ecoli)
    model_ecoli = remove_molybdate_genes(model_ecoli)
    model_ecoli.update()

    logging.info("MODEL PREPARATION COMPLETE."\
                 "\n Saving model as ../models/Guillen2021/iM1515_updated.xml")

    reframed.save_cbmodel(model_ecoli,"../models/Guillen2021/iM1515_updated.xml")

    logging.info("Building mutant models...")

    id2name_ser = map_geneid2name(model_ecoli)
    relative_abundance = load_relative_abundance_data()
    ecoli_models,name2id_38 = create_mutant_models(model_ecoli, id2name_ser, relative_abundance)
    
    logging.info("SAVING DATA.")  
    logging.info("Saving gene ID to name mapping from model in ../data/Guillen2021/generated/gene_id2name.tsv.")
    id2name_ser.to_csv("../data/Guillen2021/generated/gene_id2name.tsv",sep="\t")

    logging.info("Saving gene knockout in model for gene name ../models/Guillen2021/model_gene_KO.json. " \
                "\n NB: this accounted for effect of knockout of regulatory genes like panZ.")
    
    with open("../models/Guillen2021/model_gene_KO.json", mode='w') as writer:
        json.dump(name2id_38, writer)

    logging.info("SAVING MODELS.")
    logging.info("Saving all mutant models in ../models/Guillen2021/strains/...")
    for name,model in ecoli_models.items():
        reframed.save_cbmodel(model,"../models/Guillen2021/strains/"+model.id+".xml")
    
    logging.info("Saving top 7 most abundant strains in ../models/Guillen2021/top_7_strains/...")   
    ra_mean = relative_abundance.T.mean()
    top7_members = ra_mean[ra_mean>0.01].index

    for member in top7_members:
        model = ecoli_models[member]
        reframed.save_cbmodel(model,"../models/Guillen2021/top_7_strains/"+model.id+".xml")

if __name__ == "__main__":
    main()