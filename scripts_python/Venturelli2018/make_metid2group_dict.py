import json

metname_groups = {
    "Choline pool": {
        "choline_positive_7p50": 1
    },

    "Adenosine & Derivatives pool": {
        "adenosine_positive_6p00": 1,
        "methyl-5thioadenosine_positive_3p32": 1,
        "1-methyladenosine_positive_3p80": 1,
        "cyclic-amp_positive_9p70": 1
    },
    "Xanthine pool": {
        "xanthine_negative_7p50": 1
    },
    "Hypoxanthine pool": {
        "hypoxanthine_positive_6p50": 1,
        "inosine_negative_7p90": 1
    },
    "Adenine pool": {
        "adenine_negative_5p64": 1
    },
    "Aspartate pool": {
        "l-aspartate_positive_14p79": 1
    },
    "Glycine betaine pool": {
        "glycine_betaine_positive_7p75": 1
    },
    "Guanosine pool": {
        "guanosine_positive_10p00": 1,
        "2deoxyguanosine_positive_8p40": 1,
        "guanosine_2_3cyclic_monophosphate_positive_13p00": 1
    },
    "Guanine pool": {
        "guanine_positive_9p00": 1
    },
    "Uridine pool": {
        "uridine_negative_6p12": 1,
        "2deoxyuridine_negative_4p15": 1
    },
    "Cytidine pool": {
        "cytidine_positive_8p80": 1,
        "cytidine_5monophosphate_negative_15p65": 1,
        "2deoxycytidine_negative_6p89": 1
    },
    "Cytosine pool": {
        "cytosine_positive_7p60": 1
    },
    "Thymidine pool": {
        "thymidine_negative_3p33": 1
    },
    "Uracil pool": {
        "uracil_negative_3p95": 1,
        "5-hydroxymethyluracil_negative_5p86": 1
    },
    "Uric acid pool": {
        "uric_acid_negative_12p53": 1,
        "orotate_negative_10p17": 1
    },
    "Carnitine pool": {
        "l-carnitine_positive_10p70": 1
    },
    "Serine pool": {
        "l-serine_positive_13p90": 1
    },
    "Leucine pool": {
        "l-leucine_negative_8p29": 1,
        "histidylleucine_positive_9p41": 0.5,
        "l-leucyl-lalanine_positive_7p70": 0.5
    },
    "Isoleucine pool": {
        "l-isoleucine_negative_8p78": 1,
        "l-isoleucyl-l-proline_positive_8p10": 0.5,
        "gamma-glutamylisoleucine_positive_13p22": 0.5,
        "2-ketobutyric_acid_negative_5p14": 1
    },
    "Histidine pool": {
        "l-histidine_positive_13p09": 1,
        "histidylleucine_positive_9p41": 0.5,
        "acetyl-l-histidine_positive_10p80": 1
    },
    "Ectoine pool": {
        "ectoine_positive_10p80": 1
    },
    "Alanine pool": {
        "l-alanine_positive_12p21": 1,
        "arg-ala_positive_18p10": 0.5,
        "gamma-glutamylalanine_positive_16p20": 0.5,
        "l-leucyl-lalanine_positive_7p70": 0.5,
        "val-ala_positive_8p98": 0.5
    },
    "Lysine pool": {
        "l-lysine_positive_23p71": 1,
        "l-saccharopine_positive_15p90": 1,
        "dl-pipecolic_acid_positive_23p80": 1,
        "3-hydroxy-3-methyl-glutaric_acid_negative_14p30": 1,
        "d-lysopine_positive_9p00": 1
    },
    "Methionine pool": {
        "l-methionine_negative_9p70": 1,
        "dl-methionine_sulfoxide_positive_12p50": 1,
        "1-aminocyclopropane-1-carboxylate_positive_12p50": 1
    },
    "Asparagine pool": {
        "l-asparagine_positive_13p70": 1
    },
    "Proline pool": {
        "l-proline_positive_9p90": 1,
        "l-isoleucyl-l-proline_positive_8p10": 0.5,
        "4-hydroxy-l-proline_positive_11p95": 1
    },
    "Ornithine pool": {
        "d-ornithine_negative_24p12": 1,
        "l-ornithine_positive_23p95": 1,
        "feruloylputrescine_positive_4p80": 1
    },
    "Arginine pool": {
        "l-arginine_positive_25p00": 1,
        "arg-ala_positive_18p10": 0.5,
        "4-guanidinobutanoate_positive_12p60": 1
    },
    "Glutamate pool": {
        "l-glutamate_positive_14p68": 1,
        "gamma-glutamylalanine_positive_16p20": 0.5,
        "pyroglutamic_acid_negative_10p45": 1,
        "5-oxo-l-proline_positive_10p45": 1,
        "gamma-glutamylisoleucine_positive_13p22": 0.5
    },
    "Glutamine pool": {
        "l-glutamine_negative_13p06": 1
    },
    "Threonine pool": {
        "l-threonine_negative_12p48": 1
    },
    "Valine pool": {
        "l-valine_positive_9p90": 1,
        "val-ala_positive_8p98": 0.5,
        "val-phe_positive_5p91": 0.5
    },
    "Tyrosine pool": {
        "l-tyrosine_negative_11p05": 1
    },
    "Taurine pool": {
        "taurine_negative_12p67": 1
    },
    "Glycine pool": {
        "glycine_negative_13p29": 1
    },
    "Phenylalanine pool": {
        "l-phenylalanine_negative_8p65": 1,
        "val-phe_positive_5p91": 0.5,
        "n-acetyl-l-phenylalanine_negative_4p58": 1
    },
    "Tryptophan pool": {
        "l-tryptophan_negative_9p40": 1,
        "n-methylanthranilic_acid_positive_3p50": 1
    },
    "Citrulline pool": {
        "l-citrulline_positive_14p20": 1
    },
    "Creatinine pool": {
        "creatinine_positive_5p88": 1
    },
    "Creatine pool": {
        "creatine_positive_12p20": 1
    },
    "Vitamin B3 pool": {
        "nicotinamide_positive_3p16": 1,
        "nicotinate_negative_7p52": 1
    },
    "Biotin precursor pool": {
        "7_8-diaminononanoate_positive_19p80": 1
    },
    "Branched-chain AA precursor pool": {
        "2-aminoisobutyrate_negative_12p98": 1,
        "dl-3-amino-isobutanoic_acid_positive_13p22": 1
    },
    
    "TCA Cycle Intermediates": {
        "fumarate_negative_16p36": 1,
        "citraconic_acid_negative_4p24": 1,
        "dl-citramalic_acid_negative_15p64": 1
    },

    "Aromatic Metabolites & Phenolics": {
        "2-coumarate_positive_11p15": 1,
        "hydroxyquinol_negative_3p56": 1,
        "2-hydroxyphenylacetic_acid_negative_3p27": 1
    },
    
    "Phospholipid breakdown pool": {
        "glycero-3-phosphocholine": 1
    },
    "Glucose pool": {
        "d-gluconic_acid_negative_13p65": 1,
        "1_6-anhydrobeta-d-glucopyranose_negative_4p45": 1
    },
    "Ribose pool": {
        "l-ribulose_negative_11p90": 1
    },
    "5-Aminovaleric acid pool": {
        "5-aminovaleric_acid_positive_13p23": 1
    }
}

with open("../data/Venturelli2018/generated/metid2metname.json") as file:
    metid2metname = json.load(file)

metname2group = {metname:group for group,metnames in metname_groups.items() for metname in metnames}
metid2group = {metid:metname2group[metname] for metid, metname in metid2metname.items()}


with open("../data/Venturelli2018/generated/metname_groups.json","w") as file:
    json.dump(metname_groups,file)

with open("../data/Venturelli2018/generated/metname2group.json","w") as file:
    json.dump(metname2group,file)

with open("../data/Venturelli2018/generated/metid2group.json","w") as file:
    json.dump(metid2group,file)