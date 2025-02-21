import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from pyvis.network import Network
import json
import pandas as pd

with open("../data/Venturelli2018/generated/metid2metname.json") as json_file:
    metid2metname = json.load(json_file)
    
# Color end points for colorFader
b1='#f8cdcf' #bright pink
b2='#fa4f56' #red

c1='#ffffff' #white
c2='#0b5394' #blue
    
## From https://stackoverflow.com/questions/25668828/how-to-create-colour-gradient-in-python
def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)






def competition_network_data(results,freq_cutoff=0,metabolomics_mets_only=True):
    """
    Calculating cooperation between pairs of species by looking at the number compounds they both consume. 
    NB! This method is very different to the method for cooperation. 
    
    Input:
    results: pd.DataFrame - steadiercom results dataframe
    cutoff: float - cutoff for frequency (default: cutoff=0)
    metabolomics_mets_only: Only include metabolites in metabolomics dataset (default True)
    
    Returns:
    network_data: pd.Series: node names and edge values for comptetition network
    """
    # Find mean of frequencies for each receiver+compound combination
    results = results[["receiver","compound","frequency"]]
    
    if metabolomics_mets_only:
        results_grouped_rec = results[(results.compound.isin(metid2metname.keys()))].groupby(["receiver","compound"]).mean()
    else:
        results_grouped_rec = results.groupby(["receiver","compound"]).mean()

    # Find overlap between communities
    count_same = results_grouped_rec["frequency"]
    count_same = count_same.unstack().fillna(0)>freq_cutoff #Only for frequency higher than 0.5
    count_same = count_same.astype(int)
    count_same = count_same.dot(count_same.T) # Multiply with the transpose to get overlap

    # count same is symmetric and are only interested in half of these. 
    network_data = count_same.mask(np.tril(np.ones(count_same.shape)).astype(np.bool)).stack() 
    return network_data
    
    
def cooperation_network_data(results,metric="frequency",freq_cutoff=0,metabolomics_mets_only=True,use_count=True):
    """
    Calculating cooperation between species by looking at mean frequencies of compounds shared between them. 
    NB! This method is very different to the method for competition. 
    
    
    Input:
    results: pd.DataFrame - steadiercom results dataframe
    cutoff: float - cutoff for frequency (default: cutoff=0)
    metabolomics_mets_only: Only include metabolites in metabolomics dataset (default True)
    
    Returns:
    network_data: pd.Series: node names and edge values for comptetition network
    """
    results_env = results[(results.donor!="environment") & (results.receiver!="environment")]
    results_env = results_env[results_env.frequency>freq_cutoff]
    
    if metabolomics_mets_only:
        results_env = results_env[(results_env.compound.isin(metid2metname.keys()))]
    
    results_env = results_env[["donor","receiver",metric]]
    
    if use_count:
    # Count number of shared compounds
        results_grouped_rec_env = results_env.groupby(["donor","receiver"]).count()
    
    else:
    # Find mean of frequencies
        results_grouped_rec_env = results_env.groupby(["donor","receiver"]).mean()

    network_data = results_grouped_rec_env[metric].copy()
    
    return network_data 

def competition_centrality_metrics(network_data,cutoff=0,SP_centrality=True):
    """
    A function for finding some centrality metrics for nodes in network. 
    
    Input:
    network data: pd.Series: Edge values
    cutoff: float - cutoff for edge value
    SP_centrality: bool - centrality meassures based on shortest path (True/False)
    """
    
    # Betweenness and closeness
    G = nx.Graph()
    G.add_nodes_from(list(set(list(network_data.index.get_level_values(0).unique())+list(network_data.index.get_level_values(1).unique()))))
    
    if SP_centrality:
        edges = [(index[0],index[1],{"weight":1/data} ) for index,data in network_data.items() if data>cutoff]
    
    else:
        edges = [(index[0],index[1],{"weight":data} ) for index,data in network_data.items() if data>cutoff]
        
    G.add_edges_from(edges)
    
    if SP_centrality:
        return pd.DataFrame({"betweeness (1/data)":nx.centrality.betweenness_centrality(G,weight="weight"),
                     "closeness (1/data)":nx.centrality.closeness_centrality(G,distance="weight")})
    else:
        return pd.DataFrame({"degree (data)":dict(G.degree(weight="weight")),
                 "eigenvector (data)":nx.centrality.eigenvector_centrality(G,weight="weight")})
    

    
def cooperation_centrality_metrics(network_data,cutoff=0,SP_centrality=True,degree_type="both",SP_type="in"):
    """
    A function for finding some centrality metrics for nodes in network. 
    
    Input:
    network data: pd.Series: Edge values
    cutoff: float - cutoff for edge value
    SP_centrality: bool - centrality meassures based on shortest path (True/False)
    """
    # Betweenness and closeness
    G = nx.DiGraph()
    G.add_nodes_from(list(set(list(network_data.index.get_level_values(0).unique())+list(network_data.index.get_level_values(1).unique()))))
    
    if SP_centrality:
        edges = [(index[0],index[1],{"weight":1/data} ) for index,data in network_data.items() if data>cutoff]
    
    else:
        edges = [(index[0],index[1],{"weight":data} ) for index,data in network_data.items() if data>cutoff]
        
    G.add_edges_from(edges)
    
    if SP_centrality:
        
        if SP_type=="in":
        
            return pd.DataFrame({"betweeness (1/data)":nx.centrality.betweenness_centrality(G,weight="weight"),
                         "closeness (1/data)":nx.centrality.closeness_centrality(G,distance="weight")})
        elif SP_type=="out":
        
            return pd.DataFrame({"betweeness (1/data)":nx.centrality.betweenness_centrality(G.reverse(),weight="weight"),
                         "closeness (1/data)":nx.centrality.closeness_centrality(G.reverse(),distance="weight")})
        else:
             raise Exception("SP_type must be 'in' or 'out'")
    else:
        if degree_type=="both":
            return pd.DataFrame({"degree (data)":dict(G.degree(weight="weight")),
                     "eigenvector (data)":nx.centrality.eigenvector_centrality(G,weight="weight")})

        elif degree_type=="in":
                return pd.DataFrame({"degree (data)":dict(G.in_degree(weight="weight")),
                         "eigenvector (data)":nx.centrality.eigenvector_centrality(G,weight="weight")})
        elif degree_type=="out":
                return pd.DataFrame({"degree (data)":dict(G.out_degree(weight="weight")),
                         "eigenvector (data)":nx.centrality.eigenvector_centrality(G,weight="weight")})
        else:
            raise Exception("degree_type must be 'in', 'out' or 'both'")

def competition_network(network_data,closeness_centrality_met,cutoff=0):
    net = Network()
    
    maximum_centrality = pd.Series(closeness_centrality_met).max()

    nodes = list(set(list(network_data.index.get_level_values(0).unique())+list(network_data.index.get_level_values(1).unique())))
    
    
    color = [colorFader(b1,b2,closeness_centrality_met[node]/maximum_centrality) for node in nodes]
    value = [(closeness_centrality_met[node]*10)**2 for node in nodes]
    net.add_nodes(nodes,color=color,value=value)

    maximum_overlap = network_data.max()
    for index,value in network_data.items():
        if value>cutoff:
            net.add_edge(index[0],index[1],width=(value**2)/10,color=colorFader(c1,c2,value/maximum_overlap))

    net.barnes_hut(gravity=-1000, central_gravity=0.1, spring_length=150, spring_strength=0.0015, damping=0.09, overlap=0.8)


    template = net.templateEnv.get_template(net.path)
    net.template = template
    return net


def cooperation_network(network_data,closeness_centrality_met,cutoff=0,use_count=True):
    net = Network(directed=True)
    
    
    maximum_centrality = closeness_centrality_met.max()
    nodes = list(set(list(network_data.index.get_level_values(0).unique())+list(network_data.index.get_level_values(1).unique())))
    
    color = [colorFader(b1,b2,closeness_centrality_met[node]/maximum_centrality) for node in nodes]
    value = [(closeness_centrality_met[node]*10)**2 for node in nodes]

    net.add_nodes(nodes,color=color,value=value)
    
    
    maximum_corr = network_data.max()
    for index,value in network_data.items():
        if value>cutoff:
            
            if use_count:
                net.add_edge(index[0],index[1],width=(value**2)/10,color=colorFader(c1,c2,value/maximum_corr))

            else:
                net.add_edge(index[0],index[1],width=(value**2)*10,color=colorFader(c1,c2,value/maximum_corr)) 

    net.barnes_hut(gravity=-1000, central_gravity=0.1, spring_length=150, spring_strength=0.0015, damping=0.09, overlap=0.8)


    template = net.templateEnv.get_template(net.path)
    net.template = template
    return net

    