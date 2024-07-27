#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri July 27 14:12:50 2024

@author: Zongtai Wu based on Han et al., 2021, original code created by woochanghwang at 18/11/2020

"""

import matplotlib.pyplot as plt
import pickle
import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
import itertools

np.random.seed(111)

def load_obj(file_addr):
    with open(file_addr+ '.pkl', 'rb') as f:
        return pickle.load(f)

def igraph_to_networkx(igraph_graph):
    edges = igraph_graph.get_edgelist()
    nx_graph = nx.Graph()
    nx_graph.add_edges_from(edges)
    # Add node attributes if any
    for attr in igraph_graph.vertex_attributes():
        nx.set_node_attributes(nx_graph, {v.index: v[attr] for v in igraph_graph.vs}, attr)
    # Rename nodes to their 'name' attribute if available
    if 'name' in igraph_graph.vs.attributes():
        mapping = {v.index: v['name'] for v in igraph_graph.vs}
        nx_graph = nx.relabel_nodes(nx_graph, mapping)
    return nx_graph

"""input adresses"""
###
eigen = pd.read_csv('./data/Centrality_RWR_results.csv')
G_igraph = load_obj('./data/network_genes')
# Convert igraph to NetworkX
G = igraph_to_networkx(G_igraph)
keyprot = pd.read_csv('./data/key_proteins.txt')
inputprot = pd.read_csv('./data/in_and_out.txt', header = None)
membraneprot = pd.read_csv('./data/network_membrane_gene_names.txt', header = None)
nucleusprot = pd.read_csv('./data/network_nucleus_gene_names.txt', header = None)

# synthesise fc
node_names = list(G.nodes)
name = pd.DataFrame(columns=['gene_symbol'])
name['gene_symbol'] = node_names
# LogFC from the four studies
df1 = pd.read_csv('/Users/zongt1/VAP/project-main/RSV-meta/df1.csv')
df2 = pd.read_csv('/Users/zongt1/VAP/project-main/RSV-meta/df2.csv')
df3 = pd.read_csv('/Users/zongt1/VAP/project-main/RSV-meta/df3.csv')
df4 = pd.read_csv('/Users/zongt1/VAP/project-main/RSV-meta/df4.csv')

# Combine the dataframes into a single dataframe
combined_df = pd.concat([df1, df2, df3, df4])
# Calculate the average log fold change for each gene
avg_log_fold_change = combined_df.groupby('Name')['logFC'].mean()

# Map the average log fold change to the name dataframe
name['logFC'] = name['gene_symbol'].map(avg_log_fold_change)
name.to_csv('./data/fc_DPvsDN.csv', index=False, header=False)

# read fc
fc = pd.read_csv('./data/fc_DPvsDN.csv', header = None)


"""output adresses"""
###
address_subnetwork_plot = './visualisation/subnetwork_plot.tiff'
address_subsubnetwork_plot = './visualisation/subsubnetwork_plot.tiff'
address_network_plot = './visualisation/network_plot.tiff'
address_subsubnetworkloc_plot = './visualisation/subsubnetworkloc_plot.tiff'
address_list_subsubnetwork = './data/gene_subsubnetwork.txt'
address_list_qmidnetwork = './data/gene_qmidnetwork.txt'
address_qmidnetwork_plot = './visualisation/qmidnetwork_plot.tiff'


#%%
# Convert igraph to NetworkX
membraneprot = list(membraneprot[0])
nucleusprot = list(nucleusprot[0])
allprot = list(dict.fromkeys(list(keyprot['Gene']) + (list(inputprot.iloc[:,0]))))
fc = fc.set_index(0)
fc_mean_dict = fc.to_dict()
fc_mean_dict = fc_mean_dict[1]
eigen = eigen.set_index('Unnamed: 0')
eigen_dict = eigen.to_dict()
eigen_dict = eigen_dict['Eigen']

#use this to use all edges that contain key/input protein ....
edges = []
for protein in allprot:
    for element in G.edges:
        if protein in element:
            edges += [element]
G_sub = nx.Graph()
G_sub.add_edges_from(edges)

# use this to ONLY use key and input proteins
G_subsub = G.subgraph(allprot).copy()

G_sub.remove_nodes_from(list(nx.isolates(G_sub)))
G_subsub.remove_nodes_from(list(nx.isolates(G_subsub)))

pd.DataFrame(list(G_subsub.nodes)).to_csv('../data/subsubnetwork_genes.txt', index = False, header=None)


#%% 
def makeNetworkPlot(G, layoutstyle, fc_dict):
    # kamada kawai layout uses as few crossing edges as possible
    # multipartite layout puts membrane proteins left, cytosol proteins middle and nucleus proteins right
    
    # set node attributes: membrane 0, cytosol 1, nucleus 2
    layers = {}
    nodes = list(G.nodes)
    for node in nodes:
        if node in membraneprot:
            layers[node] = 0
        elif node in nucleusprot:
            layers[node] = 2
        else:
            layers[node] = 1
    nx.set_node_attributes(G, layers, 'layers')
    
    # generate some colormaps
    # give input proteins black node edge
    color_map_edges = []
    for node in nodes:
        if node in list(inputprot.iloc[:,0]):
            color_map_edges.append('k')
        else:
            color_map_edges.append('w') 
    
    # define node size by eigen centrality
    eigen_centrality = []
    for node in nodes:
        if node in eigen_dict:
            eigen_centrality += [eigen_dict[node]*50000]
        else:
            eigen_centrality += [0]
    
    # define node color by log2 fold change expression
    fc_cmap = []
    for node in nodes:
        if node in fc_dict:
            if np.isnan(fc_dict[node]) == True:
                fc_cmap += [0]
            else:
                fc_cmap += [fc_dict[node]]
        else:
            fc_cmap += [0]

    # make non-symmetrical colorscale where 0 is white, >0 is red, <0 is blue
    plt.set_cmap('bwr')
    a = np.vstack([np.array(fc_cmap),np.array(fc_cmap)])
    norm = mcolors.TwoSlopeNorm(vmin=a.min(), vmax = a.max(), vcenter=0)
    fc_cmap = norm(fc_cmap)
    im = plt.imshow(a, norm = norm)
    
    # plot
    fig  = plt.figure(figsize=(40,40))
    
    # put membrane proteins more to left, nucleus more to right                                   
    pos_multipartite = nx.multipartite_layout(G,subset_key = 'layers') 
    for node in nodes:
        pos_multipartite[node][0] = pos_multipartite[node][0] + (0.001*(np.random.random_sample(1)[0] - 0.002))
       
    pos_kamadakawai = nx.layout.kamada_kawai_layout(G, pos = pos_multipartite)    
    
    if layoutstyle == 'kamadakawai': 
        pos = pos_kamadakawai
    elif layoutstyle == 'multipartite': 
        pos = pos_multipartite
    else:
        print('Error! Check layout style; can be "kamadakawai" or "multipartite"' )
        
    cbar = plt.colorbar(im, shrink = 0.5,)
    cbar.ax.tick_params(labelsize=40)
    nx.draw(G,pos=pos,node_color=fc_cmap, edge_color = 'grey',  with_labels=True, font_size = 20,node_size = eigen_centrality, edgecolors = color_map_edges)             
    
    return fig

#%%
### call makeNetworkPlot function to make + save plot, choose between layouts: 'kamadakawai' or 'multipartite'
"""G"""
fig = makeNetworkPlot(G, 'kamadakawai', fc_mean_dict)
fig.savefig(address_network_plot)

"""G_sub"""
fig = makeNetworkPlot(G_sub, 'kamadakawai', fc_mean_dict)
fig.savefig(address_subnetwork_plot)

"""G_subsub"""
fig = makeNetworkPlot(G_subsub, 'kamadakawai', fc_mean_dict)
fig.savefig(address_subsubnetwork_plot)

fig = makeNetworkPlot(G_subsub, 'multipartite', fc_mean_dict)
fig.savefig(address_subsubnetworkloc_plot)