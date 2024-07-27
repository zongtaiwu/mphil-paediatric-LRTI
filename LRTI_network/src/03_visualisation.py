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
inputprot = pd.read_csv('./data/input_genes_name.txt', header = None)
membraneprot = pd.read_csv('./data/network_membrane_gene_names.txt', header = None)
nucleusprot = pd.read_csv('./data/network_nucleus_gene_names.txt', header = None)

# synthesise fc
node_names = list(G_igraph.nodes)
name = pd.DataFrame(columns=['gene_symbol'])
name['gene_symbol'] = node_names

id_dic = pd.read_csv('./data/host_gene_counts.csv')
id_dic.rename(columns={'Unnamed: 0': 'gene_id'}, inplace=True) # used to map ids
fc_dic = pd.read_csv('./data/logfc.csv')
fc_dic.rename(columns={'Unnamed: 0': 'gene_id'}, inplace=True) # used to retrieve log2foldchange

merged_df = pd.merge(name, id_dic[['gene_id', 'gene_symbol']], on="gene_symbol", how="left")
fc = pd.merge(merged_df, fc_dic[['gene_id', 'log2FoldChange']], on="gene_id", how="left")
fc.drop(columns=['gene_id'], inplace=True)
fc.columns = [''] * len(fc.columns)
fc.to_csv('./data/fc_DPvsDN.csv', index=False, header=False)

# read fc
fc = pd.read_csv('./data/fc_DPvsDN.csv', header = None)

# map to input genes names
inputgene = pd.read_csv('./data/input_genes.txt', header = None)
inputgene.rename(columns={0: 'gene_id'}, inplace=True)
inputdf = pd.merge(inputgene, id_dic[['gene_id', 'gene_symbol']], on="gene_id", how="left")
inputdf.drop(columns=['gene_id'], inplace=True)
inputdf.to_csv("./data/input_genesdf.csv", index=False, header=False)

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
    df1 = pd.read_csv(sif_address, header = None)
    listrow = []
    for index, row in df1.iterrows():
        edge1 = list( [str(row[0])] ) 
        edge2 = list( [str(row[1])] ) 
        listrow += [tuple(edge1 + edge2)]
    df_dict = {}
    for i, element in enumerate(listrow):
        df_dict[i] = df1.iloc[i, 2]
    weight = {}
    weightlist = []
    edges = list(G.edges())
    exclude = []
    for edge in edges:
        if edge in listrow:
            a = listrow.index(edge)
            arr = df_dict[a]
            switch = tuple( list( [edge[1]] ) + list( [edge[0]] ))
            if arr == '->' :
                weight[edge] = 2
                weightlist += [tuple( (edge[0] , edge[1] , 2) )]
                exclude += [tuple(edge)]
            if arr == '<-' :            
                weight[switch] = 2
                weightlist += [tuple( (edge[1] , edge[0] , 2) )]
                exclude += [tuple(switch)]
            if arr == '-|' :
                weight[edge] = 0.1
                weightlist += [tuple( (edge[0] , edge[1] , 0.1) )]
                exclude += [tuple(edge)]
            if arr == '|-' :
                weight[switch] = 0.1
                weightlist += [tuple( (edge[1] , edge[0] , 0.1) )]
                exclude += [tuple(switch)]
            if arr == '|->' :
                weight[edge] = 2
                weight[switch] = 0.1
                weightlist += [tuple( (edge[0] , edge[1] , 2) )]
                weightlist += [tuple( (edge[1] , edge[0] , 0.1) )]
                exclude += [tuple(edge)]
                exclude += [tuple(switch)]
            if arr == '<-|' :
                weight[edge] = 0.1
                weight[switch] = 2
                weightlist += [tuple( (edge[0] , edge[1] , 0.1) )]
                weightlist += [tuple( (edge[1] , edge[0] , 2) )]
                exclude += [tuple(edge)]
                exclude += [tuple(switch)]
            if arr == '<->' :
                weight[edge] = 2
                weight[switch] = 2
                weightlist += [tuple( (edge[0] , edge[1] , 2) )]
                weightlist += [tuple( (edge[1] , edge[0] , 2) )]
                exclude += [tuple(edge)]
                exclude += [tuple(switch)]
                    
                    
    for edge in edges:
        bol = False
        switch = tuple( list( [edge[1]] ) + list( [edge[0]] ))
        for element in exclude:
            if element == edge:
                bol = True
        if bol == True:
            continue             
        else:  
            weight[edge] = 1
            weightlist += [tuple( (edge[0] , edge[1] , 1) )]
            exclude += [tuple(edge)]
            
        bol = False
        for element in exclude:
            if element == switch:
                bol = True
        if bol == True:
            continue             
        else:  
            weight[switch] = 1
            weightlist += [tuple( (edge[1] , edge[0] , 1) )]
            exclude += [tuple(switch)]
    
    GG = nx.MultiDiGraph()
    GG.add_nodes_from(G.nodes)
    GG.add_weighted_edges_from(weightlist)
    
    A = nx.convert_matrix.to_numpy_array(GG, nonedge = 'nan', weight = 'weight')
    fig, ax = plt.subplots(1,1, figsize=(15,15))
    ax.imshow(np.array(A), cmap='RdYlGn')
    ax.set_xticks(np.arange(len(G.nodes)))
    ax.set_yticks(np.arange(len(G.nodes)))
    ax.set_xticklabels(G.nodes, rotation = 90)
    ax.set_yticklabels(G.nodes)
    fig.savefig(gridplot_address)
    
    return GG

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