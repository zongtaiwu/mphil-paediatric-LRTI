#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri July 27 14:12:50 2024

@author: Zongtai Wu based on Han et al., 2021, original code created by woochanghwang at 18/11/2020

Running time: depends on input, m * n = x combinations (start proteins: m, end proteins: n), 
    ca. 4 min for x = 1,800 on MacbookPro
    i.e. roughly x/500 = y minutes 

"""
import itertools
import pickle

import igraph as ig
import matplotlib.pyplot as plt
import mygene
import networkx as nx
import numpy as np
import pandas as pd
from copy import deepcopy

#%%
### make sure the working folder is the src folder (in Spyder top right) ###

""" input addresses """
### ENTER PATH/NAME OF INPUT: TEXT FILE OF START PROTEIN LIST AND END PROTEIN LIST FOR NETWORK RECONSTRUCTION (ensembl protein ID) ###
start_addr = "./data/input_proteins.txt"
end_addr = "./data/input_proteins.txt"
### INSERT PATH/NAME OF INPUT: EDGELIST FROM STRING DATABASE ###
string_addr = "./databases/9606.protein.links.v11.0.noscores.threshold400.txt"
### ENTER PATH/NAME OF INPUT: LIST OF PROTEIN ID TO GENE, DICTIONARY ###
# based on Jensen Lab dictionary which is based on alias file from STRING database
address_dictionary_protein_to_gene= './databases/list_protein_to_gene.txt'
### ENTER PATH/NAME OF INPUT: COMPARTMENT DATABASES, LIST OF GENES EXPRESSING PROTEINS THAT ARE IN MEMBRANE/NUCLEUS ###
address_compartment_membrane_genes = './databases/COMPARTMENT_membrane_gene_names.txt'
address_compartment_nucleus_genes = './databases/COMPARTMENT_nucleus_gene_names.txt'

""" output addresses """
### ENTER PATH/NAME OF OUTPUT: NEW NETWORK WITH PROTEIN IDs AS NODES###
address_network_genes = "./data/network_genes" #dont add '.pkl'
### ENTER PATH/NAME OF OUTPUT: NETWORK WITH GENE NAMES AS NODES (AS OPPOSED TO PROTEIN IDs) ###
address_network_proteins = "./data/network_proteins" #dont add '.pkl'
### ENTER PATH/NAME OF OUTPUT: LIST OF NODES AS GENE/PROTEIN NAMES ###
address_network_nodes_genes = './data/network_nodes_genes.txt'
address_network_nodes_proteins = './data/network_nodes_proteins.txt'
### ENTER PATH/NAME OF OUTPUT: LIST OF EDGES AS GENE/PROTEIN NAMES ###
address_network_edges_genes = './data/network_edges_genes.txt'
address_network_edges_proteins = './data/network_edges_proteins.txt'
### ENTER PATH/NAME FOR OUTPUT: list of network genes in membrane, list of network genes in nucleus ###
address_network_membrane_genes = './data/network_membrane_gene_names.txt'
address_network_nucleus_genes = './data/network_nucleus_gene_names.txt'
### ENTER PATH/NAME FOR OUTPUT: table of network parameters ###
address_centrality_results= "./data/Centrality_RWR_results.csv"
### ENTER PATH/NAME FOR OUTPUT: CIRCULAR PLOT OF NETWORK ###
# this is a very simple plot, just so the network can be looked at after running the code #
address_network_plot_genes = './visualisation/network_visualisation_genes.tiff'
address_network_plot_proteins = './visualisation/network_visualisation_proteins.tiff'

#%% FUNCTIONS

# def load_obj(file_addr):
#     with open(file_addr+ '.pkl', 'rb') as f:
#         return pickle.load(f)


def save_obj(obj, file_addr ):
    """function to save pickle file"""
    with open(file_addr + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def make_shortest_path_sources_to_targets(
        pair,
        graph: ig.Graph,
):
    """
    Function: Find shortest paths on STRING between all start/end protein combinations
    :input: pair of start/end proteins
    :return: edgelist of all shortest paths between start/end protein pair if existent
    ~ takes 4 min for an initial list of 43 proteins, i.e. 43x43 = ~1,850 combinations
    """
    all_path_list = []
    try:
        # exclude self-loops
        if pair[0] != pair[1]:
            for p in graph.get_all_shortest_paths(
                v=graph.vs.find(pair[0]).index,
                to=graph.vs.find(pair[1]).index,
                mode='all'
            ):
                p = [graph.vs[idx]['name'] for idx in p]
                pairs_from_p = pairwise(p)
                all_path_list += pairs_from_p
    except:
        pass

    # return new edge list
    return all_path_list


def mapping_genes_id_to(prot_list,id_from, id_to, species):
    mg = mygene.MyGeneInfo()
    return mg.querymany(prot_list, scopes=id_from, fields= id_to, species=species, as_dataframe=True, verbose = True)


def calc_network_centrality_RWR(network: ig.Graph, start_list, end_list, result_save_dir):
    '''
    Function: calculate eigenvector centrality, degree centrality
    :input: network, membrane list as start list, nucleus list as end list, address
    for output, thresholds for centralities and RWR
    :return: csv file containing centrality and RWR-values for each node,
    txt-files containing above threshold genes
    '''
    network_nodes = list(network.vs['name'])

    # eigenvector centrality
    try:
        eigenvector_centrality = network.eigenvector_centrality(directed=False, scale=True)
        eigenvector_centrality_dict = dict(zip(network_nodes, eigenvector_centrality))
    except:
        eigenvector_centrality_dict = dict()
        for node in network_nodes:
            eigenvector_centrality_dict[node] = 0.0
            
    # degree centrality
    try:
        s = 1.0 / (len(network_nodes) - 1.0)    # Normalisation
        degree_centrality = np.array(network.degree(network_nodes, loops=False)) * s
        degree_centrality_dict = dict(zip(network_nodes, degree_centrality))
    except:
        degree_centrality_dict = dict()
        for node in network_nodes:
            degree_centrality_dict[node] = 0.0

    # betweeness centrality
    try:
        s = 2/((len(network_nodes)-1)*(len(network_nodes)-2))
        between_centrality = np.array(network.betweenness(directed=False)) * s
        between_centrality_dict = dict(zip(network_nodes, between_centrality))
    except:
        between_centrality_dict = dict()
        for node in network_nodes:
            between_centrality_dict[node] = 0.0
            
    # betweeness centrality subset: starts at membrane bound proteins ends in nuclear proteins
    try:
        s = 2 / ((len(network_nodes) - 1) * (len(network_nodes) - 2))
        betweensub_centrality = np.array(network.betweenness(
            directed=False,
            sources=start_list,
            targets=end_list,
        )) * s
        betweensub_centrality_dict = dict(zip(network_nodes, betweensub_centrality))
    except:
        betweensub_centrality_dict = dict()
        for node in network_nodes:
            betweensub_centrality_dict[node] = 0.0

    # edge betweenness centrality
    # set non-zero edge betweenness centrality as edge weight
    network.es['weight'] = np.array(network.edge_betweenness(directed=False)) + 1e-6

    # pagerank 
    # returns dictionary of nodes with PageRank as value
    PR_score = network.pagerank(directed=False, weights='weight')
    PR_score = dict(zip(network_nodes, PR_score))

    # personalised pagerank: jumps back to membrane proteins when jumping
    # in COVID-paper: start_genes_for_PR = start_list #{'COVID19':1}, here:
    PRsub_score = network.personalized_pagerank(directed=False, reset_vertices=start_list, weights='weight')
    PRsub_score = dict(zip(network_nodes, PRsub_score))

    # export results as csv
    network_property_df = pd.DataFrame(
        columns=['Eigen', 'Degree', 'Between', 'Between Sub', 'RWR', 'RWR Sub'])
    for node in network_nodes:
        network_property_df.loc[node] = [eigenvector_centrality_dict[node], degree_centrality_dict[node],
                                         between_centrality_dict[node], betweensub_centrality_dict[node], 
                                         PR_score[node], PRsub_score[node]]

    network_property_df.to_csv(result_save_dir)


def make_circular_plot_network(nodes, G_gene, address_network_plot, start_list):
    '''
    Function: make a circular plot of the network, save as tiff in visualisation folder
    :input: network, list of nodes
    :return: none
    '''
    color_map = []
    for node in nodes:
        if node in start_list:
            color_map.append('lightcoral')
        else:
            color_map.append('lightcyan')   
            
    plt.subplots(figsize=(30,30))
    
    n=len(nodes)
    angle = []
    angle_dict = {}
    for i, node in zip(range(n),nodes):
        theta = 2.0*np.pi*i/n
        angle.append((np.cos(theta),np.sin(theta)))
        angle_dict[node] = theta
    pos = {}
    for node_i, node in enumerate(nodes):
        pos[node] = angle[node_i]
        
    description = nx.draw_networkx_labels(G_gene, pos, font_size=20 )
    for node, t in description.items():
        t.set_rotation(angle_dict[node]*360.0/(2.0*np.pi))
                                          
    nx.draw_circular(G_gene, node_color=color_map, with_labels=False, font_size=18, font_color ='k', edge_color = 'grey')
    plt.savefig(address_network_plot)

if __name__ == '__main__':
    # %% MAIN
    # NETWORK RECONSTRUCTION
    print("Start program. Program will tell you when finished.")

    # open start / end protein list for network reconstruction
    start_df = pd.read_csv(start_addr, sep='\t', header = None)
    start_list_prot = start_df.iloc[:,0].to_list()
    end_df = pd.read_csv(end_addr, sep='\t', header = None)
    end_list_prot = end_df.iloc[:,0].to_list()

    print("Create network:")

    with open(string_addr) as network_f:
        string_network_edges = [x.strip().split(',') for x in network_f.readlines()]
    print('open STRING list done')

    # STRING edgelist to network
    # using iGraph
    nodes = list(set(
        [each for edge in string_network_edges for each in edge]
    ))
    graph = ig.Graph(directed=False)
    graph.add_vertices(nodes)
    graph.add_edges(string_network_edges)
    graph = graph.simplify()
    print('make STRING network done')

    # Prune start/end lists
    start_list_prot = [each for each in start_list_prot if each in nodes]
    end_list_prot = [each for each in end_list_prot if each in nodes]

    ppi_deg_pair_list = list(itertools.product(start_list_prot, end_list_prot))
    # TODO: The list above holds twice the pairs we need, ie, both: A -> B and B -> A. Keep unique (undirected) pairs.
    #  There's definitely a quicker way but this works good enough for now...
    ppi_deg_pair_list_clean = []  # len(ppi_deg_pair_list_clean) = [len(ppi_deg_pair_list) + len(start_list_prot)] / 2
    for each in ppi_deg_pair_list:
        if each not in ppi_deg_pair_list_clean \
                and (each[1], each[0]) not in ppi_deg_pair_list_clean:
            ppi_deg_pair_list_clean.append(each)
    print('start/end protein combinations list done')

    shortest_paths_result = []
    for i in range(len(ppi_deg_pair_list_clean)):
        if np.mod(i+1, 100) == 0:
            print(f'Pair #{i+1} of #{len(ppi_deg_pair_list_clean)}')
        shortest_paths_result += make_shortest_path_sources_to_targets(ppi_deg_pair_list_clean[i], graph)
    shortest_paths_result = list(set(shortest_paths_result))

    # make new graph with hidden layers
    nodes_sp = list(set(
        [each for edge in shortest_paths_result for each in edge]
    ))
    G_prot = ig.Graph(directed=False)
    G_prot.add_vertices(nodes_sp)
    G_prot.add_edges(shortest_paths_result)
    G_prot = G_prot.simplify()

    # save new network with protein IDs as nodes
    save_obj(G_prot, address_network_proteins)

    #%% MAIN
    # TRANSLATION PROTEIN TO GENE

    print('\nRenaming network nodes (Protein ID to gene symbol).')

    # translate  protein IDs using local dictionary from JENSEN lab
    # open dictionary (ensembl protein IDs to gene symbols, JENSEN lab)
    network_dict = pd.read_csv(address_dictionary_protein_to_gene,  delimiter=',', header=None)
    # Keep unique only.
    network_dict = dict(sorted(network_dict.values.tolist()))

    # rename nodes in network
    G_gene = deepcopy(G_prot)

    failed_to_find = False
    for each in G_gene.vs:
        if network_dict.get(each['name']) is None:
            failed_to_find = True
        each['name'] = network_dict.get(each['name'], each['name'])   # Keep original ENSEMBL ID if not found.

    if failed_to_find:
        print('Not all nodes renamed, or more than one protein corresponds to the same gene.')

    # also translate the original input lists of proteins for later use
    start_list_gene = []
    for protein in start_list_prot:
        if protein in network_dict:
            start_list_gene += [network_dict[protein]]
        else:
            start_list_gene += [protein]
    end_list_gene = []
    for protein in end_list_prot:
        if protein in network_dict:
            end_list_gene += [network_dict[protein]]
        else:
            end_list_gene += [protein]

    # save list of nodes and network with gene names
    network_nodes_genes = pd.DataFrame(G_gene.vs['name'])
    network_nodes_genes.to_csv(address_network_nodes_genes, header=None, index=False)
    network_nodes_proteins = pd.DataFrame(G_prot.vs['name'])
    network_nodes_proteins.to_csv(address_network_nodes_proteins, header=None, index=False)

    edgelist_tmp = G_gene.get_edgelist()
    edgelist = []
    for each in edgelist_tmp:
        edgelist.append((
            G_gene.vs[each[0]]['name'],
            G_gene.vs[each[1]]['name'],
        ))
    network_edges_genes = pd.DataFrame(edgelist)
    network_edges_genes.to_csv(address_network_edges_genes, header=None, index=False)

    edgelist_tmp = G_prot.get_edgelist()
    edgelist = []
    for each in edgelist_tmp:
        edgelist.append((
            G_prot.vs[each[0]]['name'],
            G_prot.vs[each[1]]['name'],
        ))
    network_edges_proteins = pd.DataFrame(edgelist)
    network_edges_proteins.to_csv(address_network_edges_proteins, header=None, index=False)
    save_obj(G_gene, address_network_genes)

    #%% MAIN
    # PROTEIN LOCALISATION

    print('Finding proteins that are localised in nucleus/membrane.')
    # open pre-edited lists of genes in membrane and in nucleus (based on COMPARTMENT database on 30th July 2021)
    df_comp_mem = pd.read_csv(address_compartment_membrane_genes, delimiter='\t', header=None)
    df_comp_nuc = pd.read_csv(address_compartment_nucleus_genes, delimiter='\t', header=None)
    list_comp_mem = df_comp_mem.iloc[:,0].to_list()
    list_comp_nuc = df_comp_nuc.iloc[:,0].to_list()

    # convert network nodes to list
    nodes = list(G_gene.vs['name'])

    # make list of genes (=nodes) that are in membrane/nucleus
    network_membrane_gene_names = []
    network_nucleus_gene_names = []

    for node in nodes:
        if node in list_comp_mem:
            network_membrane_gene_names += [node]
        if node in list_comp_nuc:
            network_nucleus_gene_names += [node]

    # save lists as txt
    network_membrane_gene_names_list = pd.DataFrame(list(network_membrane_gene_names))
    network_nucleus_gene_names_list = pd.DataFrame(list(network_nucleus_gene_names))
    network_membrane_gene_names_list.to_csv(address_network_membrane_genes, header=None, index=False)
    network_nucleus_gene_names_list.to_csv(address_network_nucleus_genes, header=None, index=False)

    membrane_list = list(network_membrane_gene_names)
    nucleus_list = list(network_nucleus_gene_names)

    #%% MAIN
    # ANALYSE NETWORK

    # give lists, network, output address to function
    print("Analyzing pathway network using multiple centrality methods.")
    calc_network_centrality_RWR(G_gene, membrane_list, nucleus_list, address_centrality_results)

    #%% MAIN
    # VISUALISE NETWORK

    print("Generating circular plot of network.")
    # Transform to NX graphs now, as the computationally heavy part is over.
    G_gene_nx = G_gene.to_networkx(create_using=nx.Graph, vertex_attr_hashable='name')
    G_prot_nx = G_prot.to_networkx(create_using=nx.Graph, vertex_attr_hashable='name')

    make_circular_plot_network(list(G_prot_nx.nodes), G_prot_nx, address_network_plot_proteins, start_list_prot)
    make_circular_plot_network(list(G_gene_nx.nodes), G_gene_nx, address_network_plot_genes, start_list_gene)
    print("Program finished.")
