#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri July 27 14:12:50 2024

@author: Zongtai Wu based on Han et al., 2021, original code created by woochanghwang at 18/11/2020

Running time: depends network size an number of permutations

"""
import copy

#import pickle5 as pickle
import networkx as nx
import random
import time
import multiprocessing as mp
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import igraph as ig

MAX_NUMBER_OF_TRIAL = 10
FINITE_INFINITY = 999999

#%%
### make sure the working folder is the src folder (in Spyder top right) ###

""" variables """
### ENTER NUMBER OF PERMUTATIONS ###
perm_no = 1000
### ENTER NUMBER OF PROCESSES FOR POOL PROCESSING (KEEP BELOW 5 ON A NORMAL LAPTOP) ###
pro_no = 8
### ENTER THRESHOLD P-VALUE TO SELECT KEY PROTEINS ###
p_val = 0.01

""" input addresses """
### ENTER PATH/NAME INPUT: NETWORK ###
address_network = "./data/network_genes"
### ENTER PATH/NAME INPUT: CENTRALITY RESULTS ###
address_centrality = "./data/Centrality_RWR_results.csv"
### ENTER PATH/NAME INPUT: LISTS OF GENES IN MEMBRANE/NUCLEUS ###
address_membrane_nodes = './data/network_membrane_gene_names.txt'
address_nucleus_nodes = './data/network_nucleus_gene_names.txt'

""" output addresses """
### ENTER PATH/NAME OUTPUT: PERMUTATION RESULTS AS DICTIONARY .pkl ###
address_permutation_results = "./data/network_permutation3" # dont add '.pkl'
### ENTER PATH/NAME OUTPUT: p-VALUE RESULTS AS .csv ###
address_centrality_pvalue = "./data/Centrality_RWR_result_pvalue.csv"
### ENTER PATH/NAME OUTPUT: p-VALUE RESULTS AS .txt ###
address_key_proteins = "./data/key_proteins.txt"
### ENTER PATH/NAME OUTPUT: EDGE LIST OF KEY PROTEINS IN NETWORK .txt ###
address_network_key_proteins = './data/network_edges_key_proteins.txt'
### ENTER PATH/NAME OUTPUT: CIRCULAR PLOT OF KEY PROTEIN NETWORK .tiff ###
address_network_key_proteins_plot = './visualisation/network_key_proteins.tiff'

""" randomization_types:
    - "random": add same number of edges randomly between nodes of original graph
    - "preserve_topology": keep edges, shuffle nodes of original graph
    - "preserve_topology_and_node_degree": keep edges, shuffle nodes of original graph with the nodes of same degree
    - "preserve_degree_distribution": remove an edge between two random nodes with degrees k, l then add to two nodes 
        with degrees k-1 & l-1, then shuffle nodes
    - "preserve_degree_distribution_and_node_degree": remove 2 random edges between a-b and c-d where degree(a)=degree(c) 
        and degree(b)=degree(d) then add 2 edges between a-d and b-c, then shuffle nodes with the same degree
    - "erdos_renyi": creates a graph where edges are redistributed based on erdos renyi random model.
    - "barabasi_albert": creates a graph where edges are redistributed based on barabasi albert model (preferential attachment).
    
    ### To edit the randomization type, scroll to line 417 and enter type. Currently type is set to "preserve_degree_distribution", which is used in the paper.
"""
#%%

def load_obj(file_addr):
    with open(file_addr+ '.pkl', 'rb') as f:
        return pickle.load(f)

def save_obj(obj, file_addr ):
    with open(file_addr + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def randomize_graph(graph, randomization_type, allow_self_edges=False):
    """
    Creates a random network from given network as a networkx graph
    randomization_type:
        - "random": add same number of edges randomly between nodes of original graph
        - "preserve_topology": keep edges, shuffle nodes of original graph
        - "preserve_topology_and_node_degree": keep edges, shuffle nodes of original graph with the nodes of same degree
        - "preserve_degree_distribution": remove an edge between two random nodes with degrees k, l then add to two nodes with degrees k-1 & l-1, then shuffle nodes
        - "preserve_degree_distribution_and_node_degree": remove 2 random edges between a-b and c-d where degree(a)=degree(c) and degree(b)=degree(d) then add 2 edges between a-d and b-c, then shuffle nodes with the same degree
    - "erdos_renyi": creates a graph where edges are redistributed based on erdos renyi random model.
    - "barabasi_albert": creates a graph where edges are redistributed based on barabasi albert model (preferential attachment).
    """

    debug = False

    n_node = graph.number_of_nodes()
    n_edge = graph.number_of_edges()

    if randomization_type == "same_degree_sequence":
        # Takes ages to find a suitable conformation for large graphs
        sequence = graph.degree().values()
        new_graph = None
        while new_graph is None:
            new_graph = nx.random_degree_sequence_graph(sequence)
        return new_graph

    if randomization_type == "graph_tool_correlated":
        try:
            import graph_tool
        except:
            raise ValueError("Graph tool package not installed")
            return
        new_graph = graph.copy()
        graph_tool.generation.random_rewire(new_graph, model='uncorrelated', n_iter=1, edge_sweep=True,
                                            parallel_edges=False, self_loops=False, vertex_corr=None,
                                            block_membership=None, alias=True, cache_probs=True, persist=False,
                                            ret_fail=False, verbose=False)
        return new_graph

    if randomization_type == "erdos_renyi":
        # raise Exception("Work in progress")
        p = float(2 * n_edge) / (n_node * n_node - 2 * n_node)
        # Chooses each of the possible [n(n-1)]/2 edges with probability p
        new_graph = nx.erdos_renyi_graph(n_node, p)
        mapping = dict(zip(new_graph.nodes(), graph.nodes()))
        new_graph = nx.relabel_nodes(new_graph, mapping)
        available_edges = graph.edges()

        # Map graph from random model to new graph
        for edge in new_graph.edges():
            if len(available_edges) > 0:
                edge_org = available_edges.pop()
                if debug:
                    print ("From random:", (edge[0], edge[1]))
                new_graph.add_edge(edge[0], edge[1], graph.get_edge_data(edge_org[0], edge_org[1]))
            # If the random model added too many edges
            else:
                if debug:
                    print ("Removing:", edge)
                new_graph.remove_edge(edge[0], edge[1])

        # If the random model failed to add enough edges
        nodes = new_graph.nodes()
        for edge_org in available_edges:
            source_id = random.choice(nodes)
            target_id = random.choice(nodes)
            while new_graph.has_edge(source_id, target_id) or (not allow_self_edges and source_id == target_id):
                source_id = random.choice(nodes)
                target_id = random.choice(nodes)
            if debug:
                print ("Adding:", (source_id, target_id))
            new_graph.add_edge(source_id, target_id, graph.get_edge_data(edge_org[0], edge_org[1]))
        return new_graph

    if randomization_type == "barabasi_albert":
        # raise Exception("Work in progress")
        if n_edge >= n_node:
            # A graph of n nodes is grown by attaching new nodes each with m edges that are preferentially attached to existing nodes with high degree
            new_graph = nx.barabasi_albert_graph(n_node, n_edge / n_node)
            mapping = dict(zip(new_graph.nodes(), graph.nodes()))
            new_graph = nx.relabel_nodes(new_graph, mapping)
        else:
            new_graph = nx.create_empty_copy(graph)

        available_edges = graph.edges()
        degree_map = dict(nx.degree(new_graph))
        nodes = new_graph.nodes()

        # Map graph from random model to new graph
        for edge in new_graph.edges():
            if len(available_edges) > 0:
                edge_org = available_edges.pop()
                if debug:
                    print ("From random:", (edge[0], edge[1]))
                new_graph.add_edge(edge[0], edge[1], graph.get_edge_data(edge_org[0], edge_org[1]))
            # If the random model added too many edges
            else:
                nodes_to_select = [id for id, d in degree_map.items() for j in range(d + 1)]
                source_id = random.choice(nodes())
                target_id = random.choice(nodes_to_select)
                if debug:
                    print ("Removing:", (source_id, target_id))
                new_graph.remove_edge(source_id, target_id)
                degree_map[source_id] -= 1
                degree_map[target_id] -= 1

            # If the random model failed to add enough edges
        for edge_org in available_edges:
            nodes_to_select = [id for id, d in degree_map.items() for j in range(d + 1)]
            source_id = random.choice(nodes)
            target_id = random.choice(nodes_to_select)
            while new_graph.has_edge(source_id, target_id) or (not allow_self_edges and source_id == target_id):
                source_id = random.choice(nodes)
                target_id = random.choice(nodes_to_select)
            if debug:
                print ("Adding:", (source_id, target_id))
            new_graph.add_edge(source_id, target_id, graph.get_edge_data(edge_org[0], edge_org[1]))
            degree_map[source_id] += 1
            degree_map[target_id] += 1

        return new_graph

    new_graph = nx.create_empty_copy(graph)
    # new_graph.add_nodes_from(graph.nodes())

    if randomization_type == "random":
        nodes = new_graph.nodes()
        for edge in graph.edges():
            source_id = random.choice(nodes)
            target_id = random.choice(nodes)
            while new_graph.has_edge(source_id, target_id) or (not allow_self_edges and source_id == target_id):
                source_id = random.choice(nodes)
                target_id = random.choice(nodes)
            new_graph.add_edge(source_id, target_id, graph.get_edge_data(edge[0], edge[1]))

    elif randomization_type == "preserve_topology":  # shuffle_nodes
        nodes = list(graph.nodes())
        random_nodes = list(graph.nodes())
        random.shuffle(random_nodes)
        equivalences = dict([(nodes[i], random_nodes[i]) for i in range(len(nodes))])
        new_graph.add_edges_from([(equivalences[current_edge[0]], equivalences[current_edge[1]],
                                   graph.get_edge_data(current_edge[0], current_edge[1])) for current_edge in
                                  graph.edges()])

    elif randomization_type == "preserve_topology_and_node_degree":  # shuffle_nodes_within_same_degree
        nodes_by_degree = dict((degree, []) for u, degree in graph.degree())  # .values()
        graph_degree = dict(graph.degree())
        [nodes_by_degree[graph_degree[node]].append(node) for node in graph_degree]
        equivalences = {}
        for current_degree in nodes_by_degree.keys():
            nodes = nodes_by_degree[current_degree]
            random_nodes = list(nodes)
            random.shuffle(random_nodes)
            equivalences.update(dict([(nodes[i], random_nodes[i]) for i in range(len(nodes))]))
        new_graph.add_edges_from([(equivalences[current_edge[0]], equivalences[current_edge[1]],
                                   graph.get_edge_data(current_edge[0], current_edge[1])) for current_edge in
                                  graph.edges()])

    elif randomization_type == "preserve_degree_distribution":
        for current_node1, current_node2 in graph.edges():
            new_graph.add_edge(current_node1, current_node2, *graph.get_edge_data(current_node1, current_node2))
        # max_degree = sorted(zip(*list(graph.degree()))[1])[-1]  # .values()
        degree_sequence = sorted([d for n, d in graph.degree()], reverse=True)
        max_degree = max(degree_sequence)
        # print "Degree sequence", degree_sequence
        nodes_by_degree = dict((degree, {}) for degree in range(max_degree + 1))
        graph_degree = dict(graph.degree())
        [nodes_by_degree[graph_degree[node]].setdefault(node) for node in graph_degree]
        n_perturbation = random.randint(int(2 * n_edge / 3), n_edge)  # Perturb at least 66% of the edges
        for i in range(n_perturbation):
            n_trial = 0
            while True:
                n_trial += 1
                if n_trial > MAX_NUMBER_OF_TRIAL:
                    if debug:
                        print("Warning: Max number of trials exceeded in perturbation ", i)
                    break
                source_id = random.choice(list(new_graph.nodes()))
                source_degree = new_graph.degree(source_id)
                while source_degree < 1:
                    source_id = random.choice(list(new_graph.nodes()))
                    source_degree = new_graph.degree(source_id)
                target_id = random.choice(list(new_graph.neighbors(source_id)))
                target_degree = new_graph.degree(target_id)
                del nodes_by_degree[source_degree][source_id]
                nodes_by_degree[source_degree - 1].setdefault(source_id)
                if target_id == source_id:
                    target_degree -= 1
                del nodes_by_degree[target_degree][target_id]
                nodes_by_degree[target_degree - 1].setdefault(target_id)
                ## not very important to check for cases where new_source = source (v.v. for targets)
                new_target_id = random.choice(list(nodes_by_degree[target_degree - 1].keys()))
                if source_id == target_id:
                    new_source_id = new_target_id
                else:
                    new_source_id = random.choice(list(nodes_by_degree[source_degree - 1].keys()))
                if debug:
                    print(source_id, target_id, " / ", new_source_id, new_target_id)
                    print(source_degree, target_degree)
                ## check if going to add an existing edge or self edge
                if new_graph.has_edge(new_source_id, new_target_id) or (
                        not allow_self_edges and new_source_id == new_target_id):
                    del nodes_by_degree[target_degree - 1][target_id]
                    nodes_by_degree[target_degree].setdefault(target_id)
                    del nodes_by_degree[source_degree - 1][source_id]
                    nodes_by_degree[source_degree].setdefault(source_id)
                    continue
                if debug:
                    print("rm %s %s" % (source_id, target_id))
                edge_data = new_graph.get_edge_data(source_id, target_id)
                new_graph.remove_edge(source_id, target_id)
                if debug:
                    print("add %s %s" % (new_source_id, new_target_id))
                new_graph.add_edge(new_source_id, new_target_id, *edge_data)
                del nodes_by_degree[target_degree - 1][new_target_id]
                nodes_by_degree[target_degree].setdefault(new_target_id)
                if new_source_id == new_target_id and source_id != target_id:
                    source_degree += 1
                del nodes_by_degree[source_degree - 1][new_source_id]
                nodes_by_degree[source_degree].setdefault(new_source_id)
                break
        randomize_graph(new_graph, "preserve_topology")

    elif randomization_type == "preserve_degree_distribution_and_node_degree":
        ## add edges as well
        for current_node1, current_node2 in graph.edges():
            edge_data = graph.get_edge_data(current_node1, current_node2)
            # new_graph.add_edge(current_node1, current_node2, graph.get_edge_data(current_node1, current_node2))
            new_graph.add_edge(current_node1, current_node2)

        nodes_by_degree = dict((degree, {}) for u, degree in graph.degree())
        graph_degree = dict(graph.degree())
        [nodes_by_degree[graph_degree[node]].setdefault(node) for node in graph_degree]

        # if n_edge < MIN_NUMBER_OF_PERTURBATION:
        #    n_perturbation = random.randint(1, n_edge)
        # else:
        #    n_perturbation = random.randint(MIN_NUMBER_OF_PERTURBATION, n_edge)

        n_perturbation = random.randint(int(n_edge / 2), n_edge)
        for i in range(n_perturbation):
            # nodes =  list(new_graph.nodes())
            source_id = random.choice(list(new_graph.nodes()))
            # source_id = random.choice(nodes)
            # print("sourc_id", source_id)
            source_degree = new_graph.degree(source_id)
            ## find a node for which another node with the same degree exists
            # available_neighbors = []
            n_trial = 0
            while True:  # (len(nodes_by_degree[source_degree]) < 2 or len(available_neighbors) < 1):
                n_trial += 1
                if n_trial > MAX_NUMBER_OF_TRIAL:
                    if debug:
                        print("Warning: Max number of trials exceeded in perturbation ", i)
                    break
                source_id = random.choice(list(new_graph.nodes()))
                source_degree = new_graph.degree(source_id)
                if len(nodes_by_degree[source_degree]) < 2:
                    continue
                available_neighbors = []
                ## find a neighbor for which another node with the same degree exists
                # for neighbor_id in new_graph.neighbors_iter(source_id):   Networkx 1.x
                for neighbor_id in new_graph.neighbors(source_id):
                    if source_degree == new_graph.degree(neighbor_id):
                        if len(nodes_by_degree[new_graph.degree(neighbor_id)]) > 2:
                            available_neighbors.append(neighbor_id)
                    else:
                        if len(nodes_by_degree[new_graph.degree(neighbor_id)]) > 1:
                            available_neighbors.append(neighbor_id)
                if len(available_neighbors) < 1:
                    continue
                target_id = random.choice(available_neighbors)
                target_degree = new_graph.degree(target_id)
                ## select a new source node with different id
                n_trial2 = 0
                inner_break = False
                while True:
                    n_trial2 += 1
                    if n_trial2 > MAX_NUMBER_OF_TRIAL:
                        if debug:
                            print("Warning: Max number of trials exceeded in perturbation ", i)
                        inner_break = True
                        break
                    new_source_id = random.choice(list(nodes_by_degree[source_degree].keys()))
                    while new_source_id == source_id:
                        new_source_id = random.choice(list(nodes_by_degree[source_degree].keys()))
                    new_available_neighbors = []
                    ## find a neighbor as new target node for which id is different from target and has an id equivalent to target
                    for neighbor_id in new_graph.neighbors(new_source_id):
                        if target_degree == new_graph.degree(neighbor_id):
                            new_available_neighbors.append(neighbor_id)
                    if len(new_available_neighbors) < 1:
                        continue
                    new_target_id = random.choice(list(new_available_neighbors))
                    if len(new_available_neighbors) > 1:
                        while new_target_id == target_id:
                            new_target_id = random.choice(list(new_available_neighbors))
                            # print new_available_neighbors, new_target_id
                    else:
                        new_target_id = new_available_neighbors[0]
                    break
                if inner_break:
                    break
                if debug:
                    print(source_id, target_id, " / ", new_source_id, new_target_id)
                if source_id == new_target_id or new_source_id == target_id:
                    continue
                if new_graph.has_edge(source_id, new_target_id) or new_graph.has_edge(new_source_id, target_id):
                    continue
                if debug:
                    print("rm %d %d" % (source_id, target_id))
                    print("rm %d %d" % (new_source_id, new_target_id))
                edge_data_1 = new_graph.get_edge_data(source_id, target_id)
                edge_data_2 = new_graph.get_edge_data(new_source_id, new_target_id)
                new_graph.remove_edge(source_id, target_id)
                new_graph.remove_edge(new_source_id, new_target_id)
                if debug:
                    print("add %d %d" % (source_id, new_target_id))
                    print("add %d %d" % (new_source_id, target_id))
                # new_graph.add_edge(source_id, new_target_id, edge_data_1)
                # new_graph.add_edge(new_source_id, target_id, edge_data_2)
                new_graph.add_edge(source_id, new_target_id)
                new_graph.add_edge(new_source_id, target_id)

    else:
        raise Exception("Unknown randomization type %s" % randomization_type)

    return new_graph

def centrality_dict_to_list(centriality_result_dict):

    centrality_result_list = []
    for key,value in centriality_result_dict.items():
        centrality_result_list.append([key, value])

    return centrality_result_list

def analysis_centrality_RWR(input):
    original_G: ig.Graph
    original_G, process_num, no_membrane_nodes, no_nucleus_nodes = input
    start_time = time.time()

    # Rewire network
    num_edges = original_G.ecount()

    # Make copy, permutation happens in-place!
    network = copy.deepcopy(original_G)
    network.rewire(n=int(num_edges * 1.25), mode='simple')

    network_nodes = list(network.vs['name'])
    random_start_genes = random.choices(network_nodes, k=no_membrane_nodes) # in paper: number of DIPs
    random_target_genes = random.choices(network_nodes, k=no_nucleus_nodes) # in paper: number of DEPs

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
        s = 1.0 / (len(network_nodes) - 1.0)  # Normalisation
        degree_centrality = np.array(network.degree(network_nodes, loops=False)) * s
        degree_centrality_dict = dict(zip(network_nodes, degree_centrality))
    except:
        degree_centrality_dict = dict()
        for node in network_nodes:
            degree_centrality_dict[node] = 0.0

    # betweeness centrality
    try:
        s = 2 / ((len(network_nodes) - 1) * (len(network_nodes) - 2))
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
            sources=random_start_genes,
            targets=random_target_genes,
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
    PRsub_score = network.personalized_pagerank(directed=False, reset_vertices=random_start_genes, weights='weight')
    PRsub_score = dict(zip(network_nodes, PRsub_score))

    print("{} times training Runtime: {:2f} Seconds".format(process_num, ((time.time() - start_time) / 3600)))
    #print(eigenvector_centrality_dict,between_centrality_dict)
    return {'eigen' : eigenvector_centrality_dict, 'between':between_centrality_dict, 
            'betweensub': betweensub_centrality_dict, 'degree': degree_centrality_dict, 
            'rwr': PR_score, 'rwrsub': PRsub_score}


def get_permutation_result_df(permutation_result, key_method):
    # 'eigen': eigenvector_centrality_dict
    # 'between': between_centrality_dict,
    # 'degree': degree_centrality_dict,
    # 'rwr': rwr_dict

    score_list = []
    for permutation_process in permutation_result:
        network_analisys_score_dict = permutation_process[key_method]
        #print(network_analisys_score_dict)
        score_list += list(network_analisys_score_dict.values())
        #print(score_list)
        
    permutation_df= pd.DataFrame(score_list, columns=[key_method])

    return permutation_df


def get_key_proteins(centrality_pvalue_addr,key_protein_addr, p_val):

    # dont include 'betweensub_pvalue', 'rwrsub_pvalue' here:
    p_value_list = ['eigen_pvalue', 'degree_pvalue', 'between_pvalue', 'rwr_pvalue']

    df_centrality_pvalue = pd.read_csv(centrality_pvalue_addr, index_col=0)
    key_protein_list = []

    for na_method in p_value_list:
        key_genes = df_centrality_pvalue[df_centrality_pvalue[na_method]<= p_val].index.to_list()
        # print(na_method, len(key_genes))
        key_protein_list += key_genes

    key_protein_list = list(set(key_protein_list))
    # print(len(key_protein_list))

    key_protein_df = pd.DataFrame(key_protein_list, columns=['Gene'])
    key_protein_df.to_csv(key_protein_addr, index=False)
    return key_protein_list

def make_circular_plot_network(key_protein_list, G, address_network_plot):
    '''
    Function: make a circular plot of the network, save as tiff in visualisation folder
    :input: network, list of nodes
    :return: none
    '''
    G_gene = nx.Graph()
    G_gene.add_edges_from(G.edges(key_protein_list))
    nodes = list(G_gene.nodes)
            
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
                                          
    nx.draw_circular(G_gene,with_labels=False, font_size=18, font_color ='k', edge_color = 'grey')
    plt.savefig(address_network_key_proteins_plot)

def main():
    '''
    :return:
    '''
    # load network
    original_G = load_obj(address_network)
    # print(list(original_G.es))

    # make list with values 0 to permutation number (e.g. 0 to 999)
    permutation_num = list(range(perm_no))
    # make 1000 network duplicates saved in a list
    graph_list = [original_G] * len(permutation_num)
    
    no_membrane_nodes = [len(pd.read_csv(address_membrane_nodes, sep = ',', header = None))] * len(permutation_num)
    no_nucleus_nodes = [len(pd.read_csv(address_nucleus_nodes, sep = ',', header = None))] * len(permutation_num)

    pool = mp.Pool(processes=pro_no)
    print("Start permutation test")
    # fetch multiple network parameters for each permutation
    permutation_result_all = pool.map(analysis_centrality_RWR, zip(graph_list, permutation_num, no_membrane_nodes, no_nucleus_nodes))
    # permutation_result_all: [{parameter: {'gene': number, 'gene': number,..}, parameter: {'gene': number,..}}, {same for second permutation}, {third}, etc}]
    pool.close()
    # save as dictionary as -pkl file
    save_obj(permutation_result_all, address_permutation_results)
    print("Finish permutation test")

    ##########
    # Calculate empirical p-value from permutation test
    ##########
    
    print("Calculate p-value")
    # fetch each parameter as dataframe
    permutation_eigen_result_all = get_permutation_result_df(permutation_result_all, 'eigen')
    permutation_degree_result_all = get_permutation_result_df(permutation_result_all, 'degree')
    permutation_between_result_all = get_permutation_result_df(permutation_result_all, 'between')
    permutation_betweensub_result_all = get_permutation_result_df(permutation_result_all, 'betweensub')
    permutation_rwr_result_all = get_permutation_result_df(permutation_result_all, 'rwr')
    permutation_rwrsub_result_all = get_permutation_result_df(permutation_result_all, 'rwrsub')
    
    df_centrality_result = pd.read_csv(address_centrality, index_col=0)    
    
    # merge into one dataframe
    df_final_result = pd.concat(
        [permutation_eigen_result_all, permutation_degree_result_all, permutation_between_result_all,
         permutation_betweensub_result_all, permutation_rwr_result_all, permutation_rwrsub_result_all], axis=1)

    # calculate p-value
    # create new column 'eigen_pvalue', 
    df_centrality_result['eigen_pvalue'] = df_centrality_result.apply(
        lambda row: len(df_final_result[df_final_result.eigen > row['Eigen']]) / len(
            df_final_result.index),axis=1)
    
    df_centrality_result['degree_pvalue'] = df_centrality_result.apply(
        lambda row: len(df_final_result[df_final_result.degree > row['Degree']]) / len(
            df_final_result.index), axis=1)

    df_centrality_result['between_pvalue'] = df_centrality_result.apply(
        lambda row: len(df_final_result[df_final_result.between > row['Between']]) / len(
            df_final_result.index), axis=1)
    
    df_centrality_result['betweensub_pvalue'] = df_centrality_result.apply(
        lambda row: len(df_final_result[df_final_result.betweensub > row['Between Sub']]) / len(
        df_final_result.index), axis=1)

    df_centrality_result['rwr_pvalue'] = df_centrality_result.apply(
        lambda row: len(df_final_result[df_final_result.rwr > row['RWR']]) / len(
            df_final_result.index), axis=1)
    
    df_centrality_result['rwrsub_pvalue'] = df_centrality_result.apply(
        lambda row: len(df_final_result[df_final_result.rwrsub > row['RWR Sub']]) / len(
            df_final_result.index), axis=1)

    df_centrality_result.to_csv(address_centrality_pvalue)

    ################
    # Get Key proteins
    ################
    
    print( "Get key proteins, p-value {}".format(p_val))
    key_protein_list = get_key_proteins(address_centrality_pvalue, address_key_proteins, p_val)

    # Turn to networkX graph and continue!
    original_G = original_G.to_networkx(create_using=nx.Graph, vertex_attr_hashable='name')
    pd.DataFrame(list(original_G.edges(key_protein_list))).to_csv(address_network_key_proteins, index = False, header = None)
    make_circular_plot_network(key_protein_list, original_G, address_network_key_proteins_plot)
    print("Program finished.")


if __name__ == '__main__':
    main()
